!=====================================================================
! Program: tov_ml_cli.f90
!          Neutron-star sequences for every EOS in an HDF5 file
!          (replaces tov_ml.f90 of the ML_TOV project)
!
!          ./tov_ml_cli.x -i <EOS file> -o <NS file>
!                         [--rho-min 0.13] [--rho-max 1.02] [--nrho 50]
!                         [--tol 1e-8]
!
!          Input:  HDF5 file from eos_v3_cli.py / eos_v3_cli.x (or eos_v3.py)
!                  eos    (N, 3, npts): energy density (g cm^-3), pressure
!                         (dyn cm^-2), baryon number density (cm^-3)
!                  eospar (N, ncol), optional: EOS parameters
!          The number of EOSs N and of table points npts are read from the
!          file.
!
!          Output: HDF5 file with
!                  ns     (N, 3, nrho): M (Msun), R (km), lambda (1e36 g cm^2 s^2)
!                         -- the dataset of tov_ml.f90, used by the notebooks
!                  k2, I (1e45 g cm^2), beta: (N, nrho)
!                  rhoc   (nrho): central baryon densities (fm^-3)
!                  eospar (N, ncol): copied from the input, if present
!          A star whose integration does not reach the surface is NaN.
!
!          Solver: tovSolve (https://github.com/pkrastev/tovSolve, commit
!          b5bba9e), pseudo-enthalpy formalism (Lindblom 1992): libtov.f90,
!          libtov_h.f90 in this directory (changes marked "esymML").
!=====================================================================
program tov_ml_cli
  use constants
  use tov_control
  use hdf5
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

  ! --- command line ---
  character(len=1024) :: infile, outfile
  real(8)             :: rho_start, rho_end
  integer(4)          :: nsteps

  ! --- input ---
  integer(HID_T)       :: fin, dset_in, space_in, mem_space
  integer(HSIZE_T)     :: dims(3), maxdims(3), pdims(2), pmaxdims(2)
  integer(HSIZE_T)     :: offset(3), count(3)
  integer              :: rank, prank, err
  logical              :: have_par
  integer(8)           :: neos
  integer(4)           :: npts, ncol
  real(8), allocatable :: tab(:,:)          ! (npts, 3): one EOS
  real(8), allocatable :: par(:,:)          ! (ncol, neos)

  ! --- output ---
  real(8), allocatable :: ns(:,:,:)         ! (nsteps, 3, neos): M, R, lambda
  real(8), allocatable :: k2a(:,:), ia(:,:), ba(:,:)  ! (nsteps, neos)
  real(8), allocatable :: rhoc_grid(:)      ! (nsteps), fm^-3

  ! --- work ---
  real(8)    :: rhoc, rho_tmp, rho_h
  real(8)    :: M, R, k2, lambda, i_ns, beta
  integer(8) :: ind, nprint, nfail_eos, nabove
  integer(4) :: i, nfail
  integer(8) :: t0, t1, trate

  call parse_args(infile, outfile, rho_start, rho_end, nsteps, tov_eps)
  call system_clock(t0, trate)

  ! +++ Central density grid (as in tov_ml.f90) +++
  allocate(rhoc_grid(nsteps))
  rho_h   = (rho_end - rho_start) / dfloat( nsteps - 1 )
  rho_tmp = rho_start
  do i = 1, nsteps
     rhoc_grid(i) = rho_tmp
     rho_tmp = rho_tmp + rho_h
  end do

  ! +++ Open input, detect the number of EOSs and table points +++
  call h5open_f(err)
  call h5fopen_f(trim(infile), H5F_ACC_RDONLY_F, fin, err)
  if (err /= 0) call fatal('cannot open input file '//trim(infile))
  call h5dopen_f(fin, 'eos', dset_in, err)
  if (err /= 0) call fatal('no dataset "eos" in '//trim(infile))
  call h5dget_space_f(dset_in, space_in, err)
  call h5sget_simple_extent_ndims_f(space_in, rank, err)
  if (rank /= 3) then
     if (rank == 1) call fatal('dataset "eos" in '//trim(infile)//' is empty (no EOS accepted)')
     call fatal('dataset "eos" in '//trim(infile)//' is not (N, 3, npts)')
  end if
  call h5sget_simple_extent_dims_f(space_in, dims, maxdims, err)
  ! Fortran order: dims = (npts, 3, N)
  npts = int(dims(1), 4)
  neos = dims(3)
  if (dims(2) /= 3 .or. npts < 3 .or. neos < 1) call fatal('dataset "eos" is not (N, 3, npts)')

  ! --- EOS parameters, if present ---
  call h5lexists_f(fin, 'eospar', have_par, err)
  ncol = 0
  if (have_par) then
     call read_eospar()
  end if

  write(6,'(a,i0,a,i0,a,a)') 'Input: ', neos, ' EOSs with ', npts, ' points each, from ', trim(infile)
  write(6,'(a,i0,a,f0.4,a,f0.4,a,es8.1)') 'Stars: ', nsteps, ' central densities ', rho_start, &
       ' - ', rho_end, ' fm^-3, tolerance ', tov_eps
  flush(6)

  allocate(tab(npts, 3))
  allocate(ns(nsteps, 3, neos), k2a(nsteps, neos), ia(nsteps, neos), ba(nsteps, neos))

  ! --- memory space for one EOS ---
  count = (/ dims(1), 3_HSIZE_T, 1_HSIZE_T /)
  call h5screate_simple_f(2, count(1:2), mem_space, err)

  nprint    = max(1_8, neos / 20)
  nfail_eos = 0
  nabove    = 0

  ! +++ Loop over EOSs +++
  do ind = 1, neos
     offset = (/ 0_HSIZE_T, 0_HSIZE_T, int(ind - 1, HSIZE_T) /)
     call h5sselect_hyperslab_f(space_in, H5S_SELECT_SET_F, offset, count, err)
     call h5dread_f(dset_in, H5T_NATIVE_DOUBLE, tab, count(1:2), err, mem_space, space_in)
     if (err /= 0) call fatal('cannot read EOS from '//trim(infile))

     call set_eos(npts, tab(:,1), tab(:,2), tab(:,3))
     if (rho_end*fm3cm3 > tab(npts,3)) nabove = nabove + 1

     ! +++ Loop over central density +++
     nfail = 0
     do i = 1, nsteps
        rhoc = rhoc_grid(i) * fm3cm3
        call solve_tov_h_star(rhoc, M, R, k2, lambda, i_ns, beta)
        if (ieee_is_nan(M)) nfail = nfail + 1
        ns(i, 1, ind) = M
        ns(i, 2, ind) = R
        ns(i, 3, ind) = lambda
        k2a(i, ind)   = k2
        ia(i, ind)    = i_ns
        ba(i, ind)    = beta
     end do
     if (nfail > 0) nfail_eos = nfail_eos + 1

     if (mod(ind, nprint) == 0) then
        call system_clock(t1)
        write(6,'(a,i0,a,i0,a,f0.1,a)') 'EOS ', ind, '/', neos, '  (', dble(t1 - t0)/dble(trate), ' s)'
        flush(6)
     end if
  end do

  call h5sclose_f(mem_space, err)
  call h5sclose_f(space_in, err)
  call h5dclose_f(dset_in, err)
  call h5fclose_f(fin, err)

  call write_output(trim(outfile))
  call h5close_f(err)

  call system_clock(t1)
  if (nabove > 0) write(6,'(a,i0,a)') 'Warning: for ', nabove, &
       ' EOSs --rho-max is above the table; the EOS is clamped at its last point there'
  if (nfail_eos > 0) write(6,'(a,i0,a)') 'Warning: ', nfail_eos, &
       ' EOSs have stars that did not reach the surface (stored as NaN)'
  write(6,'(a,i0,a,a,a,f0.1,a)') 'Done: ', neos, ' NS sequences written to ', trim(outfile), &
       '  (', dble(t1 - t0)/dble(trate), ' s)'

contains

  subroutine read_eospar()
    integer(HID_T) :: dpar, spar
    call h5dopen_f(fin, 'eospar', dpar, err)
    call h5dget_space_f(dpar, spar, err)
    call h5sget_simple_extent_ndims_f(spar, prank, err)
    if (prank /= 2) then
       write(6,'(a)') 'Note: "eospar" is not (N, ncol); not copied'
       have_par = .false.
    else
       call h5sget_simple_extent_dims_f(spar, pdims, pmaxdims, err)
       if (pdims(2) /= dims(3)) then
          write(6,'(a)') 'Note: "eospar" does not have one row per EOS; not copied'
          have_par = .false.
       else
          ncol = int(pdims(1), 4)
          allocate(par(ncol, neos))
          call h5dread_f(dpar, H5T_NATIVE_DOUBLE, par, pdims, err)
       end if
    end if
    call h5sclose_f(spar, err)
    call h5dclose_f(dpar, err)
  end subroutine read_eospar

  ! --- Write output (Fortran (a, b, N) = Python (N, b, a)) ---
  subroutine write_output(fname)
    character(len=*), intent(in) :: fname
    integer(HID_T)   :: fout, dset
    integer(HSIZE_T) :: d3(3), d2(2), d1(1)
    character(len=64) :: cols
    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, fout, err)
    if (err /= 0) call fatal('cannot create output file '//fname)

    d3 = (/ int(nsteps, HSIZE_T), 3_HSIZE_T, int(neos, HSIZE_T) /)
    call write_dset(fout, 'ns', 3, d3, ns, dset)
    call attr_str(dset, 'columns', 'M [Msun], R [km], lambda [1e36 g cm^2 s^2]')
    call attr_str(dset, 'input', trim(infile))
    call attr_str(dset, 'solver', 'tovSolve b5bba9e, pseudo-enthalpy formalism')
    call attr_dbl(dset, 'rho_min', rho_start)
    call attr_dbl(dset, 'rho_max', rho_end)
    call attr_dbl(dset, 'nrho', dble(nsteps))
    call attr_dbl(dset, 'tol', tov_eps)
    call h5dclose_f(dset, err)

    d2 = (/ int(nsteps, HSIZE_T), int(neos, HSIZE_T) /)
    call write_dset(fout, 'k2', 2, d2, k2a, dset);   call h5dclose_f(dset, err)
    call write_dset(fout, 'I', 2, d2, ia, dset)
    call attr_str(dset, 'units', '1e45 g cm^2');     call h5dclose_f(dset, err)
    call write_dset(fout, 'beta', 2, d2, ba, dset);  call h5dclose_f(dset, err)

    d1 = int(nsteps, HSIZE_T)
    call write_dset(fout, 'rhoc', 1, d1, rhoc_grid, dset)
    call attr_str(dset, 'units', 'fm^-3');           call h5dclose_f(dset, err)

    if (have_par) then
       d2 = (/ int(ncol, HSIZE_T), int(neos, HSIZE_T) /)
       call write_dset(fout, 'eospar', 2, d2, par, dset)
       select case (ncol)
       case (2);  cols = 'L, Ksym'
       case (4);  cols = 'J0, L, Ksym, Jsym'
       case default; cols = 'unknown'
       end select
       call attr_str(dset, 'columns', trim(cols))
       call h5dclose_f(dset, err)
    end if

    call h5fclose_f(fout, err)
  end subroutine write_output

  subroutine write_dset(fid, name, rnk, d, buf, dset)
    integer(HID_T),   intent(in)  :: fid
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: rnk
    integer(HSIZE_T), intent(in)  :: d(rnk)
    real(8),          intent(in)  :: buf(*)
    integer(HID_T),   intent(out) :: dset
    integer(HID_T) :: sp
    call h5screate_simple_f(rnk, d, sp, err)
    call h5dcreate_f(fid, name, H5T_IEEE_F64LE, sp, dset, err)
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, buf(1:product(d)), (/ product(d) /), err)
    call h5sclose_f(sp, err)
    if (err /= 0) call fatal('cannot write dataset '//name)
  end subroutine write_dset

  subroutine attr_str(obj, name, val)
    integer(HID_T),   intent(in) :: obj
    character(len=*), intent(in) :: name, val
    integer(HID_T)   :: sp, tp, at
    integer(HSIZE_T) :: d1(1) = (/ 1_HSIZE_T /)
    call h5screate_f(H5S_SCALAR_F, sp, err)
    call h5tcopy_f(H5T_FORTRAN_S1, tp, err)
    call h5tset_size_f(tp, int(len(val), SIZE_T), err)
    call h5acreate_f(obj, name, tp, sp, at, err)
    call h5awrite_f(at, tp, val, d1, err)
    call h5aclose_f(at, err)
    call h5tclose_f(tp, err)
    call h5sclose_f(sp, err)
  end subroutine attr_str

  subroutine attr_dbl(obj, name, val)
    integer(HID_T),   intent(in) :: obj
    character(len=*), intent(in) :: name
    real(8),          intent(in) :: val
    integer(HID_T)   :: sp, at
    integer(HSIZE_T) :: d1(1) = (/ 1_HSIZE_T /)
    call h5screate_f(H5S_SCALAR_F, sp, err)
    call h5acreate_f(obj, name, H5T_IEEE_F64LE, sp, at, err)
    call h5awrite_f(at, H5T_NATIVE_DOUBLE, val, d1, err)
    call h5aclose_f(at, err)
    call h5sclose_f(sp, err)
  end subroutine attr_dbl

  subroutine fatal(msg)
    character(len=*), intent(in) :: msg
    write(6,'(a,a)') 'tov_ml_cli.x: error: ', msg
    stop 1
  end subroutine fatal

  ! --- Command line: -i, -o, --rho-min, --rho-max, --nrho, --tol, -h ---
  subroutine parse_args(infile, outfile, rho_start, rho_end, nsteps, tol)
    character(len=*), intent(out) :: infile, outfile
    real(8),          intent(out) :: rho_start, rho_end, tol
    integer(4),       intent(out) :: nsteps
    character(len=1024) :: arg, key, val
    integer :: ia, na, p, ios
    infile    = 'eos_set_test.h5'
    outfile   = 'ns_set_test.h5'
    rho_start = 0.13d0
    rho_end   = 1.02d0
    nsteps    = 50
    tol       = 1.0d-8
    na = command_argument_count()
    ia = 1
    do while (ia <= na)
       call get_command_argument(ia, arg)
       p = index(arg, '=')
       if (arg(1:2) == '--' .and. p > 0) then
          key = arg(1:p-1)
          val = arg(p+1:)
       else
          key = arg
          val = ''
          if (key /= '-h' .and. key /= '--help') then
             ia = ia + 1
             if (ia > na) call usage_error('argument '//trim(key)//': expected one argument')
             call get_command_argument(ia, val)
          end if
       end if
       ios = 0
       select case (trim(key))
       case ('-h', '--help')
          call usage(.true.)
          stop
       case ('-i', '--input')
          infile = val
       case ('-o', '--output')
          outfile = val
       case ('--rho-min')
          read(val, *, iostat=ios) rho_start
       case ('--rho-max')
          read(val, *, iostat=ios) rho_end
       case ('--nrho')
          read(val, *, iostat=ios) nsteps
       case ('--tol')
          read(val, *, iostat=ios) tol
       case default
          call usage_error('unrecognized arguments: '//trim(arg))
       end select
       if (ios /= 0) call usage_error('argument '//trim(key)//': invalid value: '''//trim(val)//'''')
       ia = ia + 1
    end do
    if (nsteps < 2) call usage_error('--nrho must be at least 2')
    if (.not. (rho_start > 0.0d0 .and. rho_end > rho_start)) &
         call usage_error('need 0 < --rho-min < --rho-max')
    if (.not. (tol > 0.0d0)) call usage_error('--tol must be positive')
  end subroutine parse_args

  subroutine usage(full)
    logical, intent(in) :: full
    write(6,'(a)') 'usage: tov_ml_cli.x [-h] [-i INPUT] [-o OUTPUT] [--rho-min RHO_MIN]'
    write(6,'(a)') '                    [--rho-max RHO_MAX] [--nrho NRHO] [--tol TOL]'
    if (.not. full) return
    write(6,'(a)') ''
    write(6,'(a)') 'Neutron-star sequences (M, R, lambda, k2, I, beta) for every EOS in an HDF5 file'
    write(6,'(a)') ''
    write(6,'(a)') 'options:'
    write(6,'(a)') '  -h, --help            show this help message and exit'
    write(6,'(a)') '  -i INPUT, --input INPUT'
    write(6,'(a)') '                        EOS file, dataset "eos" (N, 3, npts) (default: eos_set_test.h5)'
    write(6,'(a)') '  -o OUTPUT, --output OUTPUT'
    write(6,'(a)') '                        output file (default: ns_set_test.h5)'
    write(6,'(a)') '  --rho-min RHO_MIN     lowest central baryon density in fm^-3 (default: 0.13)'
    write(6,'(a)') '  --rho-max RHO_MAX     highest central baryon density in fm^-3 (default: 1.02)'
    write(6,'(a)') '  --nrho NRHO           number of central densities (default: 50)'
    write(6,'(a)') '  --tol TOL             relative tolerance of the ODE solver (default: 1e-8)'
  end subroutine usage

  subroutine usage_error(msg)
    character(len=*), intent(in) :: msg
    call usage(.false.)
    write(6,'(a,a)') 'tov_ml_cli.x: error: ', msg
    stop 2
  end subroutine usage_error

end program tov_ml_cli
