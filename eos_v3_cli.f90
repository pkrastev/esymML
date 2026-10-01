!=====================================================================
! Program: eos_v3_cli.f90
!          Fortran port of eos_v3.py / eos_v3_cli.py
!          Generate parametric EOSs via Monte-Carlo sampling
!          We vary J0, L, Ksym, and Jsym
!
!          ./eos_v3_cli.x -n <number of instances to try> -o <output file>
!                         [-s <random seed>] [--paper]
!
!          --paper: setup of the Galaxies 2022 paper: J0 = Jsym = 0 (only L
!          and Ksym are sampled), output tables up to 10**0.015 fm^-3, no check
!          on the pressure at the top of the table (skip 4), eospar = [L, Ksym].
!
!          Reproduces eos_v3.py bit for bit (with --paper: eos_v3.py changed
!          in the same way): same random numbers for the
!          same seed, same accepted EOSs, identical tables. The
!          reference is eos_v3.py run with NumPy's AVX-512 math disabled
!          (NPY_DISABLE_CPU_FEATURES), i.e. with NumPy using the C math
!          library. To get there:
!          - pow, log, log10 and strtod are glibc's (looked up at run time,
!            so Intel's math library is never used for them)
!          - random numbers: CPython's Mersenne Twister and random.seed()
!          - NumPy's linspace and interp, SciPy's interp1d(kind='quadratic')
!            and scipy.optimize.root(method='hybr') (MINPACK, in
!            minpack_hybrd.f) are reproduced operation by operation
!          - the band solve uses the same LAPACK (dgbsv) as SciPy
!          - compile with strict floating point: no FMA, no flush-to-zero
!=====================================================================

!=====================================================================
! MODULE EOS3_LIBM: glibc math functions and strtod
!=====================================================================
module eos3_libm
  use iso_c_binding
  implicit none

  abstract interface
     function fun2_t(x, y) bind(C)
       import :: c_double
       real(c_double), value :: x, y
       real(c_double)        :: fun2_t
     end function fun2_t
     function fun1_t(x) bind(C)
       import :: c_double
       real(c_double), value :: x
       real(c_double)        :: fun1_t
     end function fun1_t
     function strtod_t(str, endptr) bind(C)
       import :: c_char, c_ptr, c_double
       character(kind=c_char), dimension(*), intent(in) :: str
       type(c_ptr), value                               :: endptr
       real(c_double)                                   :: strtod_t
     end function strtod_t
  end interface

  interface
     function c_dlopen(filename, flag) bind(C, name='dlopen')
       import :: c_ptr, c_char, c_int
       character(kind=c_char), dimension(*), intent(in) :: filename
       integer(c_int), value                            :: flag
       type(c_ptr)                                      :: c_dlopen
     end function c_dlopen
     function c_dlsym(handle, symbol) bind(C, name='dlsym')
       import :: c_ptr, c_funptr, c_char
       type(c_ptr), value                               :: handle
       character(kind=c_char), dimension(*), intent(in) :: symbol
       type(c_funptr)                                   :: c_dlsym
     end function c_dlsym
  end interface

  procedure(fun2_t),   pointer :: pw       => null()  ! pow(x, y)
  procedure(fun1_t),   pointer :: ln       => null()  ! log(x)
  procedure(fun1_t),   pointer :: lg10     => null()  ! log10(x)
  procedure(strtod_t), pointer :: c_strtod => null()  ! strtod(str, NULL)

contains

  subroutine libm_init()
    integer(c_int), parameter :: RTLD_NOW = 2
    type(c_ptr) :: hm, hc
    hm = c_dlopen('libm.so.6'//c_null_char, RTLD_NOW)
    hc = c_dlopen('libc.so.6'//c_null_char, RTLD_NOW)
    if (.not. c_associated(hm) .or. .not. c_associated(hc)) then
       write(6,*) 'Error: cannot open libm.so.6 / libc.so.6'
       stop 1
    end if
    call c_f_procpointer(c_dlsym(hm, 'pow'//c_null_char),    pw)
    call c_f_procpointer(c_dlsym(hm, 'log'//c_null_char),    ln)
    call c_f_procpointer(c_dlsym(hm, 'log10'//c_null_char),  lg10)
    call c_f_procpointer(c_dlsym(hc, 'strtod'//c_null_char), c_strtod)
    if (.not. (associated(pw) .and. associated(ln) .and. associated(lg10) &
         .and. associated(c_strtod))) then
       write(6,*) 'Error: cannot find pow/log/log10/strtod in glibc'
       stop 1
    end if
  end subroutine libm_init

  ! --- Parse a decimal number (correctly rounded, like Python/NumPy) ---
  function str2dble(s) result(v)
    character(len=*), intent(in) :: s
    real(8) :: v
    v = c_strtod(trim(s)//c_null_char, c_null_ptr)
  end function str2dble

end module eos3_libm

!=====================================================================
! MODULE EOS3_RNG: CPython's random module (Mersenne Twister MT19937)
! random.seed(int), random.random() and random.uniform(a, b)
!=====================================================================
module eos3_rng
  implicit none
  private
  public :: rng_seed, rng_seed_urandom, rng_uniform

  integer,    parameter :: NN = 624, MM = 397
  integer(8), parameter :: MASK32   = int(z'FFFFFFFF', 8)
  integer(8), parameter :: MATRIX_A = int(z'9908B0DF', 8)
  integer(8), parameter :: UPPER    = int(z'80000000', 8)
  integer(8), parameter :: LOWER    = int(z'7FFFFFFF', 8)
  integer(8) :: mt(0:NN-1)
  integer    :: mti = NN + 1

contains

  subroutine init_genrand(s)
    integer(8), intent(in) :: s
    integer :: i
    mt(0) = iand(s, MASK32)
    do i = 1, NN-1
       mt(i) = iand(1812433253_8*ieor(mt(i-1), ishft(mt(i-1), -30)) + int(i, 8), MASK32)
    end do
    mti = NN
  end subroutine init_genrand

  subroutine init_by_array(key, klen)
    integer,    intent(in) :: klen
    integer(8), intent(in) :: key(0:klen-1)
    integer :: i, j, k
    call init_genrand(19650218_8)
    i = 1
    j = 0
    k = max(NN, klen)
    do while (k > 0)
       mt(i) = iand(ieor(mt(i), iand(ieor(mt(i-1), ishft(mt(i-1), -30))*1664525_8, MASK32)) &
                    + key(j) + int(j, 8), MASK32)
       i = i + 1
       j = j + 1
       if (i >= NN) then
          mt(0) = mt(NN-1)
          i = 1
       end if
       if (j >= klen) j = 0
       k = k - 1
    end do
    k = NN - 1
    do while (k > 0)
       mt(i) = iand(ieor(mt(i), iand(ieor(mt(i-1), ishft(mt(i-1), -30))*1566083941_8, MASK32)) &
                    - int(i, 8), MASK32)
       i = i + 1
       if (i >= NN) then
          mt(0) = mt(NN-1)
          i = 1
       end if
       k = k - 1
    end do
    mt(0) = UPPER
  end subroutine init_by_array

  ! --- random.seed(n) for an integer n ---
  subroutine rng_seed(n)
    integer(8), intent(in) :: n
    integer(8) :: a, key(0:1)
    integer    :: klen
    ! abs(n) split into 32-bit chunks, least significant first
    ! (the most negative 64-bit integer is not supported)
    a = abs(n)
    key(0) = iand(a, MASK32)
    key(1) = ishft(a, -32)
    klen = 1
    if (key(1) /= 0) klen = 2
    call init_by_array(key, klen)
  end subroutine rng_seed

  ! --- random.seed(None): seed from /dev/urandom ---
  subroutine rng_seed_urandom()
    integer(4) :: buf(4)
    integer(8) :: key(0:3)
    integer    :: u, i, ios
    open(newunit=u, file='/dev/urandom', access='stream', form='unformatted', &
         action='read', iostat=ios)
    if (ios /= 0) then
       write(6,*) 'Error: cannot read /dev/urandom'
       stop 1
    end if
    read(u) buf
    close(u)
    do i = 0, 3
       key(i) = iand(int(buf(i+1), 8), MASK32)
    end do
    call init_by_array(key, 4)
  end subroutine rng_seed_urandom

  function genrand_uint32() result(y)
    integer(8) :: y
    integer    :: kk
    integer(8) :: mag01(0:1)
    mag01(0) = 0_8
    mag01(1) = MATRIX_A
    if (mti >= NN) then
       if (mti == NN + 1) call init_genrand(5489_8)
       do kk = 0, NN-MM-1
          y = ior(iand(mt(kk), UPPER), iand(mt(kk+1), LOWER))
          mt(kk) = ieor(ieor(mt(kk+MM), ishft(y, -1)), mag01(iand(y, 1_8)))
       end do
       do kk = NN-MM, NN-2
          y = ior(iand(mt(kk), UPPER), iand(mt(kk+1), LOWER))
          mt(kk) = ieor(ieor(mt(kk+(MM-NN)), ishft(y, -1)), mag01(iand(y, 1_8)))
       end do
       y = ior(iand(mt(NN-1), UPPER), iand(mt(0), LOWER))
       mt(NN-1) = ieor(ieor(mt(MM-1), ishft(y, -1)), mag01(iand(y, 1_8)))
       mti = 0
    end if
    y = mt(mti)
    mti = mti + 1
    y = ieor(y, ishft(y, -11))
    y = ieor(y, iand(ishft(y, 7),  int(z'9D2C5680', 8)))
    y = ieor(y, iand(ishft(y, 15), int(z'EFC60000', 8)))
    y = ieor(y, ishft(y, -18))
  end function genrand_uint32

  ! --- random.random(): 53-bit float in [0, 1) ---
  function rng_random() result(r)
    real(8)    :: r
    integer(8) :: a, b
    a = ishft(genrand_uint32(), -5)
    b = ishft(genrand_uint32(), -6)
    r = (dble(a)*67108864.0d0 + dble(b))*(1.0d0/9007199254740992.0d0)
  end function rng_random

  ! --- random.uniform(a, b) = a + (b - a)*random() ---
  function rng_uniform(a, b) result(r)
    real(8), intent(in) :: a, b
    real(8) :: r
    r = a + (b - a)*rng_random()
  end function rng_uniform

end module eos3_rng

!=====================================================================
! MODULE EOS3_NUMPY: numpy.linspace and numpy.interp
!=====================================================================
module eos3_numpy
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

contains

  ! --- np.linspace(start, stop, num) ---
  subroutine linspace(start, stop, num, y)
    real(8), intent(in)  :: start, stop
    integer, intent(in)  :: num
    real(8), intent(out) :: y(num)
    real(8) :: step, tmp
    integer :: i
    step = (stop - start) / dble(num - 1)
    do i = 1, num
       tmp  = dble(i-1)*step
       y(i) = tmp + start
    end do
    y(num) = stop
  end subroutine linspace

  ! --- binary_search_with_guess() from numpy/core/src/multiarray/compiled_base.c
  !     0-based result: -1 if key < arr(1), len if key > arr(len) ---
  function bsearch_guess(key, arr, len, guess_in) result(res)
    real(8), intent(in) :: key
    integer, intent(in) :: len, guess_in
    real(8), intent(in) :: arr(0:len-1)
    integer :: res, imin, imax, imid, guess, i
    integer, parameter :: LIKELY = 8
    imin = 0
    imax = len
    guess = guess_in
    if (key > arr(len-1)) then
       res = len
       return
    else if (key < arr(0)) then
       res = -1
       return
    end if
    if (len <= 4) then
       i = 1
       do while (i < len)
          if (.not. (key >= arr(i))) exit
          i = i + 1
       end do
       res = i - 1
       return
    end if
    if (guess > len - 3) guess = len - 3
    if (guess < 1) guess = 1
    if (key < arr(guess)) then
       if (key < arr(guess-1)) then
          imax = guess - 1
          if (guess > LIKELY) then
             if (key >= arr(guess-LIKELY)) imin = guess - LIKELY
          end if
       else
          res = guess - 1
          return
       end if
    else
       if (key < arr(guess+1)) then
          res = guess
          return
       else
          if (key < arr(guess+2)) then
             res = guess + 1
             return
          else
             imin = guess + 2
             if (guess < len - LIKELY - 1) then
                if (key < arr(guess+LIKELY)) imax = guess + LIKELY
             end if
          end if
       end if
    end if
    do while (imin < imax)
       imid = imin + ishft(imax - imin, -1)
       if (key >= arr(imid)) then
          imin = imid + 1
       else
          imax = imid
       end if
    end do
    res = imin - 1
  end function bsearch_guess

  ! --- np.interp(x, xp, fp) ---
  subroutine interp(nx, x, nxp, xp, fp, res)
    integer, intent(in)  :: nx, nxp
    real(8), intent(in)  :: x(nx), xp(0:nxp-1), fp(0:nxp-1)
    real(8), intent(out) :: res(nx)
    real(8) :: lval, rval, xv, slope, r
    integer :: i, j
    lval = fp(0)
    rval = fp(nxp-1)
    if (nxp == 1) then
       do i = 1, nx
          if (x(i) < xp(0)) then
             res(i) = lval
          else if (x(i) > xp(0)) then
             res(i) = rval
          else
             res(i) = fp(0)
          end if
       end do
       return
    end if
    j = 0
    do i = 1, nx
       xv = x(i)
       if (ieee_is_nan(xv)) then
          res(i) = xv
          cycle
       end if
       j = bsearch_guess(xv, xp, nxp, j)
       if (j == -1) then
          res(i) = lval
       else if (j == nxp) then
          res(i) = rval
       else if (j == nxp - 1) then
          res(i) = fp(j)
       else if (xp(j) == xv) then
          res(i) = fp(j)
       else
          slope = (fp(j+1) - fp(j)) / (xp(j+1) - xp(j))
          r = slope*(xv - xp(j)) + fp(j)
          if (ieee_is_nan(r)) then
             r = slope*(xv - xp(j+1)) + fp(j+1)
             if (ieee_is_nan(r) .and. fp(j) == fp(j+1)) r = fp(j)
          end if
          res(i) = r
       end if
    end do
  end subroutine interp

end module eos3_numpy

!=====================================================================
! MODULE EOS3_SPLINE: scipy.interpolate.interp1d(x, y, kind='quadratic')
! = make_interp_spline(x, y, k=2) + BSpline evaluation
!=====================================================================
module eos3_spline
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  implicit none

  integer, parameter :: KSPL = 2

  type spline_t
     integer              :: n
     real(8), allocatable :: t(:)   ! knots, t(0:n+2)
     real(8), allocatable :: c(:)   ! coefficients, c(0:n-1)
     logical              :: nan_y  ! y contains NaN -> interp1d returns NaN
  end type spline_t

contains

  ! --- find_interval() from scipy/interpolate/_bspl.pyx (0-based) ---
  function find_interval(t, nt, k, xval, prev_l, extrapolate) result(res)
    integer, intent(in) :: nt, k, prev_l
    real(8), intent(in) :: t(0:nt-1), xval
    logical, intent(in) :: extrapolate
    integer :: res, l, n
    real(8) :: tb, te
    n  = nt - k - 1
    tb = t(k)
    te = t(n)
    if (ieee_is_nan(xval)) then
       res = -1
       return
    end if
    if (((xval < tb) .or. (xval > te)) .and. .not. extrapolate) then
       res = -1
       return
    end if
    if (k < prev_l .and. prev_l < n) then
       l = prev_l
    else
       l = k
    end if
    do while (l /= k)
       if (.not. (xval < t(l))) exit
       l = l - 1
    end do
    l = l + 1
    do while (l /= n)
       if (.not. (xval >= t(l))) exit
       l = l + 1
    end do
    res = l - 1
  end function find_interval

  ! --- _deBoor_D() from scipy/interpolate/src/__fitpack.h with m = 0 ---
  subroutine deboor(t, x, k, ell, res)
    integer, intent(in)    :: k, ell
    real(8), intent(in)    :: t(0:*), x
    real(8), intent(inout) :: res(0:2*k+1)   ! h = res(0:k), hh = res(k+1:2k+1)
    integer :: j, n, ind
    real(8) :: xb, xa, w
    res(0) = 1.0d0
    do j = 1, k
       res(k+1:k+j) = res(0:j-1)
       res(0) = 0.0d0
       do n = 1, j
          ind = ell + n
          xb  = t(ind)
          xa  = t(ind - j)
          if (xb == xa) then
             res(n) = 0.0d0
             cycle
          end if
          w = res(k+1 + n-1)/(xb - xa)
          res(n-1) = res(n-1) + w*(xb - x)
          res(n)   = w*(x - xa)
       end do
    end do
  end subroutine deboor

  ! --- interp1d(x, y, kind='quadratic') ---
  subroutine spline_build(n, x, y, s)
    integer,        intent(in)  :: n
    real(8),        intent(in)  :: x(0:n-1), y(0:n-1)
    type(spline_t), intent(out) :: s
    integer, parameter :: kl = KSPL, ku = KSPL, ldab = 2*KSPL + KSPL + 1
    real(8) :: ab(ldab, n), wrk(0:2*KSPL+1)
    integer :: ipiv(n), info, j, a, left, clmn, m
    s%n = n
    allocate(s%t(0:n+2), s%c(0:n-1))
    ! --- knots: k = 2 case of make_interp_spline ---
    s%t(0:2) = x(0)
    do m = 1, n-3
       s%t(m+2) = (x(m+1) + x(m)) / 2.0d0
    end do
    s%t(n:n+2) = x(n-1)
    ! --- y with NaN: interp1d fits ones and returns NaN everywhere ---
    s%nan_y = any(ieee_is_nan(y))
    if (s%nan_y) then
       s%c = 1.0d0
    else
       s%c = y
    end if
    ! --- collocation matrix in LAPACK band storage (_bspl._colloc) ---
    ab = 0.0d0
    left = KSPL
    do j = 0, n-1
       left = find_interval(s%t, n+3, KSPL, x(j), left, .false.)
       call deboor(s%t, x(j), KSPL, left, wrk)
       do a = 0, KSPL
          clmn = left - KSPL + a
          ab(kl + ku + j - clmn + 1, clmn + 1) = wrk(a)
       end do
    end do
    ! --- solve (same LAPACK gbsv as SciPy) ---
    call dgbsv(n, kl, ku, 1, ab, ldab, ipiv, s%c, n, info)
    if (info /= 0) then
       write(6,*) 'Error: collocation matrix is singular, info =', info
       stop 1
    end if
  end subroutine spline_build

  ! --- evaluate the spline at xp(1:m) (BSpline.__call__) ---
  subroutine spline_eval(s, m, xp, out)
    type(spline_t), intent(in)  :: s
    integer,        intent(in)  :: m
    real(8),        intent(in)  :: xp(m)
    real(8),        intent(out) :: out(m)
    real(8) :: work(0:2*KSPL+1), v
    integer :: ip, a, interval
    if (s%nan_y) then
       out = ieee_value(1.0d0, ieee_quiet_nan)
       return
    end if
    interval = KSPL
    do ip = 1, m
       interval = find_interval(s%t, s%n+3, KSPL, xp(ip), interval, .true.)
       if (interval < 0) then
          out(ip) = ieee_value(1.0d0, ieee_quiet_nan)
          cycle
       end if
       call deboor(s%t, xp(ip), KSPL, interval, work)
       v = 0.0d0
       do a = 0, KSPL
          v = v + s%c(interval + a - KSPL)*work(a)
       end do
       out(ip) = v
    end do
  end subroutine spline_eval

  ! --- scipy.misc.derivative(f, x0, dx, n=1, order=5) ---
  subroutine spline_deriv5(s, m, x0, dx, df)
    type(spline_t), intent(in)  :: s
    integer,        intent(in)  :: m
    real(8),        intent(in)  :: x0(m), dx
    real(8),        intent(out) :: df(m)
    integer, parameter :: iw(0:4) = (/ 1, -8, 0, 8, -1 /)
    real(8) :: xk(m), fk(m), w, sh
    integer :: k
    df = 0.0d0
    do k = 0, 4
       w  = dble(iw(k)) / 12.0d0
       sh = dble(k - 2)*dx
       xk = x0 + sh
       call spline_eval(s, m, xk, fk)
       df = df + w*fk
    end do
    df = df / dx
  end subroutine spline_deriv5

end module eos3_spline

!=====================================================================
! MODULE EOS3_PHYS: constants, grids, EOS model and beta equilibrium
!=====================================================================
module eos3_phys
  use eos3_libm
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  ! --- (i) physical / mathematical constants ---
  real(8), parameter :: mn   = 939.5656328d0            ! neutron mass
  real(8), parameter :: mp   = 938.272328d0             ! proton mass
  real(8), parameter :: mm   = 105.658389d0             ! muon mass
  real(8), parameter :: pi   = 3.1415926535897932384d0  ! pi
  real(8), parameter :: rho0 = 0.16d0                   ! density of normal nuclear matter

  ! --- (ii) conversion factors (set in phys_init) ---
  real(8) :: ev, fm3cm3, mevfm3gcm3, dycm2mevfm3, ergcm3mevfm3, hbar_c

  ! --- EOS parameters (sampling ranges) ---
  real(8), parameter :: J0_low   = -800.0d0, J0_high   = 400.0d0
  real(8), parameter :: L_low    =   30.6d0, L_high    =  86.8d0
  real(8), parameter :: Ksym_low = -400.0d0, Ksym_high = 100.0d0
  real(8), parameter :: Jsym_low = -200.0d0, Jsym_high = 800.0d0

  ! --- constants of the beta-equilibrium equations ---
  real(8) :: dm, mm2, hbc2, third, two3, half

  ! --- density grids ---
  integer, parameter :: NG = 500
  real(8) :: rhox(NG)            ! grid for beta-equilibrium
  real(8) :: xx0(NG)             ! grid for pressure / energy density
  integer :: kk                  ! number of xx0 points <= 0.09
  real(8) :: zzx(NG), lzzx(NG)   ! grid for the output tables (and log10)
  real(8) :: rlog(NG)            ! log10 of the grid for the speed of sound
  real(8) :: tx(NG), kfex(NG)    ! 3 pi^2 rho, hbar_c (3 pi^2 rho)^(1/3)

  ! --- current density for the residual functions ---
  real(8) :: g_c4, g_t, g_kfe

  ! --- paper setup (--paper): tables up to 10**0.015 fm^-3, no skip 4 ---
  logical :: paper = .false.

contains

  subroutine phys_init()
    integer :: i
    ev           = 1.60217733d-12
    fm3cm3       = 1.0d+39
    mevfm3gcm3   = (1.0d0/1.7827d0)*1.0d-12
    dycm2mevfm3  = (197.33d0/3.1616d0)*1.0d-35
    ergcm3mevfm3 = (ev*1.0d6)*fm3cm3
    hbar_c       = 197.33d0

    dm    = mn - mp
    mm2   = pw(mm, 2.0d0)
    hbc2  = pw(hbar_c, 2.0d0)
    third = 1.0d0/3.0d0
    two3  = 2.0d0/3.0d0
    half  = 0.5d0

    call linspace_(0.02d0, 1.35d0, rhox)
    call linspace_(0.04d0, 1.34d0, xx0)
    kk = count(xx0 <= 0.09d0)
    if (paper) then
       call log_grid(-14.2d0, 0.015d0, zzx)
    else
       call log_grid(-14.2d0, 1.121d0, zzx)
    end if
    call log_grid(-14.0d0, 0.09d0, rlog)
    do i = 1, NG
       rlog(i) = lg10(rlog(i))
       lzzx(i) = lg10(zzx(i))
       tx(i)   = ((3.0d0*pi)*pi)*rhox(i)
       kfex(i) = hbar_c*pw(((3.0d0*pi)*pi)*rhox(i), third)
    end do
  contains
    subroutine linspace_(a, b, y)
      use eos3_numpy, only: linspace
      real(8), intent(in)  :: a, b
      real(8), intent(out) :: y(NG)
      call linspace(a, b, NG, y)
    end subroutine linspace_
  end subroutine phys_init

  ! --- points equally spaced in log10 between 10**dlo and 10**dhi ---
  subroutine log_grid(dlo, dhi, zzx_out)
    real(8), intent(in)  :: dlo, dhi
    real(8), intent(out) :: zzx_out(NG)
    real(8) :: dstp, zz
    integer :: j
    dstp = (dhi - dlo) / dble(NG - 1)
    do j = 0, NG-1
       zz = dlo + dble(j)*dstp
       zzx_out(j+1) = pw(10.0d0, zz)
    end do
  end subroutine log_grid

  ! --- Energy of symmetric nuclear matter ---
  subroutine e_0(rho, J_0, E)
    real(8), intent(in)  :: rho(NG), J_0
    real(8), intent(out) :: E(NG)
    real(8), parameter :: E_sat = -15.8d0, K_0 = 240.0d0, rho_0 = 0.16d0
    real(8) :: x, a2, a3
    integer :: i
    a2 = (1.0d0/2.0d0)*K_0
    a3 = (1.0d0/6.0d0)*J_0
    do i = 1, NG
       x    = (rho(i) - rho_0) / (3.0d0*rho_0)
       E(i) = E_sat + a2*(x*x) + a3*pw(x, 3.0d0)
    end do
  end subroutine e_0

  ! --- Symmetry energy ---
  subroutine e_sym(rho, L, K_sym, J_sym, E)
    real(8), intent(in)  :: rho(NG), L, K_sym, J_sym
    real(8), intent(out) :: E(NG)
    real(8), parameter :: E_sym0 = 31.7d0, rho_0 = 0.16d0
    real(8) :: x, a2, a3
    integer :: i
    a2 = (1.0d0/2.0d0)*K_sym
    a3 = (1.0d0/6.0d0)*J_sym
    do i = 1, NG
       x    = (rho(i) - rho_0) / (3.0d0*rho_0)
       E(i) = E_sym0 + L*x + a2*(x*x) + a3*pw(x, 3.0d0)
    end do
  end subroutine e_sym

  ! --- Residuals of the beta-equilibrium equations (frac_e / frac_emu) ---
  !     iflag = -1 where eos_v3.py gets a NumPy warning (turned into an
  !     error there): negative base of a fractional power or a non-finite
  !     result. hybrd then stops and the EOS is rejected.
  subroutine fcn_beta(n, x, fvec, iflag)
    integer, intent(in)    :: n
    real(8), intent(inout) :: x(n)
    real(8), intent(out)   :: fvec(n)
    integer, intent(inout) :: iflag
    real(8) :: f1, f2, mue, tx3
    if (n == 3) then
       ! --- electrons only ---
       if (x(3) < 0.0d0) then
          iflag = -1
          return
       end if
       f1 = g_c4*(1.0d0 - 2.0d0*x(1)) + dm - g_kfe*pw(x(3), third)
       if (.not. ieee_is_finite(f1)) then
          iflag = -1
          return
       end if
       fvec(1) = f1
       fvec(2) = x(1) - x(3)
       fvec(3) = 1.0d0 - x(1) - x(2)
    else
       ! --- electrons plus muons ---
       tx3 = g_t*x(4)
       if (x(3) < 0.0d0 .or. tx3 < 0.0d0) then
          iflag = -1
          return
       end if
       mue = g_kfe*pw(x(3), third)
       f1  = g_c4*(1.0d0 - 2.0d0*x(1)) + dm - mue
       f2  = pw(mm2 + hbc2*pw(tx3, two3), half) - mue
       if (.not. ieee_is_finite(f1 + f2)) then
          iflag = -1
          return
       end if
       fvec(1) = f1
       fvec(2) = f2
       fvec(3) = x(1) - x(3) - x(4)
       fvec(4) = 1.0d0 - x(1) - x(2)
    end if
  end subroutine fcn_beta

  ! --- scipy.optimize.root(fun, x, method='hybr') with default options ---
  subroutine root_hybr(n, x, ok)
    integer, intent(in)    :: n
    real(8), intent(inout) :: x(n)
    logical, intent(out)   :: ok
    real(8) :: fvec(n), diag(n), fjac(n,n), r(n*(n+1)/2), qtf(n)
    real(8) :: wa1(n), wa2(n), wa3(n), wa4(n)
    integer :: info, nfev
    external hybrd
    call hybrd(fcn_beta, n, x, fvec, 1.49012d-08, 200*(n+1), n-1, n-1, &
               epsilon(1.0d0), diag, 1, 100.0d0, 0, info, nfev, fjac, n, &
               r, n*(n+1)/2, qtf, wa1, wa2, wa3, wa4)
    ok = (info >= 0)
  end subroutine root_hybr

  ! --- particle fractions at every density ---
  subroutine solve_beta(esym, ypr, yn, ye, ymu, ok)
    real(8), intent(in)  :: esym(NG)
    real(8), intent(out) :: ypr(NG), yn(NG), ye(NG), ymu(NG)
    logical, intent(out) :: ok
    real(8) :: x3(3), x4(4)
    integer :: i3
    ypr = 0.0d0
    yn  = 0.0d0
    ye  = 0.0d0
    ymu = 0.0d0
    do i3 = 1, NG
       g_c4  = 4.0d0*esym(i3)
       g_t   = tx(i3)
       g_kfe = kfex(i3)
       if (rhox(i3) <= 0.16d0) then
          x3 = (/ 1.0d-2, 1.0d0, 1.0d-2 /)
          call root_hybr(3, x3, ok)
          if (.not. ok) return
          ypr(i3) = x3(1)
          yn(i3)  = x3(2)
          ye(i3)  = x3(3)
          ymu(i3) = 0.0d0
       else
          x4 = (/ 1.0d-2, 1.0d0, 1.0d-2, 1.0d-2 /)
          call root_hybr(4, x4, ok)
          if (.not. ok) return
          ypr(i3) = x4(1)
          yn(i3)  = x4(2)
          ye(i3)  = x4(3)
          ymu(i3) = x4(4)
       end if
    end do
    ok = .true.
  end subroutine solve_beta

end module eos3_phys

!=====================================================================
! MODULE EOS3_MODEL: one EOS table from (J0, L, Ksym, Jsym)
!=====================================================================
module eos3_model
  use eos3_libm
  use eos3_numpy
  use eos3_spline
  use eos3_phys
  implicit none

  ! --- status codes (rejection reasons as in eos_v3_cli.py) ---
  integer, parameter :: ST_OK = 0, ST_SKIP1 = 1, ST_SKIP2 = 2, ST_SKIP3 = 3, &
                        ST_SKIP4 = 4, ST_SKIP5 = 5, ST_SOLVER = 6
  character(len=14), parameter :: st_name(0:6) = (/ 'ok            ', &
       'skip 1        ', 'skip 2        ', 'skip 3        ', 'skip 4        ', &
       'skip 5        ', 'solver failure' /)

  ! --- crust (APR EOS) below 0.075 fm^-3 ---
  integer              :: ncr
  real(8), allocatable :: rho_apr(:), p_apr(:), eden_apr(:)

contains

  ! --- np.loadtxt("EOS_P_A18dvUIX.txt"), first three columns ---
  subroutine load_crust(fname)
    character(len=*), intent(in) :: fname
    character(len=4096) :: line
    character(len=64)   :: tok(3)
    real(8), allocatable :: c0(:), c1(:), c2(:)
    integer :: u, ios, nrow, i, jj
    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(6,*) 'Error: ', trim(fname), ' not found.'
       stop 1
    end if
    ! --- count data rows ---
    nrow = 0
    do
       read(u, '(a)', iostat=ios) line
       if (ios /= 0) exit
       if (is_data(line)) nrow = nrow + 1
    end do
    allocate(c0(nrow), c1(nrow), c2(nrow))
    rewind(u)
    i = 0
    do
       read(u, '(a)', iostat=ios) line
       if (ios /= 0) exit
       if (.not. is_data(line)) cycle
       i = i + 1
       call first_tokens(line, tok)
       c0(i) = str2dble(tok(1))
       c1(i) = str2dble(tok(2))
       c2(i) = str2dble(tok(3))
    end do
    close(u)
    ! --- eden_apr = 10**col0*mevfm3gcm3, p_apr = 10**col1/ergcm3mevfm3 ---
    jj = 0
    do i = 1, nrow
       if (c2(i)/fm3cm3 < 0.075d0) jj = jj + 1
    end do
    ncr = jj
    allocate(rho_apr(ncr), p_apr(ncr), eden_apr(ncr))
    do i = 1, ncr
       eden_apr(i) = pw(10.0d0, c0(i))*mevfm3gcm3
       p_apr(i)    = pw(10.0d0, c1(i))/ergcm3mevfm3
       rho_apr(i)  = c2(i)/fm3cm3
    end do
  contains
    logical function is_data(s)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: a
      a = adjustl(s)
      is_data = (len_trim(a) > 0) .and. (a(1:1) /= '#')
    end function is_data
    subroutine first_tokens(s, t)
      character(len=*),  intent(in)  :: s
      character(len=64), intent(out) :: t(3)
      integer :: p, q, k
      p = 1
      do k = 1, 3
         do while (p <= len(s))
            if (s(p:p) /= ' ' .and. s(p:p) /= char(9)) exit
            p = p + 1
         end do
         q = p
         do while (q <= len(s))
            if (s(q:q) == ' ' .or. s(q:q) == char(9)) exit
            q = q + 1
         end do
         if (q <= p) then
            write(6,*) 'Error: fewer than 3 columns in ', trim(fname)
            stop 1
         end if
         t(k) = s(p:q-1)
         p = q
      end do
    end subroutine first_tokens
  end subroutine load_crust

  ! --- Build one EOS table (compute_eos in eos_v3_cli.py) ---
  subroutine compute_eos(J0, L, Ksym, Jsym, status, eos)
    real(8), intent(in)  :: J0, L, Ksym, Jsym
    integer, intent(out) :: status
    real(8), intent(out) :: eos(NG, 3)
    real(8) :: e0x(NG), esym(NG), ypr(NG), yn(NG), ye(NG), ymu(NG)
    real(8) :: sq, e_beta, eel_beta, kfmu, mumu, emu_beta, aa, bb
    real(8) :: etot(NG), eden_tmp(NG), df(NG), pz(NG), edenz(NG)
    real(8) :: rhoxx(ncr+NG), pxx(ncr+NG), edenxx(ncr+NG)
    real(8) :: eden(NG), p(NG), lp(NG), le(NG), rp(NG), rm(NG)
    real(8) :: pp(NG), pm(NG), edp(NG), edm(NG), cs2(NG)
    real(8) :: c3, c4pi, cmu4, chb3
    type(spline_t) :: f, f2
    integer :: i, nxx
    logical :: ok

    eos = 0.0d0
    call e_0(rhox, J0, e0x)
    call e_sym(rhox, L, Ksym, Jsym, esym)

    if ((e0x(NG) > e0x(NG) + esym(NG)) .or. e0x(NG) < 0.0d0 .or. esym(NG) < 0.0d0) then
       status = ST_SKIP1
       return
    end if

    ! --- Particle concentrations ---
    call solve_beta(esym, ypr, yn, ye, ymu, ok)
    if (.not. ok) then
       status = ST_SOLVER
       return
    end if

    if (yn(NG) == yn(NG-1) .and. yn(NG-1) == 1.0d0) then
       status = ST_SKIP2
       return
    end if

    ! --- beta-stable matter (total energy) ---
    c3   = (3.0d0*pi)*pi
    c4pi = (4.0d0*pi)*pi
    cmu4 = 0.5d0*pw(mm, 4.0d0)
    chb3 = c4pi*pw(hbar_c, 3.0d0)
    do i = 1, NG
       sq       = (1.0d0 - 2.0d0*ypr(i))*(1.0d0 - 2.0d0*ypr(i))
       e_beta   = e0x(i) + esym(i)*sq + mn*yn(i) + mp*ypr(i)
       eel_beta = ((hbar_c*pw(c3*ye(i), 4.0d0/3.0d0))/c4pi)*pw(rhox(i), third)
       kfmu     = hbar_c*pw((c3*rhox(i))*ymu(i), third)
       mumu     = sqrt(mm*mm + (hbar_c*hbar_c)*pw((c3*rhox(i))*ymu(i), two3))
       aa       = (mumu*kfmu)*(mumu*mumu - (0.5d0*mm)*mm)
       bb       = cmu4*ln((mumu + kfmu)/mm)
       emu_beta = ((aa - bb)/rhox(i)) / chb3
       etot(i)  = e_beta + eel_beta + emu_beta          ! Total energy
    end do

    ! --- Total pressure ---
    call spline_build(NG, rhox, etot, f)
    call spline_deriv5(f, NG, xx0, 1.0d-6, df)
    pz = (xx0*xx0)*df

    ! --- Energy density ---
    eden_tmp = rhox*etot
    call spline_build(NG, rhox, eden_tmp, f2)
    call spline_eval(f2, NG, xx0, edenz)

    ! --- crust + (xx0 > 0.09) ---
    nxx = NG - kk
    rhoxx(1:ncr)          = rho_apr
    rhoxx(ncr+1:ncr+nxx)  = xx0(kk+1:NG)
    pxx(1:ncr)            = p_apr
    pxx(ncr+1:ncr+nxx)    = pz(kk+1:NG)
    edenxx(1:ncr)         = eden_apr
    edenxx(ncr+1:ncr+nxx) = edenz(kk+1:NG)

    ! --- Interpolated values ---
    call interp(NG, zzx, ncr+nxx, rhoxx, edenxx, eden)
    call interp(NG, zzx, ncr+nxx, rhoxx, pxx, p)

    if (any(p < 0.0d0)) then
       status = ST_SKIP3
       return
    end if

    if (.not. paper .and. p(NG) < 1200.0d0) then
       status = ST_SKIP4
       return
    end if

    ! --- Speed of sound ---
    do i = 1, NG
       rp(i) = rlog(i) + 1.0d-6
       rm(i) = rlog(i) - 1.0d-6
       lp(i) = lg10(p(i))
       le(i) = lg10(eden(i))
    end do
    call interp(NG, rp, NG, lzzx, lp, pp)
    call interp(NG, rm, NG, lzzx, lp, pm)
    call interp(NG, rp, NG, lzzx, le, edp)
    call interp(NG, rm, NG, lzzx, le, edm)
    do i = 1, NG
       pp(i)  = pw(10.0d0, pp(i))
       pm(i)  = pw(10.0d0, pm(i))
       edp(i) = pw(10.0d0, edp(i))
       edm(i) = pw(10.0d0, edm(i))
    end do

    ! Simple 2-point derivative approxiamtion...........................
    cs2 = (pp - pm) / (edp - edm)

    if (any(cs2 < 0.0d0)) then
       status = ST_SKIP5
       return
    end if

    ! --- Final tables ---
    eos(:,1) = eden/mevfm3gcm3   ! Energy density
    eos(:,2) = p/dycm2mevfm3     ! Pressure
    eos(:,3) = zzx*fm3cm3        ! Number density
    status = ST_OK
  end subroutine compute_eos

end module eos3_model

!=====================================================================
! MAIN PROGRAM
!=====================================================================
program eos_v3_cli
  use eos3_libm
  use eos3_rng
  use eos3_phys
  use eos3_model
  use hdf5
  implicit none

  integer(8)          :: ntry, seed, i, nprint, t0, t1, trate
  logical             :: have_seed
  character(len=1024) :: outfile
  real(8)             :: J0, L, Ksym, Jsym
  real(8), allocatable :: pars(:,:)          ! (4, ntry): J0, L, Ksym, Jsym
  real(8), allocatable :: eos_all(:,:,:)     ! (NG, 3, cap)
  real(8), allocatable :: par_all(:,:)       ! (4, cap)
  real(8)             :: eos(NG, 3)
  integer             :: status, nacc, cap, k
  integer             :: counts(0:6)
  character(len=512)  :: rej
  character(len=20)   :: num

  call parse_args(ntry, outfile, seed, have_seed, paper)
  call system_clock(t0, trate)

  call libm_init()
  call phys_init()
  call load_crust('EOS_P_A18dvUIX.txt')

  ! --- Monte Carlo sampling of the EOS parameters ---
  ! (same draw order as eos_v3.py: J0, Jsym, L, Ksym)
  if (have_seed) then
     call rng_seed(seed)
  else
     call rng_seed_urandom()
  end if
  allocate(pars(4, ntry))
  do i = 1, ntry
     if (paper) then
        ! --- only L and Ksym are sampled (in this order) ---
        J0   = 0.0d0
        Jsym = 0.0d0
     else
        J0   = rng_uniform(J0_low, J0_high)
        Jsym = rng_uniform(Jsym_low, Jsym_high)
     end if
     L    = rng_uniform(L_low, L_high)
     Ksym = rng_uniform(Ksym_low, Ksym_high)
     pars(:, i) = (/ J0, L, Ksym, Jsym /)
  end do

  if (paper) then
     write(6,'(a,i0,a)') 'Trying ', ntry, ' instances (paper setup)'
  else
     write(6,'(a,i0,a)') 'Trying ', ntry, ' instances'
  end if
  flush(6)

  cap  = 1024
  nacc = 0
  allocate(eos_all(NG, 3, cap), par_all(4, cap))
  counts = 0
  nprint = max(1_8, ntry / 20)

  do i = 1, ntry
     call compute_eos(pars(1,i), pars(2,i), pars(3,i), pars(4,i), status, eos)
     counts(status) = counts(status) + 1
     if (status == ST_OK) then
        if (nacc == cap) call grow()
        nacc = nacc + 1
        eos_all(:, :, nacc) = eos
        par_all(:, nacc)    = pars(:, i)
     end if
     if (mod(i, nprint) == 0) then
        call system_clock(t1)
        write(6,'(a,i0,a,i0,a,i0,a,f0.1,a)') 'iter: ', i, '/', ntry, '  accepted: ', &
             nacc, '  (', dble(t1 - t0)/dble(trate), ' s)'
        flush(6)
     end if
  end do
  ! --- END OF LOOP OVER NUMBER OF SAMPLES ---

  call write_h5(trim(outfile))

  rej = ''
  do k = 1, 6
     if (counts(k) > 0) then
        if (len_trim(rej) > 0) rej = trim(rej)//','
        write(num, '(i0)') counts(k)
        rej = trim(rej)//' '//trim(st_name(k))//': '//trim(num)
     end if
  end do
  call system_clock(t1)
  write(6,'(a,a)') 'Rejected:', trim(rej)
  write(6,'(a,i0,a,i0,a,a,a,f0.1,a)') 'Accepted ', nacc, ' of ', ntry, &
       ' instances; written to ', trim(outfile), '  (', dble(t1 - t0)/dble(trate), ' s)'

contains

  subroutine grow()
    real(8), allocatable :: e2(:,:,:), p2(:,:)
    allocate(e2(NG, 3, 2*cap), p2(4, 2*cap))
    e2(:, :, 1:cap) = eos_all
    p2(:, 1:cap)    = par_all
    call move_alloc(e2, eos_all)
    call move_alloc(p2, par_all)
    cap = 2*cap
  end subroutine grow

  ! --- Write data (same datasets as eos_v3.py: eos (N,3,500), eospar (N,4);
  !     with --paper eospar is (N,2) = [L, Ksym] as in the paper files) ---
  subroutine write_h5(fname)
    character(len=*), intent(in) :: fname
    integer(HID_T)   :: file_id, space_id, dset_id
    integer(HSIZE_T) :: d3(3), d2(2), d1(1)
    integer          :: err
    real(8), allocatable :: par2(:,:)
    call h5open_f(err)
    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, err)
    if (err /= 0) then
       write(6,*) 'Error: cannot create ', fname
       stop 1
    end if
    if (nacc > 0) then
       d3 = (/ int(NG, HSIZE_T), 3_HSIZE_T, int(nacc, HSIZE_T) /)
       call h5screate_simple_f(3, d3, space_id, err)
       call h5dcreate_f(file_id, 'eos', H5T_IEEE_F64LE, space_id, dset_id, err)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eos_all(:, :, 1:nacc), d3, err)
       call h5dclose_f(dset_id, err)
       call h5sclose_f(space_id, err)
       if (paper) then
          d2 = (/ 2_HSIZE_T, int(nacc, HSIZE_T) /)
          allocate(par2(2, nacc))
          par2 = par_all(2:3, 1:nacc)
       else
          d2 = (/ 4_HSIZE_T, int(nacc, HSIZE_T) /)
          allocate(par2(4, nacc))
          par2 = par_all(:, 1:nacc)
       end if
       call h5screate_simple_f(2, d2, space_id, err)
       call h5dcreate_f(file_id, 'eospar', H5T_IEEE_F64LE, space_id, dset_id, err)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, par2, d2, err)
       call h5dclose_f(dset_id, err)
       call h5sclose_f(space_id, err)
    else
       ! --- no EOS accepted: empty 1-D datasets, as np.array([]) in Python ---
       d1 = 0_HSIZE_T
       call h5screate_simple_f(1, d1, space_id, err)
       call h5dcreate_f(file_id, 'eos', H5T_IEEE_F64LE, space_id, dset_id, err)
       call h5dclose_f(dset_id, err)
       call h5dcreate_f(file_id, 'eospar', H5T_IEEE_F64LE, space_id, dset_id, err)
       call h5dclose_f(dset_id, err)
       call h5sclose_f(space_id, err)
    end if
    call h5fclose_f(file_id, err)
    call h5close_f(err)
  end subroutine write_h5

  ! --- Command line: -n/--ntry, -o/--output, -s/--seed, -h/--help ---
  subroutine parse_args(ntry, outfile, seed, have_seed, paper_mode)
    integer(8),       intent(out) :: ntry, seed
    character(len=*), intent(out) :: outfile
    logical,          intent(out) :: have_seed, paper_mode
    character(len=1024) :: arg, key, val
    integer :: ia, na, p, ios
    ntry      = 100000
    outfile   = 'eos_set_test.h5'
    seed      = 0
    have_seed = .false.
    paper_mode = .false.
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
          if (key /= '-h' .and. key /= '--help' .and. key /= '--paper') then
             ia = ia + 1
             if (ia > na) call usage_error('argument '//trim(key)//': expected one argument')
             call get_command_argument(ia, val)
          end if
       end if
       select case (trim(key))
       case ('-h', '--help')
          call usage(.true.)
          stop
       case ('--paper')
          paper_mode = .true.
       case ('-n', '--ntry')
          read(val, *, iostat=ios) ntry
          if (ios /= 0) call usage_error('argument -n/--ntry: invalid int value: '''//trim(val)//'''')
       case ('-o', '--output')
          outfile = val
       case ('-s', '--seed')
          read(val, *, iostat=ios) seed
          if (ios /= 0) call usage_error('argument -s/--seed: invalid int value: '''//trim(val)//'''')
          have_seed = .true.
       case default
          call usage_error('unrecognized arguments: '//trim(arg))
       end select
       ia = ia + 1
    end do
    if (ntry < 1) call usage_error('--ntry must be a positive integer')
  end subroutine parse_args

  subroutine usage(full)
    logical, intent(in) :: full
    write(6,'(a)') 'usage: eos_v3_cli.x [-h] [-n NTRY] [-o OUTPUT] [--paper] [-s SEED]'
    if (.not. full) return
    write(6,'(a)') ''
    write(6,'(a)') 'Generate parametric EOSs via Monte-Carlo sampling'
    write(6,'(a)') ''
    write(6,'(a)') 'options:'
    write(6,'(a)') '  -h, --help            show this help message and exit'
    write(6,'(a)') '  -n NTRY, --ntry NTRY  number of instances (Monte-Carlo samples) to try'
    write(6,'(a)') '                        (default: 100000)'
    write(6,'(a)') '  -o OUTPUT, --output OUTPUT'
    write(6,'(a)') '                        name of the output HDF5 file (default:'
    write(6,'(a)') '                        eos_set_test.h5)'
    write(6,'(a)') '  --paper               paper setup: J0 = Jsym = 0, tables up to 10**0.015'
    write(6,'(a)') '                        fm^-3, no skip 4, eospar = [L, Ksym]'
    write(6,'(a)') '  -s SEED, --seed SEED  random seed, for reproducible runs (default: none)'
  end subroutine usage

  subroutine usage_error(msg)
    character(len=*), intent(in) :: msg
    call usage(.false.)
    write(6,'(a,a)') 'eos_v3_cli.x: error: ', msg
    stop 2
  end subroutine usage_error

end program eos_v3_cli
