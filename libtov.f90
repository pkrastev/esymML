!=====================================================================
! Modules
!.....................................................................
module constants
  implicit none
  ! Mathematical constants............................................
  real(8), parameter :: pi       = 3.1415926535897932384d0
  real(8), parameter :: eulercon = 0.577215664901532861d0   
  real(8), parameter :: a2rad    = pi/180.0d0
  real(8), parameter :: rad2a    = 180.0d0/pi
  ! Physical constants (cgs)..........................................
  real(8), parameter :: c       = 2.99792458d10
  real(8), parameter :: g       = 6.67259d-8
  real(8), parameter :: h       = 6.6260755d-27 
  real(8), parameter :: hbar    = 0.5d0 * h / pi
  real(8), parameter :: qe      = 4.8032068d-10   
  real(8), parameter :: avo     = 6.0221367d23 
  real(8), parameter :: kerg    = 1.380658d-16
  real(8), parameter :: kev     = 8.617385d-5 
  real(8), parameter :: amu     = 1.6605402d-24
  real(8), parameter :: mn      = 1.6749286d-24   
  real(8), parameter :: mp      = 1.6726231d-24   
  real(8), parameter :: me      = 9.1093897d-28   
  real(8), parameter :: rbohr   = hbar * hbar / ( me * qe * qe )
  real(8), parameter :: fine    = qe * qe / (hbar * c )  
  real(8), parameter :: hion    = 13.605698140d0  
  real(8), parameter :: ev2erg  = 1.602d-12
  real(8), parameter :: ssol    = 5.67051d-5
  real(8), parameter :: asol    = 4.0d0 * ssol / c 
  real(8), parameter :: weinlam = h * c / ( kerg * 4.965114232d0 ) 
  real(8), parameter :: weinfre = 2.821439372d0 * kerg / h   
  real(8), parameter :: rhonuc  = 2.342d14 
  ! Astrophysical constants...........................................
  real(8), parameter :: r_sun   = 6.95997d10
  real(8), parameter :: m_sun   = 1.9892d33
  real(8), parameter :: lsol    = 3.8268d33 
  real(8), parameter :: mearth  = 5.9764d27 
  real(8), parameter :: rearth  = 6.37d8 
  real(8), parameter :: ly      = 9.460528d17   
  real(8), parameter :: pc      = 3.261633d0 * ly 
  real(8), parameter :: au      = 1.495978921d13
  real(8), parameter :: secyer  = 3.1558149984d7
  ! Conversion factors................................................
  real(8), parameter :: fm3cm3      = 1.0d+39
  real(8), parameter :: mevfm3gcm3  = ( 1.0d0/1.7827d0 )*1.0d-12
  real(8), parameter :: dycm2mevfm3 = ( 197.3269788d0/3.1616d0 )*1.0d-35
  real(8), parameter :: hbar_c      = 197.3269788d0
  real(8), parameter :: pG          = (g / (c*c*c*c))*1.0d4  ! Eden, P in [1/m2]
  real(8), parameter :: mG          = (g/(c*c))*1.0d-2       ! Mass in [m]
end module constants

! MODULE EOS..........................................................
module eos
  implicit none
  integer(4)           :: np       ! Number of rows in EOS file
  real(8), allocatable :: eray(:)  ! Energy density
  real(8), allocatable :: pray(:)  ! Pressure
  real(8), allocatable :: hray(:)  ! Enthalpy
  real(8), allocatable :: xnray(:) ! Number density
end module eos

! Termination energy density and preassure............................
module terminate
  implicit none
  real(8) :: eden_term
  real(8) :: p_term
end module terminate

! esymML: solver control (tolerance, failure flag)......................
module tov_control
  implicit none
  real(8) :: tov_eps  = 1.0d-8   ! Relative tolerance of the enthalpy solver
  logical :: tov_fail = .false.  ! Set if an integration does not reach the surface
end module tov_control

!=====================================================================
! SUBROUTINE STAR_OUTPUT: Compute k2, lambda, I, beta from the surface
!                         values and write one row of results
!
! Input:
! M    -- gravitational mass (g)
! R    -- radius (cm)
! yr   -- y at the surface
! fr   -- f at the surface (Eq. (12) from arxiv: 1810.10992)
! rhoc -- central number density (1/cm3)
!=====================================================================
subroutine star_output(M, R, yr, fr, rhoc)
  use constants
  implicit none
  real(8) :: M
  real(8) :: R
  real(8) :: yr
  real(8) :: fr
  real(8) :: rhoc
  real(8) :: beta
  real(8) :: k2
  real(8) :: lambda
  real(8) :: i_ns ! Moment of inertia I

  call star_values(M, R, yr, fr, beta, k2, lambda, i_ns)

  write(6,'(7(2x,f11.6))') &
        M/m_sun,           & ! Mass, M
        R/1.0d5,           & ! Radius, R
        k2,                & ! Love number, k2
        lambda/1.0d36,     & ! Tidal deformability, lambda
        i_ns/1.0d45,       & ! Moment of inertia, I
        beta,              & ! Compactness, beta
        rhoc/fm3cm3          ! Central number density, rho_c

  return
end subroutine star_output

!=====================================================================
! esymML: SUBROUTINE STAR_VALUES: beta, k2, lambda and I of a star
!         (the computation from STAR_OUTPUT, without printing)
!
! Input:  M (g), R (cm), yr = y(R), fr = f(R)
! Output: beta, k2 (dimensionless), lambda (cm2 g s2), i_ns (g cm2)
!=====================================================================
subroutine star_values(M, R, yr, fr, beta, k2, lambda, i_ns)
  use constants
  implicit none
  real(8) :: M
  real(8) :: R
  real(8) :: yr
  real(8) :: fr
  real(8) :: beta
  real(8) :: k2
  real(8) :: lambda
  real(8) :: i_ns ! Moment of inertia I

  beta = ( g*M ) / ( R*c*c )                          ! Compactness parameter (dimentionless)
  call solve_k2( yr, beta, k2 )                       ! k2 (dimentionless)
  lambda = (2.0d0 / (3.0d0*g))*k2*(R**5)              ! lambda (cm2 g s2)
  i_ns = (((R**3)*fr)/(6.0d0+2.0d0*fr)) * (c*c/g)     ! Moment of inertia

  return
end subroutine star_values

!=====================================================================
! Load EOS file
! Format: first line is the number of rows, then one row per point
!         {Energy_Density, P, H, Number_Density} in CGS units
!         (after RNS)
! The file is read only when it differs from the one already loaded
!=====================================================================
subroutine load_eos(eos_file)
  use eos
  implicit none
  integer(4)        :: i
  integer(4)        :: iu       ! Unit number
  integer(4)        :: ios      ! I/O status
  real(8)           :: p        ! Pressure
  real(8)           :: eden     ! Energy density
  real(8)           :: h0       ! Enthalpy
  real(8)           :: n0       ! Number density
  character(len=30) :: eos_file ! EOS file
  character(len=30), save :: loaded_file = ''

  if ( eos_file == loaded_file ) return

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Read EOS file
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  open(newunit=iu, file=eos_file, status='old', action='read', iostat=ios)
  if ( ios /= 0 ) then
     write(0,*) 'ERROR load_eos: cannot open EOS file ', trim(eos_file)
     stop 1
  end if
  read(iu,*,iostat=ios) np
  if ( ios /= 0 .or. np < 3 ) then
     write(0,*) 'ERROR load_eos: bad number of rows in ', trim(eos_file)
     stop 1
  end if

  if ( allocated(pray) )  deallocate( pray )
  if ( allocated(eray) )  deallocate( eray )
  if ( allocated(hray) )  deallocate( hray )
  if ( allocated(xnray) ) deallocate( xnray )
  allocate( pray(np), eray(np), hray(np), xnray(np) )

  do i = 1, np
     read(iu,*,iostat=ios) eden, p, h0, n0
     if ( ios /= 0 ) then
        write(0,*) 'ERROR load_eos: cannot read row ', i, ' of ', trim(eos_file)
        stop 1
     end if
     eray(i)  = log10(eden)
     pray(i)  = log10(p)
     hray(i)  = log10(h0)
     xnray(i) = log10(n0)
  end do
  close(iu)

  loaded_file = eos_file

  return
end subroutine load_eos

!=====================================================================
! esymML: SUBROUTINE SET_EOS: EOS table from arrays (instead of a file)
!
! Input: n rows of energy density/c^2 (g cm^-3), pressure (dyn cm^-2)
!        and baryon number density (cm^-3), as in the RNS file format.
!        The enthalpy column is not used by the solvers (set to 0).
!=====================================================================
subroutine set_eos(n, eden, p, xn)
  use eos
  implicit none
  integer(4) :: n
  real(8)    :: eden(n)
  real(8)    :: p(n)
  real(8)    :: xn(n)
  integer(4) :: i

  if ( allocated(pray) ) then
     if ( size(pray) /= n ) deallocate( pray, eray, hray, xnray )
  end if
  if ( .not. allocated(pray) ) allocate( pray(n), eray(n), hray(n), xnray(n) )
  np = n

  do i = 1, np
     eray(i)  = log10(eden(i))
     pray(i)  = log10(p(i))
     hray(i)  = 0.0d0
     xnray(i) = log10(xn(i))
  end do

  return
end subroutine set_eos

!=====================================================================
! SUBROUTINE EOS_INTERP: Energy density and speed of sound at pressure p
!
! Input:
! p      -- pressure (dyn/cm2)
!
! Output:
! eden   -- energy density (g/cm3), as EOSINV(p, pray, eray, np, eden)
! cs2    -- speed of sound squared, cs2 = dP/dE with E = eden*c^2
!           (dimensionless)
!
! Comments:
! Same 3-point Lagrange interpolation in log10 space as LAGINT.
! cs2 = (P/E) / (dlogE/dlogP), from the derivative of the interpolant
!=====================================================================
subroutine eos_interp(p, eden, cs2)
  use constants
  use eos
  implicit none
  integer(4), parameter :: n = 3 ! Number of interpolation points
  integer(4) :: i
  integer(4) :: j
  integer(4) :: k
  real(8)    :: p
  real(8)    :: eden
  real(8)    :: cs2
  real(8)    :: x
  real(8)    :: x0
  real(8)    :: x1
  real(8)    :: x2
  real(8)    :: w(n)  ! Lagrange weights
  real(8)    :: dw(n) ! Derivatives of Lagrange weights
  real(8)    :: dle

  if ( p > 0.0d0 ) then
     x = log10(p)
  else
     x = pray(1)
  end if

  ! a binary (bisectional) search to find i so that xi(i) < x < xi(i+1)
  i = 1
  j = np
  do while ( j > i+1 )
     k = ( i + j ) / 2
     if ( x < pray(k) ) then
        j = k
     else
        i = k
     end if
  end do
  if ( i + n > np ) i = np - n + 1

  ! if x is ouside the table take a boundary value, as in LAGINT
  x  = min(max(x, pray(1)), pray(np))

  x0 = pray(i)
  x1 = pray(i+1)
  x2 = pray(i+2)

  w(1)  = (x-x1)*(x-x2) / ((x0-x1)*(x0-x2))
  w(2)  = (x-x0)*(x-x2) / ((x1-x0)*(x1-x2))
  w(3)  = (x-x0)*(x-x1) / ((x2-x0)*(x2-x1))
  dw(1) = ((x-x1)+(x-x2)) / ((x0-x1)*(x0-x2))
  dw(2) = ((x-x0)+(x-x2)) / ((x1-x0)*(x1-x2))
  dw(3) = ((x-x0)+(x-x1)) / ((x2-x0)*(x2-x1))

  eden = 10.0d0**dot_product(w, eray(i:i+2))

  dle = dot_product(dw, eray(i:i+2))
  cs2 = ( 10.0d0**x / (eden*c*c) ) / dle

  return
end subroutine eos_interp

!=====================================================================
! SUBROUTINE EOSINV: Interface to LAGINT
!                    Returns pressure, energy density, or 
!                    number density
!=====================================================================
subroutine eosinv(xx, xi, yi, ni, yy)
  implicit none
  integer(4), parameter :: order = 2 ! Order of interpolation
  integer(4)            :: ni        ! Size of arrays xi() and yi()
  real(8)               :: xx        ! Abscissa at which the interpolation is to be evaluated
  real(8)               :: yy        ! Interpolated value
  real(8)               :: xi(ni)    ! Arrays of data abscissas
  real(8)               :: yi(ni)    ! Arrays of data ordinates
  real(8), external     :: lagint
  real(8)               :: xl
  real(8)               :: yl

  if ( xx <= 0.0d0 ) then
     yy = 10.0d0**yi(1)
     return
  end if

  xl = log10(xx)
  yl = lagint(xl, xi, yi, ni, order+1)
  yy = 10.0d0**yl
  return
end subroutine eosinv
 
!=====================================================================
! FUNCTION LAGINT: Lagrange interpolation
!
! Input:
! xx    -- abscissa at which the interpolation is to be evaluated
! xi()  -- arrays of data abscissas
! yi()  -- arrays of data ordinates
! ni    -- size of the arrays xi() and yi()
! n     -- number of points for interpolation (order of interp. = n-1)
!
! Output:
! lagint- interpolated value
!
! Comments:
! if ( n > ni ) n = ni
! Program works for both equally and unequally spaced xi()
!=====================================================================
function lagint(xx, xi, yi, ni, n)
  implicit none
  integer(4) :: i
  integer(4) :: j
  integer(4) :: k
  integer(4) :: js
  integer(4) :: jl
  integer(4) :: ni
  integer(4) :: n
  integer(4) :: nn      ! Number of points used, min(n, ni)
  real(8)    :: xi(ni)
  real(8)    :: yi(ni)
  real(8)    :: lambda
  real(8)    :: y
  real(8)    :: lagint
  real(8)    :: xx

  ! check order of interpolation (n may be an expression: do not modify)
  nn = min(n, ni)

  ! if x is ouside the xi(1)-xi(ni) interval take a boundary value
  if ( xx <= xi(1) ) then
     lagint = yi(1)
     return
  end if
  if ( xx >= xi(ni) ) then
     lagint = yi(ni)
     return
  end if

  ! a binary (bisectional) search to find i so that xi(i) < x < xi(i+1)
  i = 1
  j = ni
  do while ( j > i+1 )
     k = ( i + j ) / 2
     if ( xx < xi(k) ) then
        j = k
     else
        i = k
     end if
  end do

  ! shift i that will correspond to n-th order of interpolation
  ! the search point will be in the middle in x_i, x_i+1, x_i+2 ...
  i = i + 1 - nn/2

  ! check boundaries: if i is ouside of the range [1, ... n] -> shift i
  if ( i < 1 ) i = 1
  if ( i + nn > ni ) i = ni - nn + 1

  ! Lagrange interpolation
  y = 0.0d0
  do js =i, i + nn - 1
     lambda = 1.0d0
     do jl = i, i + nn - 1
        if( jl /= js ) lambda = lambda*(xx-xi(jl))/(xi(js)-xi(jl))
     end do
     y = y + yi(js)*lambda
  end do
  lagint = y
  return
end function lagint

!=====================================================================
! SUBROUTINE SOLVE_K2: Evaluate Love number k2
!=====================================================================
subroutine solve_k2(yr, beta, k2)
  implicit none
  real(8) :: beta
  real(8) :: yr
  real(8) :: k2
  real(8) :: num
  real(8) :: den
  num = (8.0d0/5.0d0)*(beta**5)*((1.0d0-2.0d0*beta)**2)                      &
       *(2.0d0-yr+2.0d0*beta*(yr-1.0d0))
  den =  2.0d0*beta*(6.0d0-3.0d0*yr+3.0d0*beta*(5.0d0*yr-8.0d0))+            &
         4.0d0*beta*beta*beta*(13.0d0-11.0d0*yr+beta*(3.0d0*yr-2.0d0)+2.0d0* &
         beta*beta*(1.0d0+yr))+                                              &
         3.0d0*((1.0d0-2.0d0*beta)**2)*(2.0d0-yr+2.0d0*beta*(yr-1.0d0))*     &
         log(1.0d0-2.0d0*beta)
  k2 = num / den
  return
end subroutine solve_k2

!=====================================================================
! SUBROUTINE RKDP: Adaptive Dormand-Prince 5(4) Runge-Kutta step
!
! Input:
! y()    -- dependent variables at x
! dydx() -- derivatives at x
! n      -- number of equations
! x      -- independent variable
! htry   -- step size to attempt
! eps    -- required accuracy
! yscal()-- scaling vector for the error
! derivs -- subroutine derivs(x, y, dydx)
!
! Output:
! y()    -- dependent variables at x + hdid (5th order solution)
! dydx() -- derivatives at x + hdid (FSAL: reuse for the next step)
! x      -- x + hdid
! hdid   -- step size actually taken
! hnext  -- estimated step size for the next step
!
! Comments:
! Coefficients from J.R. Dormand & P.J. Prince,
! J. Comp. Appl. Math. 6, 19 (1980)
! The last stage is evaluated at the new point, so the driver does
! not need to call derivs() between steps
!=====================================================================
subroutine rkdp(y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs)
  use tov_control
  implicit none
  integer(4)         :: n
  real(8)            :: x
  real(8)            :: htry
  real(8)            :: eps
  real(8)            :: hdid
  real(8)            :: hnext
  real(8)            :: h
  real(8)            :: xnew
  real(8)            :: errmax
  real(8)            :: y(n)
  real(8)            :: dydx(n)
  real(8)            :: yscal(n)
  real(8)            :: k2(n), k3(n), k4(n), k5(n), k6(n), k7(n)
  real(8)            :: ytemp(n)
  real(8)            :: yerr(n)
  real(8), parameter :: c2  = 1.0d0/5.0d0,         c3  = 3.0d0/10.0d0,     &
                        c4  = 4.0d0/5.0d0,         c5  = 8.0d0/9.0d0
  real(8), parameter :: a21 = 1.0d0/5.0d0
  real(8), parameter :: a31 = 3.0d0/40.0d0,        a32 = 9.0d0/40.0d0
  real(8), parameter :: a41 = 44.0d0/45.0d0,       a42 = -56.0d0/15.0d0,   &
                        a43 = 32.0d0/9.0d0
  real(8), parameter :: a51 = 19372.0d0/6561.0d0,  a52 = -25360.0d0/2187.0d0, &
                        a53 = 64448.0d0/6561.0d0,  a54 = -212.0d0/729.0d0
  real(8), parameter :: a61 = 9017.0d0/3168.0d0,   a62 = -355.0d0/33.0d0,  &
                        a63 = 46732.0d0/5247.0d0,  a64 = 49.0d0/176.0d0,   &
                        a65 = -5103.0d0/18656.0d0
  ! 5th order weights (also the 7th stage row)
  real(8), parameter :: b1  = 35.0d0/384.0d0,      b3  = 500.0d0/1113.0d0, &
                        b4  = 125.0d0/192.0d0,     b5  = -2187.0d0/6784.0d0, &
                        b6  = 11.0d0/84.0d0
  ! Error weights: 5th order minus 4th order
  real(8), parameter :: e1  = 71.0d0/57600.0d0,    e3  = -71.0d0/16695.0d0, &
                        e4  = 71.0d0/1920.0d0,     e5  = -17253.0d0/339200.0d0, &
                        e6  = 22.0d0/525.0d0,      e7  = -1.0d0/40.0d0
  ! Step size control
  real(8), parameter :: safety = 0.9d0
  real(8), parameter :: pgrow  = -0.2d0
  real(8), parameter :: pshrnk = -0.25d0
  real(8), parameter :: errcon = 1.89d-4           ! (5/safety)**(1/pgrow)
  external derivs

  h = htry

  do
     ytemp = y + h*a21*dydx
     call derivs(x+c2*h, ytemp, k2)
     ytemp = y + h*(a31*dydx + a32*k2)
     call derivs(x+c3*h, ytemp, k3)
     ytemp = y + h*(a41*dydx + a42*k2 + a43*k3)
     call derivs(x+c4*h, ytemp, k4)
     ytemp = y + h*(a51*dydx + a52*k2 + a53*k3 + a54*k4)
     call derivs(x+c5*h, ytemp, k5)
     ytemp = y + h*(a61*dydx + a62*k2 + a63*k3 + a64*k4 + a65*k5)
     xnew  = x + h
     call derivs(xnew, ytemp, k6)
     ytemp = y + h*(b1*dydx + b3*k3 + b4*k4 + b5*k5 + b6*k6)
     call derivs(xnew, ytemp, k7)

     yerr   = h*(e1*dydx + e3*k3 + e4*k4 + e5*k5 + e6*k6 + e7*k7)
     errmax = maxval(abs(yerr/yscal)) / eps

     if ( errmax <= 1.0d0 ) exit

     ! Reject: shrink the step, but by no more than a factor of 10
     h = sign(max(abs(safety*h*(errmax**pshrnk)), 0.1d0*abs(h)), h)
     if ( x + h == x ) then
        ! esymML: flag the star instead of stopping the program
        write(0,*) 'WARNING rkdp: stepsize underflow at x =', x
        tov_fail = .true.
        hdid  = 0.0d0
        hnext = 0.0d0
        return
     end if
  end do

  ! Accept: grow the step, but by no more than a factor of 5
  if ( errmax > errcon ) then
     hnext = safety*h*(errmax**pgrow)
  else
     hnext = 5.0d0*h
  end if

  hdid = h
  x    = xnew
  y    = ytemp
  dydx = k7

  return
end subroutine rkdp
