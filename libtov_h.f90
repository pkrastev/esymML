!=====================================================================
! TOV solver in the enthalpy formalism
!
! L. Lindblom, ApJ 398, 569 (1992)
!
! The independent variable is the pseudo-enthalpy
!    h = int_0^P dP' / ( e(P') + P' )    ( c = G = 1 )
! measured from the surface, so the star is integrated from h = h_c
! at the center to exactly h = 0 at the surface, at P = p_term (second
! row of the EOS table).
!
! P is carried as a dependent variable ( dP/dh = e + P ), and e(P),
! cs2(P) come from EOS_INTERP. The enthalpy column of the EOS file is
! not used.
!
! Requires the modules and routines in libtov.f90
!=====================================================================

!=====================================================================
! SUBROUTINE SOLVE_TOV_H: Compute spherical NS model
!                         {M, R, k2, lambda, I, beta, rhoc}
!=====================================================================
subroutine solve_tov_h(file_in, rhoc)
  use constants
  implicit none
  real(8)           :: rhoc         ! Central number density
  real(8)           :: y(5)         ! r, m, P, y, f at the surface
  character(len=30) :: file_in      ! Input file

  ! Read EOS file.....................................................
  call load_eos(file_in)

  ! Integrate from the center to the surface..........................
  call tov_h_integrate(rhoc, y)

  ! Write out {M, R, k2, lambda, I, beta, rhoc}......................
  call star_output( y(2)/mG, y(1)*1.0d2, y(4), y(5), rhoc )

  return
end subroutine solve_tov_h

!=====================================================================
! esymML: SUBROUTINE SOLVE_TOV_H_STAR: one star for the EOS already
!         loaded (LOAD_EOS or SET_EOS), results returned, not printed
!
! Input:  rhoc -- central baryon number density (cm^-3)
! Output: M (Msun), R (km), k2, lambda (1e36 g cm2 s2), I (1e45 g cm2),
!         beta -- the columns printed by SOLVE_TOV_H.
!         All NaN if the integration did not reach the surface.
!=====================================================================
subroutine solve_tov_h_star(rhoc, M, R, k2, lambda, i_ns, beta)
  use constants
  use tov_control
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  real(8) :: rhoc
  real(8) :: M, R, k2, lambda, i_ns, beta
  real(8) :: y(5)
  real(8) :: m_g, r_cm   ! M (g), R (cm)

  tov_fail = .false.
  call tov_h_integrate(rhoc, y)

  if ( tov_fail ) then
     M = ieee_value(1.0d0, ieee_quiet_nan)
     R = M; k2 = M; lambda = M; i_ns = M; beta = M
     return
  end if

  m_g  = y(2)/mG
  r_cm = y(1)*1.0d2
  call star_values(m_g, r_cm, y(4), y(5), beta, k2, lambda, i_ns)
  M      = m_g/m_sun
  R      = r_cm/1.0d5
  lambda = lambda/1.0d36
  i_ns   = i_ns/1.0d45

  return
end subroutine solve_tov_h_star

!=====================================================================
! esymML: SUBROUTINE TOV_H_INTEGRATE: integrate one star (the body of
!         the original SOLVE_TOV_H) for the EOS already loaded
!
! Input:  rhoc -- central baryon number density (cm^-3)
! Output: y()  -- r (m), m (m), P (1/m2), y, f at the surface
!=====================================================================
subroutine tov_h_integrate(rhoc, y_out)
  use eos
  use constants
  use terminate
  use tov_control
  implicit none
  real(8)               :: y_out(5)

  integer(4), parameter :: nvar = 5     ! r, m, P, y, f
  integer(4)            :: nok
  integer(4)            :: nbad
  real(8)               :: rhoc         ! Central number density
  real(8)               :: p_c          ! Central pressure (dyn/cm2)
  real(8)               :: edenc        ! Central energy density (g/cm3)
  real(8)               :: cs2c         ! Central speed of sound squared
  real(8)               :: e_g          ! Central energy density (1/m2)
  real(8)               :: p_g          ! Central pressure (1/m2)
  real(8)               :: dedh         ! de/dh at the center (1/m2)
  real(8)               :: h_c          ! Central enthalpy
  real(8)               :: dh           ! Offset of the starting point
  real(8)               :: h0           ! Starting enthalpy
  real(8)               :: r0           ! Starting radius (m)
  real(8)               :: m0           ! Starting mass (m)
  real(8)               :: y(nvar)
  real(8)               :: eps
  real(8)               :: h_try
  real(8)               :: h_min
  real(8), external     :: delta_h

  external tov_h_derivs, rkdp

  ! Get central pressure, energy density and speed of sound...........
  call eosinv(rhoc, xnray, pray, np, p_c)
  call eos_interp(p_c, edenc, cs2c)

  ! Get termination density and pressure (surface, h = 0).............
  eden_term = 10.d0**eray(2)
  call eosinv(eden_term, eray, pray, np, p_term)

  ! Central enthalpy..................................................
  h_c = delta_h(p_term, p_c)

  ! Central values in geometric units ( c = G = 1, lengths in m )....
  e_g  = edenc*c*c*pG
  p_g  = p_c*pG
  dedh = (e_g + p_g)/cs2c

  ! Initial conditions at h0 = h_c - dh, from the series expansion
  ! around the center (Lindblom 1992, Eqs. (7) and (8))...............
  dh = 1.0d-7*h_c
  h0 = h_c - dh
  r0 = sqrt( 3.0d0*dh / (2.0d0*pi*(e_g+3.0d0*p_g)) )                     &
       * ( 1.0d0 - 0.25d0*(e_g-3.0d0*p_g-0.6d0*dedh)*dh/(e_g+3.0d0*p_g) )
  m0 = (4.0d0/3.0d0)*pi*e_g*(r0**3) * ( 1.0d0 - 0.6d0*dedh*dh/e_g )

  y(1) = r0                                      ! Radius (m)
  y(2) = m0                                      ! Gravitational mass (m)
  y(3) = p_g - (e_g+p_g)*dh                      ! Pressure (1/m2)
  y(4) = 2.0d0                                   ! y (for k2)
  y(5) = (16.0d0*pi/5.0d0)*(e_g+p_g)*(r0**2)     ! f (for I)

  ! Set parameters for integration of ODEs............................
  eps   = tov_eps                   ! esymML: default 1.0d-8, as before
  h_try = dh
  h_min = 1.0d-14*h_c

  ! Integrate from the center (h0) to the surface (h = 0).............
  call hint(y, nvar, h0, 0.0d0, eps, h_try, h_min, nok, nbad, &
       tov_h_derivs, rkdp)

  y_out = y

  return
end subroutine tov_h_integrate

!=====================================================================
! FUNCTION DELTA_H: Enthalpy difference between pressures p1 < p2
!
!    delta_h = int_p1^p2 dP / ( e(P) c^2 + P )    (dimensionless)
!
! with e(P) from EOS_INTERP. The interpolant is smooth between table
! points, so 5-point Gauss-Legendre in log10(P) on each table
! interval is accurate to round-off
!=====================================================================
function delta_h(p1, p2)
  use constants
  use eos
  implicit none
  integer(4) :: j
  integer(4) :: k
  real(8)    :: p1
  real(8)    :: p2
  real(8)    :: delta_h
  real(8)    :: xa
  real(8)    :: xb
  real(8)    :: a
  real(8)    :: b
  real(8)    :: s
  real(8)    :: x
  real(8)    :: p
  real(8)    :: eden
  real(8)    :: cs2
  real(8), parameter :: ln10 = 2.302585092994045684d0
  real(8), parameter :: tg(5) = (/ -0.9061798459386640d0, -0.5384693101056831d0, &
                                   0.0d0, 0.5384693101056831d0, 0.9061798459386640d0 /)
  real(8), parameter :: wg(5) = (/ 0.2369268850561891d0, 0.4786286704993665d0, &
                                   0.5688888888888889d0, 0.4786286704993665d0, &
                                   0.2369268850561891d0 /)

  xa = log10(p1)
  xb = log10(p2)
  s  = 0.0d0

  ! Split [xa, xb] at the table points (and past the end of the table)
  a = xa
  do k = 1, np + 1
     if ( k <= np ) then
        if ( pray(k) <= a ) cycle
        b = min(pray(k), xb)
     else
        b = xb
     end if

     do j = 1, 5
        x = 0.5d0*(a+b) + 0.5d0*(b-a)*tg(j)
        p = 10.0d0**x
        call eos_interp(p, eden, cs2)
        s = s + 0.5d0*(b-a)*wg(j) * p*ln10 / (eden*c*c + p)
     end do

     a = b
     if ( a >= xb ) exit
  end do

  delta_h = s
  return
end function delta_h

!=====================================================================
! SUBROUTINE TOV_H_DERIVS: TOV, y and f equations with the enthalpy h
!                          as the independent variable
!                          ( c = G = 1, lengths in m )
!=====================================================================
subroutine tov_h_derivs(hent, y, dydh)
  use constants
  implicit none
  real(8) :: hent ! Enthalpy (not used explicitly; h is Planck's constant)
  real(8) :: y(5)
  real(8) :: dydh(5)
  real(8) :: r
  real(8) :: m
  real(8) :: p
  real(8) :: yy
  real(8) :: ff
  real(8) :: eden
  real(8) :: ed
  real(8) :: cs2
  real(8) :: drdh
  real(8) :: r2, r3, r4, m2
  real(8) :: F
  real(8) :: Q
  real(8) :: L

  !...................................................................
  ! y(1): radius
  ! y(2): gravitational mass
  ! y(3): pressure
  ! y(4): y = r H'/H for the l = 2 tidal perturbation
  ! y(5): f = d ln(omega) / d ln(r) for slow rotation
  ! The equations do not depend on the enthalpy explicitly
  !...................................................................

  r  = y(1)
  m  = y(2)
  p  = y(3)
  yy = y(4)
  ff = y(5)

  ! Get energy density and speed of sound.............................
  call eos_interp(p/pG, eden, cs2)
  ed = eden*c*c*pG

  r2 = r*r
  r3 = r2*r
  r4 = r2*r2
  m2 = m*m
  L  = 1.0d0 / ( 1.0d0 - 2.0d0*m/r )

  ! Radius, mass and pressure.........................................
  drdh    = -r*(r-2.0d0*m) / ( m + 4.0d0*pi*r3*p )
  dydh(1) = drdh
  dydh(2) = (4.0d0*pi)*r2*ed*drdh
  dydh(3) = ed + p

  ! y: dy/dh = (dy/dr)(dr/dh)........................................
  F = ( 1.0d0 - (4.0d0*pi*r2)*( ed - p ) ) * L
  Q = (4.0d0*pi)*( 5.0d0*ed + 9.0d0*p + ((ed+p)/cs2) ) * L &
       -(6.0d0/r2)* L - (4.0d0*m2/r4)                      &
       *( ( 1.0d0 + 4.0d0*pi*r3*p/m )**2 ) * (L*L)
  dydh(4) = ( -(yy*yy/r) -(yy*F/r) -(r*Q) ) * drdh

  ! f: df/dh = (df/dr)(dr/dh)........................................
  dydh(5) = ( -(ff/r)*(ff+3.0d0)                        &
       + (4.0d0+ff)*(4.0d0*pi*r2)*(ed+p)/(r-2.0d0*m) ) * drdh

  return
end subroutine tov_h_derivs

!=====================================================================
! SUBROUTINE HINT: Integrate y() from x1 to x2 with an adaptive
!                  stepper that returns dydx at the new point (RKDP)
!
! Input:
! y()    -- dependent variables at x1
! nvar   -- number of equations
! x1, x2 -- integration interval (x2 < x1 allowed)
! eps    -- required accuracy
! h1     -- first step size to try (magnitude)
! hmin   -- smallest allowed step size (magnitude)
! derivs -- subroutine derivs(x, y, dydx)
! steper -- stepper with the interface of RKDP
!
! Output:
! y()    -- dependent variables at x2
! nok    -- number of good steps
! nbad   -- number of steps with a reduced step size
!=====================================================================
subroutine hint(y, nvar, x1, x2, eps, h1, hmin, nok, nbad, derivs, steper)
  use tov_control
  implicit none
  integer(4)            :: nvar
  integer(4)            :: nok
  integer(4)            :: nbad
  integer(4)            :: nstp
  integer(4), parameter :: maxstp = 100000
  real(8)               :: y(nvar)
  real(8)               :: x1
  real(8)               :: x2
  real(8)               :: eps
  real(8)               :: h1
  real(8)               :: hmin
  real(8)               :: x
  real(8)               :: h
  real(8)               :: hdid
  real(8)               :: hnext
  real(8)               :: yscal(nvar)
  real(8)               :: dydx(nvar)
  real(8), parameter    :: tiny = 1.0d-300 ! P at the surface is ~1e-31 (1/m2)

  external derivs, steper

  x    = x1
  h    = sign(h1, x2-x1)
  nok  = 0
  nbad = 0

  call derivs(x, y, dydx)

  do nstp = 1, maxstp
     yscal = abs(y) + tiny

     ! Do not step past x2............................................
     if ( (x+h-x2)*(x+h-x1) > 0.0d0 ) h = x2 - x

     call steper(y, dydx, nvar, x, h, eps, yscal, hdid, hnext, derivs)
     if ( tov_fail ) return                 ! esymML: step size underflow
     if ( hdid == h ) then
        nok = nok + 1
     else
        nbad = nbad + 1
     end if

     if ( (x-x2)*(x2-x1) >= 0.0d0 ) return

     h = hnext
     if ( abs(h) < hmin ) then
        write(0,*) 'WARNING hint: step size below hmin before the surface at x =', x
        tov_fail = .true.                   ! esymML
        return
     end if
  end do

  write(0,*) 'WARNING hint: too many steps before the surface at x =', x
  tov_fail = .true.                         ! esymML
  return
end subroutine hint
