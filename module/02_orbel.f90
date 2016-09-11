module orbel
! Module for routines in orbel directory
use swift
implicit none

contains
!!
!!
  function orbel_eget(e, m) result(ea)
  !----------------------------------------------------------------------
  !                    ORBEL_EGET.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
  !
  !             Input:
  !                           e ==> eccentricity anomaly. (real scalar)
  !                           m ==> mean anomaly. (real scalar)
  !             Returns:
  !                  orbel_eget ==>  eccentric anomaly. (real scalar)
  !
  !     ALGORITHM: Quartic convergence from Danby
  !     REMARKS: For results very near roundoff,  give it M between
  !           0 and 2*pi. One can condition M before calling EGET
  !           by calling my double precision function MOD2PI(M).
  !           This is not done within the routine to speed it up
  !           and because it works fine even for large M.
  !     AUTHOR: M. Duncan
  !     DATE WRITTEN: May 7,  1992.
  !     REVISIONS: May 21,  1992.  Now have it go through EXACTLY two iterations
  !                with the premise that it will only be called if
  !	         we have an ellipse with e between 0.15 and 0.8
  implicit none

  ! Inputs Only:
  real(rk) :: e, m

  ! Output Only:
  real(rk) :: ea

  ! Internals:
  real(rk) :: x, sm, cm, sx, cx
  real(rk) :: es, ec, f, fp, fpp, fppp, dx

  !----
  ! Executable code

  ! Function to solve Kepler's eqn for E (here called
  ! x) for given e and M. returns value of x.
  ! MAY 21 : FOR e < 0.18 use ESOLMD for speed and sufficient accuracy
  ! MAY 21 : FOR e > 0.8 use EHIE - this one may not converge fast enough.
  call orbel_scget(m, sm, cm)

  ! begin with a guess accurate to order ecc**3
  x = m + e*sm*(1.0_rk + e*(cm + e*(1.0_rk - 1.5_rk*sm**2)))

  ! Go through one iteration for improved estimate
  call orbel_scget(x, sx, cx)
  es = e*sx
  ec = e*cx
  f = x - es - m
  fp = 1.0_rk - ec
  fpp = es
  fppp = ec
  dx = -f/fp
  dx = -f/(fp + 0.5_rk*dx*fpp)
  dx = -f/(fp + 0.5_rk*dx*fpp + dx*dx*fppp/6.0_rk)
  ea = x + dx

  ! Do another iteration.
  ! For m between 0 and 2*pi this seems to be enough to
  ! get near roundoff error for eccentricities between 0 and 0.8
  x = ea
  call orbel_scget(x, sx, cx)
  es = e*sx
  ec = e*cx
  f = x - es - m
  fp = 1.0_rk - ec
  fpp = es
  fppp = ec
  dx = -f/fp
  dx = -f/(fp + 0.5_rk*dx*fpp)
  dx = -f/(fp + 0.5_rk*dx*fpp + dx*dx*fppp/6.0_rk)
  ea = x + dx

  return
  end function orbel_eget
!!
  function orbel_ehie(e, m) result(ea)
  !----------------------------------------------------------------------
  !                    ORBEL_EHIE.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
  !
  !             Input:
  !                           e ==> eccentricity anomaly. (real scalar)
  !                           m ==> mean anomaly. (real scalar)
  !             Returns:
  !              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
  !
  !     ALGORITHM: Use Danby's quartic for 3 iterations.
  !                Eqn. is f(x) = x - e*sin(x+M). Note  that
  !	         E = x + M. First guess is very good for e near 1.
  !	         Need to first get M between 0. and PI and use
  !		 symmetry to return right answer if M between PI and 2PI
  !     REMARKS: Modifies M so that both E and M are in range (0, TWOPI)
  !     AUTHOR: M. Duncan
  !     DATE WRITTEN: May 25, 1992.
  !     REVISIONS:
  implicit none

  ! Inputs Only:
  real(rk) :: e, m

  ! Output Only:
  real(rk) :: ea

  ! Internals:
  integer(ik), parameter :: NMAX = 3
  integer(ik) :: iflag, nper, niter
  real(rk) :: dx, x, sa, ca, esa, eca, f, fp

  !----
  ! Executable code

  ! In this section,  bring M into the range (0, TWOPI) and if
  ! the result is greater than PI,  solve for (TWOPI - M).
  iflag = 0
  nper = m/TWOPI
  m = m - nper*TWOPI
  if(m < 0.0_rk) m = m + TWOPI

  if(m > PI) then

    m = TWOPI - m
    iflag = 1

  end if

  ! Make a first guess that works well for e near 1.
  x = (6.0_rk*m)**(1.0_rk/3.0_rk) - m

  ! Iteration loop
  do niter = 1, NMAX

    call orbel_scget(x + m, sa, ca)
    esa = e*sa
    eca = e*ca
    f = x - esa
    fp = 1.0_rk - eca
    dx = -f/fp
    dx = -f/(fp + 0.5_rk*dx*esa)
    dx = -f/(fp + 0.5_rk*dx*(esa + eca*dx/3.0_rk))
    x = x + dx

  end do

  ea = m + x

  if(iflag == 1) then

    ea = TWOPI - ea
    m = TWOPI - m

  end if

  return
  end function orbel_ehie
!!
  function orbel_ehybrid(e, m) result(ea)
  !----------------------------------------------------------------------
  !                    ORBEL_EHYBRID.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
  !
  !             Input:
  !                           e ==> eccentricity anomaly. (real scalar)
  !                           m ==> mean anomaly. (real scalar)
  !             Returns:
  !              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
  !
  !     ALGORITHM: For e < 0.18 uses fast routine ESOLMD
  !	         For larger e but less than 0.8,  uses EGET
  !	         For e > 0.8 uses EHIE
  !     REMARKS: Only EHIE brings M and E into range (0, TWOPI)
  !     AUTHOR: M. Duncan
  !     DATE WRITTEN: May 25, 1992.
  !     REVISIONS: 2/26/93 hfl
  implicit none

  ! Inputs Only:
  real(rk) :: e, m

  ! Output Only:
  real(rk) :: ea

  !----
  ! Executable code

  if(e < 0.18_rk) then

    ea = orbel_esolmd(e, m)

  else

    if(e <= 0.8_rk) then

      ea = orbel_eget(e, m)

    else

      ea = orbel_ehie(e, m)

    end if

  end if

  return
  end function orbel_ehybrid
!!
  subroutine orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm, x, y, z, vx, vy, vz)
  !----------------------------------------------------------------------------
  !                          ORBEL_EL2XV.F
  !----------------------------------------------------------------------------
  !     PURPOSE: To compute cartesian positions and velocities given
  !               central mass,  ialpha ( = +1 for hyp.,  0 for para. and
  !               -1 for ellipse),  and orbital elements.
  !       input:
  !            gm       ==> G times central mass (real scalar)
  !	     ialpha   ==> conic section type ( see PURPOSE,  integer(ik) :: scalar)
  !	     a        ==> semi-major axis or pericentric distance if a parabola
  !                          (real scalar)
  !            e        ==> eccentricity (real scalar)
  !            inc      ==> inclination  (real scalar)
  !            capom    ==> longitude of ascending node (real scalar)
  !	     omega    ==> argument of perihelion (real scalar)
  !	     capm     ==> mean anomoly(real scalar)
  !
  !       Output:
  !            x, y, z    ==>  position of object (real scalars)
  !            vx, vy, vz ==>  velocity of object (real scalars)
  !
  !     ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."
  !     REMARKS: All angles are in RADIANS
  !
  !     AUTHOR:  M. Duncan.
  !     DATE WRITTEN:  May 11,  1992.
  !     REVISIONS: May 26 - now use better Kepler solver for ellipses
  !                 and hyperbolae called EHYBRID.F and FHYBRID.F
  implicit none

  ! Inputs Only:
  integer(ik) :: ialpha
  real(rk) :: gm, a, e, inc, capom, omega, capm

  ! Outputs:
  real(rk) :: x, y, z, vx, vy, vz

  ! Internals:
  real(rk) :: cape, capf, zpara, em1
  real(rk) :: sp, cp, so, co, si, ci
  real(rk) :: d11, d12, d13, d21, d22, d23
  real(rk) :: scap, ccap, shcap, chcap
  real(rk) :: sqe, sqgma, xfac1, xfac2, ri, vfac1, vfac2

  !----
  ! Executable code

  if(e < 0.0_rk) then

    write(*,'(a)') 'ERROR in orbel_el2xv: e < 0.0, setting e = 0.0'
    e = 0.0_rk

  end if

  ! check for inconsistencies between ialpha and e
  em1 = e - 1.0_rk
  if(((ialpha == 0) .and. (abs(em1) > TINY_NUMBER)) .or. ((ialpha < 0) .and. (e > 1.0_rk)) .or. ((ialpha > 0) .and. (e < 1.0_rk))) then

    write(*,'(a)')         ' ERROR in orbel_el2xv: ialpha and e are inconsistent'
    write(*,'(a,i3)')      '                       ialpha = ', ialpha
    write(*,'(a,1pe14.6)') '                            e = ', e

  end if

  ! Generate rotation matrices (on p. 42 of Fitzpatrick)
  call orbel_scget(omega, sp, cp)
  call orbel_scget(capom, so, co)
  call orbel_scget(inc, si, ci)
  d11 = cp*co - sp*so*ci
  d12 = cp*so + sp*co*ci
  d13 = sp*si
  d21 = -sp*co - cp*so*ci
  d22 = -sp*so + cp*co*ci
  d23 = cp*si

  ! Get the other quantities depending on orbit type (i.e. IALPHA)
  if(ialpha == -1) then

    cape = orbel_ehybrid(e, capm)
    call orbel_scget(cape, scap, ccap)
    sqe = sqrt(1.0_rk - e*e)
    sqgma = sqrt(gm*a)
    xfac1 = a*(ccap - e)
    xfac2 = a*sqe*scap
    ri = 1.0_rk/(a*(1.0_rk - e*ccap))
    vfac1 = -ri*sqgma*scap
    vfac2 = ri*sqgma*sqe*ccap

  end if

  if(ialpha == 1) then

    capf = orbel_fhybrid(e, capm)
    call orbel_schget(capf, shcap, chcap)
    sqe = sqrt(e*e - 1.0_rk)
    sqgma = sqrt(gm*a)
    xfac1 = a*(e - chcap)
    xfac2 = a*sqe*shcap
    ri = 1.0_rk/(a*(e*chcap - 1.0_rk))
    vfac1 = -ri*sqgma*shcap
    vfac2 = ri*sqgma*sqe*chcap

  end if

  if(ialpha == 0) then

    zpara = orbel_zget(capm)
    sqgma = sqrt(2.0_rk*gm*a)
    xfac1 = a*(1.0_rk - zpara*zpara)
    xfac2 = 2.0_rk*a*zpara
    ri = 1.0_rk/(a*(1.0_rk + zpara*zpara))
    vfac1 = -ri*sqgma*zpara
    vfac2 = ri*sqgma

  end if

  x =  d11*xfac1 + d21*xfac2
  y =  d12*xfac1 + d22*xfac2
  z =  d13*xfac1 + d23*xfac2
  vx = d11*vfac1 + d21*vfac2
  vy = d12*vfac1 + d22*vfac2
  vz = d13*vfac1 + d23*vfac2

  return
  end subroutine orbel_el2xv
!!
  function orbel_esolmd(e, m) result(ea)
  !----------------------------------------------------------------------
  !                    ORBEL_ESOLMD.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
  !
  !             Input:
  !                           e ==> eccentricity anomaly. (real scalar)
  !                           m ==> mean anomaly. (real scalar)
  !             Returns:
  !                orbel_esolmd ==>  eccentric anomaly. (real scalar)
  !
  !     ALGORITHM: Some sort of quartic convergence from Wisdom.
  !     REMARKS: ONLY GOOD FOR SMALL ECCENTRICITY SINCE IT ONLY
  !         ITERATES ONCE. (GOOD FOR PLANET CALCS.)
  !      	  ALSO DOES NOT PUT M OR E BETWEEN 0. AND 2*PI
  !     INCLUDES: needs SCGET.F
  !     AUTHOR: M. Duncan
  !     DATE WRITTEN: May 7,  1992.
  !     REVISIONS: 2/26/93 hfl
  implicit none

  ! Inputs Only:
  real(rk) :: e, m

  ! Output Only:
  real(rk) :: ea

  ! Internals:
  real(rk) :: x, sm, cm, sx, cx
  real(rk) :: es, ec, f, fp, fpp, fppp, dx

  !----
  ! Executable code

  ! Function to solve Kepler's eqn for E (here called x) for given e and M. returns value of x.
  call orbel_scget(m, sm, cm)
  x = m + e*sm*(1.0_rk + e*(cm + e*(1.0_rk - 1.5_rk*sm*sm)))

  call orbel_scget(x, sx, cx)
  es = e*sx
  ec = e*cx
  f = x - es - m
  fp = 1.0_rk - ec
  fpp = es
  fppp = ec
  dx = -f/fp
  dx = -f/(fp + 0.5_rk*dx*fpp)
  dx = -f/(fp + 0.5_rk*dx*fpp + dx*dx*fppp/6.0_rk)

  ea = x + dx

  return
  end function orbel_esolmd
!!
  function orbel_fget(e, capn) result(ea)
  !----------------------------------------------------------------------
  !                    ORBEL_FGET.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
  !
  !             Input:
  !                           e ==> eccentricity anomaly. (real scalar)
  !                        capn ==> hyperbola mean anomaly. (real scalar)
  !             Returns:
  !                  orbel_fget ==>  eccentric anomaly. (real scalar)
  !
  !     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
  !           Cel. Mech. ".  Quartic convergence from Danby's book.
  !     REMARKS:
  !     AUTHOR: M. Duncan
  !     DATE WRITTEN: May 11, 1992.
  !     REVISIONS: 2/26/93 hfl
  implicit none

  ! Inputs Only:
  real(rk) :: e, capn

  ! Output Only:
  real(rk) :: ea

  ! Internals:
  integer(ik), parameter :: IMAX = 10
  integer(ik) :: i
  real(rk) :: tmp, x, shx, chx
  real(rk) :: esh, ech, f, fp, fpp, fppp, dx

  !----
  ! Executable code

  ! Function to solve "Kepler's eqn" for F (here called x) for given e and CAPN

  ! begin with a guess proposed by Danby
  if(capn < 0.0_rk) then

    tmp = -2.0_rk*capn/e + 1.8_rk
    x = -log(tmp)

  else

    tmp = 2.0_rk*capn/e + 1.8_rk
    x = log(tmp)

  end if

  ea = x

  do i = 1, IMAX

    call orbel_schget(x, shx, chx)
    esh = e*shx
    ech = e*chx
    f = esh - x - capn
    !write(*,*) 'i, x, f: ', i, x, f
    fp = ech - 1.0_rk
    fpp = esh
    fppp = ech
    dx = -f/fp
    dx = -f/(fp + 0.5_rk*dx*fpp)
    dx = -f/(fp + 0.5_rk*dx*fpp + dx*dx*fppp/6.0_rk)
    ea = x + dx

    ! If we have converged here there's no point in going on
    if(abs(dx) <= TINY_NUMBER) return
    x = ea

  end do

  write(*,*) 'FGET: RETURNING WITHOUT COMPLETE CONVERGENCE'

  return
  end function orbel_fget
!!
  function orbel_fhybrid(e, capn) result(ea)
  !----------------------------------------------------------------------
  !                    ORBEL_FHYBRID.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
  !
  !             Input:
  !                           e ==> eccentricity anomaly. (real scalar)
  !                        capn ==> hyperbola mean anomaly. (real scalar)
  !             Returns:
  !               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
  !
  !     ALGORITHM: For abs(N) < 0.636*ecc -0.6 ,  use FLON
  !	         For larger N,  uses FGET
  !     REMARKS:
  !     AUTHOR: M. Duncan
  !     DATE WRITTEN: May 26, 1992.
  !     REVISIONS:
  !     REVISIONS: 2/26/93 hfl
  implicit none

  ! Inputs Only:
  real(rk) :: e, capn

  ! Output Only:
  real(rk) :: ea

  ! Internals:
  real(rk) :: abn

  !----
  ! Executable code

  abn = capn
  if(capn < 0.0_rk) abn = -abn

  if(abn < 0.636_rk*e - 0.6_rk) then

    ea = orbel_flon(e, capn)

  else

    ea = orbel_fget(e, capn)

  end if

  return
  end function orbel_fhybrid
!!
  function orbel_flon(e, capn) result(ea)
  !----------------------------------------------------------------------
  !                    ORBEL_FLON.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
  !
  !             Input:
  !                           e ==> eccentricity anomaly. (real scalar)
  !                        capn ==> hyperbola mean anomaly. (real scalar)
  !             Returns:
  !                  orbel_flon ==>  eccentric anomaly. (real scalar)
  !
  !     ALGORITHM: Uses power series for N in terms of F and Newton, s method
  !     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
  !     AUTHOR: M. Duncan
  !     DATE WRITTEN: May 26, 1992.
  !     REVISIONS:
  implicit none

  ! Inputs Only:
  real(rk) :: e, capn

  ! Output Only:
  real(rk) :: ea

  ! Internals:
  integer(ik), parameter :: IMAX = 10
  real(rk), parameter :: a11 = 156.0_rk, a9 = 17160.0_rk, a7 = 1235520.0_rk, a5 = 51891840.0_rk, a3 = 1037836800.0_rk
  real(rk), parameter :: b11 = 11.0_rk*a11, b9 = 9.0_rk*a9, b7 = 7.0_rk*a7, b5 = 5.0_rk*a5, b3 = 3.0_rk*a3
  integer(ik) :: i, iflag
  real(rk) :: a, b, sq, biga, bigb
  real(rk) :: x, x2
  real(rk) :: f, fp, dx
  real(rk) :: diff
  real(rk) :: a0, a1, b1

  !----
  ! Executable code


  ! Function to solve "Kepler's eqn" for F (here called x) for given e and CAPN. Only good for smallish CAPN

  iflag = 0
  if(capn < 0.0_rk) then

    iflag = 1
    capn = -capn

  end if

  a1 = 6227020800.0_rk*(1.0_rk - 1.0_rk/e)
  a0 = -6227020800.0_rk*capn/e
  b1 = a1

  ! Set iflag nonzero if capn < 0.,  in which case solve for -capn
  ! and change the sign of the final answer for F.
  ! Begin with a reasonable guess based on solving the cubic for small F
  a = 6.0_rk*(e - 1.0_rk)/e
  b = -6.0_rk*capn/e
  sq = sqrt(0.25_rk*b*b + a*a*a/27.0_rk)
  biga = (sq - 0.5_rk*b)**(1.0_rk/3.0_rk)
  bigb = -(sq + 0.5_rk*b)**(1.0_rk/3.0_rk)
  x = biga + bigb
  !write(*,'(a,1pe14.6)') 'cubic = ', x**3 + a*x + b
  ea = x

  ! If capn is tiny (or zero) no need to go further than cubic even for e = 1.0
  if(capn >= TINY_NUMBER) then

    do i = 1, IMAX

      x2 = x*x
      f = a0 + x*(a1 + x2*(a3 + x2*(a5 + x2*(a7 + x2*(a9 + x2*(a11 + x2))))))
      fp = b1 + x2*(b3 + x2*(b5 + x2*(b7 + x2*(b9 + x2*(b11 + 13.0_rk*x2)))))
      dx = -f/fp
      !write(*,'(a)')                    'i, dx, x, f:'
      !write(*,'(1x,i3,3(2x,1p1e22.15))') i, dx, x, f
      ea = x + dx

      ! If we have converged here there's no point in going on
      if(abs(dx) <= TINY_NUMBER) exit
      x = ea

    end do

  end if

  ! Check if capn was originally negative
  if(iflag == 1) then

    ea = -ea
    capn = -capn

  end if

  ! Abnormal return here - we've gone thru the loop IMAX times without convergence
  if(i >= imax) then

    write(*,'(a)') 'FLON: RETURNING WITHOUT COMPLETE CONVERGENCE'
    diff = e*sinh(ea) - ea - capn
    write(*,'(a)') 'N, F, ecc*sinh(F) - F - N: '
    write(*,'(3(1pe14.6))') capn, ea, diff

  end if

  return
  end function orbel_flon
!!
  subroutine orbel_scget(angle, sx, cx)
  !----------------------------------------------------------------------
  !	                  ORBEL_SCGET.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Given an angle,  efficiently compute sin and cos.
  !
  !        Input:
  !             angle ==> angle in radians (real scalar)
  !
  !        Output:
  !             sx    ==>  sin(angle)  (real scalar)
  !             cx    ==>  cos(angle)  (real scalar)
  !
  !     ALGORITHM: Obvious from the code
  !     REMARKS: The HP 700 series won't return correct answers for sin
  !       and cos if the angle is bigger than 3e7. We first reduce it
  !       to the range [0, 2pi) and use the sqrt rather than cos (it's faster)
  !       BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES!
  !     AUTHOR:  M. Duncan.
  !     DATE WRITTEN:  May 6,  1992.
  !     REVISIONS:
  implicit none

  ! Inputs Only:
  real(rk) :: angle

  ! Output:
  real(rk) :: sx, cx

  ! Internals:
  integer(ik) :: nper
  real(rk) :: x

  !-----------------
  ! Executable code
  !-----------------

  nper = angle/TWOPI
  x = angle - nper*TWOPI
  if(x < 0.0_rk) x = x + TWOPI

  sx = sin(x)
  cx = sqrt(1.0_rk - sx*sx)
  if((x > PIBY2) .and. (x < PI3BY2)) cx = -cx

  return
  end subroutine orbel_scget
!!
  subroutine orbel_schget(angle, shx, chx)
  !----------------------------------------------------------------------
  !	                  ORBEL_SCHGET.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Given an angle,  efficiently compute sinh and cosh.
  !
  !        Input:
  !             angle ==> angle in radians (real scalar)
  !
  !        Output:
  !             shx    ==>  sinh(angle)  (real scalar)
  !             chx    ==>  cosh(angle)  (real scalar)
  !
  !     ALGORITHM: Obvious from the code
  !     REMARKS: Based on the routine SCGET for sine's and cosine's.
  !       We use the sqrt rather than cosh (it's faster)
  !       BE SURE THE ANGLE IS IN RADIANS AND IT CAN'T BE LARGER THAN 300
  !       OR OVERFLOWS WILL OCCUR!
  !     AUTHOR:  M. Duncan.
  !     DATE WRITTEN:  May 6,  1992.
  !     REVISIONS:
  implicit none

  ! Inputs Only:
  real(rk) :: angle

  ! Output:
  real(rk) :: shx, chx

  !-----------------
  ! Executable code
  !-----------------

  shx = sinh(angle)
  chx = sqrt(1.0_rk + shx*shx)

  return
  end subroutine orbel_schget
!!
  subroutine orbel_xv2aeq(x, y, z, vx, vy, vz, gmsum, ialpha, a, e, q)
  !----------------------------------------------------------------------
  !	               ORBEL_XV2AEQ.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Given the cartesian position and velocity of an orbit, 
  !       compute the osculating orbital elements a,  e,  and q only.
  !
  !       input:
  !            x, y, z    ==>  position of object (real scalars)
  !            vx, vy, vz ==>  velocity of object (real scalars)
  !            gmsum       ==> G*(M1+M2) (real scalar)
  !
  !       Output:
  !	     ialpha   ==> conic section type ( see PURPOSE,  integer scalar)
  !	     a        ==> semi-major axis or pericentric distance if a parabola
  !                          (real scalar)
  !            e        ==> eccentricity (real scalar)
  !            q        ==> perihelion distance (real scalar); q = a(1 - e)
  !
  !     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech."
  !     REMARKS: Based on M. Duncan's orbel_xv2el.f
  !      This routine is generally applied to study (hyperbolic) close
  !       encounters of test particles with planets.
  !     AUTHOR:  L. Dones
  !     DATE WRITTEN:  February 24,  1994
  !     REVISIONS:
  implicit none

  ! Inputs Only:
  real(rk) :: x, y, z, vx, vy, vz, gmsum

  ! Outputs
  integer(ik) :: ialpha
  real(rk) :: a, e, q

  ! Internals:
  real(rk) :: hx, hy, hz, h2, r, v2, energy, fac

  !-----------------
  ! Executable code
  !-----------------

  ! Compute the angular momentum H, and thereby the inclination INC
  hx = y*vz - z*vy
  hy = z*vx - x*vz
  hz = x*vy - y*vx
  h2 = hx*hx + hy*hy + hz*hz

  ! Compute the radius R and velocity squared V2, and the dot
  ! product RDOTV, the energy per unit mass ENERGY.
  r = sqrt(x*x + y*y + z*z)
  v2 = vx*vx + vy*vy + vz*vz
  energy = 0.5_rk*v2 - gmsum/r

  ! Determine type of conic section and label it via IALPHA
  if(abs(energy*r/gmsum) < sqrt(TINY_NUMBER)) then

    ialpha = 0

  else

    if(energy < 0.0_rk) ialpha = -1
    if(energy > 0.0_rk) ialpha = 1

  end if

  ! Depending on the conic type, determine the remaining elements

  ! ELLIPSE
  if(ialpha == -1) then

    a = -0.5_rk*gmsum/energy
    fac = 1.0_rk - h2/(gmsum*a)

    if(fac > TINY_NUMBER) then

      e = sqrt(fac)

    else

      e = 0.0_rk

    end if

    q = a*(1.0_rk - e)

  end if


  ! HYPERBOLA
  if(ialpha == 1) then

    a = 0.5_rk*gmsum/energy
    fac = h2/(gmsum*a)

    if(fac > TINY_NUMBER) then

      e = sqrt(1.0_rk + fac)

      ! have to insert minus sign in expression for q because this code
      ! takes a > 0, even for a hyperbola
      q = -a*(1.0_rk - e)

    else

      ! we only get here if a hyperbola is essentially a parabola
      ! so we calculate e accordingly to avoid singularities
      e = 1.0_rk
      q = 0.5_rk*h2/gmsum

    end if

  end if

  ! PARABOLA: (NOTE - in this case "a",  which is formally infinite,
  !            is arbitrarily set equal to the pericentric distance q).
  if(ialpha == 0) then

    a =  0.5_rk*h2/gmsum
    e = 1.0_rk
    q = a

  end if

  return
  end subroutine orbel_xv2aeq
!!
  subroutine orbel_xv2aqt(x, y, z, vx, vy, vz, gmsum, ialpha, a, q, capm, tperi)
  !-----------------------------------------------------------------------
  !	               ORBEL_XV2AQT.F
  !-----------------------------------------------------------------------
  !     PURPOSE:  Given the cartesian position and velocity of an orbit, 
  !       compute the osculating orbital elements a,  e,  and q only.
  !
  !       input:
  !            x, y, z    ==>  position of object (real scalars)
  !            vx, vy, vz ==>  velocity of object (real scalars)
  !            gmsum       ==> G*(M1+M2) (real scalar)
  !
  !       Output:
  !	     ialpha   ==> conic section type ( see PURPOSE,  integer scalar)
  !	     a        ==> semi-major axis or pericentric distance if a parabola
  !                          (real scalar)
  !            q        ==> perihelion distance (real scalar); q = a(1 - e)
  !          capm       ==> mean anomoly(real scalar)
  !          tperi      ==> time to next or last perihelion,  which ever is less
  !                         (real scalar)
  !
  !     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech."
  !     REMARKS: Based on M. Duncan's orbel_xv2el.f
  !      This routine is generally applied to study (hyperbolic) close
  !       encounters of test particles with planets.
  !     AUTHOR:  Hal Levison
  !     DATE WRITTEN:  8/7/01
  !     REMARKS: The tperi may not be correct for parabolic orbits.
  !              I Think it is OK but beware!
  !     REVISIONS:
  implicit none

  ! Inputs Only:
  real(rk) :: x, y, z, vx, vy, vz, gmsum

  ! Outputs
  integer(ik) :: ialpha
  real(rk) :: a, q, capm, tperi

  ! Internals:
  real(rk) :: hx, hy, hz, h2, r, v2, energy, fac, vdotr, cape, e
  real(rk) :: capf, tmpf, meanmo, face, w

  !----
  ! Executable code

  ! Compute the angular momentum H, and thereby the inclination INC
  hx = y*vz - z*vy
  hy = z*vx - x*vz
  hz = x*vy - y*vx
  h2 = hx*hx + hy*hy + hz*hz

  ! Compute the radius R and velocity squared V2, and the dot
  ! product RDOTV, the energy per unit mass ENERGY.
  r = sqrt(x*x + y*y + z*z)
  v2 = vx*vx + vy*vy + vz*vz
  energy = 0.5_rk*v2 - gmsum/r
  vdotr = x*vx + y*vy + z*vz

  ! Determine type of conic section and label it via IALPHA
  if(abs(energy*r/gmsum) < sqrt(TINY_NUMBER)) then

    ialpha = 0

  else

    if(energy < 0.0_rk) ialpha = -1
    if(energy > 0.0_rk) ialpha = 1

  end if

  ! Depending on the conic type,  determine the remaining elements

  ! ELLIPSE
  if(ialpha == -1) then

    a = -0.5_rk*gmsum/energy
    fac = 1.0_rk - h2/(gmsum*a)

    if(fac > TINY_NUMBER) then

      e = sqrt(fac)
      face = (a - r)/(a*e)

      ! Apr. 16/93 : watch for case where face is slightly outside unity
      if(face > 1.0_rk) then

        cape = 0.0_rk

      else

        if(face > -1.0_rk) then

          cape = acos(face)

        else

          cape = PI

        end if

      end if

      if(vdotr < 0.0_rk) cape = TWOPI - cape

    else

      e = 0.0_rk
      cape = 0.0_rk

    end if

    capm = cape - e*sin(cape)
    q = a*(1.0_rk - e)

  end if

  ! HYPERBOLA
  if(ialpha == 1) then

    a = 0.5_rk*gmsum/energy
    fac = h2/(gmsum*a)

    if(fac > TINY_NUMBER) then

      e = sqrt(1.0_rk + fac)

      ! have to insert minus sign in expression for q because this code
      ! takes a > 0,  even for a hyperbola
      q = -a*(1.0_rk - e)

      tmpf = (a + r)/(a*e)
      if(tmpf < 1.0_rk) tmpf = 1.0_rk

      capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_rk))
      if(vdotr < 0.0_rk) capf = - capf

    else

      ! we only get here if a hyperbola is essentially a parabola
      ! so we calculate e accordingly to avoid singularities
      e = 1.0_rk
      q = 0.5_rk*h2/gmsum

      tmpf = (a + r)/(a*e)
      capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_rk))

    end if

    capm = e*sinh(capf) - capf

  end if

  ! PARABOLA: (NOTE - in this case "a",  which is formally infinite,
  !         is arbitrarily set equal to the pericentric distance q).
  if(ialpha == 0) then

    a = 0.5_rk*h2/gmsum
    e = 1.0_rk
    q = a
    w = acos(2.0_rk*a/r - 1.0_rk)
    if(vdotr < 0.0_rk) w = TWOPI - w
    tmpf = tan(0.5_rk*w)
    capm = tmpf*(1.0_rk + tmpf*tmpf/3.0_rk)

  end if

  meanmo = sqrt(gmsum/a**3)
  if((capm < PI) .or. (ialpha >= 0)) then

    tperi = -1.0_rk*capm/meanmo

  else

    tperi = -1.0_rk*(capm - TWOPI)/meanmo

  end if

  return
  end subroutine orbel_xv2aqt
!!
  subroutine orbel_xv2el(x, y, z, vx, vy, vz, gmsum, ialpha, a, e, inc, capom, omega, capm)
  !----------------------------------------------------------------------
  !	               ORBEL_XV2EL.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Given the cartesian position and velocity of an orbit, 
  !       compute the osculating orbital elements.
  !
  !       input:
  !            x, y, z    ==>  position of object (real scalars)
  !            vx, vy, vz ==>  velocity of object (real scalars)
  !            gmsum       ==> G*(M1+M2) (real scalar)
  !
  !       Output:
  !	     ialpha   ==> conic section type ( see PURPOSE,  integer scalar)
  !	     a        ==> semi-major axis or pericentric distance if a parabola
  !                          (real scalar)
  !            e        ==> eccentricity (real scalar)
  !            inc      ==> inclination  (real scalar)
  !            capom    ==> longitude of ascending node (real scalar)
  !	     omega    ==> argument of perihelion (real scalar)
  !	     capm     ==> mean anomoly(real scalar)
  !
  !     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech."
  !     REMARKS:  If the inclination INC is less than TINY_NUMBER,  we
  !       arbitrarily choose the longitude of the ascending node LGNODE
  !       to be 0.0 (so the ascending node is then along the X axis).  If
  !       the  eccentricity E is less than SQRT(TINY_NUMBER),  we arbitrarily
  !       choose the argument of perihelion to be 0.
  !     AUTHOR:  M. Duncan.
  !     DATE WRITTEN:  May 8, 1992.
  !     REVISIONS: 2/4/2K
  implicit none

  ! Inputs Only:
  real(rk) :: x, y, z, vx, vy, vz, gmsum

  ! Outputs
  integer(ik) :: ialpha
  real(rk) :: a, e, inc, capom, omega, capm

  ! Internals:
  real(rk) :: hx, hy, hz, h2, h, r, v2, v, vdotr, energy, fac, face, cape, capf, tmpf
  real(rk) :: cw, sw, w, u

  !----
  ! Executable code

  ! Compute the angular momentum H, and thereby the inclination INC.
  hx = y*vz - z*vy
  hy = z*vx - x*vz
  hz = x*vy - y*vx
  h2 = hx*hx + hy*hy +hz*hz
  h = sqrt(h2)
  if(hz > h) hz = h ! Hal's fix
  inc = acos(hz/h)

  ! Compute longitude of ascending node CAPOM and the argument of latitude u.
  fac = sqrt(hx**2 + hy**2)/h

  if(fac < TINY_NUMBER) then

    capom = 0.0_rk
    u = atan2(y, x)
    if(abs(inc - PI) < 10.0_rk*TINY_NUMBER) u = -u

  else

    capom = atan2(hx, -hy)
    u = atan2(z/sin(inc), x*cos(capom) + y*sin(capom))

  end if

  if(capom < 0.0_rk) capom = capom + TWOPI
  if(u < 0.0_rk) u = u + TWOPI

  ! Compute the radius R and velocity squared V2,  and the dot
  ! product RDOTV, the energy per unit mass ENERGY.
  r = sqrt(x*x + y*y + z*z)
  v2 = vx*vx + vy*vy + vz*vz
  v = sqrt(v2)
  vdotr = x*vx + y*vy + z*vz
  energy = 0.5_rk*v2 - gmsum/r

  ! Determine type of conic section and label it via IALPHA
  if(abs(energy*r/gmsum) < sqrt(TINY_NUMBER)) then

    ialpha = 0

  else

    if(energy < 0.0_rk) ialpha = -1
    if(energy > 0.0_rk) ialpha = 1

  end if

  ! Depending on the conic type,  determine the remaining elements

  ! ELLIPSE
  if(ialpha == -1) then

    a = -0.5_rk*gmsum/energy
    fac = 1.0_rk - h2/(gmsum*a)

    if(fac > TINY_NUMBER) then

      e = sqrt(fac)
      face = (a - r)/(a*e)

      ! Apr. 16/93 : watch for case where face is slightly outside unity
      if(face > 1.0_rk) then

        cape = 0.0_rk

      else

        if(face > -1.0_rk) then

          cape = acos(face)

        else

          cape = PI

        end if

      end if

      if(vdotr < 0.0_rk) cape = TWOPI - cape
      cw = (cos(cape) - e)/(1.0_rk - e*cos(cape))
      sw = sqrt(1.0_rk - e*e)*sin(cape)/(1.0_rk - e*cos(cape))
      w = atan2(sw, cw)
      if(w < 0.0_rk) w = w + TWOPI

    else

      e = 0.0_rk
      w = u
      cape = u

    end if

    capm = cape - e*sin(cape)
    omega = u - w
    if(omega < 0.0_rk) omega = omega + TWOPI
    omega = omega - TWOPI*int(omega/TWOPI)

  end if

  ! HYPERBOLA:
  if(ialpha == 1) then

    a = 0.5_rk*gmsum/energy
    fac = h2/(gmsum*a)

    if(fac > TINY_NUMBER) then

      e = sqrt(1.0_rk + fac)
      tmpf = (a + r)/(a*e)
      if(tmpf < 1.0_rk) tmpf = 1.0_rk
      capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_rk))
      if(vdotr < 0.0_rk) capf = -capf
      cw = (e - cosh(capf))/(e*cosh(capf) - 1.0_rk)
      sw = sqrt(e*e - 1.0_rk)*sinh(capf)/(e*cosh(capf) - 1.0_rk)
      w = atan2(sw, cw)
      if(w < 0.0_rk) w = w + TWOPI

    else

      ! we only get here if a hyperbola is essentially a parabola
      ! so we calculate e and w accordingly to avoid singularities
      e = 1.0_rk
      tmpf = 0.5_rk*h2/gmsum
      w = acos(2.0_rk*tmpf/r - 1.0_rk)
      if(vdotr < 0.0_rk) w = TWOPI - w
      tmpf = (a + r)/(a*e)
      capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_rk))

    end if

    capm = e*sinh(capf) - capf
    omega = u - w
    if(omega < 0.0_rk) omega = omega + TWOPI
    omega = omega - TWOPI*int(omega/TWOPI)

  end if

  ! PARABOLA: (NOTE - in this case we use "a" to mean pericentric distance)
  if(ialpha == 0) then

    a = 0.5_rk*h2/gmsum
    e = 1.0_rk
    w = acos(2.0_rk*a/r - 1.0_rk)
    if(vdotr < 0.0_rk) w = TWOPI - w
    tmpf = tan(0.5_rk*w)
    capm = tmpf*(1.0_rk + tmpf*tmpf/3.0_rk)
    omega = u - w
    if(omega < 0.0_rk) omega = omega + TWOPI
    omega = omega - TWOPI*int(omega/TWOPI)

  end if

  return
  end subroutine orbel_xv2el
!!
  function orbel_zget(q) result(ea)
  !----------------------------------------------------------------------
  !                    ORBEL_ZGET.F
  !----------------------------------------------------------------------
  !     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola
  !          given Q (Fitz. notation.)
  !
  !             Input:
  !                           q ==>  parabola mean anomaly. (real scalar)
  !             Returns:
  !                  orbel_zget ==>  eccentric anomaly. (real scalar)
  !
  !     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
  !     REMARKS: For a parabola we can solve analytically.
  !     AUTHOR: M. Duncan
  !     DATE WRITTEN: May 11, 1992.
  !     REVISIONS: May 27 - corrected it for negative Q and use power
  !	      series for small Q.
  implicit none

  ! Inputs Only:
  real(rk) :: q

  ! Output Only:
  real(rk) :: ea

  ! Internals:
  integer(ik) :: iflag
  real(rk) :: x, tmp

  !----
  ! Executable code

  iflag = 0
  if(q < 0.0_rk) then

    iflag = 1
    q = -q

  end if

  if(q < 1.0e-3_rk) then

    ea = q*(1.0_rk - (q*q/3.0_rk)*(1.0_rk - q*q))

  else

    x = 0.5_rk*(3.0_rk*q + sqrt(9.0_rk*q**2 + 4.0_rk))
    tmp = x**(1.0_rk/3.0_rk)
    ea = tmp - 1.0_rk/tmp

  end if

  if(iflag == 1) then

    ea = -ea
    q = -q

  end if

  return
  end function orbel_zget
!!
!!
end module orbel