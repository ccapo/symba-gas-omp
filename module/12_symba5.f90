module symba5
! Module for routines in symba5 directory
!
! Remarks:  Copied from symba5.inc
! Author:  Hal Levison
! Date:    3/20/97
! Last revision:
use swift
use orbel
use util
use coord
use mvs
use obl
use anal
use helio
use discard
use rmvs
implicit none

! Maximum number of encounters
integer(ik), parameter :: NENMAX = 262144 ! must be less than 2^18 = 262144

! Ratio of the number of time steps in the adjoining shells
integer(ik), parameter :: NTENC = 3

! Scale factor for hill's sphere to take shorter time step
real(rk), parameter :: RHSCALE = 6.5_rk

! Ratio of shell radii squared
!real(rk), parameter :: rshell = 0.48075_rk ! rshell ~ ntenc^(-2/3)
!real(rk), parameter :: rshell = real(NTENC, rk)**(-2.0_rk/3.0_rk) ! rshell ~ ntenc^(-2/3)
real(rk), parameter :: RSHELL = 0.480749856769136133_rk ! rshell ~ ntenc^(-2/3)

! Maximum number of perigee passages about a massive body before merging satellite with massive body
! This applies to those particles that become bound (however temporarily) during a close encounter,
! and limits the damage that such trapped particles do to the overall performance of the code.
integer(ik), parameter :: NBOUNDMAX = 5

contains
!!
!!
  subroutine symba5_chk2(rhill, xh, yh, zh, vxh, vyh, vzh, dt, irec, icflg)
  !*************************************************************************
  !                            SYMBA5_CHK2.F
  !*************************************************************************
  ! This subroutine checks to see if there are encounters
  !
  !             Input:
  !                 rhill         ==>  Radius of hill sphere (real array)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 ip1, ip2       ==>  The two bodies to check (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 xh, yh, zh      ==>  initial position in helio coord
  !                                    (real arrays)
  !                 vxh, vyh, vzh   ==>  initial velocity in helio coord
  !                                    (real arrays)
  !                 dt            ==>  time step  (real scalor)
  !                 irec          ==>  current recursion level (int scalar)
  !             Output:
  !                 icflg         ==> ecounter?  = 1 Yes
  !                                              = 0 No (int scalar)
  !                 lvdotr        ==> = .true. if i & j are approaching
  !                                   = .false if i & j are receding
  !                                     (logical*1 scalar)
  !
  ! Remarks: Based on plh_chk.f.  Same as symba_chk.f
  ! Authors:  Hal Levison
  ! Date: 03/20/97
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: irec
  real(rk) :: xh(2), yh(2), zh(2), dt
  real(rk) :: vxh(2), vyh(2), vzh(2), rhill(2)

  ! Outputs
  integer(ik) :: icflg

  ! Internals
  real(rk) :: r2crit, r2min
  real(rk) :: xr, yr, zr, vxr, vyr, vzr
  real(rk) :: dr2, dv2, vdotr, tmin

  !-----------------
  ! Executable code
  !-----------------

  r2crit = ((rhill(1) + rhill(2))*RHSCALE*RSHELL**irec)**2

  xr = xh(2) - xh(1)
  yr = yh(2) - yh(1)
  zr = zh(2) - zh(1)
  dr2 = xr**2 + yr**2 + zr**2

  vxr = vxh(2) - vxh(1)
  vyr = vyh(2) - vyh(1)
  vzr = vzh(2) - vzh(1)
  dv2 = vxr**2 + vyr**2 + vzr**2

  vdotr = xr*vxr + yr*vyr + zr*vzr

  tmin = -vdotr/dv2

  if(vdotr > 0.0_rk) then

    if(dr2 >= r2crit) then

      icflg = 0

    else

      icflg = 1

    end if

  else

    ! We are converging, so we need to calculate the minimum separation attained in time dt
    tmin = -vdotr/dv2

    if(tmin < dt) then

      r2min = dr2 - (vdotr**2)/dv2

    else

      r2min = dr2 + 2.0_rk*vdotr*dt + dv2*dt**2

    end if

    r2min = min(r2min, dr2) ! Really make sure

    if(r2min <= r2crit) then

      icflg = 1

    else

      icflg = 0

    end if

  end if

  return
  end subroutine symba5_chk2
!!
  subroutine symba5_chk(rhill, nbod, ip1, ip2, mass, xh, yh, zh, vxh, vyh, vzh, dt, irec, ieflag, lvdotr)
  !*************************************************************************
  !                            SYMBA5_CHK.F
  !*************************************************************************
  ! This subroutine checks to see if there are encounters
  !
  !             Input:
  !                 rhill         ==>  Radius of hill sphere (real array)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 ip1, ip2       ==>  The two bodies to check (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 xh, yh, zh      ==>  initial position in helio coord
  !                                    (real arrays)
  !                 vxh, vyh, vzh   ==>  initial velocity in helio coord
  !                                    (real arrays)
  !                 dt            ==>  time step  (real scalor)
  !                 irec          ==>  current recursion level (int scalar)
  !             Output:
  !                 ieflag        ==> ecounter?  = 1 Yes
  !                                              =  0 No (integer(ik) :: scalar)
  !                 lvdotr        ==> = .true. if i & j are approaching
  !                                   = .false if i & j are receding
  !                                     (logical*1 scalar)
  !
  ! Remarks: Based on plh_chk.f.  Same as symba_chk.f
  ! Authors:  Hal Levison
  ! Date:   3/20/97
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: nbod, irec, ip1, ip2
  real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod), dt
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod), rhill(nbod)

  ! Outputs
  integer(ik) :: ieflag
  logical(lk) :: lvdotr

  ! Internals
  real(rk) :: r2crit, r2critp, rcrit
  real(rk) :: xr, yr, zr, vxr, vyr, vzr
  real(rk) :: vdotr

  !-----------------
  ! Executable code
  !-----------------

  rcrit = (rhill(ip1) + rhill(ip2))*RHSCALE*RSHELL**irec
  r2crit = rcrit*rcrit
  r2critp = -1.0_rk          ! not used here

  xr = xh(ip2) - xh(ip1)
  yr = yh(ip2) - yh(ip1)
  zr = zh(ip2) - zh(ip1)
  vxr = vxh(ip2) - vxh(ip1)
  vyr = vyh(ip2) - vyh(ip1)
  vzr = vzh(ip2) - vzh(ip1)
  call rmvs_chk_ind(xr, yr, zr, vxr, vyr, vzr, dt, r2crit, r2critp, ieflag)

  vdotr = xr*vxr + yr*vyr + zr*vzr
  lvdotr = vdotr < 0.0_rk

  !rr2 = xr**2 + yr**2 + zr**2
  !vr2 = vxr**2 + vyr**2 + vzr**2
  !massc = mass(ip1)
  !energy = 0.5_rk*vr2 - massc/sqrt(rr2)
  !lbound = (energy < 0.0_rk) .and. (rr2 < rhill(ip1)

  return
  end subroutine symba5_chk
!!
  subroutine symba5_enc_drift(nbod, ielc, ielst, ielev, irec, mass, xh, yh, zh, vxb, vyb, vzb, dt)
  !*************************************************************************
  !                        SYMBA5_ENC_DRIFT.F90
  !*************************************************************************
  ! This subroutine loops thorugh the particles and calls the danby routine
  !
  !             Input:
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 ielev         ==>  Level of particles (int array)
  !                 irec          ==>  current level of the code
  !                 mass          ==>  mass of bodies (real array)
  !                 xh, yh, zh      ==>  initial position in helio coord
  !                                    (real arrays)
  !                 vxb, vyb, vzb   ==>  initial position in bary coord
  !                                    (real arrays)
  !                 dt            ==>  time step
  !             Output:
  !                 xh, yh, zh      ==>  final position in helio coord
  !                                       (real arrays)
  !                 vxb, vyb, vzb   ==>  final position in bary coord
  !                                       (real arrays)
  !
  ! Remarks:  Based on helio_drift.f
  ! Authors:  Hal Levison
  ! Date:    1/20.97
  ! Last revision:
  implicit none

  ! Inputs Only:
  integer(ik) :: nbod, irec, ielc
  integer(ik) :: ielst(NENMAX,2), ielev(nbod)
  real(rk) :: mass(nbod), dt

  ! Inputs and Outputs:
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

  ! Internals:
  logical(lk) :: ldflag(nbod)
  integer(ik) :: i, j, ie, iflag
  real(rk) :: mu

  !----
  ! Executable code
  !-----------------

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, ldflag)
  do i = 1, nbod

    ldflag(i) = .true.

  end do
  !$OMP END PARALLEL DO

  ! Make a copy of the central body's mass
  mu = mass(1)

  ! Take a drift forward dth
  do ie = 1, ielc

    i = ielst(ie,1)
    j = ielst(ie,2)

    ! Perform a drift for particle i
    if((ielev(i) == irec) .and. (mass(i) /= 0.0_rk) .and. ldflag(i)) then

      call drift_one(mu, xh(i), yh(i), zh(i), vxb(i), vyb(i), vzb(i), dt, iflag)
      ldflag(i) = .false. ! Set so we don't return here

      if(iflag /= 0) then

        write(*,'(a,i7,a)')         'Particle:    ', i, ' is lost!!!!!!!!!'
        write(*,'(a,2(1pe14.6))')   'Mass, dt:    ', mu, dt
        write(*,'(a,3(1pe14.6))')   'Helio. pos.: ', xh(i), yh(i), zh(i)
        write(*,'(a,3(1pe14.6),/)') 'Bary. vel.:  ', vxb(i), vyb(i), vzb(i)
        call util_exit(FAILURE)

      end if

    end if

    ! Perform a drift for particle j
    if((ielev(j) == irec) .and. (mass(j) /= 0.0_rk) .and. ldflag(j)) then

      call drift_one(mu, xh(j), yh(j), zh(j), vxb(j), vyb(j), vzb(j), dt, iflag)
      ldflag(j) = .false. ! Set so we don't return here

      if(iflag /= 0) then

        write(*,'(a,i7,a)')         'Particle:    ', j, ' is lost!!!!!!!!!'
        write(*,'(a,2(1pe14.6))')   'Mass, dt:    ', mu, dt
        write(*,'(a,3(1pe14.6))')   'Helio. pos.: ', xh(j), yh(j), zh(j)
        write(*,'(a,3(1pe14.6),/)') 'Bary. vel.:  ', vxb(j), vyb(j), vzb(j)
        call util_exit(FAILURE)

      end if

    end if

  end do

  return
  end subroutine symba5_enc_drift
!!
  subroutine symba5_getacch(nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh, mtiny, ielc, ielst)
  !*************************************************************************
  !                        SYMBA5_GETACCH.F
  !*************************************************************************
  ! This subroutine calculates the acceleration on the massive particles
  ! in the HELIOCENTRIC frame.
  !             Input:
  !                 nbod        ==>  number of massive bodies (int scalor)
  !                 nbodm       ==>  Location of last massive body(int scalar)
  !                 mass        ==>  mass of bodies (real array)
  !                 j2rp2, j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
  !                                     (real scalars)
  !                 xh, yh, zh    ==>  position in heliocentric coord (real arrays)
  !                 mtiny       ==>  Small mass  (real array)
  !                ielc           ==>  number of encounters (integer*2 scalar)
  !                ielst          ==>  list of ecnounters (2D integer*2 array)
  !             Output:
  !                 axh, ayh, azh ==>  acceleration in helio coord (real arrays)
  !
  ! Remarks: Based on helio_getacch.f,  but does not include the forces of
  !          an body B on body A,  if body B and A are having an encounter.
  ! Author:  Hal Levison
  ! Date:    3/20/97
  ! Last revision: 11/22/97
  implicit none

  ! Inputs:
  integer(ik) :: nbod, nbodm
  integer(ik) :: ielst(NENMAX,2), ielc
  real(rk) :: mass(nbod), j2rp2, j4rp4, mtiny
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)

  ! Outputs:
  real(rk) :: axh(nbod), ayh(nbod), azh(nbod)

  ! Internals:
  integer(ik) :: i, j, k
  real(rk) :: aoblx(nbod), aobly(nbod), aoblz(nbod), aoblx0, aobly0, aoblz0
  real(rk) :: ir3h(nbod), irh(nbod)
  real(rk) :: dx, dy, dz, dr2, idr32, faci, facj, xhi, yhi, zhi, massi
  real(rk) :: daxj, dayj, dazj
  !real(rk) :: daxerr, dayerr, dazerr

  !-----------------
  ! Executable code
  !-----------------

  ! Zero things
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, axh, ayh, azh)
  do i = 1, nbod

    axh(i) = 0.0_rk
    ayh(i) = 0.0_rk
    azh(i) = 0.0_rk

  end do
  !$OMP END PARALLEL DO

  ! Compute the acceleration between gravitating bodies, the acceleration on non-gravitating
  ! bodies due to the gravitating bodies, and their back reaction on the gravitating bodies.
  ! No mutual gravitational interaction is computed for bodies nbodm + 1 to nbod.
  !
  ! N.B. The parallel approach outlined below is not reccomended for a large number of
  ! gravitating bodies (i.e. nbodm >~ 1000), otherwise the overhead for creating a destroying
  ! the parallel region in the inner loop will start to become a factor.
  do i = 2, nbodm

    ! Initialize the sum of the acceleration from bodies j
    daxj = 0.0_rk !; daxerr = 0.0_rk
    dayj = 0.0_rk !; dayerr = 0.0_rk
    dazj = 0.0_rk !; dazerr = 0.0_rk

    ! Extract a copy of the information for particle i, and make a copy available for each thread
    xhi = xh(i); yhi = yh(i); xhi = xh(i)
    massi = mass(i)

    !FIRSTPRIVATE(i, xhi, yhi, zhi, massi, daxerr, dayerr, dazerr)
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(i, xhi, yhi, zhi, massi) &
    !$OMP PRIVATE(j, dx, dy, dz, dr2, idr32, faci, facj) SHARED(nbod, mass, xh, yh, zh, axh, ayh, azh) &
    !$OMP REDUCTION(+ : daxj, dayj, dazj)
    do j = i + 1, nbod

      dx = xh(j) - xh(i)
      dy = yh(j) - yh(i)
      dz = zh(j) - zh(i)
      dr2 = dx**2 + dy**2 + dz**2

      idr32 = 1.0_rk/(dr2*sqrt(dr2))
      faci = massi*idr32
      facj = mass(j)*idr32

      ! Sum the acceleration for gravitating body i due to body j (gravitating/non-gravitating)
      daxj = daxj + facj*dx
      dayj = dayj + facj*dy
      dazj = dazj + facj*dz
      !daxj = util_kahan_sum(daxj, facj*dx, daxerr)
      !dayj = util_kahan_sum(dayj, facj*dy, dayerr)
      !dazj = util_kahan_sum(dazj, facj*dz, dazerr)

      ! Update the acceleration for body j (gravitating/non-gravitating) due to gravitating body i
      axh(j) = axh(j) - faci*dx
      ayh(j) = ayh(j) - faci*dy
      azh(j) = azh(j) - faci*dz

    end do
    !$OMP END PARALLEL DO

    ! Update the acceleration for gravitating body i due to bodies j (gravitating/non-gravitating)
    axh(i) = axh(i) + daxj
    ayh(i) = ayh(i) + dayj
    azh(i) = azh(i) + dazj

  end do

  ! Now subtract off anyone in an encounter
  do k = 1, ielc

    i = ielst(k,1)
    j = ielst(k,2)

    dx = xh(j) - xh(i)
    dy = yh(j) - yh(i)
    dz = zh(j) - zh(i)
    dr2 = dx**2 + dy**2 + dz**2

    idr32 = 1.0_rk/(dr2*sqrt(dr2))
    faci = mass(i)*idr32
    facj = mass(j)*idr32

    axh(j) = axh(j) + faci*dx
    ayh(j) = ayh(j) + faci*dy
    azh(j) = azh(j) + faci*dz

    axh(i) = axh(i) - facj*dx
    ayh(i) = ayh(i) - facj*dy
    azh(i) = azh(i) - facj*dz

  end do

  ! Now do j2 and j4 stuff
  if(j2rp2 /= 0.0_rk) then

    call getacch_ir3(nbod, 2, xh, yh, zh, ir3h, irh)
    call obl_acc(nbod, mass, j2rp2, j4rp4, xh, yh, zh, irh, aoblx, aobly, aoblz)

    aoblx0 = aoblx(1); aobly0 = aobly(1); aoblz0 = aoblz(1)

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(aoblx0, aobly0, aoblz0) &
    !$OMP SHARED(nbod, mass, axh, ayh, azh, aoblx, aobly, aoblz)
    do i = 2, nbod

      if(mass(i) /= 0.0_rk) then

        axh(i) = axh(i) + aoblx(i) - aoblx0
        ayh(i) = ayh(i) + aobly(i) - aobly0
        azh(i) = azh(i) + aoblz(i) - aoblz0

      end if

    end do
    !$OMP END PARALLEL DO

  end if

  return
  end subroutine symba5_getacch
!!
  subroutine symba5_helio_drift(nbod, ielev, irec, mass, xh, yh, zh, vxb, vyb, vzb, dt)
  !*************************************************************************
  !                        SYMBA5_HELIO_DRIFT.F
  !*************************************************************************
  ! This subroutine loops thorugh the particles and calls the danby routine
  !
  !             Input:
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 ielev         ==>  Level of particles (int array)
  !                 irec          ==>  current level of the code
  !                 mass          ==>  mass of bodies (real array)
  !                 xh, yh, zh      ==>  initial position in helio coord
  !                                    (real arrays)
  !                 vxb, vyb, vzb   ==>  initial position in bary coord
  !                                    (real arrays)
  !                 dt            ==>  time step
  !             Output:
  !                 xh, yh, zh      ==>  final position in helio coord
  !                                       (real arrays)
  !                 vxb, vyb, vzb   ==>  final position in bary coord
  !                                       (real arrays)
  !
  ! Remarks:  Based on helio_drift.f
  ! Authors:  Hal Levison
  ! Date:    1/20.97
  ! Last revision:
  implicit none

  ! Inputs Only:
  integer(ik) :: nbod, irec
  integer(ik) :: ielev(nbod)
  real(rk) :: mass(nbod), dt

  ! Inputs and Outputs:
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

  ! Internals:
  logical(lk) :: lflag
  integer(ik) :: j, iflag(nbod)
  real(rk) :: mu

  !----
  ! Executable code
  !-----------------

  ! Make a copy of the central body's mass, and set its flag (which should not be used)
  mu = mass(1)
  iflag(1) = 0
  lflag = .false.

  ! Take a drift forward dth
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(irec, mu, dt) PRIVATE(j) &
  !$OMP SHARED(nbod, mass, xh, yh, zh, vxb, vyb, vzb, iflag, ielev, lflag)
  do j = 2, nbod

    iflag(j) = 0
    if((ielev(j) == irec) .and. (mass(j) /= 0.0_rk)) call drift_one(mu, xh(j), yh(j), zh(j), vxb(j), vyb(j), vzb(j), dt, iflag(j))
    if(iflag(j) /= 0) lflag = .true.

  end do
  !$OMP END PARALLEL DO

  ! If any drifts were unsuccessful, print all of them, then abort
  if(lflag) then

    do j = 2, nbod

      if(iflag(j) /= 0) then

        write(*,'(a,i7,a)')         'Particle:    ', j, ' is lost!!!!!!!!!'
        write(*,'(a,2(1pe14.6))')   'Mass, dt:    ', mu, dt
        write(*,'(a,3(1pe14.6))')   'Helio. pos.: ', xh(j), yh(j), zh(j)
        write(*,'(a,3(1pe14.6),/)') 'Bary. vel.:  ', vxb(j), vyb(j), vzb(j)

      end if

    end do

    ! Abandon ship!
    call util_exit(FAILURE)

  end if

  return
  end subroutine symba5_helio_drift
!!
  subroutine symba5_helio_getacch(iflg, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh)
  !*************************************************************************
  !                        SYMBA5_HELIO_GETACCH.F
  !*************************************************************************
  ! This subroutine calculates the acceleration on the massive particles
  ! in the HELIOCENTRIC frame.
  !             Input:
  !                 iflg        ==>  =0 calculate forces (int scalor)
  !                                  =1 don't
  !                 nbod        ==>  number of massive bodies (int scalor)
  !                 nbodm       ==>  The last massive particle
  !                                  (int scalor)
  !                 mass        ==>  mass of bodies (real array)
  !                 j2rp2, j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
  !                                     (real scalars)
  !                 xh, yh, zh    ==>  position in heliocentric coord (real arrays)
  !             Output:
  !                 axh, ayh, azh ==>  acceleration in helio coord (real arrays)
  !
  ! Remarks Based on helio_getacch.f
  ! Author:  Hal Levison
  ! Date:    9/12/99
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: nbod, nbodm, iflg
  real(rk) :: mass(nbod), j2rp2, j4rp4
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)

  ! Outputs:
  real(rk) :: axh(nbod), ayh(nbod), azh(nbod)

  ! Internals:
  integer(ik) :: i, j
  real(rk) :: aoblx(nbod), aobly(nbod), aoblz(nbod), aoblx0, aobly0, aoblz0
  real(rk), save :: axhl(NTPMAX), ayhl(NTPMAX), azhl(NTPMAX)
  real(rk) :: ir3h(nbod), irh(nbod)
  real(rk) :: dx, dy, dz, dr2, idr32, faci, facj, xhi, yhi, zhi, massi
  real(rk) :: daxj, dayj, dazj
  !real(rk) :: daxerr, dayerr, dazerr

  !----
  ! Executable code
  !-----------------

  if(iflg == 0) then

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, axhl, ayhl, azhl)
    do i = 1, nbod

      axhl(i) = 0.0_rk
      ayhl(i) = 0.0_rk
      azhl(i) = 0.0_rk

    end do
    !$OMP END PARALLEL DO

    ! Compute the acceleration between gravitating bodies, the acceleration on non-gravitating
    ! bodies due to the gravitating bodies, and their back reaction on the gravitating bodies.
    ! No mutual gravitational interaction is computed for bodies nbodm + 1 to nbod.
    !
    ! N.B. The parallel approach outlined below is not reccomended for a large number of
    ! gravitating bodies (i.e. nbodm >~ 1000), otherwise the overhead for creating a destroying
    ! the parallel region in the inner loop will start to become a factor.
    do i = 2, nbodm

      ! Initialize the sum of the acceleration from bodies j
      daxj = 0.0_rk !; daxerr = 0.0_rk
      dayj = 0.0_rk !; dayerr = 0.0_rk
      dazj = 0.0_rk !; dazerr = 0.0_rk

      ! Extract a copy of the information for particle i, and make a copy available for each thread
      xhi = xh(i); yhi = yh(i); zhi = zh(i)
      massi = mass(i)

      !FIRSTPRIVATE(i, xhi, yhi, zhi, massi, daxerr, dayerr, dazerr)
      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(i, xhi, yhi, zhi, massi) &
      !$OMP PRIVATE(j, dx, dy, dz, dr2, idr32, faci, facj) SHARED(nbod, mass, xh, yh, zh, axhl, ayhl, azhl) &
      !$OMP REDUCTION(+ : daxj, dayj, dazj)
      do j = i + 1, nbod

        dx = xh(j) - xhi
        dy = yh(j) - yhi
        dz = zh(j) - zhi
        dr2 = dx**2 + dy**2 + dz**2

        idr32 = 1.0_rk/(dr2*sqrt(dr2))
        faci = massi*idr32
        facj = mass(j)*idr32

        ! Sum the acceleration for gravitating body i due to body j (gravitating/non-gravitating)
        daxj = daxj + facj*dx
        dayj = dayj + facj*dy
        dazj = dazj + facj*dz
        !daxj = util_kahan_sum(daxj, facj*dx, daxerr)
        !dayj = util_kahan_sum(dayj, facj*dy, dayerr)
        !dazj = util_kahan_sum(dazj, facj*dz, dazerr)

        ! Update the acceleration for body j (gravitating/non-gravitating) due to gravitating body i
        axhl(j) = axhl(j) - faci*dx
        ayhl(j) = ayhl(j) - faci*dy
        azhl(j) = azhl(j) - faci*dz

      end do
      !$OMP END PARALLEL DO

      ! Update the acceleration for gravitating body i due to bodies j (gravitating/non-gravitating)
      axhl(i) = axhl(i) + daxj
      ayhl(i) = ayhl(i) + dayj
      azhl(i) = azhl(i) + dazj

    end do

  end if

  ! Now do j2 and j4 stuff
  if(j2rp2 /= 0.0_rk) then

    call getacch_ir3(nbod, 2, xh, yh, zh, ir3h, irh)
    call obl_acc(nbod, mass, j2rp2, j4rp4, xh, yh, zh, irh, aoblx, aobly, aoblz)

    aoblx0 = aoblx(1); aobly0 = aobly(1); aoblz0 = aoblz(1)

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(aoblx0, aobly0, aoblz0) &
    !$OMP SHARED(nbod, mass, axh, ayh, azh, axhl, ayhl, azhl, aoblx, aobly, aoblz)
    do i = 1, nbod

      axh(i) = axhl(i) + aoblx(i) - aoblx0
      ayh(i) = ayhl(i) + aobly(i) - aobly0
      azh(i) = azhl(i) + aoblz(i) - aoblz0

    end do
    !$OMP END PARALLEL DO

  else

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, axh, ayh, azh, axhl, ayhl, azhl)
    do i = 1, nbod

      axh(i) = axhl(i)
      ayh(i) = ayhl(i)
      azh(i) = azhl(i)

    end do
    !$OMP END PARALLEL DO

  end if

  return
  end subroutine symba5_helio_getacch
!!
  subroutine symba5_kick(nbod, mass, irec, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dt, sgn, ielc, ielst)
  !*************************************************************************
  !                             SYMBA5_KICK.F
  !*************************************************************************
  ! Do a symba5 kick
  !
  !             Input:
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 irec          ==>  recursion level  (integer scalar)
  !                 iecnt         ==>  The number of objects that each planet 
  !                                    is encountering (int*2 array)
  !                 ielev         ==>  The level that this particle should go
  !                                             (int*2 array)
  !                 rhill         ==>  Hill sphere of planet (real Scalar)
  !                 xh, yh, zh      ==>  initial position in helio coord 
  !                                    (real arrays)
  !                 vxb, vyb, vzb   ==>  initial velocity in bari coord 
  !                                    (real arrays)
  !                dt             ==>  timestep  (real scalar)
  !                sgn            ==>  add or subtract the force (real scalar)
  !                ielc           ==>  number of encounters (integer*2 scalar)
  !                ielst          ==>  list of ecnounters (2D integer*2 array)
  !            Output:
  !                 vxb, vyb, vzb   ==>  final velocity in bari coord 
  !                                    (real arrays)
  !
  ! Remarks: Uses Man Hoi's force
  ! Authors:  Hal Levison 
  ! Date:   3/20/97
  ! Last revision: 
  implicit none

  ! Inputs Only:
  integer(ik) :: nbod, irec
  integer(ik) :: ielev(nbod)
  integer(ik) :: ielst(NENMAX,2), ielc
  real(rk) :: mass(nbod), dt, rhill(nbod), sgn
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)

  ! Inputs and Outputs:
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

  ! Internals:
  integer(ik) :: i, j, irm1, irecl, ie
  real(rk) :: dx, dy, dz, fac, ris, r
  real(rk) :: ri, rr, r2, faci, facj, ir3, rim1

  !----
  ! Executable code 
  !-----------------

  irm1 = irec - 1
  if(sgn < 0.0_rk) then

    irecl = irec - 1

  else

    irecl = irec

  end if

  ! calculate the accelerations
  do ie = 1, ielc

    i = ielst(ie,1)
    j = ielst(ie,2)

    if((ielev(i) >= irm1) .and. (ielev(j) >= irm1)) then

      ri = ((rhill(i) + rhill(j))*RHSCALE*RSHELL**irecl)**2
      rim1 = ri*RSHELL**2

      dx = xh(j) - xh(i)
      dy = yh(j) - yh(i)
      dz = zh(j) - zh(i)
      r2 = dx**2 + dy**2 + dz**2
      ir3 = 1.0_rk/(r2*sqrt(r2))

      if(r2 < rim1) then

        fac = 0.0_rk

      else if(r2 < ri) then

        ris = sqrt(ri)
        r = sqrt(r2)
        rr = (ris - r)/(ris*(1.0_rk - RSHELL))
        fac = (1.0_rk - 3.0_rk*rr**2 + 2.0_rk*rr**3)*ir3

      else

        fac = ir3

      end if

      ! apply the kick for body i
      facj = sgn*mass(j)*fac*dt
      vxb(i) = vxb(i) + facj*dx
      vyb(i) = vyb(i) + facj*dy
      vzb(i) = vzb(i) + facj*dz

      ! apply the kick for body j
      faci = sgn*mass(i)*fac*dt
      vxb(j) = vxb(j) - faci*dx
      vyb(j) = vyb(j) - faci*dy
      vzb(j) = vzb(j) - faci*dz

    end if

  end do

  return
  end subroutine symba5_kick
!!
  subroutine symba5_merge(t, dt, nbod, nbodm, ip1, ip2, mass, xh, yh, zh, vxb, vyb, vzb, &
             ireci, lvdotrold, ibound, rpl, mergelst, mergecnt, rhill, eoff, ielc, ielst)
  !*************************************************************************
  !                            SYMBA5_MERGE.F
  !*************************************************************************
  ! This subroutine checks to see if there are encounters
  !
  !             Input:
  !                 t             ==>  current time (real scalar)
  !                 dt            ==>  time step (real scalar)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 ip1, ip2       ==>  The two bodies to check (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 xh, yh, zh      ==>  initial position in helio coord 
  !                                    (real arrays)
  !                 vxb, vyb, vzb   ==>  initial velocity in helio coord 
  !                                    (real arrays)
  !                 ireci         ==>  current recursion level (int scalar)
  !                 lvdotrold     ==>  old radial velocity test
  !                                   = .true. if i & j are approaching
  !                                   = .false if i & j are receding
  !                                     (logical*1 scalar)
  !                 iecnt         ==>  The number of objects that each planet 
  !                                    is encountering (int*2 array)
  !                 rpl           ==>  physical size of a planet.
  !                                    (real array)
  !             mergelst          ==>  list of mergers (int array)
  !             mergecnt          ==>  count of mergers (int array)
  !             rhill             ==>  Hill sphere of planet (real Scalar)
  !             eoff              ==>  Energy offset (real scalar)
  !                ielc           ==>  number of encounters (integer*2 scalar)
  !                ielst          ==>  list of ecnounters (2D integer*2 array)
  !
  !             Output:  Changed only if a Megrer happens
  !                 mass          ==>  mass of bodies (real array)
  !                 xh, yh, zh      ==>  initial position in helio coord 
  !                                    (real arrays)
  !                 vxb, vyb, vzb   ==>  initial velocity in helio coord 
  !                                    (real arrays)
  !                 iecnt         ==>  The number of objects that each planet 
  !                                    is encountering (int*2 array)
  !                 rpl           ==>  physical size of a planet.
  !                                    (real array)
  !             mergelst          ==>  list of mergers (int array)
  !             mergecnt          ==>  count of mergers (int array)
  !             rhill             ==>  Hill sphere of planet (real Scalar)
  !             eoff              ==>  Energy offset (real scalar)
  !                ielc           ==>  number of encounters (integer*2 scalar)
  !                ielst          ==>  list of ecnounters (2D integer*2 array)
  !
  ! Remarks: 
  ! Authors:  Hal Levison
  ! Date:   1/2/97
  ! Last revision: 1/24/97
  implicit none

  ! Inputs:
  logical(lk) :: lvdotrold
  integer(ik) :: nbod, nbodm, ireci, ip1, ip2
  real(rk) :: t, dt

  ! Inputs and Outputs:
  integer(ik) :: mergelst(NENMAX,2), mergecnt
  integer(ik) :: ielst(NENMAX,2), ielc, ibound(nbod)
  real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod), eoff
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod), rpl(nbod), rhill(nbod)

  ! Internals
  logical(lk) :: lswitch, lcross
  integer(ik) :: ialpha, ip1l, ip2l, itmp
  real(rk) :: xr, yr, zr, vxr, vyr, vzr, vdotr, tcross2, dt2
  real(rk) :: rlim, rlim2, rr2, vr2, gmsum, energy, a, e, q

  !-----------------!
  ! Executable Code !
  !-----------------!

  xr = xh(ip2) - xh(ip1)
  yr = yh(ip2) - yh(ip1)
  zr = zh(ip2) - zh(ip1)
  rr2 = xr**2 + yr**2 + zr**2

  rlim = rpl(ip1) + rpl(ip2)
  rlim2 = rlim**2
  if(rlim == 0.0_rk) return ! <=== NOTE!!!!

  vxr = vxb(ip2) - vxb(ip1)
  vyr = vyb(ip2) - vyb(ip1)
  vzr = vzb(ip2) - vzb(ip1)
  vdotr = xr*vxr + yr*vyr + zr*vzr
  vr2 = vxr**2 + vyr**2 + vzr**2

  tcross2 = rr2/vr2
  dt2 = dt**2

  lswitch = (lvdotrold .and. (vdotr > 0.0_rk))
  lcross = (tcross2 <= dt2)

  gmsum = mass(ip1) + mass(ip2)
  call orbel_xv2aeq(xr, yr, zr, vxr, vyr, vzr, gmsum, ialpha, a, e, q)

  ! Sort indices in decending order by mass
  ip1l = ip1
  ip2l = ip2
  if(mass(ip2l) > mass(ip1l)) then

    itmp = ip1l
    ip1l = ip2l
    ip2l = itmp

  end if

  ! If particle ip2 is bound to ip1 (ialpha < 0), orbit smaller than 0.1*R_hill of ip1 and it has come in for perigee passage,
  ! then increment ibound, otherwise set to zero.
  if((ialpha < 0) .and. (a < rhill(ip1l)/10.0_rk) .and. lswitch) then

    ibound(ip2l) = ibound(ip2l) + 1

  else

    ibound(ip2l) = 0

  end if

  ! If the sum of the radii is larger than the distance between ip1 and ip2, then merge
  if(rlim2 >= rr2) then

    !ip1l = ip1
    !ip2l = ip2
    call discard_mass_merge5(t, nbod, nbodm, ip1l, ip2l, mass, xh, yh, zh, vxb, vyb, vzb, rpl, eoff, ielc, ielst)
    mergecnt = mergecnt + 1
    mergelst(mergecnt,1) = ip1l
    mergelst(mergecnt,2) = ip2l
    rhill(ip2l) = 0.0_rk
    ibound(ip2l) = 0

    call util_hills1(mass(1), mass(ip1l), xh(ip1l), yh(ip1l), zh(ip1l), vxb(ip1l), vyb(ip1l), vzb(ip1l), rhill(ip1l))
    return ! <== NOTE!!!!

  end if

  ! If the sign of vdotr has changed from approaching to receding (i.e. from vdotr < 0.0 to vdotr > 0.0)
  ! and if crossing time is less than the current time step, and periagee is less than rlim, then merge
  if(lswitch .and. lcross .and. (q < rlim)) then

    !ip1l = ip1
    !ip2l = ip2
    call discard_mass_merge5(t, nbod, nbodm, ip1l, ip2l, mass, xh, yh, zh, vxb, vyb, vzb, rpl, eoff, ielc, ielst)
    mergecnt = mergecnt + 1
    mergelst(mergecnt,1) = ip1l
    mergelst(mergecnt,2) = ip2l
    rhill(ip2l) = 0.0_rk
    ibound(ip2l) = 0

    call util_hills1(mass(1), mass(ip1l), xh(ip1l), yh(ip1l), zh(ip1l), vxb(ip1l), vyb(ip1l), vzb(ip1l), rhill(ip1l))
    return ! <== NOTE!!!!

  end if

  ! If particle ip2 has been flagged to be bound to particle ip1 more than four times (four perigee passages),
  ! consider ip2 to be in a satellite of ip1 and argue that we can merge ip2 with ip1 due to the presence of gas
  ! drag, which is not implemented in these recursive routines.  This can speed up the code since we avoid
  ! performing 10^5+ time steps to accurately resolve such orbits, but this may accidentily remove some
  ! particles that are only captured temporarily.
  if(ibound(ip2l) >= NBOUNDMAX) then

    !write(*,'(a,i6,a,i9,a)') " Particle: ", ip2, " has been bound ", ibound(ip2), " times"
    !write(*,'(2(1x,1pe22.15))') sqrt(rr2), a
    
    !ip1l = ip1
    !ip2l = ip2
    call discard_mass_merge5(t, nbod, nbodm, ip1l, ip2l, mass, xh, yh, zh, vxb, vyb, vzb, rpl, eoff, ielc, ielst)
    mergecnt = mergecnt + 1
    mergelst(mergecnt,1) = ip1l
    mergelst(mergecnt,2) = ip2l
    rhill(ip2l) = 0.0_rk
    ibound(ip2l) = 0

    call util_hills1(mass(1), mass(ip1l), xh(ip1l), yh(ip1l), zh(ip1l), vxb(ip1l), vyb(ip1l), vzb(ip1l), rhill(ip1l))
    return ! <== NOTE!!!!

  end if

  return
  end subroutine symba5_merge
!!
  subroutine symba5_nbodm(nbod, mass, mtiny, nbodm)
  !*************************************************************************
  !                            SYMBA5_NBODM.F
  !*************************************************************************
  ! Returns the location of the last massive body in the list
  !
  !             Input:
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 mtiny         ==>  Small mass  (real array)
  !             Output:
  !                 nbodm         ==>  location of the last massive body
  !                                    (int scalar)
  !
  ! Remarks:  If all the objects are massive,  then nbodm=nbod-1 so that
  !           the do loops will have the correct limits.
  ! Authors:  Hal Levison
  ! Date:    3/20/97
  ! Last revision: 1/29/06
  implicit none

  ! Inputs Only:
  integer(ik) :: nbod
  real(rk) :: mass(nbod),  mtiny

  ! Outputs only
  integer(ik) :: nbodm

  ! Internals
  integer(ik) :: i
  real(rk) :: mtiny_private

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! This assumes that the particle identifiers are already sorted in order of decreasing mass
  if(mass(nbod) > mtiny) then

    nbodm = nbod - 1

    write(*,'(a,i6,a)') " Out of ", nbod, " bodies, all are massive."

  else

    ! Initialize the value of nbodm
    nbodm = 1

    ! Make a private copy of mtiny available for each thread
    mtiny_private = mtiny

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) &
    !$OMP FIRSTPRIVATE(mtiny_private) SHARED(nbod, mass) REDUCTION(MAX : nbodm)
    do i = 2, nbod - 1

      if(mass(i) > mtiny_private) nbodm = i

    end do
    !$OMP END PARALLEL DO

    write(*,'(2(a,i6),a)') " Out of ", nbod, " bodies, ", nbodm, " are massive."

  end if

  return
  end subroutine symba5_nbodm
!!
  subroutine symba5_step_helio(i1st, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, dt)
  !*************************************************************************
  !                            SYMBA5_STEP_HELIO.F
  !*************************************************************************
  ! This subroutine takes a step in helio coord.
  ! Does a KICK than a DRIFT than a KICK.
  ! ONLY DOES MASSIVE PARTICLES
  !
  !             Input:
  !                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 j2rp2, j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
  !                                     (real scalars)
  !                 xh, yh, zh      ==>  initial position in helio coord 
  !                                    (real arrays)
  !                 vxh, vyh, vzh   ==>  initial velocity in helio coord 
  !                                    (real arrays)
  !                 dt            ==>  time step
  !             Output:
  !                 xh, yh, zh      ==>  final position in helio coord 
  !                                       (real arrays)
  !                 vxh, vyh, vzh   ==>  final velocity in helio coord 
  !                                       (real arrays)
  ! Remarks: Based on helio_step_pl.f but does not pass the intermediate
  !          positions and velocities back for the TP to use.
  ! Authors:  Hal Levison 
  ! Date:    3/20/97
  ! Last revision: 12/13/00
  implicit none

  ! Inputs Only: 
  integer(ik) :: nbod, i1st, nbodm
  real(rk) :: mass(nbod), dt, j2rp2, j4rp4

  ! Inputs and Outputs:
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

  ! Internals:
  integer(ik) :: i1stloc
  real(rk) :: dth, msys
  real(rk) :: axh(nbod), ayh(nbod), azh(nbod)
  real(rk), save :: vxb(NTPMAX), vyb(NTPMAX), vzb(NTPMAX)
  real(rk) :: ptxb, ptyb, ptzb            ! Not used here
  real(rk) :: ptxe, ptye, ptze

  !----
  ! Executable code
  !-----------------

  dth = 0.5_rk*dt

  i1stloc = i1st
  if(i1st == 0) then

    ! Convert vel to bery to jacobi coords
    call coord_vh2b(nbod, mass, vxh, vyh, vzh, vxb, vyb, vzb, msys)
    i1st = 1              ! turn this off

  end if

  ! Do the linear drift due to momentum of the Sun
  call helio_lindrift(nbod, mass, vxb, vyb, vzb, dth, xh, yh, zh, ptxb, ptyb, ptzb)

  ! Get the accelerations in helio frame, if first time step
  call symba5_helio_getacch(i1stloc, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh)

  ! Apply a heliocentric kick for a half dt 
  call kickvh(nbod, vxb, vyb, vzb, axh, ayh, azh, dth)

  ! Drift in helio coords for the full step 
  call helio_drift(nbod, mass, xh, yh, zh, vxb, vyb, vzb, dt)

  ! Get the accelerations in helio frame, regardless if it is the first step
  call symba5_helio_getacch(0, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh)

  ! Apply a heliocentric kick for a half dt 
  call kickvh(nbod, vxb, vyb, vzb, axh, ayh, azh, dth)

  ! Do the linear drift due to momentum of the Sun
  call helio_lindrift(nbod, mass, vxb, vyb, vzb, dth, xh, yh, zh, ptxe, ptye, ptze)

  ! convert back to helio velocities
  call coord_vb2h(nbod, mass, vxb, vyb, vzb, vxh, vyh, vzh)

  return
  end subroutine symba5_step_helio
!!
  subroutine symba5_step_interp(time, ielev, nbod, nbodm, mass, rhill, j2rp2, j4rp4, rpl, &
             xh, yh, zh, vxh, vyh, vzh, dt, mergelst, mergecnt, eoff, ielc, ielst, mtiny, ibound)
  !*************************************************************************
  !                            SYMBA5_STEP_INTERP.F
  !*************************************************************************
  !
  !             Input:
  !                 time          ==> Current time (real scalar)
  !                 ielev         ==>  The level that this particle should go
  !                                             (int*2 array)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 nbodm         ==>  Location of last massive body(int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 rhill         ==>  Radius of hill sphere (real array)
  !                 j2rp2, j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
  !                                     (real scalars)
  !                 xh, yh, zh      ==>  initial position in helio coord 
  !                                    (real arrays)
  !                 vxh, vyh, vzh   ==>  initial velocity in helio coord 
  !                                    (real arrays)
  !                 dt            ==>  time step
  !                 rpl           ==>  physical size of a planet.
  !                                    (real array)
  !                 eoff          ==>  Energy offset (real scalar)
  !                ielc           ==>  number of encounters (integer*2 scalar)
  !                ielst          ==>  list of ecnounters (2D integer*2 array)
  !                mtiny          ==>  Small mass  (real array)
  !             Output:
  !                 xh, yh, zh      ==>  final position in helio coord 
  !                                       (real arrays)
  !                 vxh, vyh, vzh   ==>  final velocity in helio coord 
  !                                       (real arrays)
  !                 rpl           ==>  Recalculated physical size of a planet.
  !                                    if merger happened (real array)
  !                 nbod          ==>  Recalculated number of massive bodies 
  !                                    if merger happened (int scalar)
  !                 nbodm         ==>  Location of last massive body(int scalar)
  !                 mass          ==>  Recalculated mass of bodies 
  !                                    if merger happened (real array)
  !                 mergelst      ==>  list of mergers (int array)
  !                 mergecnt      ==>  count of mergers (int array)
  !                 eoff          ==>  Energy offset (real scalar)
  ! Remarks: 
  ! Authors:  Hal Levison 
  ! Date:    11/21/96
  ! Last revision: 5/13/99
  implicit none

  ! Inputs Only:
  integer(ik) :: ielst(NENMAX,2), ielc
  integer(ik) :: ielev(nbod)
  real(rk) :: mass(nbod), dt, j2rp2, j4rp4, time, mtiny

  ! Inputs and Outputs:
  integer(ik) :: nbod, nbodm, ibound(nbod)
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)
  real(rk) :: rpl(nbod), eoff
  real(rk) :: rhill(nbod)

  ! Outputs
  integer(ik) :: mergelst(NENMAX,2), mergecnt

  ! Internals:
  logical(lk) :: lvdotr(NTPMAX) ! Used by symba_step_recur
  integer(ik) :: i, irec
  real(rk) :: dth, msys
  real(rk), save :: axh(NTPMAX), ayh(NTPMAX), azh(NTPMAX)
  real(rk), save :: vxb(NTPMAX), vyb(NTPMAX), vzb(NTPMAX)
  real(rk) :: ptxb, ptyb, ptzb            ! Not used here
  real(rk) :: ptxe, ptye, ptze

  !----
  ! Executable code
  !-----------------

  ! Half the time step
  dth = 0.5_rk*dt

  ! Initialize recursion index and merge counter
  irec = 0
  mergecnt = 0

  ! Convert vel to bary to jacobi coords
  call coord_vh2b(nbod, mass, vxh, vyh, vzh, vxb, vyb, vzb, msys)
  !write(100,'(a,3(1x,1pe22.15))') 'vxh(1:3) = ', vxh(1:3)
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  ! Do the linear drift due to momentum of the Sun
  call helio_lindrift(nbod, mass, vxb, vyb, vzb, dth, xh, yh, zh, ptxb, ptyb, ptzb)
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  ! Get the accelerations in helio frame. For each object
  ! only include those guys that it is not encountering with. 
  call symba5_getacch(nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh, mtiny, ielc, ielst)
  !write(100,'(a,3(1x,1pe22.15))') 'axh(1:3) = ', axh(1:3)

  ! Apply a heliocentric kick for a half dt 
  call kickvh(nbod, vxb, vyb, vzb, axh, ayh, azh, dth)
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  ! Do a recursion step for full dt
  call symba5_helio_drift(nbod, ielev, -1, mass, xh, yh, zh, vxb, vyb, vzb, dt)
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  call symba5_step_recur(time, nbod, nbodm, mass, irec, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, &
       rpl, mergelst, mergecnt, dt, eoff, lvdotr, ibound, ielc, ielst) 
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  ! Get the accelerations in helio frame. For each object
  ! only include those guys that it is not encountering with. 
  call symba5_getacch(nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh, mtiny, ielc, ielst)
  !write(100,'(a,3(1x,1pe22.15))') 'axh(1:3) = ', axh(1:3)

  ! Apply a heliocentric kick for a half dt 
  call kickvh(nbod, vxb, vyb, vzb, axh, ayh, azh, dth)
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  ! Do the linear drift due to momentum of the Sun
  call helio_lindrift(nbod, mass, vxb, vyb, vzb, dth, xh, yh, zh, ptxe, ptye, ptze)
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  ! convert back to helio velocities
  call coord_vb2h(nbod, mass, vxb, vyb, vzb, vxh, vyh, vzh)
  !write(100,'(a,3(1x,1pe22.15))') 'vxh(1:3) = ', vxh(1:3)
  !write(100,*) '--------'

  return
  end subroutine symba5_step_interp
!!
  subroutine symba5_step_pl(i1st, time, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, dt, &
             lclose, rpl, isenc, mergelst, mergecnt, iecnt, eoff, rhill, mtiny, ibound)
  !----------------------------------------------------------------------------
  !				SYMBA5_STEP_PL.F90
  !----------------------------------------------------------------------------
  !
  !             Input:
  !                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
  !                 time          ==>  current time (real scalar)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 nbodm         ==>  location of the last massie body
  !                                    (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
  !                                     (real scalars)
  !                 xh,yh,zh      ==>  initial position in helio coord
  !                                    (real arrays)
  !                 vxh,vyh,vzh   ==>  initial velocity in helio coord
  !                                    (real arrays)
  !                 dt            ==>  time step
  !                 lclose        ==> .true. --> check for close encounters
  !                                      (logical*2 scalar)
  !                 rpl           ==>  physical size of a planet.
  !                                    (real array)
  !                 eoff          ==>  Energy offset (real scalar)
  !                 rhill         ==>  size of planet's hills sphere
  !                                    (real array)
  !                 mtiny         ==>  Small mass  (real array)
  !             Output:
  !                 xh,yh,zh      ==>  final position in helio coord
  !                                       (real arrays)
  !                 vxh,vyh,vzh   ==>  final velocity in helio coord
  !                                       (real arrays)
  !                 rpl           ==>  Recalculated physical size of a planet.
  !                                    if merger happened (real array)
  !                 nbod          ==>  Recalculated number of massive bodies
  !                                    if merger happened (int scalar)
  !                 mass          ==>  Recalculated mass of bodies
  !                                    if merger happened (real array)
  !                 isenc         ==>  0 --> No encounter during last dt
  !                                    1 --> There was encounters
  !                                     (integer scalar)
  !                 mergelst      ==>  list of mergers (int array)
  !                 mergecnt      ==>  count of mergers (int array)
  !                 iecnt         ==>  Number of encounters (int*2 array)
  !                 eoff          ==>  Energy offset (real scalar)
  !                 rhill         ==>  size of planet's hills sphere
  !                                    (real array)
  !
  ! Remarks: Based on symba2_step_pl.f
  ! Authors:  Hal Levison
  ! Date:    11/27/97
  ! Last revision:
  !$ use omp_lib
  implicit none

  ! Inputs Only: 
  logical(lk) :: lclose
  integer(ik) :: nbod, i1st, nbodm
  real(rk) :: mass(nbod), dt, time, j2rp2, j4rp4, mtiny

  ! Inputs and Outputs:
  integer(ik) :: ibound(nbod)
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)
  real(rk) :: rpl(nbod), rhill(nbod), eoff

  ! Outputs only
  integer(ik) :: isenc
  integer(ik) :: iecnt(nbod), ielev(nbod)
  integer(ik) :: mergelst(NENMAX,2), mergecnt

  ! Internals
  logical(lk) :: lvdotr            ! Not used in the routine
  integer(ik) :: i, j, k, ieflg, irec
  integer(ik) :: ielst(NENMAX,2), ielc
  integer(ik) :: id, ielc_shared(nthreads), istart(nthreads), ielst_shared(nbod)
  real(rk) :: rhillij(2), xhij(2), yhij(2), zhij(2), vxhij(2), vyhij(2), vzhij(2)

  !----
  ! Executable code
  !-----------------

  ! Initialize relevant variables
  isenc = 0
  ielc = 0
  irec = 0

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, iecnt, ielev)
  do i = 1, nbod

    iecnt(i) = 0
    ielev(i) = -1

  end do
  !$OMP END PARALLEL DO

  ! Check for encounters
  if(lclose) then

    do i = 2, nbodm

      ! Extract the information for body i, and make a copy available for each thread
      rhillij(1) = rhill(i)
      xhij(1) = xh(i); yhij(1) = yh(i); zhij(1) = zh(i)
      vxhij(1) = vxh(i); vyhij(1) = vyh(i); vzhij(1) = vzh(i)

      !$OMP PARALLEL DEFAULT(NONE) PRIVATE(j, id, ieflg) &
      !$OMP FIRSTPRIVATE(i, nthreads, irec, xhij, yhij, zhij, vxhij, vyhij, vzhij, rhillij, dt) &
      !$OMP SHARED(nbod, rhill, xh, yh, zh, vxh, vyh, vzh, istart, ielc_shared, ielst_shared, iecnt, ielev)

      id = 1                              ! Set thread identifier for *serial* case
      !$ id = omp_get_thread_num() + 1    ! Set thread identifier for *parallel* case
      ielc_shared(id) = 0                 ! Initialize encounter counter for each thread
      istart(id) = (id - 1)*nbod/nthreads ! Set the starting position in ielst_shared for each thread

      !$OMP DO SCHEDULE(STATIC)
      do j = i + 1, nbod

        rhillij(2) = rhill(j)
        xhij(2) = xh(j); yhij(2) = yh(j); zhij(2) = zh(j)
        vxhij(2) = vxh(j); vyhij(2) = vyh(j); vzhij(2) = vzh(j)
        call symba5_chk2(rhillij, xhij, yhij, zhij, vxhij, vyhij, vzhij, dt, irec, ieflg)

        if(ieflg /= 0) then

          iecnt(j) = iecnt(j) + 1
          ielev(j) = 0
          ielc_shared(id) = ielc_shared(id) + 1 ! Update counter *before* storing index j since istart starts at 0
          ielst_shared(istart(id) + ielc_shared(id)) = j

        end if

      end do
      !$OMP END DO NOWAIT

      !$OMP END PARALLEL

      ! If there any encounters, append to the global encounter list for body i
      if(sum(ielc_shared) > 0) then

        ! Set global encounter flag
        isenc = 1

        ! Set the recursion level to 0 for body i
        ielev(i) = 0

        ! Search through each thread's local encounter list
        do k = 1, nthreads

          ! If an encounter was found for thread k, append ielc_shared(k) entries to ielst
          if(ielc_shared(k) > 0) then

            do j = 1, ielc_shared(k)

              iecnt(i) = iecnt(i) + 1
              ielst(ielc + j, 1) = i
              ielst(ielc + j, 2) = ielst_shared(istart(k) + j)

            end do

            ! Update the global encounter counter
            ielc = ielc + ielc_shared(k)

            if(ielc > NENMAX) then

              write(*,'(a)') 'ERROR: Encounter matrix is filled.'
              write(*,'(a)') 'STOPPING'
              call util_exit(FAILURE)

            end if

          end if

        end do

      end if

    end do

  end if

  ! do a step
  if(isenc == 0) then

    call symba5_step_helio(i1st, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, dt)
    mergecnt = 0

  else

    call symba5_step_interp(time, ielev, nbod, nbodm, mass, rhill, j2rp2, j4rp4, rpl, xh, yh, zh, vxh, vyh, vzh, &
         dt, mergelst, mergecnt, eoff, ielc, ielst, mtiny, ibound)
    i1st = 0

  end if

  ! Print number of encounters and mergers found in this time step
  !write(100,*) time, ielc, mergecnt
  !do i = 1, ielc
  !  write(100,'(2i9)') ielst(i,:)
  !end do

  return
  end subroutine symba5_step_pl
!!
  recursive subroutine symba5_step_recur(t, nbod, nbodm, mass, ireci, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, &
                       rpl, mergelst, mergecnt, dt0, eoff, lvdotr, ibound, ielc, ielst)
  !------------------------------------------------------------------------------
  !				SYMBA5_STEP_RECUR.F90
  !------------------------------------------------------------------------------
  !
  !             Input:
  !                 t             ==>  time (real Scalar)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 nbodm         ==>  Location of last massive body(int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 ireci         ==>  Input recursion level  (integer scalar)
  !                 ilevl         ==>  largest recursion level used
  !                                    (integer array)
  !                 ielev         ==>  The level that this particle should go
  !                                             (int*2 array)
  !                 j2rp2, j4rp4   ==>  J2*radii_pl^2 and J4*radii_pl^4
  !                                     (real scalars)
  !                 rhill         ==>  Hill sphere of planet (real Scalar)
  !                 xh, yh, zh      ==>  initial position in helio coord
  !                                    (real arrays)
  !                 vxb, vyb, vzb   ==>  initial velocity in bari coord
  !                                    (real arrays)
  !                dt0            ==>  Global timestep  (real scalar)
  !                rpl            ==>  physical size of a planet.
  !                                    (real array)
  !                eoff           ==>  Energy offset (real scalar)
  !                ielc           ==>  number of encounters (integer*2 scalar)
  !                ielst          ==>  list of ecnounters (2D integer*2 array)
  !             Output:
  !                 xh, yh, zh      ==>  final position in helio coord
  !                                       (real arrays)
  !                 vxb, vyb, vzb   ==>  final velocity in bari coord
  !                                       (real arrays)
  !             mergelst          ==>  list of mergers (int array)
  !             mergecnt          ==>  count of mergers (int array)
  !                 rpl           ==>  Recalculated physical size of a planet.
  !                                    if merger happened (real array)
  !                 mass          ==>  Recalculated mass of bodies
  !                                    if merger happened (real array)
  !                eoff           ==>  Energy offset (real scalar)
  !                lvdotr         ==> vdotr relative flag
  !                                   = .true. if i, j are receding
  !                                   = .false is approaching
  !                                     (2D logical*1 array)
  !                                   Used internally,  but only need 1 copy.
  !
  ! Remarks: If a merger occurs,  does not change nbod and puts the mass
  !          of one of the particles to zero.
  ! Authors:  Hal Levison
  ! Date:    3/20/97
  ! Last revision: 5/13/99
  implicit none

  ! Inputs Only:
  integer(ik) :: nbod, ireci, nbodm
  integer(ik) :: ielev(nbod)
  integer(ik) :: ielst(NENMAX,2), ielc
  real(rk) :: mass(nbod), dt0, rhill(nbod), t

  ! Inputs and Outputs:
  logical(lk) :: lvdotr(NTPMAX)
  integer(ik) :: ibound(nbod)
  integer(ik) :: mergelst(NENMAX,2), mergecnt
  real(rk) :: xh(nbod), yh(nbod), zh(nbod), eoff
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod), rpl(nbod)

  ! Internals:
  integer(ik) :: i, j, ie
  integer(ik) :: icflg, it, irecp, ieflg
  real(rk) :: dtl, dth, sgn

  !----
  ! Executable code

  dtl = dt0/real(NTENC**ireci, rk)
  dth = 0.5_rk*dtl

  if(dtl/dt0 <= TINY_NUMBER) then

    write(*,'(a)') 'Warning in SYMBA_STEP_RECUR:'
    write(*,'(a)') ' Local timestep too small '
    write(*,'(a)') ' Roundoff will be important!!!! '
    call util_exit(FAILURE)

  end if

  irecp = ireci + 1

  if(ireci == 0) then

    ! Do we need to go deeper?
    icflg = 0
    do ie = 1, ielc

      i = ielst(ie,1)
      j = ielst(ie,2)

      if((ielev(i) >= ireci) .and. (ielev(j) >= ireci)) then

        call symba5_chk(rhill, nbod, i, j, mass, xh, yh, zh, vxb, vyb, vzb, dtl, irecp, ieflg, lvdotr(ie))

        if(ieflg /= 0) then

          icflg = 1
          ielev(i) = irecp
          ielev(j) = irecp

        end if

      end if

    end do

    sgn = 1.0_rk
    call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

    call symba5_enc_drift(nbod, ielc, ielst, ielev, ireci, mass, xh, yh, zh, vxb, vyb, vzb, dtl)

    if(icflg /= 0) call symba5_step_recur(t, nbod, nbodm, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, &
                        rpl, mergelst, mergecnt, dt0, eoff, lvdotr, ibound, ielc, ielst)

    sgn = 1.0_rk
    call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

    ! look for mergers
    do ie = 1, ielc

      i = ielst(ie,1)
      j = ielst(ie,2)

      if((ielev(i) >= ireci) .and. (ielev(j) >= ireci)) then

        call symba5_merge(t, dtl, nbod, nbodm, i, j, mass, xh, yh, zh, vxb, vyb, vzb, &
             ireci, lvdotr(ie), ibound, rpl, mergelst, mergecnt, rhill, eoff, ielc, ielst)

      end if

      if(ielev(i) == irecp) ielev(i) = ireci
      if(ielev(j) == irecp) ielev(j) = ireci

    end do

  else

    do it = 1, NTENC

      ! Do we need to go deeper?
      icflg = 0
      do ie = 1, ielc

        i = ielst(ie,1)
        j = ielst(ie,2)

        if((ielev(i) >= ireci) .and. (ielev(j) >= ireci)) then

          call symba5_chk(rhill, nbod, i, j, mass, xh, yh, zh, vxb, vyb, vzb, dtl, irecp, ieflg, lvdotr(ie))

          if(ieflg /= 0) then

            icflg = 1
            ielev(i) = irecp
            ielev(j) = irecp

          end if

        end if

      end do

      sgn = 1.0_rk
      call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

      sgn = -1.0_rk
      call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

      call symba5_enc_drift(nbod, ielc, ielst, ielev, ireci, mass, xh, yh, zh, vxb, vyb, vzb, dtl)

      if(icflg /= 0) call symba5_step_recur(t, nbod, nbodm, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, &
                          rpl, mergelst, mergecnt, dt0, eoff, lvdotr, ibound, ielc, ielst)

      sgn = 1.0_rk
      call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

      sgn = -1.0_rk
      call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

      ! look for mergers
      do ie = 1, ielc

        i = ielst(ie,1)
        j = ielst(ie,2)

        if((ielev(i) >= ireci) .and. (ielev(j) >= ireci)) then

          call symba5_merge(t, dtl, nbod, nbodm, i, j, mass, xh, yh, zh, vxb, vyb, vzb, &
               ireci, lvdotr(ie), ibound, rpl, mergelst, mergecnt, rhill, eoff, ielc, ielst)

        end if

        if(ielev(i) == irecp) ielev(i) = ireci
        if(ielev(j) == irecp) ielev(j) = ireci

      end do

    end do

  end if

  return
  end subroutine symba5_step_recur
!!
!!
end module symba5