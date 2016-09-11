module helio
! Module for routines in the helio directory
use swift
use util
use mvs
implicit none

contains
!!
!!
  subroutine helio_drift(nbod, mass, xh, yh, zh, vxb, vyb, vzb, dt)
  !-------------------------------------------------------------------------
  !                        HELIO_DRIFT.F
  !-------------------------------------------------------------------------
  ! This subroutine loops thorugh the particles and calls the danby routine
  !
  !             Input:
  !                 nbod          ==>  number of massive bodies (int scalar)
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
  ! Remarks:  Based on drift.f
  ! Authors:  Hal Levison
  ! Date:    11/14/96
  ! Last revision: 1/8/97  for symba
  implicit none

  ! Inputs Only:
  integer(ik) :: nbod
  real(rk) :: mass(nbod), dt

  ! Inputs and Outputs:
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

  ! Internals:
  logical(lk) :: lflag
  integer(ik) :: i, iflag(nbod)
  real(rk) :: mu

  !----
  ! Executable code

  ! Make a copy of the central body's mass, and set its flag (which should not be used)
  mu = mass(1)
  iflag(1) = 0
  lflag = .false.

  ! Take a drift forward dth
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(mu, dt) &
  !$OMP SHARED(nbod, mass, xh, yh, zh, vxb, vyb, vzb, iflag, lflag)
  do i = 2, nbod

    iflag(i) = 0
    if(mass(i) /= 0.0_rk) call drift_one(mu, xh(i), yh(i), zh(i), vxb(i), vyb(i), vzb(i), dt, iflag(i))
    if(iflag(i) /= 0) lflag = .true.

  end do
  !$OMP END PARALLEL DO

  ! If any drifts were flagged as unsuccessful, print them all
  if(lflag) then

    do i = 2, nbod

      if(iflag(i) /= 0) then

        write(*,'(a,i7,a)')         'Particle:    ', i, ' is lost!!!!!!!!!'
        write(*,'(a,2(1pe14.6))')   'Mass, dt:    ', mass(i), dt
        write(*,'(a,3(1pe14.6))')   'Helio. pos.: ', xh(i), yh(i), zh(i)
        write(*,'(a,3(1pe14.6),/)') 'Bary. vel.:  ', vxb(i), vyb(i), vzb(i)

      end if

    end do

    ! Abandon Ship!
    call util_exit(FAILURE)

  end if

  return
  end subroutine helio_drift
!!
  subroutine helio_lindrift(nbod, mass, vxb, vyb, vzb, dt, xh, yh, zh, px, py, pz)
  !-------------------------------------------------------------------------
  !                            HELIO_LINDRIFT.F
  !-------------------------------------------------------------------------
  ! This subroutine takes a linear drift due to mometum of Sun
  !
  !             Input:
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 vxb, vyb, vzb ==>  velocity in bary coord
  !                                    (real arrays)
  !                 dt            ==>  time step
  !                 xh, yh, zh    ==>  initial position in helio coord
  !                                    (real arrays)
  !             Output:
  !                 xh, yh, zh    ==>  final position in helio coord
  !                                    (real arrays)
  !                 px, py, pz    ==>  momentum of sun: tp's need this
  !                                    (real scalars)
  !
  ! Remarks: Bases on Martin's code h2.f
  ! Authors:  Hal Levison
  ! Date:    11/14/96
  ! Last revision: 1/8/97
  implicit none

  ! Inputs Only:
  integer(ik) :: nbod
  real(rk) :: mass(nbod), dt
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

  ! Inputs and Outputs:
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)

  ! Outputs Only:
  real(rk) :: px, py, pz

  ! Internals:
  integer(ik) :: i
  real(rk) :: mstar
  !real(rk) :: pxerr, pyerr, pzerr
  real(rk) :: pxtmp, pytmp, pztmp

  !----
  ! Executable code

  ! Initialize the temporary variables and their error variables
  mstar = mass(1)
  px = 0.0_rk !; pxerr = 0.0_rk
  py = 0.0_rk !; pyerr = 0.0_rk
  pz = 0.0_rk !; pzerr = 0.0_rk

  !FIRSTPRIVATE(mu, dt, pxerr, pyerr, pzerr)
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, pxtmp, pytmp, pztmp) FIRSTPRIVATE(mstar, dt) &
  !$OMP SHARED(nbod, mass, xh, yh, zh, vxb, vyb, vzb, px, py, pz)

  !$OMP DO SCHEDULE(STATIC) REDUCTION(+ : px, py, pz)
  do i = 2, nbod

    px = px + mass(i)*vxb(i)
    py = py + mass(i)*vyb(i)
    pz = pz + mass(i)*vzb(i)
    !px = util_kahan_sum(px, mass(i)*vxb(i), pxerr)
    !py = util_kahan_sum(py, mass(i)*vyb(i), pyerr)
    !pz = util_kahan_sum(pz, mass(i)*vzb(i), pzerr)

  end do

  ! Make a private copy of the momentum vector
  pxtmp = px/mstar
  pytmp = py/mstar
  pztmp = pz/mstar

  !$OMP DO SCHEDULE(STATIC)
  do i = 2, nbod

    if(mass(i) /= 0.0_rk) then

      xh(i) = xh(i) + pxtmp*dt
      yh(i) = yh(i) + pytmp*dt
      zh(i) = zh(i) + pztmp*dt

    end if

  end do
  !$OMP END DO NOWAIT

  !$OMP END PARALLEL

  return
  end subroutine helio_lindrift
!!
!!
end module helio