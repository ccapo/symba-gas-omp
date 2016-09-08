subroutine helio_drift(nbod, mass, xh, yh, zh, vxb, vyb, vzb, dt)
!*************************************************************************
!                        HELIO_DRIFT.F
!*************************************************************************
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
use module_swift
use module_interfaces, except_this_one => helio_drift
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
