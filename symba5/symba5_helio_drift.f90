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
use module_swift
use module_symba5
use module_interfaces, except_this_one => symba5_helio_drift
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
