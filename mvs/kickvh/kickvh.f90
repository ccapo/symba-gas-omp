subroutine kickvh(nbod, vxh, vyh, vzh, axh, ayh, azh, dt)
!*************************************************************************
!                        KICKVH.F
!*************************************************************************
! To kick the velocity components vxh(*) by axh(*)*dt
!
!             Input:
!                 nbod          ==>  number of bodies (int scalar)
!                 vxh, vyh, vzh   ==>  initial velocity in helio coord
!                                    (real arrays)
!                 axh, ayh, azh   ==>  acceleration in helio coord
!                                    (real arrays)
!                 dt            ==>  time step
!             Output:
!                 vxh, vyh, vzh   ==>  final velocity in helio coord
!                                    (real arrays)
!
!     ALGORITHM: Obvious
!     REMARKS:  Only alters particles 2 thru nbod since Sun is #1
!
!     AUTHOR:  M. Duncan.
!     DATE WRITTEN:  Feb. 2,  1993.
!     REVISIONS: 2/18/93   HFL
use module_swift
implicit none

! Inputs Only:
integer(ik) :: nbod
real(rk) :: axh(nbod), ayh(nbod), azh(nbod)
real(rk) :: dt

! Inputs and Output:
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

! Internals:
integer(ik) :: i

!----
! Executable code

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(dt) &
!$OMP SHARED(nbod, vxh, vyh, vzh, axh, ayh, azh)
do i = 2, nbod

  vxh(i) = vxh(i) + axh(i)*dt
  vyh(i) = vyh(i) + ayh(i)*dt
  vzh(i) = vzh(i) + azh(i)*dt

end do
!$OMP END PARALLEL DO

return
end subroutine kickvh