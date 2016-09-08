subroutine drift_one(mu, x, y, z, vx, vy, vz, dt, iflg)
!*************************************************************************
!                        DRIFT_ONE.F
!*************************************************************************
! This subroutine does the danby-type drift for one particle,  using
! appropriate vbles and redoing a drift if the accuracy is too poor
! (as flagged by the integer iflg).
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mu            ==>  mass of central body (real scalar)
!                 x, y, z         ==>  initial position in jacobi coord
!                                    (real scalar)
!                 vx, vy, vz      ==>  initial position in jacobi coord
!                                    (real scalar)
!                 dt            ==>  time step
!             Output:
!                 x, y, z         ==>  final position in jacobi coord
!                                       (real scalars)
!                 vx, vy, vz      ==>  final position in jacobi coord
!                                       (real scalars)
!                 iflg          ==>  integer (zero for successful step)
!
! Authors:  Hal Levison & Martin Duncan
! Date:    2/10/93
! Last revision: 2/10/93
use module_swift
use module_interfaces, except_this_one => drift_one
implicit none

! Inputs Only:
real(rk) :: mu, dt

! Inputs and Outputs:
real(rk) :: x, y, z
real(rk) :: vx, vy, vz

! Output
integer(ik) :: iflg

! Internals:
integer(ik) :: i
real(rk) :: dttmp

!----
! Executable code

call drift_dan(mu, x, y, z, vx, vy, vz, dt, iflg)

if(iflg /= 0) then

  dttmp = 0.1_rk*dt

  do i = 1, 10

    call drift_dan(mu, x, y, z, vx, vy, vz, dttmp, iflg)
    if(iflg /= 0) return ! Abandon all hope ye who enter here

  end do

end if

return
end subroutine drift_one
