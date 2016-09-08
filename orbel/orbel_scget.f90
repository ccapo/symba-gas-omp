subroutine orbel_scget(angle, sx, cx)
!**********************************************************************
!	                  ORBEL_SCGET.F
!**********************************************************************
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
use module_swift
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