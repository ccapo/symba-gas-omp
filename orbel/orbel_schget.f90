subroutine orbel_schget(angle, shx, chx)
!**********************************************************************
!	                  ORBEL_SCHGET.F
!**********************************************************************
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
use module_swift
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