subroutine drift_kepu_stumpff(x, c0, c1, c2, c3)
!*************************************************************************
!                        DRIFT_KEPU_STUMPFF.F
!*************************************************************************
! subroutine for the calculation of stumpff functions
! see Danby p.172  equations 6.9.15
!
!             Input:
!                 x             ==>  argument
!             Output:
!                 c0, c1, c2, c3   ==>  c's from p171-172
!                                       (real scalors)
! Author:  Hal Levison
! Date:    2/3/93
! Last revision: 2/3/93
use module_swift
implicit none

! Inputs:
real(rk) :: x

! Outputs:
real(rk) :: c0, c1, c2, c3

! Internals:
real(rk), parameter :: xm = 0.1_rk
integer(ik) :: n, i

!----
! Executable code

n = 0
do while(abs(x) >= xm)

  n = n + 1
  x = x/4.0_rk

end do

c2 = (1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x/182.0_rk)/132.0_rk)/90.0_rk)/56.0_rk)/30.0_rk) &
     /12.0_rk)/2.0_rk
c3 = (1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x/210.0_rk)/156.0_rk)/110.0_rk)/72.0_rk)/42.0_rk) &
     /20.0_rk)/6.0_rk
c1 = 1.0_rk - x*c3
c0 = 1.0_rk - x*c2

if(n /= 0) then

  do i = n, 1, -1

    c3 = (c2 + c0*c3)/4.0_rk
    c2 = (c1**2)/2.0_rk
    c1 = c0*c1
    c0 = 2.0_rk*c0**2 - 1.0_rk
    x = 4.0_rk*x

  end do

end if

return
end subroutine drift_kepu_stumpff