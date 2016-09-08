function orbel_zget(q) result(ea)
!**********************************************************************
!                    ORBEL_ZGET.F
!**********************************************************************
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
use module_swift
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