function orbel_ehybrid(e, m) result(ea)
!**********************************************************************
!                    ORBEL_EHYBRID.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           m ==> mean anomaly. (real scalar)
!             Returns:
!              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: For e < 0.18 uses fast routine ESOLMD
!	         For larger e but less than 0.8,  uses EGET
!	         For e > 0.8 uses EHIE
!     REMARKS: Only EHIE brings M and E into range (0, TWOPI)
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 25, 1992.
!     REVISIONS: 2/26/93 hfl
use module_swift
use module_interfaces, except_this_one => orbel_ehybrid
implicit none

! Inputs Only:
real(rk) :: e, m

! Output Only:
real(rk) :: ea

!----
! Executable code

if(e < 0.18_rk) then

  ea = orbel_esolmd(e, m)

else

  if(e <= 0.8_rk) then

    ea = orbel_eget(e, m)

  else

    ea = orbel_ehie(e, m)

  end if

end if

return
end function orbel_ehybrid
