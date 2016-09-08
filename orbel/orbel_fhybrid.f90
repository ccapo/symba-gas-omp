function orbel_fhybrid(e, capn) result(ea)
!**********************************************************************
!                    ORBEL_FHYBRID.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                        capn ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: For abs(N) < 0.636*ecc -0.6 ,  use FLON
!	         For larger N,  uses FGET
!     REMARKS:
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 26, 1992.
!     REVISIONS:
!     REVISIONS: 2/26/93 hfl
use module_swift
use module_interfaces, except_this_one => orbel_fhybrid
implicit none

! Inputs Only:
real(rk) :: e, capn

! Output Only:
real(rk) :: ea

! Internals:
real(rk) :: abn

!----
! Executable code

abn = capn
if(capn < 0.0_rk) abn = -abn

if(abn < 0.636_rk*e - 0.6_rk) then

  ea = orbel_flon(e, capn)

else

  ea = orbel_fget(e, capn)

end if

return
end function orbel_fhybrid
