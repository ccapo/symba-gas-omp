function orbel_esolmd(e, m) result(ea)
!**********************************************************************
!                    ORBEL_ESOLMD.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           m ==> mean anomaly. (real scalar)
!             Returns:
!                orbel_esolmd ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Some sort of quartic convergence from Wisdom.
!     REMARKS: ONLY GOOD FOR SMALL ECCENTRICITY SINCE IT ONLY
!         ITERATES ONCE. (GOOD FOR PLANET CALCS.)
!      	  ALSO DOES NOT PUT M OR E BETWEEN 0. AND 2*PI
!     INCLUDES: needs SCGET.F
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 7,  1992.
!     REVISIONS: 2/26/93 hfl
use module_swift
use module_interfaces, except_this_one => orbel_esolmd
implicit none

! Inputs Only:
real(rk) :: e, m

! Output Only:
real(rk) :: ea

! Internals:
real(rk) :: x, sm, cm, sx, cx
real(rk) :: es, ec, f, fp, fpp, fppp, dx

!----
! Executable code

! Function to solve Kepler's eqn for E (here called x) for given e and M. returns value of x.
call orbel_scget(m, sm, cm)
x = m + e*sm*(1.0_rk + e*(cm + e*(1.0_rk - 1.5_rk*sm*sm)))

call orbel_scget(x, sx, cx)
es = e*sx
ec = e*cx
f = x - es - m
fp = 1.0_rk - ec
fpp = es
fppp = ec
dx = -f/fp
dx = -f/(fp + 0.5_rk*dx*fpp)
dx = -f/(fp + 0.5_rk*dx*fpp + dx*dx*fppp/6.0_rk)

ea = x + dx

return
end function orbel_esolmd
