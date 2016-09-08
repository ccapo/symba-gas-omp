function orbel_eget(e, m) result(ea)
!**********************************************************************
!                    ORBEL_EGET.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           m ==> mean anomaly. (real scalar)
!             Returns:
!                  orbel_eget ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Quartic convergence from Danby
!     REMARKS: For results very near roundoff,  give it M between
!           0 and 2*pi. One can condition M before calling EGET
!           by calling my double precision function MOD2PI(M).
!           This is not done within the routine to speed it up
!           and because it works fine even for large M.
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 7,  1992.
!     REVISIONS: May 21,  1992.  Now have it go through EXACTLY two iterations
!                with the premise that it will only be called if
!	         we have an ellipse with e between 0.15 and 0.8
use module_swift
use module_interfaces, except_this_one => orbel_eget
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

! Function to solve Kepler's eqn for E (here called
! x) for given e and M. returns value of x.
! MAY 21 : FOR e < 0.18 use ESOLMD for speed and sufficient accuracy
! MAY 21 : FOR e > 0.8 use EHIE - this one may not converge fast enough.
call orbel_scget(m, sm, cm)

! begin with a guess accurate to order ecc**3
x = m + e*sm*(1.0_rk + e*(cm + e*(1.0_rk - 1.5_rk*sm**2)))

! Go through one iteration for improved estimate
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

! Do another iteration.
! For m between 0 and 2*pi this seems to be enough to
! get near roundoff error for eccentricities between 0 and 0.8
x = ea
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
end function orbel_eget
