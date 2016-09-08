function orbel_ehie(e, m) result(ea)
!**********************************************************************
!                    ORBEL_EHIE.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           m ==> mean anomaly. (real scalar)
!             Returns:
!              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Use Danby's quartic for 3 iterations.
!                Eqn. is f(x) = x - e*sin(x+M). Note  that
!	         E = x + M. First guess is very good for e near 1.
!	         Need to first get M between 0. and PI and use
!		 symmetry to return right answer if M between PI and 2PI
!     REMARKS: Modifies M so that both E and M are in range (0, TWOPI)
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 25, 1992.
!     REVISIONS:
use module_swift
use module_interfaces, except_this_one => orbel_ehie
implicit none

! Inputs Only:
real(rk) :: e, m

! Output Only:
real(rk) :: ea

! Internals:
integer(ik), parameter :: NMAX = 3
integer(ik) :: iflag, nper, niter
real(rk) :: dx, x, sa, ca, esa, eca, f, fp

!----
! Executable code

! In this section,  bring M into the range (0, TWOPI) and if
! the result is greater than PI,  solve for (TWOPI - M).
iflag = 0
nper = m/TWOPI
m = m - nper*TWOPI
if(m < 0.0_rk) m = m + TWOPI

if(m > PI) then

  m = TWOPI - m
  iflag = 1

end if

! Make a first guess that works well for e near 1.
x = (6.0_rk*m)**(1.0_rk/3.0_rk) - m

! Iteration loop
do niter = 1, NMAX

  call orbel_scget(x + m, sa, ca)
  esa = e*sa
  eca = e*ca
  f = x - esa
  fp = 1.0_rk - eca
  dx = -f/fp
  dx = -f/(fp + 0.5_rk*dx*esa)
  dx = -f/(fp + 0.5_rk*dx*(esa + eca*dx/3.0_rk))
  x = x + dx

end do

ea = m + x

if(iflag == 1) then

  ea = TWOPI - ea
  m = TWOPI - m

end if

return
end function orbel_ehie
