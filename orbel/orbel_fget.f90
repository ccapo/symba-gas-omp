function orbel_fget(e, capn) result(ea)
!**********************************************************************
!                    ORBEL_FGET.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                        capn ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!                  orbel_fget ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
!           Cel. Mech. ".  Quartic convergence from Danby's book.
!     REMARKS:
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 11, 1992.
!     REVISIONS: 2/26/93 hfl
use module_swift
use module_interfaces, except_this_one => orbel_fget
implicit none

! Inputs Only:
real(rk) :: e, capn

! Output Only:
real(rk) :: ea

! Internals:
integer(ik), parameter :: IMAX = 10
integer(ik) :: i
real(rk) :: tmp, x, shx, chx
real(rk) :: esh, ech, f, fp, fpp, fppp, dx

!----
! Executable code

! Function to solve "Kepler's eqn" for F (here called x) for given e and CAPN

! begin with a guess proposed by Danby
if(capn < 0.0_rk) then

  tmp = -2.0_rk*capn/e + 1.8_rk
  x = -log(tmp)

else

  tmp = 2.0_rk*capn/e + 1.8_rk
  x = log(tmp)

end if

ea = x

do i = 1, IMAX

  call orbel_schget(x, shx, chx)
  esh = e*shx
  ech = e*chx
  f = esh - x - capn
  !write(*,*) 'i, x, f: ', i, x, f
  fp = ech - 1.0_rk
  fpp = esh
  fppp = ech
  dx = -f/fp
  dx = -f/(fp + 0.5_rk*dx*fpp)
  dx = -f/(fp + 0.5_rk*dx*fpp + dx*dx*fppp/6.0_rk)
  ea = x + dx

  ! If we have converged here there's no point in going on
  if(abs(dx) <= TINY) return
  x = ea

end do

write(*,*) 'FGET: RETURNING WITHOUT COMPLETE CONVERGENCE'

return
end function orbel_fget
