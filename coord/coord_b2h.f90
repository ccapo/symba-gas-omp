subroutine coord_b2h(nbod, mass, xb, yb, zb, vxb, vyb, vzb,  xh, yh, zh, vxh, vyh, vzh)
!***********************************************************************
!	                    COORD_B2H.F
!***********************************************************************
!     PURPOSE: Converts from Barycentric to Helio coords.
!     ARGUMENTS:  Input is
!                    nbod ==> number of bodies (must be less than NBMAX)
!                             (integer)
!	             mass(*) ==>  masses (real array)
!                                 NOT USED BUT INCLUDED IN ORDER TO HAVE
!                                 SYMMETRY IN SUBROUTINE CALLS
!		     xb(*), yb(*), zb(*) ==> Barycentric particle coords
!                                          (real array)
!		     vxb(*), vyb(*), vzb(*) ==> Barycentric particle velocities
!                                             (real array)
!                 Returned are
!                    xh(*), yh(*), zh(*) ==> Helio particle positions
!                                          (real array)
!                    vxh(*), vyh(*), vzh(*) ==> Helio particle velocities
!                                            (real array)
!
!     ALGORITHM: Obvious
!     REMARKS:  Can of course use this to get coords. relative to any body.
!              by changing the one subtracted off.
!
!     Authors:  Martin Duncan
!     WRITTEN:  Jan 27/93
!     REVISIONS: 2/17/95  HFL
use module_swift
implicit none

! Inputs:
integer(ik) :: nbod
real(rk) :: mass(nbod)
real(rk) :: xb(nbod), yb(nbod), zb(nbod)
real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

! Outputs:
real(rk) :: xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

! Internals:
integer(ik) :: i
real(rk) :: xb0, yb0, zb0, vxb0, vyb0, vzb0

!----
! Executable code

! Store the barycentric position and velocity of the central body
xb0 = xb(1)
yb0 = yb(1)
zb0 = zb(1)
vxb0 = vxb(1)
vyb0 = vyb(1)
vzb0 = vzb(1)

! Compute the heliocentric position and velocity of all the bodies
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) &
!$OMP FIRSTPRIVATE(xb0, yb0, zb0, vxb0, vyb0, vzb0) &
!$OMP SHARED(nbod, xb, yb, zb, vxb, vyb, vzb, xh, yh, zh, vxh, vyh, vzh)
do i = 1, nbod

  xh(i) = xb(i) - xb0
  yh(i) = yb(i) - yb0
  zh(i) = zb(i) - zb0
  vxh(i) = vxb(i) - vxb0
  vyh(i) = vyb(i) - vyb0
  vzh(i) = vzb(i) - vzb0

end do
!$OMP END PARALLEL DO

return
end subroutine coord_b2h