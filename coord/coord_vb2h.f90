subroutine coord_vb2h(nbod, mass, vxb, vyb, vzb, vxh, vyh, vzh)
!***********************************************************************
!	                    COORD_VB2H.F
!***********************************************************************
!     PURPOSE: Converts from Barycentric to Helio coords.
!               Velocity only
!     ARGUMENTS:  Input is
!                    nbod ==> number of bodies (must be less than NBMAX)
!                             (integer)
!	             mass(*) ==>  masses (real array)
!                                 NOT USED BUT INCLUDED IN ORDER TO HAVE
!                                 SYMMETRY IN SUBROUTINE CALLS
!		     vxb(*), vyb(*), vzb(*) ==> Barycentric particle velocities
!                                             (real array)
!                 Returned are
!                    vxh(*), vyh(*), vzh(*) ==> Helio particle velocities
!                                            (real array)
!
!     ALGORITHM: Obvious
!     Authors:  Hal Levison
!     WRITTEN:  11/14/96
!     REVISIONS: 11/21/96
use module_swift
use module_interfaces, except_this_one => coord_vb2h
implicit none

! Inputs:
integer(ik) :: nbod
real(rk) :: mass(nbod)
real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

! Outputs:
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

! Internals:
integer(ik) :: i
real(rk) :: vx, vy, vz
!real(rk) :: vxerr, vyerr, vzerr
real(rk) :: mstar, vxtmp, vytmp, vztmp

!-----------------
! Executable code
!-----------------

! Initialize temporary variables and their error variables
mstar = mass(1)
vx = -mass(2)*vxb(2) !; vxerr = 0.0_rk
vy = -mass(2)*vyb(2) !; vyerr = 0.0_rk
vz = -mass(2)*vzb(2) !; vzerr = 0.0_rk

!FIRSTPRIVATE(mstar, vxerr, vyerr, vzerr)
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, vxtmp, vytmp, vztmp) FIRSTPRIVATE(mstar) &
!$OMP SHARED(nbod, mass, vxh, vyh, vzh, vxb, vyb, vzb, vx, vy, vz)

!$OMP DO SCHEDULE(STATIC) REDUCTION(+ : vx, vy, vz)
do i = 3, nbod

  vx = vx - mass(i)*vxb(i)
  vy = vy - mass(i)*vyb(i)
  vz = vz - mass(i)*vzb(i)
  !vx = util_kahan_sum(vx, -mass(i)*vxb(i), vxerr)
  !vy = util_kahan_sum(vy, -mass(i)*vyb(i), vyerr)
  !vz = util_kahan_sum(vz, -mass(i)*vzb(i), vzerr)

end do
!$OMP END DO

vxtmp = vx/mstar
vytmp = vy/mstar
vztmp = vz/mstar

!$OMP DO SCHEDULE(STATIC)
do i = 2, nbod

  vxh(i) = vxb(i) - vxtmp
  vyh(i) = vyb(i) - vytmp
  vzh(i) = vzb(i) - vztmp

end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

vxh(1) = 0.0_rk
vyh(1) = 0.0_rk
vzh(1) = 0.0_rk

!mstar = mass(1)
!vxb(1) = vx/mstar
!vyb(1) = vy/mstar
!vzb(1) = vz/mstar

return
end subroutine coord_vb2h
