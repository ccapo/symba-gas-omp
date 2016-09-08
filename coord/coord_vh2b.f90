 subroutine coord_vh2b(nbod, mass, vxh, vyh, vzh, vxb, vyb, vzb, msys)
!***********************************************************************
!	                    COORD_VH2B.F
!***********************************************************************
!     PURPOSE: Converts from Heliocentric to Barycentric coords.
!              Velocity only
!     ARGUMENTS:  Input is
!                    nbod ==> number of bodies (must be less than NBMAX)
!                             (integer)
!	             mass(*) ==>  masses (real array)
!		     vxh(*), vyh(*), vzh(*) ==> heliocentric particle velocities
!                                             (real array)
!                 Returned are
!                    vxb(*), vyb(*), vzb(*) ==> bary. particle velocities
!                                            (real array)
!                    msys              ==>  Total mass of of system
!                                            (real scalar)
!     Authors:  Hal Levison
!     ALGORITHM: Obvious
!     WRITTEN:  11/14/96
!     REVISIONS: 11/21/96
use module_swift
use module_interfaces, except_this_one => coord_vh2b
implicit none

! Inputs:
integer(ik) :: nbod
real(rk) :: mass(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

! Outputs:
real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

! Internals:
integer(ik) :: i
real(rk) :: msys, vx, vy, vz
!real(rk) :: merr, vxerr, vyerr, vzerr
real(rk) :: vxtmp, vytmp, vztmp

!----
! Executable code

msys = mass(1) !; merr = 0.0_rk
vx = 0.0_rk !; vxerr = 0.0_rk
vy = 0.0_rk !; vyerr = 0.0_rk
vz = 0.0_rk !; vzerr = 0.0_rk

!FIRSTPRIVATE(merr, vxerr, vyerr, vzerr)
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, vxtmp, vytmp, vztmp) &
!$OMP SHARED(nbod, mass, vxh, vyh, vzh, vxb, vyb, vzb, msys, vx, vy, vz)

!$OMP DO SCHEDULE(STATIC) REDUCTION(+ : msys, vx, vy, vz)
do i = 2, nbod

  msys = msys + mass(i)
  vx = vx + mass(i)*vxh(i)
  vy = vy + mass(i)*vyh(i)
  vz = vz + mass(i)*vzh(i)
  !msys = util_kahan_sum(msys, mass(i), merr)
  !vx = util_kahan_sum(vx, mass(i)*vxh(i), vxerr)
  !vy = util_kahan_sum(vy, mass(i)*vyh(i), vyerr)
  !vz = util_kahan_sum(vz, mass(i)*vzh(i), vzerr)

end do
!$OMP END DO

vxtmp = -vx/msys
vytmp = -vy/msys
vztmp = -vz/msys

!$OMP DO SCHEDULE(STATIC)
do i = 2, nbod

  vxb(i) = vxh(i) + vxtmp
  vyb(i) = vyh(i) + vytmp
  vzb(i) = vzh(i) + vztmp

end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

vxb(1) = -vx/msys
vyb(1) = -vy/msys
vzb(1) = -vz/msys

return
end subroutine coord_vh2b
