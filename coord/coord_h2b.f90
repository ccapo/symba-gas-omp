subroutine coord_h2b(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xb, yb, zb, vxb, vyb, vzb, msys)
!***********************************************************************
!	                    COORD_H2B.F
!***********************************************************************
!     PURPOSE: Converts from Heliocentric to Barycentric coords.
!     ARGUMENTS:  Input is
!                    nbod ==> number of bodies (must be less than NBMAX)
!                             (integer)
!	             mass(*) ==>  masses (real array)
!		     xh(*), yh(*), zh(*) ==> heliocentric particle coords
!                                          (real array)
!		     vxh(*), vyh(*), vzh(*) ==> heliocentric particle velocities
!                                             (real array)
!                 Returned are
!                    xb(*), yb(*), zb(*) ==> bary. particle positions
!                                          (real array)
!                    vxb(*), vyb(*), vzb(*) ==> bary. particle velocities
!                                            (real array)
!                    msys              ==>  Total mass of of system
!                                            (real scalar)
!     Authors:  Martin Duncan
!     ALGORITHM: Obvious
!     WRITTEN:  Jan 27/93
!     REVISIONS: 2/22/94  HFL
use module_swift
use module_interfaces, except_this_one => coord_h2b
implicit none

! Inputs:
integer(ik) :: nbod
real(rk) :: mass(nbod)
real(rk) :: xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

! Outputs:
real(rk) :: xb(nbod), yb(nbod), zb(nbod)
real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

! Internals:
integer(ik) :: i
real(rk) :: msys, x, y, z, vx, vy, vz
!real(rk) :: merr, xerr, yerr, zerr, vxerr, vyerr, vzerr
real(rk) :: xtmp, ytmp, ztmp, vxtmp, vytmp, vztmp

!----
! Executable code

! Initialize the total system mass, and the temporary variables and their error variables
msys = 0.0_rk !; merr = 0.0_rk
x = 0.0_rk !; xerr = 0.0_rk
y = 0.0_rk !; yerr = 0.0_rk
z = 0.0_rk !; zerr = 0.0_rk
vx = 0.0_rk !; vxerr = 0.0_rk
vy = 0.0_rk !; vyerr = 0.0_rk
vz = 0.0_rk !; vzerr = 0.0_rk

!FIRSTPRIVATE(merr, xerr, yerr, zerr, vxerr, vyerr, vzerr)
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, xtmp, ytmp, ztmp, vxtmp, vytmp, vztmp) &
!$OMP SHARED(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xb, yb, zb, vxb, vyb, vzb, msys, x, y, z, vx, vy, vz)

! Compute the total mass of the system, along with the mass-weighted sum of positions and velocities
!$OMP DO SCHEDULE(STATIC) REDUCTION(+ : msys, x, y, z, vx, vy, vz)
do i = 1, nbod

  msys = msys + mass(i)
  x = x + mass(i)*xh(i)
  y = y + mass(i)*yh(i)
  z = z + mass(i)*zh(i)
  vx = vx + mass(i)*vxh(i)
  vy = vy + mass(i)*vyh(i)
  vz = vz + mass(i)*vzh(i)
  !msys = util_kahan_sum(msys, mass(i), merr)
  !x = util_kahan_sum(x, mass(i)*xh(i), xerr)
  !y = util_kahan_sum(y, mass(i)*yh(i), yerr)
  !z = util_kahan_sum(z, mass(i)*zh(i), zerr)
  !vx = util_kahan_sum(vx, mass(i)*vxh(i), vxerr)
  !vy = util_kahan_sum(vy, mass(i)*vyh(i), vyerr)
  !vz = util_kahan_sum(vz, mass(i)*vzh(i), vzerr)

end do
!$OMP END DO

xtmp = x/msys
ytmp = y/msys
ztmp = z/msys
vxtmp = vx/msys
vytmp = vy/msys
vztmp = vz/msys

!$OMP DO SCHEDULE(STATIC)
do i = 1, nbod

  xb(i) = xh(i) - xtmp
  yb(i) = yh(i) - ytmp
  zb(i) = zh(i) - ztmp
  vxb(i) = vxh(i) - vxtmp
  vyb(i) = vyh(i) - vytmp
  vzb(i) = vzh(i) - vztmp

end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

return
end subroutine coord_h2b
