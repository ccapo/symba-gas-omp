subroutine util_hills(nbod, mass, xh, yh, zh, vxh, vyh, vzh, r2hill)
!-------------------------------------------------------------------------
!    UTIL_HILLS.F90
!-------------------------------------------------------------------------
! This subroutine calculates the hill's sphere for the planets
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh, yh, zh      ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxh, vyh, vzh   ==>  initial velocity in helio coord 
!                                    (real arrays)
!             Output:
!                  r2hill       ==>  the SQUARE of the planet's hill's sphere 
!                                    (real array)
!
!
! Remarks: 
! Authors:  Hal Levison 
! Date:    2/19/93
! Last revision: 1/6/97
!-------------------------------------------------------------------------
use module_swift
use module_interfaces, except_this_one => util_hills
implicit none

! Inputs: 
integer(ik) :: nbod
real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

! Outputs
real(rk) :: r2hill(nbod)

! Internals
integer(ik) :: i
real(rk) :: mass0, mu, energy, ap, rhill, r, v2

!-----------------!
! Executable Code !
!-----------------!

r2hill(1) = 0.0_rk
mass0 = mass(1)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, mu, r, v2, energy, ap, rhill) &
!$OMP FIRSTPRIVATE(mass0) SHARED(nbod, mass, xh, yh, zh, vxh, vyh, vzh, r2hill)
do i = 2, nbod

  if(mass(i) /= 0.0_rk) then

    mu = mass0*mass(i)/(mass0 + mass(i))
    r = sqrt(xh(i)*xh(i) + yh(i)*yh(i) + zh(i)*zh(i))
    v2 = vxh(i)*vxh(i) + vyh(i)*vyh(i) + vzh(i)*vzh(i)
    energy = 0.5_rk*mu*v2 - mass0*mass(i)/r
    ap = -mass0*mass(i)/(2.0_rk*energy)
    rhill = ap*(mu/(3.0_rk*mass0))**(1.0_rk/3.0_rk)
    r2hill(i) = rhill*rhill

  else

    r2hill(i) = 0.0_rk

  end if

end do
!$OMP END PARALLEL DO

return
end subroutine util_hills
