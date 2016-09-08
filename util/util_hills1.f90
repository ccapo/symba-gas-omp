subroutine util_hills1(mstar, mpl, xh, yh, zh, vxh, vyh, vzh, rhill)
!-------------------------------------------------------------------------
!    UTIL_HILLS1.F90
!-------------------------------------------------------------------------
! This subroutine calculates the hill's sphere for the planets
!
!             Input:
!                 mstar          ==>  mass of sun (real scalar)
!                 mpl           ==>  mass of sun (real scalar)
!                 xh, yh, zh      ==>  position of pl in helio coord 
!                                    (real scalars)
!                 vxh, vyh, vzh   ==>  velocity of pl in helio coord 
!                                    (real scalars)
!             Output:
!                  rhill        ==>  the radius of planet's hill's sphere 
!                                    (real scalar)
!
!
! Remarks: Based on util_hill
! Authors:  Hal Levison 
! Date:    1/8/97
! Last revision:
!-------------------------------------------------------------------------
use module_swift
use module_interfaces, except_this_one => util_hills1
implicit none

! Inputs:
real(rk) :: mstar, mpl, xh, yh, zh
real(rk) :: vxh, vyh, vzh

! Outputs
real(rk) :: rhill

! Internals
real(rk) :: mu, energy, ap, r, v2

!-----------------!
! Executable code !
!-----------------!

mu = mstar*mpl/(mstar + mpl)
r = sqrt(xh*xh + yh*yh + zh*zh)
v2 = vxh*vxh + vyh*vyh + vzh*vzh
energy = 0.5_rk*mu*v2 - mstar*mpl/r
ap = -mstar*mpl/(2.0_rk*energy)
rhill = ap*(mu/(3.0_rk*mstar))**(1.0_rk/3.0_rk)

return
end subroutine util_hills1
