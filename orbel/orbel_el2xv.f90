subroutine orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm, x, y, z, vx, vy, vz)
!****************************************************************************
!                          ORBEL_EL2XV.F
!****************************************************************************
!     PURPOSE: To compute cartesian positions and velocities given
!               central mass,  ialpha ( = +1 for hyp.,  0 for para. and
!               -1 for ellipse),  and orbital elements.
!       input:
!            gm       ==> G times central mass (real scalar)
!	     ialpha   ==> conic section type ( see PURPOSE,  integer(ik) :: scalar)
!	     a        ==> semi-major axis or pericentric distance if a parabola
!                          (real scalar)
!            e        ==> eccentricity (real scalar)
!            inc      ==> inclination  (real scalar)
!            capom    ==> longitude of ascending node (real scalar)
!	     omega    ==> argument of perihelion (real scalar)
!	     capm     ==> mean anomoly(real scalar)
!
!       Output:
!            x, y, z    ==>  position of object (real scalars)
!            vx, vy, vz ==>  velocity of object (real scalars)
!
!     ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."
!     REMARKS: All angles are in RADIANS
!
!     AUTHOR:  M. Duncan.
!     DATE WRITTEN:  May 11,  1992.
!     REVISIONS: May 26 - now use better Kepler solver for ellipses
!                 and hyperbolae called EHYBRID.F and FHYBRID.F
use module_swift
use module_interfaces, except_this_one => orbel_el2xv
implicit none

! Inputs Only:
integer(ik) :: ialpha
real(rk) :: gm, a, e, inc, capom, omega, capm

! Outputs:
real(rk) :: x, y, z, vx, vy, vz

! Internals:
real(rk) :: cape, capf, zpara, em1
real(rk) :: sp, cp, so, co, si, ci
real(rk) :: d11, d12, d13, d21, d22, d23
real(rk) :: scap, ccap, shcap, chcap
real(rk) :: sqe, sqgma, xfac1, xfac2, ri, vfac1, vfac2

!----
! Executable code

if(e < 0.0_rk) then

  write(*,'(a)') 'ERROR in orbel_el2xv: e < 0.0, setting e = 0.0'
  e = 0.0_rk

end if

! check for inconsistencies between ialpha and e
em1 = e - 1.0_rk
if(((ialpha == 0) .and. (abs(em1) > TINY)) .or. ((ialpha < 0) .and. (e > 1.0_rk)) .or. ((ialpha > 0) .and. (e < 1.0_rk))) then

  write(*,'(a)')         ' ERROR in orbel_el2xv: ialpha and e are inconsistent'
  write(*,'(a,i3)')      '                       ialpha = ', ialpha
  write(*,'(a,1pe14.6)') '                            e = ', e

end if

! Generate rotation matrices (on p. 42 of Fitzpatrick)
call orbel_scget(omega, sp, cp)
call orbel_scget(capom, so, co)
call orbel_scget(inc, si, ci)
d11 = cp*co - sp*so*ci
d12 = cp*so + sp*co*ci
d13 = sp*si
d21 = -sp*co - cp*so*ci
d22 = -sp*so + cp*co*ci
d23 = cp*si

! Get the other quantities depending on orbit type (i.e. IALPHA)
if(ialpha == -1) then

  cape = orbel_ehybrid(e, capm)
  call orbel_scget(cape, scap, ccap)
  sqe = sqrt(1.0_rk - e*e)
  sqgma = sqrt(gm*a)
  xfac1 = a*(ccap - e)
  xfac2 = a*sqe*scap
  ri = 1.0_rk/(a*(1.0_rk - e*ccap))
  vfac1 = -ri*sqgma*scap
  vfac2 = ri*sqgma*sqe*ccap

end if

if(ialpha == 1) then

  capf = orbel_fhybrid(e, capm)
  call orbel_schget(capf, shcap, chcap)
  sqe = sqrt(e*e - 1.0_rk)
  sqgma = sqrt(gm*a)
  xfac1 = a*(e - chcap)
  xfac2 = a*sqe*shcap
  ri = 1.0_rk/(a*(e*chcap - 1.0_rk))
  vfac1 = -ri*sqgma*shcap
  vfac2 = ri*sqgma*sqe*chcap

end if

if(ialpha == 0) then

  zpara = orbel_zget(capm)
  sqgma = sqrt(2.0_rk*gm*a)
  xfac1 = a*(1.0_rk - zpara*zpara)
  xfac2 = 2.0_rk*a*zpara
  ri = 1.0_rk/(a*(1.0_rk + zpara*zpara))
  vfac1 = -ri*sqgma*zpara
  vfac2 = ri*sqgma

end if

x =  d11*xfac1 + d21*xfac2
y =  d12*xfac1 + d22*xfac2
z =  d13*xfac1 + d23*xfac2
vx = d11*vfac1 + d21*vfac2
vy = d12*vfac1 + d22*vfac2
vz = d13*vfac1 + d23*vfac2

return
end subroutine orbel_el2xv
