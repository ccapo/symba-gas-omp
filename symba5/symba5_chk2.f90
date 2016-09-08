subroutine symba5_chk2(rhill, xh, yh, zh, vxh, vyh, vzh, dt, irec, icflg)
!*************************************************************************
!                            SYMBA5_CHK2.F
!*************************************************************************
! This subroutine checks to see if there are encounters
!
!             Input:
!                 rhill         ==>  Radius of hill sphere (real array)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 ip1, ip2       ==>  The two bodies to check (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh, yh, zh      ==>  initial position in helio coord
!                                    (real arrays)
!                 vxh, vyh, vzh   ==>  initial velocity in helio coord
!                                    (real arrays)
!                 dt            ==>  time step  (real scalor)
!                 irec          ==>  current recursion level (int scalar)
!             Output:
!                 icflg         ==> ecounter?  = 1 Yes
!                                              = 0 No (int scalar)
!                 lvdotr        ==> = .true. if i & j are approaching
!                                   = .false if i & j are receding
!                                     (logical*1 scalar)
!
! Remarks: Based on plh_chk.f.  Same as symba_chk.f
! Authors:  Hal Levison
! Date: 03/20/97
! Last revision:
use module_swift
use module_symba5
use module_interfaces, except_this_one => symba5_chk2
implicit none

! Inputs:
integer(ik) :: irec
real(rk) :: xh(2), yh(2), zh(2), dt
real(rk) :: vxh(2), vyh(2), vzh(2), rhill(2)

! Outputs
integer(ik) :: icflg

! Internals
real(rk) :: r2crit, r2min
real(rk) :: xr, yr, zr, vxr, vyr, vzr
real(rk) :: dr2, dv2, vdotr, tmin

!-----------------
! Executable code
!-----------------

r2crit = ((rhill(1) + rhill(2))*RHSCALE*RSHELL**irec)**2

xr = xh(2) - xh(1)
yr = yh(2) - yh(1)
zr = zh(2) - zh(1)
dr2 = xr**2 + yr**2 + zr**2

vxr = vxh(2) - vxh(1)
vyr = vyh(2) - vyh(1)
vzr = vzh(2) - vzh(1)
dv2 = vxr**2 + vyr**2 + vzr**2

vdotr = xr*vxr + yr*vyr + zr*vzr

tmin = -vdotr/dv2

if(vdotr > 0.0_rk) then

  if(dr2 >= r2crit) then

    icflg = 0

  else

    icflg = 1

  end if

else

  ! We are converging, so we need to calculate the minimum separation attained in time dt
  tmin = -vdotr/dv2

  if(tmin < dt) then

    r2min = dr2 - (vdotr**2)/dv2

  else

    r2min = dr2 + 2.0_rk*vdotr*dt + dv2*dt**2

  end if

  r2min = min(r2min, dr2) ! Really make sure

  if(r2min <= r2crit) then

    icflg = 1

  else

    icflg = 0

  end if

end if

return
end subroutine symba5_chk2
