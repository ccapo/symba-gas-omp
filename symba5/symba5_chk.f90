subroutine symba5_chk(rhill, nbod, ip1, ip2, mass, xh, yh, zh, vxh, vyh, vzh, dt, irec, ieflag, lvdotr)
!*************************************************************************
!                            SYMBA5_CHK.F
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
!                 ieflag        ==> ecounter?  = 1 Yes
!                                              =  0 No (integer(ik) :: scalar)
!                 lvdotr        ==> = .true. if i & j are approaching
!                                   = .false if i & j are receding
!                                     (logical*1 scalar)
!
! Remarks: Based on plh_chk.f.  Same as symba_chk.f
! Authors:  Hal Levison
! Date:   3/20/97
! Last revision:
use module_swift
use module_symba5
use module_interfaces, except_this_one => symba5_chk
implicit none

! Inputs:
integer(ik) :: nbod, irec, ip1, ip2
real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod), dt
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod), rhill(nbod)

! Outputs
integer(ik) :: ieflag
logical(lk) :: lvdotr

! Internals
real(rk) :: r2crit, r2critp, rcrit
real(rk) :: xr, yr, zr, vxr, vyr, vzr
real(rk) :: vdotr

!-----------------
! Executable code
!-----------------

rcrit = (rhill(ip1) + rhill(ip2))*RHSCALE*RSHELL**irec
r2crit = rcrit*rcrit
r2critp = -1.0_rk          ! not used here

xr = xh(ip2) - xh(ip1)
yr = yh(ip2) - yh(ip1)
zr = zh(ip2) - zh(ip1)
vxr = vxh(ip2) - vxh(ip1)
vyr = vyh(ip2) - vyh(ip1)
vzr = vzh(ip2) - vzh(ip1)
call rmvs_chk_ind(xr, yr, zr, vxr, vyr, vzr, dt, r2crit, r2critp, ieflag)

vdotr = xr*vxr + yr*vyr + zr*vzr
lvdotr = vdotr < 0.0_rk

!rr2 = xr**2 + yr**2 + zr**2
!vr2 = vxr**2 + vyr**2 + vzr**2
!massc = mass(ip1)
!energy = 0.5_rk*vr2 - massc/sqrt(rr2)
!lbound = (energy < 0.0_rk) .and. (rr2 < rhill(ip1)

return
end subroutine symba5_chk
