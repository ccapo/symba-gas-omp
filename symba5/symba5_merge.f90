subroutine symba5_merge(t, dt, nbod, nbodm, ip1, ip2, mass, xh, yh, zh, vxb, vyb, vzb, &
           ireci, lvdotrold, ibound, rpl, mergelst, mergecnt, rhill, eoff, ielc, ielst)
!*************************************************************************
!                            SYMBA5_MERGE.F
!*************************************************************************
! This subroutine checks to see if there are encounters
!
!             Input:
!                 t             ==>  current time (real scalar)
!                 dt            ==>  time step (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 ip1, ip2       ==>  The two bodies to check (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh, yh, zh      ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxb, vyb, vzb   ==>  initial velocity in helio coord 
!                                    (real arrays)
!                 ireci         ==>  current recursion level (int scalar)
!                 lvdotrold     ==>  old radial velocity test
!                                   = .true. if i & j are approaching
!                                   = .false if i & j are receding
!                                     (logical*1 scalar)
!                 iecnt         ==>  The number of objects that each planet 
!                                    is encountering (int*2 array)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!             mergelst          ==>  list of mergers (int array)
!             mergecnt          ==>  count of mergers (int array)
!             rhill             ==>  Hill sphere of planet (real Scalar)
!             eoff              ==>  Energy offset (real scalar)
!                ielc           ==>  number of encounters (integer*2 scalar)
!                ielst          ==>  list of ecnounters (2D integer*2 array)
!
!             Output:  Changed only if a Megrer happens
!                 mass          ==>  mass of bodies (real array)
!                 xh, yh, zh      ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxb, vyb, vzb   ==>  initial velocity in helio coord 
!                                    (real arrays)
!                 iecnt         ==>  The number of objects that each planet 
!                                    is encountering (int*2 array)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!             mergelst          ==>  list of mergers (int array)
!             mergecnt          ==>  count of mergers (int array)
!             rhill             ==>  Hill sphere of planet (real Scalar)
!             eoff              ==>  Energy offset (real scalar)
!                ielc           ==>  number of encounters (integer*2 scalar)
!                ielst          ==>  list of ecnounters (2D integer*2 array)
!
! Remarks: 
! Authors:  Hal Levison
! Date:   1/2/97
! Last revision: 1/24/97
use module_swift
use module_symba5
use module_interfaces, except_this_one => symba5_merge
implicit none

! Inputs:
logical(lk) :: lvdotrold
integer(ik) :: nbod, nbodm, ireci, ip1, ip2
real(rk) :: t, dt

! Inputs and Outputs:
integer(ik) :: mergelst(NENMAX,2), mergecnt
integer(ik) :: ielst(NENMAX,2), ielc, ibound(nbod)
real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod), eoff
real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod), rpl(nbod), rhill(nbod)

! Internals
logical(lk) :: lswitch, lcross
integer(ik) :: ialpha, ip1l, ip2l, itmp
real(rk) :: xr, yr, zr, vxr, vyr, vzr, vdotr, tcross2, dt2
real(rk) :: rlim, rlim2, rr2, vr2, gmsum, energy, a, e, q

!-----------------!
! Executable Code !
!-----------------!

xr = xh(ip2) - xh(ip1)
yr = yh(ip2) - yh(ip1)
zr = zh(ip2) - zh(ip1)
rr2 = xr**2 + yr**2 + zr**2

rlim = rpl(ip1) + rpl(ip2)
rlim2 = rlim**2
if(rlim == 0.0_rk) return ! <=== NOTE!!!!

vxr = vxb(ip2) - vxb(ip1)
vyr = vyb(ip2) - vyb(ip1)
vzr = vzb(ip2) - vzb(ip1)
vdotr = xr*vxr + yr*vyr + zr*vzr
vr2 = vxr**2 + vyr**2 + vzr**2

tcross2 = rr2/vr2
dt2 = dt**2

lswitch = (lvdotrold .and. (vdotr > 0.0_rk))
lcross = (tcross2 <= dt2)

gmsum = mass(ip1) + mass(ip2)
call orbel_xv2aeq(xr, yr, zr, vxr, vyr, vzr, gmsum, ialpha, a, e, q)

! Sort indices in decending order by mass
ip1l = ip1
ip2l = ip2
if(mass(ip2l) > mass(ip1l)) then

  itmp = ip1l
  ip1l = ip2l
  ip2l = itmp

end if

! If particle ip2 is bound to ip1 (ialpha < 0), orbit smaller than 0.1*R_hill of ip1 and it has come in for perigee passage,
! then increment ibound, otherwise set to zero.
if((ialpha < 0) .and. (a < rhill(ip1l)/10.0_rk) .and. lswitch) then

  ibound(ip2l) = ibound(ip2l) + 1

else

  ibound(ip2l) = 0

end if

! If the sum of the radii is larger than the distance between ip1 and ip2, then merge
if(rlim2 >= rr2) then

  !ip1l = ip1
  !ip2l = ip2
  call discard_mass_merge5(t, nbod, nbodm, ip1l, ip2l, mass, xh, yh, zh, vxb, vyb, vzb, rpl, eoff, ielc, ielst)
  mergecnt = mergecnt + 1
  mergelst(mergecnt,1) = ip1l
  mergelst(mergecnt,2) = ip2l
  rhill(ip2l) = 0.0_rk
  ibound(ip2l) = 0

  call util_hills1(mass(1), mass(ip1l), xh(ip1l), yh(ip1l), zh(ip1l), vxb(ip1l), vyb(ip1l), vzb(ip1l), rhill(ip1l))
  return ! <== NOTE!!!!

end if

! If the sign of vdotr has changed from approaching to receding (i.e. from vdotr < 0.0 to vdotr > 0.0)
! and if crossing time is less than the current time step, and periagee is less than rlim, then merge
if(lswitch .and. lcross .and. (q < rlim)) then

  !ip1l = ip1
  !ip2l = ip2
  call discard_mass_merge5(t, nbod, nbodm, ip1l, ip2l, mass, xh, yh, zh, vxb, vyb, vzb, rpl, eoff, ielc, ielst)
  mergecnt = mergecnt + 1
  mergelst(mergecnt,1) = ip1l
  mergelst(mergecnt,2) = ip2l
  rhill(ip2l) = 0.0_rk
  ibound(ip2l) = 0

  call util_hills1(mass(1), mass(ip1l), xh(ip1l), yh(ip1l), zh(ip1l), vxb(ip1l), vyb(ip1l), vzb(ip1l), rhill(ip1l))
  return ! <== NOTE!!!!

end if

! If particle ip2 has been flagged to be bound to particle ip1 more than four times (four perigee passages),
! consider ip2 to be in a satellite of ip1 and argue that we can merge ip2 with ip1 due to the presence of gas
! drag, which is not implemented in these recursive routines.  This can speed up the code since we avoid
! performing 10^5+ time steps to accurately resolve such orbits, but this may accidentily remove some
! particles that are only captured temporarily.
if(ibound(ip2l) >= NBOUNDMAX) then

  !write(*,'(a,i6,a,i9,a)') " Particle: ", ip2, " has been bound ", ibound(ip2), " times"
  !write(*,'(2(1x,1pe22.15))') sqrt(rr2), a
  
  !ip1l = ip1
  !ip2l = ip2
  call discard_mass_merge5(t, nbod, nbodm, ip1l, ip2l, mass, xh, yh, zh, vxb, vyb, vzb, rpl, eoff, ielc, ielst)
  mergecnt = mergecnt + 1
  mergelst(mergecnt,1) = ip1l
  mergelst(mergecnt,2) = ip2l
  rhill(ip2l) = 0.0_rk
  ibound(ip2l) = 0

  call util_hills1(mass(1), mass(ip1l), xh(ip1l), yh(ip1l), zh(ip1l), vxb(ip1l), vyb(ip1l), vzb(ip1l), rhill(ip1l))
  return ! <== NOTE!!!!

end if

return
end subroutine symba5_merge
