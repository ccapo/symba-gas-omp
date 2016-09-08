subroutine symba5_step_interp(time, ielev, nbod, nbodm, mass, rhill, j2rp2, j4rp4, rpl, &
           xh, yh, zh, vxh, vyh, vzh, dt, mergelst, mergecnt, eoff, ielc, ielst, mtiny, ibound)
!*************************************************************************
!                            SYMBA5_STEP_INTERP.F
!*************************************************************************
!
!             Input:
!                 time          ==> Current time (real scalar)
!                 ielev         ==>  The level that this particle should go
!                                             (int*2 array)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 nbodm         ==>  Location of last massive body(int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 rhill         ==>  Radius of hill sphere (real array)
!                 j2rp2, j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 xh, yh, zh      ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxh, vyh, vzh   ==>  initial velocity in helio coord 
!                                    (real arrays)
!                 dt            ==>  time step
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!                 eoff          ==>  Energy offset (real scalar)
!                ielc           ==>  number of encounters (integer*2 scalar)
!                ielst          ==>  list of ecnounters (2D integer*2 array)
!                mtiny          ==>  Small mass  (real array)
!             Output:
!                 xh, yh, zh      ==>  final position in helio coord 
!                                       (real arrays)
!                 vxh, vyh, vzh   ==>  final velocity in helio coord 
!                                       (real arrays)
!                 rpl           ==>  Recalculated physical size of a planet.
!                                    if merger happened (real array)
!                 nbod          ==>  Recalculated number of massive bodies 
!                                    if merger happened (int scalar)
!                 nbodm         ==>  Location of last massive body(int scalar)
!                 mass          ==>  Recalculated mass of bodies 
!                                    if merger happened (real array)
!                 mergelst      ==>  list of mergers (int array)
!                 mergecnt      ==>  count of mergers (int array)
!                 eoff          ==>  Energy offset (real scalar)
! Remarks: 
! Authors:  Hal Levison 
! Date:    11/21/96
! Last revision: 5/13/99
use module_swift
use module_symba5
use module_interfaces, except_this_one => symba5_step_interp
implicit none

! Inputs Only:
integer(ik) :: ielst(NENMAX,2), ielc
integer(ik) :: ielev(nbod)
real(rk) :: mass(nbod), dt, j2rp2, j4rp4, time, mtiny

! Inputs and Outputs:
integer(ik) :: nbod, nbodm, ibound(nbod)
real(rk) :: xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)
real(rk) :: rpl(nbod), eoff
real(rk) :: rhill(nbod)

! Outputs
integer(ik) :: mergelst(NENMAX,2), mergecnt

! Internals:
logical(lk) :: lvdotr(NTPMAX) ! Used by symba_step_recur
integer(ik) :: i, irec
real(rk) :: dth, msys
real(rk), save :: axh(NTPMAX), ayh(NTPMAX), azh(NTPMAX)
real(rk), save :: vxb(NTPMAX), vyb(NTPMAX), vzb(NTPMAX)
real(rk) :: ptxb, ptyb, ptzb            ! Not used here
real(rk) :: ptxe, ptye, ptze

!----
! Executable code
!-----------------

! Half the time step
dth = 0.5_rk*dt

! Initialize recursion index and merge counter
irec = 0
mergecnt = 0

! Convert vel to bary to jacobi coords
call coord_vh2b(nbod, mass, vxh, vyh, vzh, vxb, vyb, vzb, msys)
!write(100,'(a,3(1x,1pe22.15))') 'vxh(1:3) = ', vxh(1:3)
!write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

! Do the linear drift due to momentum of the Sun
call helio_lindrift(nbod, mass, vxb, vyb, vzb, dth, xh, yh, zh, ptxb, ptyb, ptzb)
!write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

! Get the accelerations in helio frame. For each object
! only include those guys that it is not encountering with. 
call symba5_getacch(nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh, mtiny, ielc, ielst)
!write(100,'(a,3(1x,1pe22.15))') 'axh(1:3) = ', axh(1:3)

! Apply a heliocentric kick for a half dt 
call kickvh(nbod, vxb, vyb, vzb, axh, ayh, azh, dth)
!write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

! Do a recursion step for full dt
call symba5_helio_drift(nbod, ielev, -1, mass, xh, yh, zh, vxb, vyb, vzb, dt)
!write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

call symba5_step_recur(time, nbod, nbodm, mass, irec, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, &
     rpl, mergelst, mergecnt, dt, eoff, lvdotr, ibound, ielc, ielst) 
!write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

! Get the accelerations in helio frame. For each object
! only include those guys that it is not encountering with. 
call symba5_getacch(nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh, mtiny, ielc, ielst)
!write(100,'(a,3(1x,1pe22.15))') 'axh(1:3) = ', axh(1:3)

! Apply a heliocentric kick for a half dt 
call kickvh(nbod, vxb, vyb, vzb, axh, ayh, azh, dth)
!write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

! Do the linear drift due to momentum of the Sun
call helio_lindrift(nbod, mass, vxb, vyb, vzb, dth, xh, yh, zh, ptxe, ptye, ptze)
!write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

! convert back to helio velocities
call coord_vb2h(nbod, mass, vxb, vyb, vzb, vxh, vyh, vzh)
!write(100,'(a,3(1x,1pe22.15))') 'vxh(1:3) = ', vxh(1:3)
!write(100,*) '--------'

return
end subroutine symba5_step_interp
