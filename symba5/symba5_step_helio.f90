subroutine symba5_step_helio(i1st, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, dt)
!*************************************************************************
!                            SYMBA5_STEP_HELIO.F
!*************************************************************************
! This subroutine takes a step in helio coord.
! Does a KICK than a DRIFT than a KICK.
! ONLY DOES MASSIVE PARTICLES
!
!             Input:
!                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 j2rp2, j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 xh, yh, zh      ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxh, vyh, vzh   ==>  initial velocity in helio coord 
!                                    (real arrays)
!                 dt            ==>  time step
!             Output:
!                 xh, yh, zh      ==>  final position in helio coord 
!                                       (real arrays)
!                 vxh, vyh, vzh   ==>  final velocity in helio coord 
!                                       (real arrays)
! Remarks: Based on helio_step_pl.f but does not pass the intermediate
!          positions and velocities back for the TP to use.
! Authors:  Hal Levison 
! Date:    3/20/97
! Last revision: 12/13/00
use module_swift
use module_interfaces, except_this_one => symba5_step_helio
implicit none

! Inputs Only: 
integer(ik) :: nbod, i1st, nbodm
real(rk) :: mass(nbod), dt, j2rp2, j4rp4

! Inputs and Outputs:
real(rk) :: xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

! Internals:
integer(ik) :: i1stloc
real(rk) :: dth, msys
real(rk) :: axh(nbod), ayh(nbod), azh(nbod)
real(rk), save :: vxb(NTPMAX), vyb(NTPMAX), vzb(NTPMAX)
real(rk) :: ptxb, ptyb, ptzb            ! Not used here
real(rk) :: ptxe, ptye, ptze

!----
! Executable code
!-----------------

dth = 0.5_rk*dt

i1stloc = i1st
if(i1st == 0) then

  ! Convert vel to bery to jacobi coords
  call coord_vh2b(nbod, mass, vxh, vyh, vzh, vxb, vyb, vzb, msys)
  i1st = 1              ! turn this off

end if

! Do the linear drift due to momentum of the Sun
call helio_lindrift(nbod, mass, vxb, vyb, vzb, dth, xh, yh, zh, ptxb, ptyb, ptzb)

! Get the accelerations in helio frame, if first time step
call symba5_helio_getacch(i1stloc, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh)

! Apply a heliocentric kick for a half dt 
call kickvh(nbod, vxb, vyb, vzb, axh, ayh, azh, dth)

! Drift in helio coords for the full step 
call helio_drift(nbod, mass, xh, yh, zh, vxb, vyb, vzb, dt)

! Get the accelerations in helio frame, regardless if it is the first step
call symba5_helio_getacch(0, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh)

! Apply a heliocentric kick for a half dt 
call kickvh(nbod, vxb, vyb, vzb, axh, ayh, azh, dth)

! Do the linear drift due to momentum of the Sun
call helio_lindrift(nbod, mass, vxb, vyb, vzb, dth, xh, yh, zh, ptxe, ptye, ptze)

! convert back to helio velocities
call coord_vb2h(nbod, mass, vxb, vyb, vzb, vxh, vyh, vzh)

return
end subroutine symba5_step_helio
