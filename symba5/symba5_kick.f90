subroutine symba5_kick(nbod, mass, irec, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dt, sgn, ielc, ielst)
!*************************************************************************
!                             SYMBA5_KICK.F
!*************************************************************************
! Do a symba5 kick
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 irec          ==>  recursion level  (integer scalar)
!                 iecnt         ==>  The number of objects that each planet 
!                                    is encountering (int*2 array)
!                 ielev         ==>  The level that this particle should go
!                                             (int*2 array)
!                 rhill         ==>  Hill sphere of planet (real Scalar)
!                 xh, yh, zh      ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxb, vyb, vzb   ==>  initial velocity in bari coord 
!                                    (real arrays)
!                dt             ==>  timestep  (real scalar)
!                sgn            ==>  add or subtract the force (real scalar)
!                ielc           ==>  number of encounters (integer*2 scalar)
!                ielst          ==>  list of ecnounters (2D integer*2 array)
!            Output:
!                 vxb, vyb, vzb   ==>  final velocity in bari coord 
!                                    (real arrays)
!
! Remarks: Uses Man Hoi's force
! Authors:  Hal Levison 
! Date:   3/20/97
! Last revision: 
use module_swift
use module_symba5
implicit none

! Inputs Only:
integer(ik) :: nbod, irec
integer(ik) :: ielev(nbod)
integer(ik) :: ielst(NENMAX,2), ielc
real(rk) :: mass(nbod), dt, rhill(nbod), sgn
real(rk) :: xh(nbod), yh(nbod), zh(nbod)

! Inputs and Outputs:
real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

! Internals:
integer(ik) :: i, j, irm1, irecl, ie
real(rk) :: dx, dy, dz, fac, ris, r
real(rk) :: ri, rr, r2, faci, facj, ir3, rim1

!----
! Executable code 
!-----------------

irm1 = irec - 1
if(sgn < 0.0_rk) then

  irecl = irec - 1

else

  irecl = irec

end if

! calculate the accelerations
do ie = 1, ielc

  i = ielst(ie,1)
  j = ielst(ie,2)

  if((ielev(i) >= irm1) .and. (ielev(j) >= irm1)) then

    ri = ((rhill(i) + rhill(j))*RHSCALE*RSHELL**irecl)**2
    rim1 = ri*RSHELL**2

    dx = xh(j) - xh(i)
    dy = yh(j) - yh(i)
    dz = zh(j) - zh(i)
    r2 = dx**2 + dy**2 + dz**2
    ir3 = 1.0_rk/(r2*sqrt(r2))

    if(r2 < rim1) then

      fac = 0.0_rk

    else if(r2 < ri) then

      ris = sqrt(ri)
      r = sqrt(r2)
      rr = (ris - r)/(ris*(1.0_rk - RSHELL))
      fac = (1.0_rk - 3.0_rk*rr**2 + 2.0_rk*rr**3)*ir3

    else

      fac = ir3

    end if

    ! apply the kick for body i
    facj = sgn*mass(j)*fac*dt
    vxb(i) = vxb(i) + facj*dx
    vyb(i) = vyb(i) + facj*dy
    vzb(i) = vzb(i) + facj*dz

    ! apply the kick for body j
    faci = sgn*mass(i)*fac*dt
    vxb(j) = vxb(j) - faci*dx
    vyb(j) = vyb(j) - faci*dy
    vzb(j) = vzb(j) - faci*dz

  end if

end do

return
end subroutine symba5_kick
