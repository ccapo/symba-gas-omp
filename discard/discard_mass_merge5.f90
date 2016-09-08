subroutine discard_mass_merge5(time, nbod, nbodm, ip1, ip2, mass, xh, yh, zh, vxh, vyh, vzh, rpl, eoff, ielc, ielst)
!*************************************************************************
!                            DISCARD_MASS_MERGE5.F
!*************************************************************************
! Merge two massive bodies
!
!             Input:
!                 time          ==>  current time (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 ip1, ip2       ==>  planets to merge (real scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh, yh, zh      ==>   position in helio coord
!                                    (real arrays)
!                 vxh, vyh, vzh   ==>   pl vel in helio coord
!                                    (real arrays)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!                 eoff          ==> Amount of energy lost due to discards
!                                          (real scalar)
!                ielc           ==>  number of encounters (integer*2 scalar)
!                ielst          ==>  list of ecnounters (2D integer*2 array)
!             Output:
!                 mass          ==>  recalculated mass of bodies (real array)
!                 xh, yh, zh      ==>  recalculated position in helio coord
!                                    (real arrays)
!                 vxh, vyh, vzh   ==>  recalculated pl vel in helio coord
!                                    (real arrays)
!                 rpl           ==>  recalculated physical sizes of a planet.
!                                    (real array)
!                 eoff          ==> Updated amount of energy lost from discards
!                                          (real scalar)
!                ielc           ==>  number of encounters (integer*2 scalar)
!                ielst          ==>  list of ecnounters (2D integer*2 array)
!
! Remarks:
!
! Authors:  Hal Levison
! Date:    12/30/96
! Last revision: 1/30/97
use module_swift
use module_symba5
use module_interfaces, except_this_one => discard_mass_merge5
implicit none

! Inputs:
integer(ik) :: ip1, ip2
real(rk) :: time

! Input and Output
integer(ik) :: nbod, nbodm
integer(ik) :: ielc, ielst(NENMAX,2)
real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod), rpl(nbod)
real(rk) :: eoff

! internal
integer(ik) :: i, j, itmp
real(rk) :: mtot, m1, r1
real(rk) :: x1, y1, z1
real(rk) :: vx1, vy1, vz1
real(rk) :: m2, r2
real(rk) :: x2, y2, z2
real(rk) :: vx2, vy2, vz2
real(rk) :: j2rp2, j4rp4, ke, pot, energy1, energy2, eltot(NDIM)

!-----------------!
! Executable code !
!-----------------!

j2rp2 = 0.0_rk
j4rp4 = 0.0_rk
call anal_energy_discard5(0, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, ke, pot, energy1, eltot)

!if(mass(ip2) > mass(ip1)) then
!
!  itmp = ip1
!  ip1 = ip2
!  ip2 = itmp
!
!end if

write(*,'(a20,i6,a5,i6,a7,1pe12.5)') 'Merging particles: ', ip1,  ' and ',  ip2, ' at t =', time

x1 = xh(ip1)
y1 = yh(ip1)
z1 = zh(ip1)
vx1 = vxh(ip1)
vy1 = vyh(ip1)
vz1 = vzh(ip1)
m1 = mass(ip1)
r1 = rpl(ip1)

x2 = xh(ip2)
y2 = yh(ip2)
z2 = zh(ip2)
vx2 = vxh(ip2)
vy2 = vyh(ip2)
vz2 = vzh(ip2)
m2 = mass(ip2)
r2 = rpl(ip2)

! Note:  I am just putting these guys together here,  which is
!        clearly wrong.  I should integrate back to the time
!        of close approach.
mtot = m1 + m2
rpl(ip1) = (r1**3 + r2**3)**(1.0_rk/3.0_rk)
vxh(ip1) = (m1*vx1 + m2*vx2)/mtot
vyh(ip1) = (m1*vy1 + m2*vy2)/mtot
vzh(ip1) = (m1*vz1 + m2*vz2)/mtot
mass(ip1) = mtot

! Put in zeros for the rest the second particle
xh(ip2) = x2*1.0e10_rk   ! so danby does not fail
yh(ip2) = y2*1.0e10_rk
zh(ip2) = z2*1.0e10_rk
vxh(ip2) = 0.0_rk
vyh(ip2) = 0.0_rk
vzh(ip2) = 0.0_rk
mass(ip2) = 0.0_rk
rpl(ip2) = 0.0_rk

! Remove any encounters with ip2
j = 1
do while(j <= ielc)

  if((ielst(j,1) == ip2) .or. (ielst(j,2) == ip2)) then

    do i = j + 1, ielc

      ielst(i - 1,1) = ielst(i,1)
      ielst(i - 1,2) = ielst(i,2)

    end do

    ielc = ielc - 1

  else

    j = j + 1

  end if

end do

call anal_energy_discard5(0, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, ke, pot, energy2, eltot)
eoff = eoff + energy1 - energy2

call io_discard_merge(time, ip1, ip2, m1, r1, x1, y1, z1, vx1, vy1, vz1, m2, r2, x2, y2, z2, vx2, vy2, vz2, &
                      mass(ip1), rpl(ip1), xh(ip1), yh(ip1), zh(ip1), vxh(ip1), vyh(ip1), vzh(ip1))

return
end subroutine discard_mass_merge5
