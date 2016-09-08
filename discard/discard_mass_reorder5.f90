subroutine discard_mass_reorder5(ip, nbod, mass, xh, yh, zh, vxh, vyh, vzh, rpl, rhill, isperih, ibound)
!*************************************************************************
!                            DISCARD_MASS_REORDER5.F
!*************************************************************************
! Remove a massive body
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 ip            ==>  planets to remove (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh, yh, zh      ==>   position in helio coord
!                                    (real arrays)
!                 vxh, vyh, vzh   ==>   pl vel in helio coord
!                                    (real arrays)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!                 rhill         ==>  size of a planet's hill's sphere.
!                                    (real array)
!                 isperih       ==> heliocentric peri flags. (real array)
!             Output:
!                 ip            ==>  planets to remove (int scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh, yh, zh      ==>   position in helio coord
!                                    (real arrays)
!                 vxh, vyh, vzh   ==>   pl vel in helio coord
!                                    (real arrays)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!                 rhill         ==>  size of a planet's hill's sphere.
!                                    (real array)
!                 isperih       ==> heliocentric peri flags. (real array)
!
! Remarks:
!
! Authors:  Hal Levison
! Date:    1/2/97
! Last revision: 5/13/99
use module_swift
implicit none

! Inputs:
integer(ik) :: ip

! Input and Output
integer(ik) :: nbod
integer(ik) :: isperih(nbod), ibound(nbod)
real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod), rpl(nbod)
real(rk) :: rhill(nbod)
real(rk) :: rdrag(NTPMAX), dragc(NTPMAX)

common / sizedist / rdrag, dragc

! internal
integer(ik) :: i

!-----
! Executable code

! Move all the values from bodies [ip + 1, nbod] up one index
do i = ip, nbod - 1

  xh(i) = xh(i + 1)
  yh(i) = yh(i + 1)
  zh(i) = zh(i + 1)
  vxh(i) = vxh(i + 1)
  vyh(i) = vyh(i + 1)
  vzh(i) = vzh(i + 1)
  mass(i) = mass(i + 1)
  rpl(i) = rpl(i + 1)
  rhill(i) = rhill(i + 1)
  isperih(i) = isperih(i + 1)
  ! *** 04/24/08 -- CCC ***
  rdrag(i) = rdrag(i + 1)
  dragc(i) = dragc(i + 1)
  ibound(i) = ibound(i + 1)

end do

! Decrease the nummber of bodies by one
nbod = nbod - 1

return
end subroutine discard_mass_reorder5
