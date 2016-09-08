subroutine anal_energy_write(t, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, iu, fopenstat, eoff)
!*************************************************************************
!                          ANAL_ENERGY_WRITE.F
!*************************************************************************
! Writes the energy of the total system (massive bodies) wrt time.
!
!      Input:
!            t             ==>  current time
!            nbod          ==>  number of massive bodies (int scalar)
!            mass          ==>  mass of bodies (real array)
!            j2rp2         ==>  scaled value of j2 moment (real(rk) :: scalar)
!            j4rp4         ==>  scaled value of j4 moment (real(rk) :: scalar)
!            xh,yh,zh      ==>  current position in helio coord 
!                               (real arrays)
!            vxh,vyh,vzh   ==>  current velocity in helio coord 
!                               (real arrays)
!            iu            ==>  unit to write to (int scalar)
!            fopenstat     ==>  The status flag for the open
!                                statements of the output files.
!                                      (character*80)
!            eoff          ==> An energy offset that is added to the energy
!                                      (real(rk) :: scalar)
!
! Remarks: 
! Authors:  Hal Levison
! Date:    3/4/93
! Last revision: 12/27/96
use module_swift
use module_interfaces, except_this_one => anal_energy_write
implicit none

! Inputs:
integer(ik) :: nbod, nbodm, iu
real(rk) :: mass(nbod), t, j2rp2, j4rp4, eoff
real(rk) :: xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)
character(len = 24) :: fopenstat

! Internals
integer(ik), save :: i1st = 0
real(rk) :: energy, eltot(NDIM), ke, pot

!----
! Executable code

! Compute and print initial ke, pot, energy and ang. mom.
call anal_energy_discard5(0, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, ke, pot, energy, eltot)

energy = energy + eoff

call io_energy_write(i1st, t, energy, eltot, iu, fopenstat)

if(i1st == 0) i1st = 1

return
end subroutine anal_energy_write
