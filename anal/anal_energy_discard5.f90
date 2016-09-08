subroutine anal_energy_discard5(iflg, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, ke, pot, energy, eltot)
!-------------------------------------------------------------------------
!                          ANAL_ENERGY_DISCARD5.F90
!-------------------------------------------------------------------------
! Calculates the energy of the total system (massive bodies) wrt time.
! returns the total energy of n objects by direct pairwise summation
! G = 1., and we allow diff. masses.  Also returns square of total ang. mom.
!
!      Input:
!            iflg          ==>  use to turn on/off energy calculation
!            t             ==>  current time
!            nbod          ==>  number of massive bodies (int scalar)
!            nbodm         ==>  Location of last massive body(int scalar)
!            mass          ==>  mass of bodies (real array)
!            j2rp2         ==>  scaled value of j2 moment (real*8 scalar)
!            j4rp4         ==>  scaled value of j4 moment (real*8 scalar)
!            xh,yh,zh      ==>  current position in heliocentric coord 
!                               (real arrays)
!            vxh,vyh,vzh   ==>  current velocity in heliocentric coord 
!                               (real arrays)
!
!      Output:
!            ke            ==>  kinetic energy
!            pot           ==>  potential energy
!            energy        ==>  Total energy
!            eltot         ==>  components of total angular momentum
!                              (real array)
!
! Remarks: Based on anal_energy
! Authors:  Hal Levison
! Date:  12/16/06
! Last revision:  
use module_swift
use module_interfaces, except_this_one => anal_energy_discard5
implicit none

! Inputs: 
integer(ik) :: iflg, nbod, nbodm
real(rk) :: mass(nbod), j2rp2, j4rp4
real(rk) :: xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

! Output
real(rk) :: energy, eltot(NDIM), ke, pot

! Internals
logical, save :: leuse = .true.

!-----------------!
! Executable code !
!-----------------!

if(iflg < 0) then

  leuse = .false.
  return               !  <==== NOTE

else if(iflg > 0) then

  leuse = .true.
  return               !  <==== NOTE

end if

! iflg = 0

if(leuse) then

  call anal_energy(nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, ke, pot, energy, eltot)

else

  ke = 0.0_rk
  pot = 0.0_rk
  energy = 0.0_rk
  eltot = 0.0_rk

end if


return	
end subroutine anal_energy_discard5
