module anal
! Module for routines in the anal directory
use swift
use util
use coord
use mvs
use obl
use io
implicit none

contains
!!
!!
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
!!
  subroutine anal_energy(nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, ke, pe, energy, l)
  !-------------------------------------------------------------------------
  !                          ANAL_ENERGY.F
  !-------------------------------------------------------------------------
  ! Calculates the energy of the total system (massive bodies) wrt time.
  ! returns the total energy of n objects by direct pairwise summation
  ! G = 1., and we allow diff. masses.  Also returns square of total ang. mom.
  !
  !      Input:
  !            t             ==>  current time
  !            nbod          ==>  number of massive bodies (int scalar)
  !            mass          ==>  mass of bodies (real array)
  !            j2rp2         ==>  scaled value of j2 moment (real(rk) :: scalar)
  !            j4rp4         ==>  scaled value of j4 moment (real(rk) :: scalar)
  !            xh,yh,zh      ==>  current position in heliocentric coord
  !                               (real arrays)
  !            vxh,vyh,vzh   ==>  current velocity in heliocentric coord
  !                               (real arrays)
  !
  !      Output:
  !            ke            ==>  kinetic energy
  !            pe            ==>  potential energy
  !            energy        ==>  Total energy
  !            l             ==>  components of total angular momentum
  !                               (real array)
  !
  ! Remarks:
  ! Authors:  Martin Duncan
  ! Date:  ?
  ! Last revision:  1/24/97 HFL
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm
  real(rk) :: mass(nbod), j2rp2, j4rp4
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

  ! Output variables
  real(rk) :: energy, l(NDIM), ke, pe

  ! Internal variables
  integer(ik) :: i, j
  real(rk) :: dx, dy, dz, dr2, xbi, ybi, zbi, massi, oblpe, msys, lx, ly, lz
  !real(rk) :: kerr, perr, lxerr, lyerr, lzerr
  real(rk) :: irh(nbod), ir3h(nbod)
  real(rk) :: xb(nbod), yb(nbod), zb(nbod)
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Convert from heliocentric to barycentric coordinates
  call coord_h2b(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xb, yb, zb, vxb, vyb, vzb, msys)

  ! Initialize the kinetic and potential energy, along with the components of angular momentum
  ke = 0.5_rk*mass(nbod)*(vxb(nbod)**2 + vyb(nbod)**2 + vzb(nbod)**2)
  pe = 0.0_rk
  lx = mass(nbod)*(yb(nbod)*vzb(nbod) - zb(nbod)*vyb(nbod))
  ly = mass(nbod)*(zb(nbod)*vxb(nbod) - xb(nbod)*vzb(nbod))
  lz = mass(nbod)*(xb(nbod)*vyb(nbod) - yb(nbod)*vxb(nbod))

  ! Initialize the correction factors for the Kahan summation formula
  !kerr = 0.0_rk; perr = 0.0_rk
  !lxerr = 0.0_rk; lyerr = 0.0_rk; lzerr = 0.0_rk

  do i = 1, nbodm

    ! Components of the angular momentum
    lx = lx + mass(i)*(yb(i)*vzb(i) - zb(i)*vyb(i))
    ly = ly + mass(i)*(zb(i)*vxb(i) - xb(i)*vzb(i))
    lz = lz + mass(i)*(xb(i)*vyb(i) - yb(i)*vxb(i))
    !lx = util_kahan_sum(lx, mass(i)*(yb(i)*vzb(i) - zb(i)*vyb(i)), lxerr)
    !ly = util_kahan_sum(ly, mass(i)*(zb(i)*vxb(i) - xb(i)*vzb(i)), lyerr)
    !lz = util_kahan_sum(lz, mass(i)*(xb(i)*vyb(i) - yb(i)*vxb(i)), lzerr)

    ! Kinetic energy
    ke = ke + 0.5_rk*mass(i)*(vxb(i)**2 + vyb(i)**2 + vzb(i)**2)
    !ke = util_kahan_sum(ke, 0.5_rk*mass(i)*(vxb(i)**2 + vyb(i)**2 + vzb(i)**2), kerr)

    ! Extract the information for body i, and make a copy available for each thread
    massi = mass(i)
    xbi = xb(i)
    ybi = yb(i)
    zbi = zb(i)

    ! Potential energy
    !OMP FIRSTPRIVATE(i, massi, xbi, ybi, zbi, perr) LASTPRIVATE(perr)
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(j, dx, dy, dz, dr2) &
    !$OMP FIRSTPRIVATE(i, massi, xbi, ybi, zbi) SHARED(nbod, mass, xb, yb, zb) &
    !$OMP REDUCTION(+ : pe)
    do j = i + 1, nbod

      if((massi /= 0.0_rk) .and. (mass(j) /= 0.0_rk)) then

        dx = xb(j) - xbi
        dy = yb(j) - ybi
        dz = zb(j) - zbi
        dr2 = dx**2 + dy**2 + dz**2
        pe = pe - (massi*mass(j))/sqrt(dr2)
        !pe = util_kahan_sum(pe, -(massi*mass(j))/sqrt(dr2), perr)

      end if

    end do
    !$OMP END PARALLEL DO

  end do

  !FIRSTPRIVATE(lxerr, lyerr, lzerr, kerr)
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) &
  !$OMP SHARED(nbod, nbodm, mass, xb, yb, zb, vxb, vyb, vzb) REDUCTION(+ : lx, ly, lz, ke)
  do i = nbodm + 1, nbod - 1

    ! Components of the angular momentum
    lx = lx + mass(i)*(yb(i)*vzb(i) - zb(i)*vyb(i))
    ly = ly + mass(i)*(zb(i)*vxb(i) - xb(i)*vzb(i))
    lz = lz + mass(i)*(xb(i)*vyb(i) - yb(i)*vxb(i))
    !lx = util_kahan_sum(lx, mass(i)*(yb(i)*vzb(i) - zb(i)*vyb(i)), lxerr)
    !ly = util_kahan_sum(ly, mass(i)*(zb(i)*vxb(i) - xb(i)*vzb(i)), lyerr)
    !lz = util_kahan_sum(lz, mass(i)*(xb(i)*vyb(i) - yb(i)*vxb(i)), lzerr)

    ! Kinetic energy
    ke = ke + 0.5_rk*mass(i)*(vxb(i)**2 + vyb(i)**2 + vzb(i)**2)
    !ke = util_kahan_sum(ke, 0.5_rk*mass(i)*(vxb(i)**2 + vyb(i)**2 + vzb(i)**2), kerr)

  end do
  !$OMP END PARALLEL DO

  ! Store components of the angular momentum in a vector
  l = (/ lx, ly, lz /)

  ! If oblateness terms are present, then compute contribution to potential
  ! --- Put into a separate subroutine to avoid unnecessarily creating irh and ir3h vectors if oblateness terms are absent ---
  if(j2rp2 /= 0.0_rk) then

    call getacch_ir3(nbod, 2, xh, yh, zh, ir3h, irh)
    call obl_pot(nbod, mass, j2rp2, j4rp4, xh, yh, zh, irh, oblpe)
    pe = pe + oblpe

  end if

  ! Total energy of the system
  energy = ke + pe

  return
  end subroutine anal_energy
!!
  subroutine anal_energy_write(t, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, iu, fopenstat, eoff)
  !-------------------------------------------------------------------------
  !                          ANAL_ENERGY_WRITE.F
  !-------------------------------------------------------------------------
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
!!
!!
end module anal