subroutine obl_pot(nbod, mass, j2rp2, j4rp4, xh, yh, zh, irh, oblpot)
!***************************************************************************
!			OBL_POT.F
!*************************************************************************
! OBL_POT returns the total potential in the barycentric frame for NBOD
! particles due to the oblateness of mass(1) using
! the values of J2RP2 and J4RP4 passed into the routine.
! (J2RP2 for example is the product of
! J_2 times the square of the central body's radius)
! Here we return the potential produced
! only by the J2 and J4 terms (i.e. including
! neither the monopole nor higher order terms).
!
!             Input:
!                 nbod     ==>  number of massive bodies (incl. central one)
!                 mass(*)  ==>  masses of particles (real(rk) :: array)
!                 j2rp2    ==>  scaled value of j2 moment (real(rk) :: scalar)
!                 j4rp4    ==>  scaled value of j4 moment (real(rk) :: scalar)
!                 xh(*), yh(*), zh(*)   ==>  HELIO. positions of particles
!                                    (real(rk) :: vectors)
!                 irh(*)   ==> 1./ magnitude of radius vector (real(rk) :: vector)
!                                (passed in to save calcs.)
!
!             Output:
!                 oblpot  ==>  BARY. potential
!                                        (real(rk) :: scalar)
!
! Remarks:
! Authors:  Martin Duncan
! Date:    3/4/94
! Last revision:
use module_swift
use module_interfaces, except_this_one => obl_pot
implicit none

! Inputs Only:
integer(ik) :: nbod
real(rk) :: mass(nbod)
real(rk) :: j2rp2, j4rp4
real(rk) :: xh(nbod), yh(nbod), zh(nbod), irh(nbod)

! Output
real(rk) :: oblpot

! Internals
integer(ik) :: i
real(rk) :: rinv2, t0, t1, t2, t3
real(rk) :: p2, p4, mstar
!real(rk) :: oerr

!----
! Executable code

! Sum all the the bary terms for each "planet" due to central oblate "sun"
oblpot = 0.0_rk !; oerr = 0.0_rk
mstar = mass(1)

!FIRSTPRIVATE(mu, j2rp2, j4rp4, oerr)
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(mstar, j2rp2, j4rp4) &
!$OMP PRIVATE(i, rinv2, t0, t1, t2, t3, p2, p4) SHARED(nbod, mass, xh, yh, zh, irh) REDUCTION(+ : oblpot)
do i = 2, nbod

  ! Note that here we assume we know inverse of radius rather than calc. it from (x, y, z) to save the sqrt.
  rinv2 = irh(i)**2
  t0 = mstar*mass(i)*rinv2*irh(i)
  t1 = j2rp2
  t2 = rinv2*zh(i)**2
  t3 = j4rp4*rinv2

  p2 = 0.5_rk*(3.0_rk*t2 - 1.0_rk)
  p4 = 0.125_rk*((35.0_rk*t2 - 30.0_rk)*t2 + 3.0_rk)

  oblpot = oblpot + t0*(t1*p2 + t3*p4)
  !oblpot = util_kahan_sum(oblpot, t0*(t1*p2 + t3*p4), oerr)

end do
!$OMP END PARALLEL DO

return
end subroutine obl_pot
