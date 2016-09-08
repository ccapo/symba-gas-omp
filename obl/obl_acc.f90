subroutine obl_acc(nbod, mass, j2rp2, j4rp4, xh, yh, zh, irh, aoblx, aobly, aoblz)
!***************************************************************************
!			OBL_ACC.F
!*************************************************************************
! OBL_ACC returns the BARYCENTRIC x, y, z components of the acc. on NBOD
! particles due to the oblateness of mass(1) using
! the values of J2RP2 and J4RP4 passed into the routine.
! (J2RP2 for example is the product of
! J_2 times the square of the central body's radius)
! Here we return the net acc. produced
! only by the J2 and J4 terms (i.e. including
! neither the monopole nor higher order terms).
!
!             Input:
!                 nbod     ==>  number of massive bodies (incl. central one)
!                 mass(*)  ==>  masses of particles (real(rk) :: array)
!                 j2rp2    ==>  scaled value of j2 moment (real(rk) :: scalar)
!                 j4rp4    ==>  scaled value of j4 moment (real(rk) :: scalar)
!                                    (real(rk) :: vectors)
!                 xh(*), yh(*), zh(*)   ==>  HELIO. positions of particles
!                 irh(*)   ==> 1./ magnitude of radius vector (real(rk) :: vector)
!                                (passed in to save calcs.)
!             Output:
!               aoblx(*), aobly(*), aoblz(*)  ==>  BARY. components of accel
!                                        (real(rk) :: vectors)
!
! Remarks:  aoblx(1) (for example) contains x-component of
!           bary. acc. of central body
! Authors:  Martin Duncan
! Date:    3/4/94
! Last revision:
use module_swift
use module_interfaces, except_this_one => obl_acc
implicit none

! Inputs Only:
integer(ik) :: nbod
real(rk) :: j2rp2, j4rp4
real(rk) :: mass(nbod)
real(rk) :: xh(nbod), yh(nbod), zh(nbod), irh(nbod)

! Output
real(rk) :: aoblx(nbod), aobly(nbod), aoblz(nbod)

! Internals
integer(ik) :: i
real(rk) :: rinv2, t0, t1, t2, t3, mstar
real(rk) :: fac1, fac2, aoblx0, aobly0, aoblz0
!real(rk) :: axerr, ayerr, azerr

!----
! Executable code

! Initialize the bary. acc. of the Sun
aoblx0 = 0.0_rk !; axerr = 0.0_rk
aobly0 = 0.0_rk !; ayerr = 0.0_rk
aoblz0 = 0.0_rk !; azerr = 0.0_rk
mstar = mass(1)

!FIRSTPRIVATE(mu, imu, j2rp2, j4rp4, axerr, ayerr, azerr)
! First get the barycentric acceleration of each body due to oblate central body
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(nbod, mass, xh, yh, zh, irh, aoblx, aobly, aoblz) &
!$OMP FIRSTPRIVATE(mstar, j2rp2, j4rp4) PRIVATE(i, rinv2, t0, t1, t2, t3, fac1, fac2) &
!$OMP REDUCTION(+ : aoblx0, aobly0, aoblz0)
do i = 2, nbod

  ! Note that here we assume we know inverse of radius rather than calc. it
  ! from (x, y, z) to save the sqrt.
  rinv2 = irh(i)**2
  t0 = -mstar*irh(i)*rinv2**2
  t1 = 1.5_rk*j2rp2
  t2 = rinv2*zh(i)**2
  t3 = 1.875_rk*j4rp4*rinv2

  fac1 = t0*(t1 - t3 - (5.0_rk*t1 - (14.0_rk - 21.0_rk*t2)*t3)*t2)
  fac2 = 2.0_rk*t0*(t1 - (2.0_rk - (14.0_rk*t2/3.0_rk))*t3)

  aoblx(i) = fac1*xh(i)
  aobly(i) = fac1*yh(i)
  aoblz(i) = (fac1 + fac2)*zh(i)

  ! Compute the bary. acc. of Sun due to all the planets
  aoblx0 = aoblx0 + mass(i)*aoblx(i)/mstar
  aobly0 = aobly0 + mass(i)*aobly(i)/mstar
  aoblz0 = aoblz0 + mass(i)*aoblz(i)/mstar
  !aoblx0 = util_kahan_sum(aoblx0, mass(i)*aoblx(i)/mstar, axerr)
  !aobly0 = util_kahan_sum(aobly0, mass(i)*aobly(i)/mstar, ayerr)
  !aoblz0 = util_kahan_sum(aoblz0, mass(i)*aoblz(i)/mstar, azerr)

end do
!$OMP END PARALLEL DO

! Store the bary. acc. of Sun due to all the planets
aoblx(1) = -aoblx0
aobly(1) = -aobly0
aoblz(1) = -aoblz0

return
end subroutine obl_acc
