subroutine symba5_helio_getacch(iflg, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh)
!*************************************************************************
!                        SYMBA5_HELIO_GETACCH.F
!*************************************************************************
! This subroutine calculates the acceleration on the massive particles
! in the HELIOCENTRIC frame.
!             Input:
!                 iflg        ==>  =0 calculate forces (int scalor)
!                                  =1 don't
!                 nbod        ==>  number of massive bodies (int scalor)
!                 nbodm       ==>  The last massive particle
!                                  (int scalor)
!                 mass        ==>  mass of bodies (real array)
!                 j2rp2, j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 xh, yh, zh    ==>  position in heliocentric coord (real arrays)
!             Output:
!                 axh, ayh, azh ==>  acceleration in helio coord (real arrays)
!
! Remarks Based on helio_getacch.f
! Author:  Hal Levison
! Date:    9/12/99
! Last revision:
use module_swift
use module_interfaces, except_this_one => symba5_helio_getacch
implicit none

! Inputs:
integer(ik) :: nbod, nbodm, iflg
real(rk) :: mass(nbod), j2rp2, j4rp4
real(rk) :: xh(nbod), yh(nbod), zh(nbod)

! Outputs:
real(rk) :: axh(nbod), ayh(nbod), azh(nbod)

! Internals:
integer(ik) :: i, j
real(rk) :: aoblx(nbod), aobly(nbod), aoblz(nbod), aoblx0, aobly0, aoblz0
real(rk), save :: axhl(NTPMAX), ayhl(NTPMAX), azhl(NTPMAX)
real(rk) :: ir3h(nbod), irh(nbod)
real(rk) :: dx, dy, dz, dr2, idr32, faci, facj, xhi, yhi, zhi, massi
real(rk) :: daxj, dayj, dazj
!real(rk) :: daxerr, dayerr, dazerr

!----
! Executable code
!-----------------

if(iflg == 0) then

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, axhl, ayhl, azhl)
  do i = 1, nbod

    axhl(i) = 0.0_rk
    ayhl(i) = 0.0_rk
    azhl(i) = 0.0_rk

  end do
  !$OMP END PARALLEL DO

  ! Compute the acceleration between gravitating bodies, the acceleration on non-gravitating
  ! bodies due to the gravitating bodies, and their back reaction on the gravitating bodies.
  ! No mutual gravitational interaction is computed for bodies nbodm + 1 to nbod.
  !
  ! N.B. The parallel approach outlined below is not reccomended for a large number of
  ! gravitating bodies (i.e. nbodm >~ 1000), otherwise the overhead for creating a destroying
  ! the parallel region in the inner loop will start to become a factor.
  do i = 2, nbodm

    ! Initialize the sum of the acceleration from bodies j
    daxj = 0.0_rk !; daxerr = 0.0_rk
    dayj = 0.0_rk !; dayerr = 0.0_rk
    dazj = 0.0_rk !; dazerr = 0.0_rk

    ! Extract a copy of the information for particle i, and make a copy available for each thread
    xhi = xh(i); yhi = yh(i); zhi = zh(i)
    massi = mass(i)

    !FIRSTPRIVATE(i, xhi, yhi, zhi, massi, daxerr, dayerr, dazerr)
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(i, xhi, yhi, zhi, massi) &
    !$OMP PRIVATE(j, dx, dy, dz, dr2, idr32, faci, facj) SHARED(nbod, mass, xh, yh, zh, axhl, ayhl, azhl) &
    !$OMP REDUCTION(+ : daxj, dayj, dazj)
    do j = i + 1, nbod

      dx = xh(j) - xhi
      dy = yh(j) - yhi
      dz = zh(j) - zhi
      dr2 = dx**2 + dy**2 + dz**2

      idr32 = 1.0_rk/(dr2*sqrt(dr2))
      faci = massi*idr32
      facj = mass(j)*idr32

      ! Sum the acceleration for gravitating body i due to body j (gravitating/non-gravitating)
      daxj = daxj + facj*dx
      dayj = dayj + facj*dy
      dazj = dazj + facj*dz
      !daxj = util_kahan_sum(daxj, facj*dx, daxerr)
      !dayj = util_kahan_sum(dayj, facj*dy, dayerr)
      !dazj = util_kahan_sum(dazj, facj*dz, dazerr)

      ! Update the acceleration for body j (gravitating/non-gravitating) due to gravitating body i
      axhl(j) = axhl(j) - faci*dx
      ayhl(j) = ayhl(j) - faci*dy
      azhl(j) = azhl(j) - faci*dz

    end do
    !$OMP END PARALLEL DO

    ! Update the acceleration for gravitating body i due to bodies j (gravitating/non-gravitating)
    axhl(i) = axhl(i) + daxj
    ayhl(i) = ayhl(i) + dayj
    azhl(i) = azhl(i) + dazj

  end do

end if

! Now do j2 and j4 stuff
if(j2rp2 /= 0.0_rk) then

  call getacch_ir3(nbod, 2, xh, yh, zh, ir3h, irh)
  call obl_acc(nbod, mass, j2rp2, j4rp4, xh, yh, zh, irh, aoblx, aobly, aoblz)

  aoblx0 = aoblx(1); aobly0 = aobly(1); aoblz0 = aoblz(1)

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(aoblx0, aobly0, aoblz0) &
  !$OMP SHARED(nbod, mass, axh, ayh, azh, axhl, ayhl, azhl, aoblx, aobly, aoblz)
  do i = 1, nbod

    axh(i) = axhl(i) + aoblx(i) - aoblx0
    ayh(i) = ayhl(i) + aobly(i) - aobly0
    azh(i) = azhl(i) + aoblz(i) - aoblz0

  end do
  !$OMP END PARALLEL DO

else

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, axh, ayh, azh, axhl, ayhl, azhl)
  do i = 1, nbod

    axh(i) = axhl(i)
    ayh(i) = ayhl(i)
    azh(i) = azhl(i)

  end do
  !$OMP END PARALLEL DO

end if

return
end subroutine symba5_helio_getacch
