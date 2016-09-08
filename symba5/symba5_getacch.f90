subroutine symba5_getacch(nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, axh, ayh, azh, mtiny, ielc, ielst)
!*************************************************************************
!                        SYMBA5_GETACCH.F
!*************************************************************************
! This subroutine calculates the acceleration on the massive particles
! in the HELIOCENTRIC frame.
!             Input:
!                 nbod        ==>  number of massive bodies (int scalor)
!                 nbodm       ==>  Location of last massive body(int scalar)
!                 mass        ==>  mass of bodies (real array)
!                 j2rp2, j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 xh, yh, zh    ==>  position in heliocentric coord (real arrays)
!                 mtiny       ==>  Small mass  (real array)
!                ielc           ==>  number of encounters (integer*2 scalar)
!                ielst          ==>  list of ecnounters (2D integer*2 array)
!             Output:
!                 axh, ayh, azh ==>  acceleration in helio coord (real arrays)
!
! Remarks: Based on helio_getacch.f,  but does not include the forces of
!          an body B on body A,  if body B and A are having an encounter.
! Author:  Hal Levison
! Date:    3/20/97
! Last revision: 11/22/97
use module_swift
use module_symba5
use module_interfaces, except_this_one => symba5_getacch
implicit none

! Inputs:
integer(ik) :: nbod, nbodm
integer(ik) :: ielst(NENMAX,2), ielc
real(rk) :: mass(nbod), j2rp2, j4rp4, mtiny
real(rk) :: xh(nbod), yh(nbod), zh(nbod)

! Outputs:
real(rk) :: axh(nbod), ayh(nbod), azh(nbod)

! Internals:
integer(ik) :: i, j, k
real(rk) :: aoblx(nbod), aobly(nbod), aoblz(nbod), aoblx0, aobly0, aoblz0
real(rk) :: ir3h(nbod), irh(nbod)
real(rk) :: dx, dy, dz, dr2, idr32, faci, facj, xhi, yhi, zhi, massi
real(rk) :: daxj, dayj, dazj
!real(rk) :: daxerr, dayerr, dazerr

!-----------------
! Executable code
!-----------------

! Zero things
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, axh, ayh, azh)
do i = 1, nbod

  axh(i) = 0.0_rk
  ayh(i) = 0.0_rk
  azh(i) = 0.0_rk

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
  xhi = xh(i); yhi = yh(i); xhi = xh(i)
  massi = mass(i)

  !FIRSTPRIVATE(i, xhi, yhi, zhi, massi, daxerr, dayerr, dazerr)
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(i, xhi, yhi, zhi, massi) &
  !$OMP PRIVATE(j, dx, dy, dz, dr2, idr32, faci, facj) SHARED(nbod, mass, xh, yh, zh, axh, ayh, azh) &
  !$OMP REDUCTION(+ : daxj, dayj, dazj)
  do j = i + 1, nbod

    dx = xh(j) - xh(i)
    dy = yh(j) - yh(i)
    dz = zh(j) - zh(i)
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
    axh(j) = axh(j) - faci*dx
    ayh(j) = ayh(j) - faci*dy
    azh(j) = azh(j) - faci*dz

  end do
  !$OMP END PARALLEL DO

  ! Update the acceleration for gravitating body i due to bodies j (gravitating/non-gravitating)
  axh(i) = axh(i) + daxj
  ayh(i) = ayh(i) + dayj
  azh(i) = azh(i) + dazj

end do

! Now subtract off anyone in an encounter
do k = 1, ielc

  i = ielst(k,1)
  j = ielst(k,2)

  dx = xh(j) - xh(i)
  dy = yh(j) - yh(i)
  dz = zh(j) - zh(i)
  dr2 = dx**2 + dy**2 + dz**2

  idr32 = 1.0_rk/(dr2*sqrt(dr2))
  faci = mass(i)*idr32
  facj = mass(j)*idr32

  axh(j) = axh(j) + faci*dx
  ayh(j) = ayh(j) + faci*dy
  azh(j) = azh(j) + faci*dz

  axh(i) = axh(i) - facj*dx
  ayh(i) = ayh(i) - facj*dy
  azh(i) = azh(i) - facj*dz

end do

! Now do j2 and j4 stuff
if(j2rp2 /= 0.0_rk) then

  call getacch_ir3(nbod, 2, xh, yh, zh, ir3h, irh)
  call obl_acc(nbod, mass, j2rp2, j4rp4, xh, yh, zh, irh, aoblx, aobly, aoblz)

  aoblx0 = aoblx(1); aobly0 = aobly(1); aoblz0 = aoblz(1)

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(aoblx0, aobly0, aoblz0) &
  !$OMP SHARED(nbod, mass, axh, ayh, azh, aoblx, aobly, aoblz)
  do i = 2, nbod

    if(mass(i) /= 0.0_rk) then

      axh(i) = axh(i) + aoblx(i) - aoblx0
      ayh(i) = ayh(i) + aobly(i) - aobly0
      azh(i) = azh(i) + aoblz(i) - aoblz0

    end if

  end do
  !$OMP END PARALLEL DO

end if

return
end subroutine symba5_getacch
