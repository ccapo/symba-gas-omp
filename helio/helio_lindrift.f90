subroutine helio_lindrift(nbod, mass, vxb, vyb, vzb, dt, xh, yh, zh, px, py, pz)
!*************************************************************************
!                            HELIO_LINDRIFT.F
!*************************************************************************
! This subroutine takes a linear drift due to mometum of Sun
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 vxb, vyb, vzb ==>  velocity in bary coord
!                                    (real arrays)
!                 dt            ==>  time step
!                 xh, yh, zh    ==>  initial position in helio coord
!                                    (real arrays)
!             Output:
!                 xh, yh, zh    ==>  final position in helio coord
!                                    (real arrays)
!                 px, py, pz    ==>  momentum of sun: tp's need this
!                                    (real scalars)
!
! Remarks: Bases on Martin's code h2.f
! Authors:  Hal Levison
! Date:    11/14/96
! Last revision: 1/8/97
use module_swift
use module_interfaces, except_this_one => helio_lindrift
implicit none

! Inputs Only:
integer(ik) :: nbod
real(rk) :: mass(nbod), dt
real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

! Inputs and Outputs:
real(rk) :: xh(nbod), yh(nbod), zh(nbod)

! Outputs Only:
real(rk) :: px, py, pz

! Internals:
integer(ik) :: i
real(rk) :: mstar
!real(rk) :: pxerr, pyerr, pzerr
real(rk) :: pxtmp, pytmp, pztmp

!----
! Executable code

! Initialize the temporary variables and their error variables
mstar = mass(1)
px = 0.0_rk !; pxerr = 0.0_rk
py = 0.0_rk !; pyerr = 0.0_rk
pz = 0.0_rk !; pzerr = 0.0_rk

!FIRSTPRIVATE(mu, dt, pxerr, pyerr, pzerr)
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, pxtmp, pytmp, pztmp) FIRSTPRIVATE(mstar, dt) &
!$OMP SHARED(nbod, mass, xh, yh, zh, vxb, vyb, vzb, px, py, pz)

!$OMP DO SCHEDULE(STATIC) REDUCTION(+ : px, py, pz)
do i = 2, nbod

  px = px + mass(i)*vxb(i)
  py = py + mass(i)*vyb(i)
  pz = pz + mass(i)*vzb(i)
  !px = util_kahan_sum(px, mass(i)*vxb(i), pxerr)
  !py = util_kahan_sum(py, mass(i)*vyb(i), pyerr)
  !pz = util_kahan_sum(pz, mass(i)*vzb(i), pzerr)

end do

! Make a private copy of the momentum vector
pxtmp = px/mstar
pytmp = py/mstar
pztmp = pz/mstar

!$OMP DO SCHEDULE(STATIC)
do i = 2, nbod

  if(mass(i) /= 0.0_rk) then

    xh(i) = xh(i) + pxtmp*dt
    yh(i) = yh(i) + pytmp*dt
    zh(i) = zh(i) + pztmp*dt

  end if

end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

return
end subroutine helio_lindrift
