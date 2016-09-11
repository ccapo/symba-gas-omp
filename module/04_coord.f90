module coord
! Module for the routines contained in coord
use swift
use util
implicit none

contains
!!
!!
  subroutine coord_b2h(nbod, mass, xb, yb, zb, vxb, vyb, vzb,  xh, yh, zh, vxh, vyh, vzh)
  !-----------------------------------------------------------------------
  !	                    COORD_B2H.F
  !-----------------------------------------------------------------------
  !     PURPOSE: Converts from Barycentric to Helio coords.
  !     ARGUMENTS:  Input is
  !                    nbod ==> number of bodies (must be less than NBMAX)
  !                             (integer)
  !	             mass(*) ==>  masses (real array)
  !                                 NOT USED BUT INCLUDED IN ORDER TO HAVE
  !                                 SYMMETRY IN SUBROUTINE CALLS
  !		     xb(*), yb(*), zb(*) ==> Barycentric particle coords
  !                                          (real array)
  !		     vxb(*), vyb(*), vzb(*) ==> Barycentric particle velocities
  !                                             (real array)
  !                 Returned are
  !                    xh(*), yh(*), zh(*) ==> Helio particle positions
  !                                          (real array)
  !                    vxh(*), vyh(*), vzh(*) ==> Helio particle velocities
  !                                            (real array)
  !
  !     ALGORITHM: Obvious
  !     REMARKS:  Can of course use this to get coords. relative to any body.
  !              by changing the one subtracted off.
  !
  !     Authors:  Martin Duncan
  !     WRITTEN:  Jan 27/93
  !     REVISIONS: 2/17/95  HFL
  implicit none

  ! Inputs:
  integer(ik) :: nbod
  real(rk) :: mass(nbod)
  real(rk) :: xb(nbod), yb(nbod), zb(nbod)
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

  ! Outputs:
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

  ! Internals:
  integer(ik) :: i
  real(rk) :: xb0, yb0, zb0, vxb0, vyb0, vzb0

  !----
  ! Executable code

  ! Store the barycentric position and velocity of the central body
  xb0 = xb(1)
  yb0 = yb(1)
  zb0 = zb(1)
  vxb0 = vxb(1)
  vyb0 = vyb(1)
  vzb0 = vzb(1)

  ! Compute the heliocentric position and velocity of all the bodies
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) &
  !$OMP FIRSTPRIVATE(xb0, yb0, zb0, vxb0, vyb0, vzb0) &
  !$OMP SHARED(nbod, xb, yb, zb, vxb, vyb, vzb, xh, yh, zh, vxh, vyh, vzh)
  do i = 1, nbod

    xh(i) = xb(i) - xb0
    yh(i) = yb(i) - yb0
    zh(i) = zb(i) - zb0
    vxh(i) = vxb(i) - vxb0
    vyh(i) = vyb(i) - vyb0
    vzh(i) = vzb(i) - vzb0

  end do
  !$OMP END PARALLEL DO

  return
  end subroutine coord_b2h
!!
  subroutine coord_h2b(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xb, yb, zb, vxb, vyb, vzb, msys)
  !-----------------------------------------------------------------------
  !	                    COORD_H2B.F
  !-----------------------------------------------------------------------
  !     PURPOSE: Converts from Heliocentric to Barycentric coords.
  !     ARGUMENTS:  Input is
  !                    nbod ==> number of bodies (must be less than NBMAX)
  !                             (integer)
  !	             mass(*) ==>  masses (real array)
  !		     xh(*), yh(*), zh(*) ==> heliocentric particle coords
  !                                          (real array)
  !		     vxh(*), vyh(*), vzh(*) ==> heliocentric particle velocities
  !                                             (real array)
  !                 Returned are
  !                    xb(*), yb(*), zb(*) ==> bary. particle positions
  !                                          (real array)
  !                    vxb(*), vyb(*), vzb(*) ==> bary. particle velocities
  !                                            (real array)
  !                    msys              ==>  Total mass of of system
  !                                            (real scalar)
  !     Authors:  Martin Duncan
  !     ALGORITHM: Obvious
  !     WRITTEN:  Jan 27/93
  !     REVISIONS: 2/22/94  HFL
  implicit none

  ! Inputs:
  integer(ik) :: nbod
  real(rk) :: mass(nbod)
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

  ! Outputs:
  real(rk) :: xb(nbod), yb(nbod), zb(nbod)
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

  ! Internals:
  integer(ik) :: i
  real(rk) :: msys, x, y, z, vx, vy, vz
  !real(rk) :: merr, xerr, yerr, zerr, vxerr, vyerr, vzerr
  real(rk) :: xtmp, ytmp, ztmp, vxtmp, vytmp, vztmp

  !----
  ! Executable code

  ! Initialize the total system mass, and the temporary variables and their error variables
  msys = 0.0_rk !; merr = 0.0_rk
  x = 0.0_rk !; xerr = 0.0_rk
  y = 0.0_rk !; yerr = 0.0_rk
  z = 0.0_rk !; zerr = 0.0_rk
  vx = 0.0_rk !; vxerr = 0.0_rk
  vy = 0.0_rk !; vyerr = 0.0_rk
  vz = 0.0_rk !; vzerr = 0.0_rk

  !FIRSTPRIVATE(merr, xerr, yerr, zerr, vxerr, vyerr, vzerr)
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, xtmp, ytmp, ztmp, vxtmp, vytmp, vztmp) &
  !$OMP SHARED(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xb, yb, zb, vxb, vyb, vzb, msys, x, y, z, vx, vy, vz)

  ! Compute the total mass of the system, along with the mass-weighted sum of positions and velocities
  !$OMP DO SCHEDULE(STATIC) REDUCTION(+ : msys, x, y, z, vx, vy, vz)
  do i = 1, nbod

    msys = msys + mass(i)
    x = x + mass(i)*xh(i)
    y = y + mass(i)*yh(i)
    z = z + mass(i)*zh(i)
    vx = vx + mass(i)*vxh(i)
    vy = vy + mass(i)*vyh(i)
    vz = vz + mass(i)*vzh(i)
    !msys = util_kahan_sum(msys, mass(i), merr)
    !x = util_kahan_sum(x, mass(i)*xh(i), xerr)
    !y = util_kahan_sum(y, mass(i)*yh(i), yerr)
    !z = util_kahan_sum(z, mass(i)*zh(i), zerr)
    !vx = util_kahan_sum(vx, mass(i)*vxh(i), vxerr)
    !vy = util_kahan_sum(vy, mass(i)*vyh(i), vyerr)
    !vz = util_kahan_sum(vz, mass(i)*vzh(i), vzerr)

  end do
  !$OMP END DO

  xtmp = x/msys
  ytmp = y/msys
  ztmp = z/msys
  vxtmp = vx/msys
  vytmp = vy/msys
  vztmp = vz/msys

  !$OMP DO SCHEDULE(STATIC)
  do i = 1, nbod

    xb(i) = xh(i) - xtmp
    yb(i) = yh(i) - ytmp
    zb(i) = zh(i) - ztmp
    vxb(i) = vxh(i) - vxtmp
    vyb(i) = vyh(i) - vytmp
    vzb(i) = vzh(i) - vztmp

  end do
  !$OMP END DO NOWAIT

  !$OMP END PARALLEL

  return
  end subroutine coord_h2b
!!
  subroutine coord_vb2h(nbod, mass, vxb, vyb, vzb, vxh, vyh, vzh)
  !-----------------------------------------------------------------------
  !	                    COORD_VB2H.F
  !-----------------------------------------------------------------------
  !     PURPOSE: Converts from Barycentric to Helio coords.
  !               Velocity only
  !     ARGUMENTS:  Input is
  !                    nbod ==> number of bodies (must be less than NBMAX)
  !                             (integer)
  !	             mass(*) ==>  masses (real array)
  !                                 NOT USED BUT INCLUDED IN ORDER TO HAVE
  !                                 SYMMETRY IN SUBROUTINE CALLS
  !		     vxb(*), vyb(*), vzb(*) ==> Barycentric particle velocities
  !                                             (real array)
  !                 Returned are
  !                    vxh(*), vyh(*), vzh(*) ==> Helio particle velocities
  !                                            (real array)
  !
  !     ALGORITHM: Obvious
  !     Authors:  Hal Levison
  !     WRITTEN:  11/14/96
  !     REVISIONS: 11/21/96
  implicit none

  ! Inputs:
  integer(ik) :: nbod
  real(rk) :: mass(nbod)
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

  ! Outputs:
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

  ! Internals:
  integer(ik) :: i
  real(rk) :: vx, vy, vz
  !real(rk) :: vxerr, vyerr, vzerr
  real(rk) :: mstar, vxtmp, vytmp, vztmp

  !-----------------
  ! Executable code
  !-----------------

  ! Initialize temporary variables and their error variables
  mstar = mass(1)
  vx = -mass(2)*vxb(2) !; vxerr = 0.0_rk
  vy = -mass(2)*vyb(2) !; vyerr = 0.0_rk
  vz = -mass(2)*vzb(2) !; vzerr = 0.0_rk

  !FIRSTPRIVATE(mstar, vxerr, vyerr, vzerr)
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, vxtmp, vytmp, vztmp) FIRSTPRIVATE(mstar) &
  !$OMP SHARED(nbod, mass, vxh, vyh, vzh, vxb, vyb, vzb, vx, vy, vz)

  !$OMP DO SCHEDULE(STATIC) REDUCTION(+ : vx, vy, vz)
  do i = 3, nbod

    vx = vx - mass(i)*vxb(i)
    vy = vy - mass(i)*vyb(i)
    vz = vz - mass(i)*vzb(i)
    !vx = util_kahan_sum(vx, -mass(i)*vxb(i), vxerr)
    !vy = util_kahan_sum(vy, -mass(i)*vyb(i), vyerr)
    !vz = util_kahan_sum(vz, -mass(i)*vzb(i), vzerr)

  end do
  !$OMP END DO

  vxtmp = vx/mstar
  vytmp = vy/mstar
  vztmp = vz/mstar

  !$OMP DO SCHEDULE(STATIC)
  do i = 2, nbod

    vxh(i) = vxb(i) - vxtmp
    vyh(i) = vyb(i) - vytmp
    vzh(i) = vzb(i) - vztmp

  end do
  !$OMP END DO NOWAIT

  !$OMP END PARALLEL

  vxh(1) = 0.0_rk
  vyh(1) = 0.0_rk
  vzh(1) = 0.0_rk

  !mstar = mass(1)
  !vxb(1) = vx/mstar
  !vyb(1) = vy/mstar
  !vzb(1) = vz/mstar

  return
  end subroutine coord_vb2h
!!
  subroutine coord_vh2b(nbod, mass, vxh, vyh, vzh, vxb, vyb, vzb, msys)
  !-----------------------------------------------------------------------
  !	                    COORD_VH2B.F
  !-----------------------------------------------------------------------
  !     PURPOSE: Converts from Heliocentric to Barycentric coords.
  !              Velocity only
  !     ARGUMENTS:  Input is
  !                    nbod ==> number of bodies (must be less than NBMAX)
  !                             (integer)
  !	             mass(*) ==>  masses (real array)
  !		     vxh(*), vyh(*), vzh(*) ==> heliocentric particle velocities
  !                                             (real array)
  !                 Returned are
  !                    vxb(*), vyb(*), vzb(*) ==> bary. particle velocities
  !                                            (real array)
  !                    msys              ==>  Total mass of of system
  !                                            (real scalar)
  !     Authors:  Hal Levison
  !     ALGORITHM: Obvious
  !     WRITTEN:  11/14/96
  !     REVISIONS: 11/21/96
  implicit none

  ! Inputs:
  integer(ik) :: nbod
  real(rk) :: mass(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

  ! Outputs:
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

  ! Internals:
  integer(ik) :: i
  real(rk) :: msys, vx, vy, vz
  !real(rk) :: merr, vxerr, vyerr, vzerr
  real(rk) :: vxtmp, vytmp, vztmp

  !----
  ! Executable code

  msys = mass(1) !; merr = 0.0_rk
  vx = 0.0_rk !; vxerr = 0.0_rk
  vy = 0.0_rk !; vyerr = 0.0_rk
  vz = 0.0_rk !; vzerr = 0.0_rk

  !FIRSTPRIVATE(merr, vxerr, vyerr, vzerr)
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, vxtmp, vytmp, vztmp) &
  !$OMP SHARED(nbod, mass, vxh, vyh, vzh, vxb, vyb, vzb, msys, vx, vy, vz)

  !$OMP DO SCHEDULE(STATIC) REDUCTION(+ : msys, vx, vy, vz)
  do i = 2, nbod

    msys = msys + mass(i)
    vx = vx + mass(i)*vxh(i)
    vy = vy + mass(i)*vyh(i)
    vz = vz + mass(i)*vzh(i)
    !msys = util_kahan_sum(msys, mass(i), merr)
    !vx = util_kahan_sum(vx, mass(i)*vxh(i), vxerr)
    !vy = util_kahan_sum(vy, mass(i)*vyh(i), vyerr)
    !vz = util_kahan_sum(vz, mass(i)*vzh(i), vzerr)

  end do
  !$OMP END DO

  vxtmp = -vx/msys
  vytmp = -vy/msys
  vztmp = -vz/msys

  !$OMP DO SCHEDULE(STATIC)
  do i = 2, nbod

    vxb(i) = vxh(i) + vxtmp
    vyb(i) = vyh(i) + vytmp
    vzb(i) = vzh(i) + vztmp

  end do
  !$OMP END DO NOWAIT

  !$OMP END PARALLEL

  vxb(1) = -vx/msys
  vyb(1) = -vy/msys
  vzb(1) = -vz/msys

  return
  end subroutine coord_vh2b
!!
!!
end module coord