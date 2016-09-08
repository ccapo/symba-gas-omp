subroutine discard_massive5(time, dt, nbod, nbodm, mass, xh, yh, zh, vxh, vyh, vzh, rmin, rmax, rmaxu, qmin, lclose, rpl, rhill, &
           isenc, mergelst, mergecnt, iecnt, eoff, i1st, ibound)
!*************************************************************************
!                            DISCARD_MASSIVE5.F
!*************************************************************************
! This subroutine checks to see if a massive body should be discarded or
! merged.
!
!             Input:
!                 time          ==>  current time (real scalar)
!                 dt            ==>  time step  (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh, yh, zh      ==>   position in helio coord
!                                    (real arrays)
!                 vxh, vyh, vzh   ==>   pl vel in helio coord
!                                    (real arrays)
!                 rmin, rmax      ==>  maximum and min distance from Sun
!                                     if <0  then don't check
!                                        (real scalar)
!                 rmaxu          ==>  maximum distance from Sun in not bound
!                                     if <0  then don't check
!                                        (real scalar)
!                  qmin          ==> Smallest perihelion distance
!                                      if <0  then don't check
!                                          (real scalar)
!                 lclose        ==> .true. --> marge particles if they
!                                    get too close. Read in that
!                                    distance in io_init_pl
!                                      (logical(lk) ::*2 scalar)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!                 rhill         ==>  size of a planet's hill's sphere.
!                                    (real array)
!                 isenc         ==>  0 --> No encounter during last dt
!                                    1 --> There was encounters
!                                     (integer scalar)
!                 eoff          ==> Amount of energy lost due to discards
!                                          (real scalar)
!                 mergelst      ==>  list of mergers (int array)
!                 mergecnt      ==>  count of mergers (int array)
!                 iecnt         ==>  Number of encounters (int*2 array)
!                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
!             Output:
!                 nbod          ==>  recalculated number of massive bodies
!                                       (int scalar)
!                 mass          ==>  recalculated mass of bodies (real array)
!                 xh, yh, zh      ==>  recalculated position in helio coord
!                                    (real arrays)
!                 vxh, vyh, vzh   ==>  recalculated pl vel in helio coord
!                                    (real arrays)
!                 rpl           ==> recalculated physical sizes of a planet.
!                                    (real array)
!                 rhill         ==>  reordered size of planet's hill's sphere.
!                                    (real array)
!                 eoff          ==> Updated amount of energy lost from discards
!                                          (real scalar)
!                 i1st          ==>  set to 0 if reordered (int scalar)
!
!
! Remarks:
!
! Authors:  Hal Levison
! Date:    12/30/96
! Last revision: 5/13/99
use module_swift
use module_symba5
use module_interfaces, except_this_one => discard_massive5
implicit none

! Inputs:
logical(lk) :: lclose
integer(ik) :: mergelst(NENMAX,2), mergecnt, isenc
integer(ik) :: iecnt(nbod)
real(rk) :: time, dt
real(rk) :: rmin, rmax, rmaxu, qmin

! Input and Output
integer(ik) :: nbod, nbodm, i1st, ibound(nbod)
real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)
real(rk) :: eoff, rpl(nbod), rhill(nbod)

! internal
logical(lk) :: lrflg(nbod), lwhy
integer(ik) :: i, j, iu, iflag, i1, i2, iwhy(nbod)
integer(ik), save :: isperih(NTPMAX)
real(rk) :: xb(nbod), yb(nbod), zb(nbod)
real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)
real(rk) :: rmin2, rmax2, rmaxu2, energy
real(rk) :: ei, ef, ke, pot, eltot(NDIM), vdotr
real(rk) :: rh2, rb2, vb2, msys
character(len = 24) :: cdummy

!-----------------!
! Executable code !
!-----------------!

! set things up
lwhy = .false.
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, lrflg, iwhy)
do i = 1, nbod

  lrflg(i) = .true.
  iwhy(i) = 0

end do
!$OMP END PARALLEL DO

do i = 1, mergecnt

  i2 = mergelst(i,2)

  if(lrflg(i2)) then

    lrflg(i2) = .false.

  else

    mergelst(i,2) = -1

  end if

end do

! take care of mergers
do i = 1, mergecnt

  i1 = mergelst(i,1)
  i2 = mergelst(i,2)
  vdotr = xh(i1)*vxh(i1) + yh(i1)*vyh(i1) + zh(i1)*vzh(i1)
  if(vdotr > 0.0_rk) then

    isperih(i1) = 1

  else

    isperih(i1) = -1

  end if

  if(i2 > 0) then

    call discard_mass_reorder5(i2, nbod, mass, xh, yh, zh, vxh, vyh, vzh, rpl, rhill, isperih, ibound)
    i1st = 0

    ! If the index of other particles due for removal are larger than i2, decrease each by one
    do j = i + 1, mergecnt

      if(mergelst(j,1) > i2) mergelst(j,1) = mergelst(j,1) - 1
      if(mergelst(j,2) > i2) mergelst(j,2) = mergelst(j,2) - 1

    end do

  end if

end do

! check for position
if((rmin >= 0.0_rk) .or. (rmax >= 0.0_rk) .or. (rmaxu >= 0.0_rk)) then

  rmin2 = rmin*rmin
  rmax2 = rmax*rmax
  rmaxu2 = rmaxu*rmaxu

  call coord_h2b(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xb, yb, zb, vxb, vyb, vzb, msys)

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, rh2, rb2, vb2, energy) &
  !$OMP FIRSTPRIVATE(rmin, rmax, rmaxu, rmin2, rmax2, rmaxu2, msys, time) &
  !$OMP SHARED(nbod, xh, yh, zh, xb, yb, zb, vxb, vyb, vzb, iwhy, iecnt, lwhy)
  do i = 2, nbod

    rh2 = xh(i)**2 + yh(i)**2 + zh(i)**2
    if((rmax >= 0.0_rk) .and. (rh2 > rmax2)) then

      write(*,*) 'Particle ', i, ' too far from Sun at t = ', time
      iwhy(i) = -3
      lwhy = .true. ! Set global removal flag

    end if

    if((rmin >= 0.0_rk) .and. (rh2 < rmin2)) then

      write(*,*) 'Particle ', i, ' too close from Sun at t = ', time
      iwhy(i) = 1
      lwhy = .true. ! Set global removal flag

    end if

    if((iecnt(i) == 0) .and. (rmaxu >= 0.0_rk) .and. (iwhy(i) == 0)) then

      rb2 = xb(i)**2 + yb(i)**2 + zb(i)**2
      vb2 = vxb(i)**2 + vyb(i)**2 + vzb(i)**2
      energy = 0.5_rk*vb2 - msys/sqrt(rb2)

      if((energy > 0.0_rk) .and. (rb2 > rmaxu2)) then

        write(*,*) 'Particle ', i, ' is unbound and too far from barycenter at t = ', time
        iwhy(i) = -2
        lwhy = .true. ! Set global removal flag

      end if

    end if

  end do
  !$OMP END PARALLEL DO

end if

! check perihelion distance
if(qmin >= 0.0_rk) call discard_mass_peri(time, nbod, iecnt, mass, xh, yh, zh, vxh, vyh, vzh, qmin, iwhy, lwhy, isperih)

! If there are any particles that need to be removed
if(lwhy) then

  ! Unit number for discard output file
  iu = 40

  ! Compute the total energy of the system *before* removing particles
  call anal_energy_discard5(0, nbod, nbodm, mass, 0.0_rk, 0.0_rk, xh, yh, zh, vxh, vyh, vzh, ke, pot, ei, eltot)

  ! Loop over particles
  do i = 2, nbod

    ! If this particle has been flagged for removal
    if(iwhy(i) /= 0) then

      ! Write the relevant information to discard file
      call io_discard_mass(1, time, i, mass(i), rpl(i), xh(i), yh(i), zh(i), vxh(i), vyh(i), vzh(i), iu, iwhy(i), cdummy)

      ! Update the list of removal flags by shifting all entries below i up one
      if(i < nbod) then

        do j = i, nbod - 1

          iwhy(j) = iwhy(j + 1)

        end do

      end if

      ! Update relevant quantities by shifting all entries below i up one
      call discard_mass_reorder5(i, nbod, mass, xh, yh, zh, vxh, vyh, vzh, rpl, rhill, isperih, ibound)
      i1st = 0

    end if

  end do

  ! Compute the total energy of the system *after* removing particles, then update the energy offset
  call anal_energy_discard5(0, nbod, nbodm, mass, 0.0_rk, 0.0_rk, xh, yh, zh, vxh, vyh, vzh, ke, pot, ef, eltot)
  eoff = eoff + ei - ef

endif

return
end subroutine discard_massive5
