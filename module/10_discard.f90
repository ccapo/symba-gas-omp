module discard
! Module for routines in discard directory
use swift
use util
use coord
use anal
use io
implicit none

! Maximum number of encounters
integer(ik), private, parameter :: NENMAX = 262144 ! must be less than 2^18 = 262144

contains
!!
!!
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
!!
  subroutine discard_mass_merge5(time, nbod, nbodm, ip1, ip2, mass, xh, yh, zh, vxh, vyh, vzh, rpl, eoff, ielc, ielst)
  !*************************************************************************
  !                            DISCARD_MASS_MERGE5.F
  !*************************************************************************
  ! Merge two massive bodies
  !
  !             Input:
  !                 time          ==>  current time (real scalar)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 ip1, ip2       ==>  planets to merge (real scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 xh, yh, zh      ==>   position in helio coord
  !                                    (real arrays)
  !                 vxh, vyh, vzh   ==>   pl vel in helio coord
  !                                    (real arrays)
  !                 rpl           ==>  physical size of a planet.
  !                                    (real array)
  !                 eoff          ==> Amount of energy lost due to discards
  !                                          (real scalar)
  !                ielc           ==>  number of encounters (integer*2 scalar)
  !                ielst          ==>  list of ecnounters (2D integer*2 array)
  !             Output:
  !                 mass          ==>  recalculated mass of bodies (real array)
  !                 xh, yh, zh      ==>  recalculated position in helio coord
  !                                    (real arrays)
  !                 vxh, vyh, vzh   ==>  recalculated pl vel in helio coord
  !                                    (real arrays)
  !                 rpl           ==>  recalculated physical sizes of a planet.
  !                                    (real array)
  !                 eoff          ==> Updated amount of energy lost from discards
  !                                          (real scalar)
  !                ielc           ==>  number of encounters (integer*2 scalar)
  !                ielst          ==>  list of ecnounters (2D integer*2 array)
  !
  ! Remarks:
  !
  ! Authors:  Hal Levison
  ! Date:    12/30/96
  ! Last revision: 1/30/97
  implicit none

  ! Inputs:
  integer(ik) :: ip1, ip2
  real(rk) :: time

  ! Input and Output
  integer(ik) :: nbod, nbodm
  integer(ik) :: ielc, ielst(NENMAX,2)
  real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod), rpl(nbod)
  real(rk) :: eoff

  ! internal
  integer(ik) :: i, j, itmp
  real(rk) :: mtot, m1, r1
  real(rk) :: x1, y1, z1
  real(rk) :: vx1, vy1, vz1
  real(rk) :: m2, r2
  real(rk) :: x2, y2, z2
  real(rk) :: vx2, vy2, vz2
  real(rk) :: j2rp2, j4rp4, ke, pot, energy1, energy2, eltot(NDIM)

  !-----------------!
  ! Executable code !
  !-----------------!

  j2rp2 = 0.0_rk
  j4rp4 = 0.0_rk
  call anal_energy_discard5(0, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, ke, pot, energy1, eltot)

  !if(mass(ip2) > mass(ip1)) then
  !
  !  itmp = ip1
  !  ip1 = ip2
  !  ip2 = itmp
  !
  !end if

  write(*,'(a20,i6,a5,i6,a7,1pe12.5)') 'Merging particles: ', ip1,  ' and ',  ip2, ' at t =', time

  x1 = xh(ip1)
  y1 = yh(ip1)
  z1 = zh(ip1)
  vx1 = vxh(ip1)
  vy1 = vyh(ip1)
  vz1 = vzh(ip1)
  m1 = mass(ip1)
  r1 = rpl(ip1)

  x2 = xh(ip2)
  y2 = yh(ip2)
  z2 = zh(ip2)
  vx2 = vxh(ip2)
  vy2 = vyh(ip2)
  vz2 = vzh(ip2)
  m2 = mass(ip2)
  r2 = rpl(ip2)

  ! Note:  I am just putting these guys together here,  which is
  !        clearly wrong.  I should integrate back to the time
  !        of close approach.
  mtot = m1 + m2
  rpl(ip1) = (r1**3 + r2**3)**(1.0_rk/3.0_rk)
  vxh(ip1) = (m1*vx1 + m2*vx2)/mtot
  vyh(ip1) = (m1*vy1 + m2*vy2)/mtot
  vzh(ip1) = (m1*vz1 + m2*vz2)/mtot
  mass(ip1) = mtot

  ! Put in zeros for the rest the second particle
  xh(ip2) = x2*1.0e10_rk   ! so danby does not fail
  yh(ip2) = y2*1.0e10_rk
  zh(ip2) = z2*1.0e10_rk
  vxh(ip2) = 0.0_rk
  vyh(ip2) = 0.0_rk
  vzh(ip2) = 0.0_rk
  mass(ip2) = 0.0_rk
  rpl(ip2) = 0.0_rk

  ! Remove any encounters with ip2
  j = 1
  do while(j <= ielc)

    if((ielst(j,1) == ip2) .or. (ielst(j,2) == ip2)) then

      do i = j + 1, ielc

        ielst(i - 1,1) = ielst(i,1)
        ielst(i - 1,2) = ielst(i,2)

      end do

      ielc = ielc - 1

    else

      j = j + 1

    end if

  end do

  call anal_energy_discard5(0, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, ke, pot, energy2, eltot)
  eoff = eoff + energy1 - energy2

  call io_discard_merge(time, ip1, ip2, m1, r1, x1, y1, z1, vx1, vy1, vz1, m2, r2, x2, y2, z2, vx2, vy2, vz2, &
                        mass(ip1), rpl(ip1), xh(ip1), yh(ip1), zh(ip1), vxh(ip1), vyh(ip1), vzh(ip1))

  return
  end subroutine discard_mass_merge5
!!
  subroutine discard_mass_peri(time, nbod, iecnt, mass, xh, yh, zh, vxh, vyh, vzh, qmin, iwhy, lwhy, isperi)
  !*************************************************************************
  !                            DISCARD_MASS_PERI.F
  !*************************************************************************
  ! This subroutine checks to see if a partical should be discarded because
  ! of its perihelion distance gets too small
  !
  !             Input:
  !                 time           ==>  current time (real scalar)
  !                 nbod           ==>  number of test bodies (int scalar)
  !                 iecnt          ==>  Number of encounters (int*2 array)
  !                 mass           ==>  mass of bodies (real array)
  !                 xh, yh, zh       ==>   part position in helio coord
  !                                      (real arrays)
  !                 vxh, vyh, vzh    ==>   part vel in helio coord
  !                                      (real arrays)
  !                 qmin           ==>  Smallest perihelion distance
  !                                      (real scalar)
  !                 iwhy           ==>  status of the object
  !                                      (integer array)
  !             Output:
  !                 iwhy           ==>  status of the object
  !                                      (integer array)
  !                 isperi         ==> = 0 if tp went through peri
  !                                    =-1 if tp pre peri
  !                                    = 1 if tp post peri
  !                                         (integer array)
  ! Remarks: Based on discard_peri
  ! Authors:  Hal Levison
  ! Date:    12/30/96
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: nbod
  integer(ik) :: iecnt(nbod)
  real(rk) :: mass(nbod), time, qmin
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

  ! Input and Output
  logical(lk) :: lwhy
  integer(ik) :: iwhy(nbod)
  integer(ik) :: isperi(nbod)

  ! internal
  logical(lk), save :: lfirst = .true.
  integer(ik) :: i
  real(rk) :: peri(nbod)
  logical(lk) :: lperi(nbod)

  !-----
  ! Executable code

  if(lfirst) then ! if first time through, set things up

    call util_mass_peri(0, nbod, xh, yh, zh, vxh, vyh, vzh, mass, isperi, peri, lperi)
    lfirst = .false.

  else

    call util_mass_peri(1, nbod, xh, yh, zh, vxh, vyh, vzh, mass, isperi, peri, lperi)

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(time, qmin) &
    !$OMP SHARED(nbod, isperi, iecnt, peri, iwhy, lwhy)
    do i = 2, nbod

      if((isperi(i) == 0) .and. (iecnt(i) == 0) .and. (peri(i) <= qmin)) then

        write(*,'(a,i7,a,1pe14.6)') 'Particle ', i, ' perihelion distance too small at t = ', time
        iwhy(i) = -4
        lwhy = .true.

      end if

    end do
    !$OMP END PARALLEL DO

  end if

  return
  end subroutine discard_mass_peri
!!
  subroutine discard_mass_reorder5(ip, nbod, mass, xh, yh, zh, vxh, vyh, vzh, rpl, rhill, isperih, ibound)
  !*************************************************************************
  !                            DISCARD_MASS_REORDER5.F
  !*************************************************************************
  ! Remove a massive body
  !
  !             Input:
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 ip            ==>  planets to remove (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 xh, yh, zh      ==>   position in helio coord
  !                                    (real arrays)
  !                 vxh, vyh, vzh   ==>   pl vel in helio coord
  !                                    (real arrays)
  !                 rpl           ==>  physical size of a planet.
  !                                    (real array)
  !                 rhill         ==>  size of a planet's hill's sphere.
  !                                    (real array)
  !                 isperih       ==> heliocentric peri flags. (real array)
  !             Output:
  !                 ip            ==>  planets to remove (int scalar)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 xh, yh, zh      ==>   position in helio coord
  !                                    (real arrays)
  !                 vxh, vyh, vzh   ==>   pl vel in helio coord
  !                                    (real arrays)
  !                 rpl           ==>  physical size of a planet.
  !                                    (real array)
  !                 rhill         ==>  size of a planet's hill's sphere.
  !                                    (real array)
  !                 isperih       ==> heliocentric peri flags. (real array)
  !
  ! Remarks:
  !
  ! Authors:  Hal Levison
  ! Date:    1/2/97
  ! Last revision: 5/13/99
  implicit none

  ! Inputs:
  integer(ik) :: ip

  ! Input and Output
  integer(ik) :: nbod
  integer(ik) :: isperih(nbod), ibound(nbod)
  real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod), rpl(nbod)
  real(rk) :: rhill(nbod)
  real(rk) :: rdrag(NTPMAX), dragc(NTPMAX)

  common / sizedist / rdrag, dragc

  ! internal
  integer(ik) :: i

  !-----
  ! Executable code

  ! Move all the values from bodies [ip + 1, nbod] up one index
  do i = ip, nbod - 1

    xh(i) = xh(i + 1)
    yh(i) = yh(i + 1)
    zh(i) = zh(i + 1)
    vxh(i) = vxh(i + 1)
    vyh(i) = vyh(i + 1)
    vzh(i) = vzh(i + 1)
    mass(i) = mass(i + 1)
    rpl(i) = rpl(i + 1)
    rhill(i) = rhill(i + 1)
    isperih(i) = isperih(i + 1)
    ! *** 04/24/08 -- CCC ***
    rdrag(i) = rdrag(i + 1)
    dragc(i) = dragc(i + 1)
    ibound(i) = ibound(i + 1)

  end do

  ! Decrease the nummber of bodies by one
  nbod = nbod - 1

  return
  end subroutine discard_mass_reorder5
!!
!!
end module discard