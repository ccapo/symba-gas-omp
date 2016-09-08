subroutine symba5_enc_drift(nbod, ielc, ielst, ielev, irec, mass, xh, yh, zh, vxb, vyb, vzb, dt)
!*************************************************************************
!                        SYMBA5_ENC_DRIFT.F90
!*************************************************************************
! This subroutine loops thorugh the particles and calls the danby routine
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 ielev         ==>  Level of particles (int array)
!                 irec          ==>  current level of the code
!                 mass          ==>  mass of bodies (real array)
!                 xh, yh, zh      ==>  initial position in helio coord
!                                    (real arrays)
!                 vxb, vyb, vzb   ==>  initial position in bary coord
!                                    (real arrays)
!                 dt            ==>  time step
!             Output:
!                 xh, yh, zh      ==>  final position in helio coord
!                                       (real arrays)
!                 vxb, vyb, vzb   ==>  final position in bary coord
!                                       (real arrays)
!
! Remarks:  Based on helio_drift.f
! Authors:  Hal Levison
! Date:    1/20.97
! Last revision:
use module_swift
use module_symba5
use module_interfaces, except_this_one => symba5_enc_drift
! use omp_lib
implicit none

! Inputs Only:
integer(ik) :: nbod, irec, ielc
integer(ik) :: ielst(NENMAX,2), ielev(nbod)
real(rk) :: mass(nbod), dt

! Inputs and Outputs:
real(rk) :: xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod)

! Internals:
logical(lk) :: ldflag(nbod)
integer(ik) :: i, j, ie, iflag
real(rk) :: mu

!----
! Executable code
!-----------------

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, ldflag)
do i = 1, nbod

  ldflag(i) = .true.

end do
!$OMP END PARALLEL DO

! Make a copy of the central body's mass
mu = mass(1)

! Take a drift forward dth
do ie = 1, ielc

  i = ielst(ie,1)
  j = ielst(ie,2)

  ! Perform a drift for particle i
  if((ielev(i) == irec) .and. (mass(i) /= 0.0_rk) .and. ldflag(i)) then

    call drift_one(mu, xh(i), yh(i), zh(i), vxb(i), vyb(i), vzb(i), dt, iflag)
    ldflag(i) = .false. ! Set so we don't return here

    if(iflag /= 0) then

      write(*,'(a,i7,a)')         'Particle:    ', i, ' is lost!!!!!!!!!'
      write(*,'(a,2(1pe14.6))')   'Mass, dt:    ', mu, dt
      write(*,'(a,3(1pe14.6))')   'Helio. pos.: ', xh(i), yh(i), zh(i)
      write(*,'(a,3(1pe14.6),/)') 'Bary. vel.:  ', vxb(i), vyb(i), vzb(i)
      call util_exit(FAILURE)

    end if

  end if

  ! Perform a drift for particle j
  if((ielev(j) == irec) .and. (mass(j) /= 0.0_rk) .and. ldflag(j)) then

    call drift_one(mu, xh(j), yh(j), zh(j), vxb(j), vyb(j), vzb(j), dt, iflag)
    ldflag(j) = .false. ! Set so we don't return here

    if(iflag /= 0) then

      write(*,'(a,i7,a)')         'Particle:    ', j, ' is lost!!!!!!!!!'
      write(*,'(a,2(1pe14.6))')   'Mass, dt:    ', mu, dt
      write(*,'(a,3(1pe14.6))')   'Helio. pos.: ', xh(j), yh(j), zh(j)
      write(*,'(a,3(1pe14.6),/)') 'Bary. vel.:  ', vxb(j), vyb(j), vzb(j)
      call util_exit(FAILURE)

    end if

  end if

end do

return
end subroutine symba5_enc_drift
