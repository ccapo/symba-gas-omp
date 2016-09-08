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
use module_swift
use module_interfaces, except_this_one => discard_mass_peri
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
