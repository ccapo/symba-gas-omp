subroutine symba5_nbodm(nbod, mass, mtiny, nbodm)
!*************************************************************************
!                            SYMBA5_NBODM.F
!*************************************************************************
! Returns the location of the last massive body in the list
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 mtiny         ==>  Small mass  (real array)
!             Output:
!                 nbodm         ==>  location of the last massive body
!                                    (int scalar)
!
! Remarks:  If all the objects are massive,  then nbodm=nbod-1 so that
!           the do loops will have the correct limits.
! Authors:  Hal Levison
! Date:    3/20/97
! Last revision: 1/29/06
use module_swift
use module_symba5
implicit none

! Inputs Only:
integer(ik) :: nbod
real(rk) :: mass(nbod),  mtiny

! Outputs only
integer(ik) :: nbodm

! Internals
integer(ik) :: i
real(rk) :: mtiny_private

!-----------------!
! Executable Code !
!-----------------!

! This assumes that the particle identifiers are already sorted in order of decreasing mass
if(mass(nbod) > mtiny) then

  nbodm = nbod - 1

  write(*,'(a,i6,a)') " Out of ", nbod, " bodies, all are massive."

else

  ! Initialize the value of nbodm
  nbodm = 1

  ! Make a private copy of mtiny available for each thread
  mtiny_private = mtiny

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) &
  !$OMP FIRSTPRIVATE(mtiny_private) SHARED(nbod, mass) REDUCTION(MAX : nbodm)
  do i = 2, nbod - 1

    if(mass(i) > mtiny_private) nbodm = i

  end do
  !$OMP END PARALLEL DO

  write(*,'(2(a,i6),a)') " Out of ", nbod, " bodies, ", nbodm, " are massive."

end if

return
end subroutine symba5_nbodm
