subroutine io_discard_mass(init,time,id,m1,r1,x1,y1,z1,vx1,vy1,vz1,iu,iwhy,fopenstat)
!------------------------------------------------------------------------------------
!				IO_DISCARD_MASS.F90
!------------------------------------------------------------------------------------
! Write out information about a discarded massive body.
!
!             Input:
!                 init          ==>  initiize flag if = 0 initialize and return
!                                                     = 1 run through
!                 id            ==> particle number (int scalar)
!                 time          ==>  current time (real scalar)
!                 m1            ==>  Mass of pl (real scalar)
!                 r1            ==>  Radius of pl 2 (real scalar)
!                 x1,y1,z1      ==>  current position of pl 1 in helio coord
!                                    (real scalar)
!                 vx1,vy1,vz1   ==>  current velocity of pl 1 in helio coord
!                                    (real scalar)
!                 iu            ==> IO unit (int scalar)
!                 iwhy          ==> reason for discard (int scalar)
!                 fopenstat     ==>  The status flag for the open
!                                      statements of the output files.
!                                          (character*80)
! Remarks:
! Authors:  Hal Levison
! Date:    12/30/96
! Last revision:
use module_swift
use module_io
use module_interfaces, except_this_one => io_discard_mass
implicit none

! Inputs:
integer(ik) :: iwhy,iu,init,id
real(rk) :: time
real(rk) :: m1,r1
real(rk) :: x1,y1,z1
real(rk) :: vx1,vy1,vz1
character(len = *) :: fopenstat

! Internals
integer(ik) :: ierr

!-----------------
! Executable code
!-----------------

if(init == 0) then

  ! Replaced fopenstat with 'unknown' for initialization call
  !call io_open(iu,'discard_mass.dat',fopenstat,'FORMATTED',ierr)
  call io_open(iu,'discard_mass.dat','unknown','FORMATTED',ierr)

  ! If there was an error and fopenstat='append' then
  ! try to open as **unknown** (CCC -- 09/05/07)
  if(ierr /= 0) then

    if((fopenstat == 'append') .or. (fopenstat == 'APPEND')) then

      call io_open(iu,'discard_mass.dat','unknown','FORMATTED',ierr)

    end if

  end if

  if(ierr /= 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_discard_mass: '
    write(*,'(a)') '  Could not open discard output file'
    call util_exit(FAILURE)

  end if

  return ! <== NOTE!!!!

else

  call io_open(iu,'discard_mass.dat','append','FORMATTED',ierr)

end if

write(iu,'(1x,1pe22.15,1x,i4)') time,iwhy
write(iu,'("-1",1x,i6,1x,2(1x,1pe22.15))') id,m1,r1
write(iu,'(3(1x,1pe22.15))') x1,y1,z1
write(iu,'(3(1x,1pe22.15))') vx1,vy1,vz1

close(unit = iu)

return
end subroutine io_discard_mass
