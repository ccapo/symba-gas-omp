subroutine io_write_radius_r(time,nbod,radius,oname,iu,fopenstat)
!-------------------------------------------------------------------------
!			IO_WRITE_RADIUS_R.F90
!-------------------------------------------------------------------------
! Write out radiuses
!
! Input:
!  time       ==>  current time (real scalar)
!  nbod       ==>  number of massive bodies (int scalar)
!  radius     ==>  radii of bodies (real array)
!  oname      ==>  output file name (character string)
!  iu         ==>  unit number to write to
!  fopenstat  ==>  The status flag for the open
!                  statements of the output files.
!                  (character*80)
!
! Remarks: Based on io_write_frame
! Authors:  Hal Levison
! Date:    1/9/97
! Last revision: 11/2/99
use module_swift
use module_io
use module_interfaces, except_this_one => io_write_radius_r
implicit none

! Inputs:
integer(ik) :: nbod,iu
real(rk) :: radius(nbod),time
character(len = *) :: oname,fopenstat

! Internals
integer(ik) :: i,ierr
integer(ik), save :: i1st = 0 ! i1st = 0 first time through; i1st = 1 after
integer(integer2) :: nbod2
real(real4) :: radius4(NTPMAX), ttmp

!-----------------
! Executable code
!-----------------

! If first time through open file
if(i1st == 0) then

  call io_open(iu,'radius.' // oname,fopenstat,'UNFORMATTED',ierr)

  if(ierr /= 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_radius_r:'
    write(*,'(a)') '  Could not open binary output file'
    call util_exit(FAILURE)

  end if

  i1st = 1

else

  call io_open(iu,'radius.' // oname,'append','UNFORMATTED',ierr)

end if

do i = 1, nbod

  radius4(i) = radius(i)

end do
nbod2 = nbod
ttmp = time

write(iu) ttmp,nbod2
write(iu) (radius4(i), i = 1, nbod2)

close(unit = iu)

return
end subroutine io_write_radius_r
