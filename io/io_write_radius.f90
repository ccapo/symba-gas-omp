subroutine io_write_radius(time,nbod,radius,oname,iu,fopenstat)
!-------------------------------------------------------------------------
!				IO_WRITE_RADIUS.F90
!-------------------------------------------------------------------------
! NEW VERSION OF THIS, USES FXDR
! Write out radiii
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
! Date:    11/2/99
! Last revision:
use module_swift
use module_io
use module_fxdr
use module_interfaces, except_this_one => io_write_radius
implicit none

! Inputs:
integer(ik) :: nbod,iu
real(rk) :: radius(nbod),time
character(len = *) :: oname,fopenstat

! Internals
integer(ik) :: ierr,i
integer(ik), save :: i1st = 0 ! i1st = 0 first time through; i1st = 1 after
real(real4) :: radius4(NTPMAX), ttmp

!-----------------
! Executable code
!-----------------

! If first time through open file
if(i1st == 0) then

  if((fopenstat == 'append') .or. (fopenstat == 'APPEND')) then

    call io_open_fxdr('radius.' // oname, 'a', .true., iu, ierr)

  else

    if((fopenstat == 'new') .or. (fopenstat == 'NEW')) then

      call io_open_fxdr('radius.' // oname, 'w', .true., iu, ierr)

      if(iu >= 0) then

        write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
        write(*,'(a)') '  Binary output file exists'
        call util_exit(FAILURE)

      end if

    end if

    call io_open_fxdr('radius.' // oname, 'w', .true., iu, ierr)

  end if

  if(ierr < 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
    write(*,'(a)') '  Could not open binary output file'
    call util_exit(FAILURE)

  end if

  i1st = 1

else

  call io_open_fxdr('radius.' // oname, 'a', .true., iu, ierr)

  if(ierr < 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
    write(*,'(a)') '  Could not open binary output file with append'
    call util_exit(FAILURE)

  end if

end if

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, radius, radius4)
do i = 1, nbod

  radius4(i) = radius(i)

end do
!$OMP END PARALLEL DO

ttmp = time

ierr = ixdrreal( iu, ttmp )
!ierr = ixdrdouble( iu, time )
if(ierr < 0) then

  write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
  write(*,'(a)') '  Could not write time'
  call util_exit(FAILURE)

end if

ierr = ixdrint( iu, nbod )
if(ierr < 0) then

  write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
  write(*,'(a)') '  Could not write nbod'
  call util_exit(FAILURE)

end if

ierr = ixdrrmat( iu, nbod, radius4 )
!ierr = ixdrdmat( iu, nbod, radius )
if(ierr < 0) then

  write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
  write(*,'(a)') '  Could not write radii'
  call util_exit(FAILURE)

end if

ierr = ixdrclose(iu)
if(ierr < 0) then

  write(*,'(a)') ' SWIFT ERROR: in io_write_frame: '
  write(*,'(a)') '  Could not close radius output file'
  call util_exit(FAILURE)

end if

return
end subroutine io_write_radius
