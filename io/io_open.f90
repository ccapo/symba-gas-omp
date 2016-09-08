subroutine io_open(iu,fname,fopenstat,format,ierr)
!-------------------------------------------------------------------------
!				IO_OPEN.F90
!-------------------------------------------------------------------------
! open files
!
!             Input:
!                 iu              ==>  unit number (integer scalar)
!                 fname           ==>  file name (character*80)
!                 fopenstat       ==>  The status flag for the open
!                                      statements of the output files.
!                                          (character*80)
!                 format          ==>  format string (character*80)
!             Output:
!                 ierr            ==>  output from iostat
!
! Remarks:
! Authors:  Hal Levison
! Date:    3/3/94
! Last revision: 1/30/98
use module_swift
use module_io
implicit none

! Inputs:
integer(ik) :: iu
character(len = *) :: fname, fopenstat, format

! Outputs:
integer(ik) :: ierr

!-----------------!
! Executable code !
!-----------------!

if((fopenstat == 'append') .or. (fopenstat == 'APPEND')) then

  open(unit = iu, file = fname, status = 'old', POSITION = 'append', FORM = format, IOSTAT = ierr)

  if(ierr /= 0) then

    write(*,*) 'Warning:  Could not open ', fname, ' with position = append.'
    write(*,*) '          Will open as status = new'

    open(unit = iu, file = fname, status = 'new', form = format, iostat = ierr)

  end if

else

  open(unit = iu, file = fname, status = fopenstat, form = format, iostat = ierr)

end if


return
end subroutine io_open
