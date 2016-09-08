subroutine io_open_fxdr(fname,fopenstat,lflg,iu,ierr)
!-------------------------------------------------------------------------
!				IO_OPEN_FXDR.F90
!-------------------------------------------------------------------------
! open files using the fxdr routines
!
! Input:
!   fname           ==>  file name (character*80)
!   fopenstat       ==>  The status flag for the open
!                        statements of the output files.
!                        (character*1)
!   lflg            ==>  if .TRUE., then FXDR routines return even if
!                        there is an I/O error.  If .FALSE., then routines
!                        halt on I/O error (logical scalar)
!
! Output:
!   iu              ==>  unit number (integer scalar)
!   ierr            ==>  output from iostat
!
! Remarks:
! Authors:  Hal Levison
! Date:    11/3/99
! Last revision:
use module_swift
use module_io
use module_fxdr
implicit none

! Inputs:
logical(lk) :: lflg
integer(ik) :: iu
character(len = *) :: fname
character(len = 1) :: fopenstat

! Outputs:
integer(ik) :: ierr

! Internals:

!-----------------!
! Executable code !
!-----------------!

iu = initxdr(fname, fopenstat, lflg)

if(iu > 0) then

  ierr = 0

else

  ierr = iu

end if

return
end subroutine io_open_fxdr
