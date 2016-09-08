function io_read_hdr_r(iu,time,nbod,nleft) result(read_hdr_r)
!-------------------------------------------------------------------------
!			IO_READ_HDR_R.F90
!-------------------------------------------------------------------------
! Read in header part of the real*4 file
!
!             Input:
!                 iu           ==> unit number to write to
!             Output:
!                 time         ==> current time (real scalar)
!                 nbod         ==> number of massive bodies (int scalar)
!                 nleft        ==> number of active tp (int scalar)
!
!             Returns:
!               io_read_hdr_r  ==> = 0 read ok
!                                   != 0 read failed is set to iostat variable
! Remarks:
! Authors:  Hal Levison
! Date:    2/22/94
! Last revision:
use module_swift
use module_io
implicit none

! Inputs:
integer(ik) :: iu

! Output
integer(ik) :: nbod,nleft,read_hdr_r
real(rk) :: time

! Internals
integer(ik) :: ierr
integer(integer2) :: nleft2,nbod2
real(real4) :: ttmp

!-----------------
! Executable code
!-----------------

read(iu,iostat=ierr) ttmp,nbod2,nleft2
read_hdr_r = ierr
if(ierr /= 0) return

nbod = nbod2
nleft = nleft2
time = ttmp

return
end function io_read_hdr_r