function io_read_hdr(iu,time,nbod,nleft) result(read_hdr)
!-------------------------------------------------------------------------
!				IO_READ_HDR.F90
!-------------------------------------------------------------------------
! NEW VERSION OF THIS, USES FXDR
!
!             Input:
!                 iu            ==> unit number to write to
!             Output:
!                 time          ==> current time (real scalar)
!                 nbod          ==> number of massive bodies (int scalar)
!                 nleft         ==> number of active tp (int scalar)
!
!             Returns:
!               io_read_hdr     ==> = 0 read ok
!                                   !=0 read failed, set to iostat variable
! Remarks:
! Authors:  Hal Levison
! Date:    11/2/99
! Last revision:
use module_swift
use module_io
use module_fxdr
implicit none

! Inputs:
integer(ik) :: iu

! Output
integer(ik) :: nbod,nleft,read_hdr
real(rk) :: time

! Internals
integer(ik) :: nn(2)
integer(ik) :: ierr
real(real4) :: ttmp

!-----------------
! Executable code
!-----------------

ierr = ixdrreal( iu, ttmp )
!ierr = ixdrdouble( iu, time )
read_hdr = ierr

if(ierr /= 0) return

ierr = ixdrimat( iu, 2, nn )
read_hdr = ierr

if(ierr /= 0) return

nbod = nn(1)
nleft = nn(2)
time = ttmp

return
end function io_read_hdr
