subroutine io_write_hdr_r(iu,time,nbod,ntp,istat)
!-------------------------------------------------------------------------
!				IO_WRITE_HDR_R.F90
!-------------------------------------------------------------------------
! Write out header part of the real*4 binary file
!
! Input:
!  iu     ==>  unit number to write to
!  time   ==>  current time (real scalar)
!  nbod   ==>  number of massive bodies (int scalar)
!  ntp    ==>  number of massive bodies (int scalar)
!  istat  ==>  status of the test paricles
!
! Remarks:
! Authors:  Hal Levison
! Date:    2/22/94
! Last revision:
use module_swift
use module_io
implicit none

! Inputs:
integer(ik) :: nbod,ntp,istat(ntp+1,NSTAT),iu
real(rk) :: time

! Internals
integer(ik) :: i
integer(integer2) :: nleft,nbod2
real(real4) :: ttmp

!-----------------
! Executable code
!-----------------

! Calculate number of remaining test particles
nleft = 0
do i = 1, ntp

  if(istat(i,1) == 0) nleft = nleft + 1

end do

nbod2 = nbod
ttmp = time

write(iu) ttmp,nbod2,nleft

return
end subroutine io_write_hdr_r
