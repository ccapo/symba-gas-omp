subroutine io_write_hdr(iu,time,nbod,ntp,istat)
!-------------------------------------------------------------------------
!			IO_WRITE_HDR.F90
!-------------------------------------------------------------------------
! NEW VERSION OF THIS, USES FXDR
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
! Date:    11/2/99
! Last revision:
use module_swift
use module_io
use module_fxdr
use module_interfaces, only: util_exit
implicit none

! Inputs:
integer(ik) :: nbod,ntp,istat(ntp+1,NSTAT),iu
real(rk) :: time

! Internals
integer(ik) :: i,ierr
integer(ik) :: nleft,nn(2)
real(real4) :: ttmp
!real(rk) :: ttmp

!-----------------
! Executable code
!-----------------

! Calculate number of remaining test particles
nleft = 0
do i = 1, ntp

  if(istat(i,1) == 0) nleft = nleft + 1

end do

ttmp = time

ierr = ixdrreal( iu, ttmp )
!ierr = ixdrdouble( iu, ttmp )
if(ierr > 0) then

  write(*,'(a)') ' SWIFT ERROR: in io_write_hdr: '
  write(*,'(a)') '  Could not write time'
  call util_exit(1)

end if

nn(1) = nbod
nn(2) = nleft
ierr = ixdrimat( iu, 2, nn )
if(ierr > 0) then

  write(*,'(a)') ' SWIFT ERROR: in io_write_hdr: '
  write(*,'(a)') '  Could not write nbod and nleft'
  call util_exit(1)

end if

return
end subroutine io_write_hdr
