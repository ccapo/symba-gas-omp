function io_read_radius_r(time,nbod,radius,iu) result(read_radius_r)
!-------------------------------------------------------------------------
!			IO_READ_RADIUS_R.F90
!-------------------------------------------------------------------------
! Read in the radius file.
!
! Output:
!  time  ==>  current time (real scalar)
!  nbod  ==>  number of radiusive bodies (int scalar)
!  radius  ==>  radius of bodies (real array)
!  iu    ==>  unit number to read from
!
! Returns:
!  read_radius  ==>  = 0 read ok
!                 != 0 read failed is set to iostat variable
!
! Remarks: Based on io_read_frame
! Authors:  Hal Levison
! Date:    1/9/97
! Last revision: 11/2/99
use module_swift
use module_io
implicit none

! Inputs:
integer(ik) :: iu

! Outputs
integer(ik) :: nbod, read_radius_r
real(rk) :: radius(nbod),time

! Internals
integer(ik) :: i,ierr
integer(integer2) :: nbod2
real(real4) :: radius4(NTPMAX),ttmp

!-----------------
! Executable code
!-----------------

read(iu, iostat = ierr) ttmp,nbod2
read_radius_r = ierr
if(ierr /= 0) return

read(iu, iostat = ierr) (radius4(i), i = 1, nbod2)
read_radius_r = ierr
if(ierr /= 0) return

do i = 1, nbod2

  radius(i) = radius4(i)

end do
nbod = nbod2
time = ttmp

return
end function io_read_radius_r
