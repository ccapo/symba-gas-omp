function io_read_radius(time,nbod,radius,iu) result(read_radius)
!-------------------------------------------------------------------------
!				IO_READ_RADIUS.F90
!-------------------------------------------------------------------------
! NEW VERSION OF THIS, USES FXDR
! Read in the radius file.
!
! Output:
!  time  ==>  current time (real scalar)
!  nbod  ==>  number of radiusive bodies (int scalar)
!  radius  ==>  radius of bodies (real array)
!  iu    ==>  unit number to read from
!
! Returns:
!  io_read_radius ==> = 0 read ok
!                  != 0 read failed is set to iostat variable
!
! Remarks: Based on io_read_frame
! Authors:  Hal Levison
! Date:    11/2/99
! Last revision:
use module_swift
use module_io
use module_fxdr
implicit none

! Inputs:
integer(ik) :: iu

! Outputs
integer(ik) :: nbod, read_radius
real(rk) :: radius(nbod),time

! Internals
integer(ik) :: i,ierr
real(real4) :: radius4(NTPMAX), ttmp

!-----------------
! Executable code
!-----------------

ierr = ixdrreal( iu, ttmp )
!ierr = ixdrdouble( iu, time )
read_radius = ierr
if(ierr /= 0) return
time = ttmp

ierr = ixdrint( iu, nbod )
read_radius = ierr
if(ierr /= 0) return

ierr = ixdrrmat( iu, nbod, radius4 )
!ierr = ixdrdmat( iu, nbod, radius )
read_radius = ierr
if(ierr /= 0) return

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, radius, radius4)
do i = 1, nbod

  radius(i) = radius4(i)

end do
!$OMP END PARALLEL DO

return
end function io_read_radius
