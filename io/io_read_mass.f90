function io_read_mass(time,nbod,mass,iu) result(read_mass)
!-------------------------------------------------------------------------
!				IO_READ_MASS.F90
!-------------------------------------------------------------------------
! NEW VERSION OF THIS, USES FXDR
! Read in the mass file.
!
! Output:
!  time  ==>  current time (real scalar)
!  nbod  ==>  number of massive bodies (int scalar)
!  mass  ==>  mass of bodies (real array)
!  iu    ==>  unit number to read from
!
! Returns:
!  io_read_mass ==> = 0 read ok
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
integer(ik) :: nbod, read_mass
real(rk) :: mass(nbod),time

! Internals
integer(ik) :: i,ierr
real(real4) :: mass4(NTPMAX), ttmp

!-----------------
! Executable code
!-----------------

ierr = ixdrreal( iu, ttmp )
!ierr = ixdrdouble( iu, time )
read_mass = ierr
if(ierr /= 0) return
time = ttmp

ierr = ixdrint( iu, nbod )
read_mass = ierr
if(ierr /= 0) return

ierr = ixdrrmat( iu, nbod, mass4 )
!ierr = ixdrdmat( iu, nbod, mass )
read_mass = ierr
if(ierr /= 0) return

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, mass, mass4)
do i = 1, nbod

  mass(i) = mass4(i)

end do
!$OMP END PARALLEL DO

return
end function io_read_mass
