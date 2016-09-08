function io_read_mass_r(time,nbod,mass,iu) result(read_mass_r)
!-------------------------------------------------------------------------
!			IO_READ_MASS_R.F90
!-------------------------------------------------------------------------
! Read in the mass file.
!
! Output:
!  time  ==>  current time (real scalar)
!  nbod  ==>  number of massive bodies (int scalar)
!  mass  ==>  mass of bodies (real array)
!  iu    ==>  unit number to read from
!
! Returns:
!  read_mass  ==>  = 0 read ok
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
integer(ik) :: nbod, read_mass_r
real(rk) :: mass(nbod),time

! Internals
integer(ik) :: i,ierr
integer(integer2) :: nbod2
real(real4) :: mass4(NTPMAX),ttmp

!-----------------
! Executable code
!-----------------

read(iu, iostat = ierr) ttmp,nbod2
read_mass_r = ierr
if(ierr /= 0) return

read(iu, iostat = ierr) (mass4(i), i = 1, nbod2)
read_mass_r = ierr
if(ierr /= 0) return

do i = 1, nbod2

  mass(i) = mass4(i)

end do
nbod = nbod2
time = ttmp

return
end function io_read_mass_r