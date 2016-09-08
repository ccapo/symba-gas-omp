function io_read_line_r(iu,id,a,e,inc,capom,omega,capm) result(read_line_r)
!-------------------------------------------------------------------------
!			IO_READ_LINE_R.F90
!-------------------------------------------------------------------------
! Read one line from real*4 binary file.
!
! Input:
!   iu     ==> unit number to write to
!
! Output:
!   a      ==> semi-major axis or pericentric distance if a parabola
!              (real scalar)
!   e      ==> eccentricity (real scalar)
!   inc    ==> inclination  (real scalar)
!   capom  ==> longitude of ascending node (real scalar)
!   omega  ==> argument of perihelion (real scalar)
!   capm   ==> mean anomoly(real scalar)
!
! Returns:
!   read_line_r ==> = 0 read ok
!                  != 0 read failed is set to iostat variable
!
! Remarks:
! Authors:  Hal Levison
! Date:    2/22/94
! Last revision:
use module_swift
use module_io
implicit none

! Inputs:
integer(ik) :: iu

! Output:
integer(ik) :: id, read_line_r
real(rk) :: a,e,inc,capom,omega,capm

! Internals
integer(ik) :: ierr
integer(integer2) :: id2
real(real4) :: a4,e4,inc4,capom4,omega4,capm4

!-----------------
! Executable code
!-----------------

read(iu, iostat = ierr) id2, a4, e4, inc4, capom4, omega4, capm4
read_line_r = ierr
if(ierr /= 0) return

id = id2
a = a4
e = e4
inc = inc4
capom = capom4
capm = capm4
omega = omega4

return
end function io_read_line_r