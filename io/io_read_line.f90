function io_read_line(iu,id,a,e,inc,capom,omega,capm) result(read_line)
!-------------------------------------------------------------------------
!				IO_READ_LINE.F90
!-------------------------------------------------------------------------
! NEW VERSION OF THIS, USES FXDR
!
!      Input:
!            iu       ==> unit number to write to
!      Output:
!            id       ==> Particle id number (int scalar)
!            a        ==> semi-major axis or pericentric distance if a parabola
!                          (real scalar)
!            e        ==> eccentricity (real scalar)
!            inc      ==> inclination  (real scalar)
!            capom    ==> longitude of ascending node (real scalar)
!            omega    ==> argument of perihelion (real scalar)
!            capm     ==> mean anomoly(real scalar)
!       Returns:
!      io_read_line    ==>   =0 read ok
!                           !=0 read failed is set to iostat variable
!
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

! Output:
integer(ik) :: id, read_line
real(rk) :: a,e,inc,capom,omega,capm

! Internals
integer(ik) :: ierr
real(real4) :: orbel(6)
!real(rk) :: orbel(6)

!-----------------
! Executable code
!-----------------

ierr = ixdrint(iu, id)
read_line = ierr
if(ierr /= 0) return

ierr = ixdrrmat( iu, 6, orbel )
!ierr = ixdrdmat( iu, 6, orbel )
read_line = ierr
if(ierr /= 0) return

a = orbel(1)
e = orbel(2)
inc = orbel(3)
capom = orbel(4)
omega = orbel(5)
capm = orbel(6)

return
end function io_read_line
