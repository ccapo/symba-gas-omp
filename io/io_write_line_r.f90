subroutine io_write_line_r(iu,id,a,e,inc,capom,omega,capm)
!-------------------------------------------------------------------------
!			IO_WRITE_LINE_R.F90
!-------------------------------------------------------------------------
! Write out one line to real*4 binary file.
!
! Input:
!  iu       ==> unit number to write to
!  a        ==> semi-major axis or pericentric distance if a parabola
!               real scalar)
!  e        ==> eccentricity (real scalar)
!  inc      ==> inclination  (real scalar)
!  capom    ==> longitude of ascending node (real scalar)
!  omega    ==> argument of perihelion (real scalar)
!  capm     ==> mean anomoly(real scalar)
!
! Remarks:
! Authors:  Hal Levison
! Date:    2/22/94
! Last revision:
use module_swift
use module_io
implicit none

! Inputs:
integer(ik) :: iu,id
real(rk) :: a,e,inc,capom,omega,capm

! Internals
integer(integer2) :: id2
real(real4) :: a4,e4,inc4,capom4,omega4,capm4

!-----------------
! Executable code
!-----------------

id2 = id

a4 = a
e4 = e
inc4 = inc
capom4 = capom
capm4 = capm
omega4 = omega

write(iu) id2,a4,e4,inc4,capom4,omega4,capm4

return
end subroutine io_write_line_r