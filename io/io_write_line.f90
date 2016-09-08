subroutine io_write_line(iu,id,a,e,inc,capom,omega,capm)
!-------------------------------------------------------------------------
!			IO_WRITE_LINE.F90
!-------------------------------------------------------------------------
! NEW VERSION OF THIS, USES FXDR
!
! Input:
!  iu       ==> unit number to write to
!  a        ==> semi-major axis or pericentric distance if a parabola
!               (real scalar)
!  e        ==> eccentricity (real scalar)
!  inc      ==> inclination  (real scalar)
!  capom    ==> longitude of ascending node (real scalar)
!  omega    ==> argument of perihelion (real scalar)
!  capm     ==> mean anomoly(real scalar)
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
integer(ik) :: iu,id
real(rk) :: a,e,inc,capom,omega,capm

! Internals
integer(ik) :: ierr
real(real4) :: orbel(6)
!real(rk) :: orbel(6)

!-----------------
! Executable code
!-----------------

ierr = ixdrint(iu, id)
if(ierr > 0) then

  write(*,'(a)') ' SWIFT ERROR: in io_write_line: '
  write(*,'(a)') '  Could not write id'
  call util_exit(1)

end if

orbel(1) = a
orbel(2) = e
orbel(3) = inc
orbel(4) = capom
orbel(5) = omega
orbel(6) = capm

ierr = ixdrrmat( iu, 6, orbel )
!ierr = ixdrdmat( iu, 6, orbel )
if(ierr.gt.0) then

  write(*,'(a)') ' SWIFT ERROR: in io_write_line: '
  write(*,'(a)') '  Could not write orbit elements'
  call util_exit(1)

end if

return
end subroutine io_write_line
