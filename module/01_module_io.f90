module module_io
! Module for io routines
!
! Author:  Hal Levison
! Date:    2/21/94
! Last revision: 2/24/94
use module_swift
implicit none

! Number of bytes in iflgchk
integer(ik), parameter :: IO_NBITS = 6

! Bit 0 set ==> write big binary data file
! Bit 1 set ==> write real*4 binary file rather than int*2: ignored if bit0 == F
! Bit 2 set ==> calc energy of system wrt time
! Bit 3 set ==> calc jacobi of the test particles
! Bit 4 set ==> check if particles are removed
! Bit 5 set ==> include J2 and J4 terms
! Bit 6 set ==> symba5

end module module_io