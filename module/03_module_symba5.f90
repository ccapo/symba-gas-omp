module module_symba5
! Module for the SyMBA5 subroutines
!
! Remarks:  Copied from symba5.inc
! Author:  Hal Levison
! Date:    3/20/97
! Last revision:
use module_swift
implicit none

! Maximum number of encounters
integer(ik), parameter :: NENMAX = 262144 ! must be less than 2^18 = 262144

! Ratio of the number of time steps in the adjoining shells
integer(ik), parameter :: NTENC = 3

! Scale factor for hill's sphere to take shorter time step
real(rk), parameter :: RHSCALE = 6.5_rk

! Ratio of shell radii squared
!real(rk), parameter :: rshell = 0.48075_rk ! rshell ~ ntenc^(-2/3)
!real(rk), parameter :: rshell = real(NTENC, rk)**(-2.0_rk/3.0_rk) ! rshell ~ ntenc^(-2/3)
real(rk), parameter :: RSHELL = 0.480749856769136133_rk ! rshell ~ ntenc^(-2/3)

! Maximum number of perigee passages about a massive body before merging satellite with massive body
! This applies to those particles that become bound (however temporarily) during a close encounter,
! and limits the damage that such trapped particles do to the overall performance of the code.
integer(ik), parameter :: NBOUNDMAX = 5

end module module_symba5
