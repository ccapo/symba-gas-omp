module module_rmvs
! Module for the encounter subroutines
!
! Author:  Hal Levison
! Date:    3/7/93
! Last revision:
use module_swift
implicit none

! Scale factor for hill's sphere to take shorter time step
real(rk), parameter :: RHSCALE = 3.5_rk

! Scale factor for hill's sphere to go to planet centric coord.
real(rk), parameter :: RHPSCALE = 1.0_rk

! Ratio of the number of time steps in the encounter (helocentric) vs normal
integer(ik), parameter :: NTENC = 10

! Ratio of the number of time steps in the planetcentric encounter vs heliocentric
integer(ik), parameter :: NTPHENC = 3

! Ratio of the number of time steps in the encounter (planetcentric) vs normal
integer(ik), parameter :: NTPENC = NTENC*NTPHENC

end module module_rmvs