function util_power_law(p, xmin, xmax) result(x)
!------------------------------------------------------------------
!    UTIL_POWER_LAW.F90
!------------------------------------------------------------------
! Returns a random number drawn from a power-law distribution
! with exponent p, where xmin and xmax are the minimum and
! maximum range for the random number.
!
! Input:  p    ==> Exponent of power-law distribution
!         xmin ==> Minimum value for random number
!         xmax ==> Maximum value for random number
!
! Output: x    ==> Random number drawn from power-law distrubution
!
! By: Chris Capobianco
! Date: 02/04/09
!------------------------------------------------------------------
use module_swift
use module_interfaces, only: util_randomu
implicit none

! Input variables
real(rk) :: p, xmin, xmax

! Output variable
real(rk) :: x

! Internal variables
real(rk) :: p1, ip1, xi

p1 = 1.0_rk - p
ip1 = 1.0_rk/p1
xi = util_randomu()

if(p /= 1.0_rk) then

  x = (xmin**p1 + xi*(xmax**p1 - xmin**p1))**ip1

else

  x = xmin*(xmax/xmin)**xi

end if

return
end function util_power_law
