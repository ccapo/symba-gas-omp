function util_rayleigh(rms, xmin, xmax) result(x)
!--------------------------------------------------------------
!    UTIL_RAYLEIGH.F90
!--------------------------------------------------------------
! Returns a random number drawn from a rayleigh distribution
! with dispersion rms, where xmin and xmax are the minimum and
! maximum range for the random number.
!--------------------------------------------------------------
use module_swift
use module_interfaces, only: util_randomu
implicit none

! Input variables
real(rk) :: rms, xmin, xmax

! Output variable
real(rk) :: x

! Internal variables
real(rk) :: fmin, fmax, xi

fmin = -exp(-0.5_rk*(xmin/rms)**2)
fmax = -exp(-0.5_rk*(xmax/rms)**2)
xi = util_randomu()

x = rms*sqrt(-2.0_rk*log((fmin - fmax)*xi - fmin))

return
end function util_rayleigh
