function util_randomu() result(ran)
!--------------------------------------------------------------
!    UTIL_RANDOMU.F90
!--------------------------------------------------------------
! Generates a uniform random number between (0.0,1.0) using
! Intel Fortran 90 random number generator
!--------------------------------------------------------------
use module_swift
implicit none

! Output variable
real(rk) :: ran

call random_number(ran)

return
end function util_randomu
