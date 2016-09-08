subroutine getacch_ir3(nbod, istart, x, y, z, ir3, ir)
!*************************************************************************
!                        GETACCH_IR3.F
!*************************************************************************
! Calculate r^-3 for an array of particles
!             Input:
!                 nbod     ==>  number of massive bodies (int scalor)
!                istart    ==>  body to start with (int scalor)
!                 x, y, z    ==>  positions (real arrays)
!             Output:
!                 ir3       ==>  r^-3  (real array)
!                 ir        ==>  r^-1  (real array)
!
! Author:  Hal Levison  
! Date:    2/2/93
! Last revision: 2/24/94
use module_swift
implicit none

! Inputs: 
integer(ik) :: nbod, istart
real(rk) :: x(nbod), y(nbod), z(nbod)

! Outputs:
real(rk) :: ir3(nbod)
real(rk) :: ir(nbod)

! Internals:
integer(ik) :: i
real(rk) :: r2

!----
! Executable code

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, r2) SHARED(nbod, istart, x, y, z, ir, ir3)
do i = istart, nbod

  r2 = x(i)**2 + y(i)**2 + z(i)**2
  ir(i) = 1.0_rk/sqrt(r2)
  ir3(i) = ir(i)/r2

end do
!$OMP END PARALLEL DO

return
end subroutine getacch_ir3