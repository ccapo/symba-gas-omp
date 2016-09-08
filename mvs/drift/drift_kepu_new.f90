subroutine drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg)
!*************************************************************************
!                        DRIFT_KEPU_NEW.F
!*************************************************************************
! subroutine for solving kepler's equation in universal variables.
! using NEWTON'S METHOD
!
!             Input:
!                 s             ==>  inital value of universal variable
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 s             ==>  final value of universal variable
!                 fp            ==>  f' from p170
!                                       (real scalors)
!                 c1, c2, c3      ==>  c's from p171-172
!                                       (real scalors)
!                 iflg          ==>  =0 if converged; !=0 if not
!
! Author:  Hal Levison
! Date:    2/3/93
! Last revision: 4/21/93
use module_swift
use module_interfaces, except_this_one => drift_kepu_new
implicit none

! Inputs:
real(rk) :: s, dt, r0, mu, alpha, u

! Outputs:
integer(ik) :: iflg
real(rk) :: fp, c1, c2, c3

! Internals:
integer(ik) :: nc
real(rk) :: x, c0, ds
real(rk) :: f, fpp, fppp, fdt

!----
! Executable code

do nc = 0, 6

  x = s*s*alpha
  call drift_kepu_stumpff(x, c0, c1, c2, c3)
  c1 = c1*s
  c2 = c2*s*s
  c3 = c3*s*s*s
  f = r0*c1 + u*c2 + mu*c3 - dt
  fp = r0*c0 + u*c1 + mu*c2
  fpp = (mu - r0*alpha)*c1 + u*c0
  fppp = (mu - r0*alpha)*c0 - u*alpha*c1
  ds = - f/fp
  ds = - f/(fp + 0.5_rk*ds*fpp)
  ds = -f/(fp + 0.5_rk*ds*fpp + ds*ds*fppp/6.0_rk)
  s = s + ds
  fdt = f/dt

  ! quartic convergence
  if(fdt*fdt < DANBYB*DANBYB) then

    iflg = 0 ! Newton's method succeeded
    exit

  else

    iflg = 1 ! Newton's method failed

  end if

end do

return
end subroutine drift_kepu_new
