subroutine drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg)
!*************************************************************************
!                        DRIFT_KEPU_LAG.F
!*************************************************************************
! subroutine for solving kepler's equation in universal variables.
! using LAGUERRE'S METHOD
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
!                 iflgn          ==>  =0 if converged; !=0 if not
!
! Author:  Hal Levison
! Date:    2/3/93
! Last revision: 4/21/93
use module_swift
use module_interfaces, except_this_one => drift_kepu_lag
implicit none

! Inputs:
real(rk) :: s, dt, r0, mu, alpha, u

! Outputs:
integer(ik) :: iflg
real(rk) :: fp, c1, c2, c3

! Internals:
!integer(ik), parameter :: NTMP = NLAG2 + 1
real(rk), parameter :: ln = 5.0_rk
integer(ik) :: nc, ncmax
real(rk) :: x, fpp, ds, c0, f
real(rk) :: fdt

!----
! Executable code

! To get close approach needed to take lots of iterations if alpha<0
if(alpha < 0.0_rk) then

  ncmax = NLAG2

else

  ncmax = NLAG1

end if

! start laguere's method
do nc = 0, ncmax

  x = s*s*alpha
  call drift_kepu_stumpff(x, c0, c1, c2, c3)
  c1 = c1*s
  c2 = c2*s*s
  c3 = c3*s*s*s
  f = r0*c1 + u*c2 + mu*c3 - dt
  fp = r0*c0 + u*c1 + mu*c2
  fpp = (mu - r0*alpha)*c1 + u*c0
  ds = -ln*f/(fp + sign(1.0_rk, fp)*sqrt(abs((ln - 1.0_rk)*(ln - 1.0_rk)*fp*fp - (ln - 1.0_rk)*ln*f*fpp)))
  s = s + ds

  fdt = f/dt

  ! quartic convergence
  if(fdt*fdt < DANBYB*DANBYB) then

    iflg = 0 ! Laguerre's method succeeded
    exit

  else

    iflg = 2 ! Laguerre's method failed

  end if

end do

return
end subroutine drift_kepu_lag
