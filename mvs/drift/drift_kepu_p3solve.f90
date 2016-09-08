subroutine drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflg)
!*************************************************************************
!                        DRIFT_KEPU_P3SOLVE.F
!*************************************************************************
! Returns the real root of cubic often found in solving kepler
! problem in universal variables.
!
!             Input:
!                 dt            ==>  time step (real scalar)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalar)
!                 mu            ==>  Reduced mass of system (real scalar)
!                 alpha         ==>  Twice the binding energy (real scalar)
!                 u             ==>  Vel. dot radial vector (real scalar)
!             Output:
!                 s             ==>  solution of cubic eqn for the
!                                    universal variable
!                 iflg          ==>  success flag ( = 0 if O.K.) (integer)
!
! Author:  Martin Duncan
! Date:    March 12/93
! Last revision: March 12/93
use module_swift
implicit none

! Inputs:
real(rk) :: dt, r0, mu, alpha, u

! Outputs:
integer(ik) :: iflg
real(rk) :: s

! Internals:
real(rk) :: denom, a0, a1, a2, q, r, sq2, sq, p1, p2

!----
! Executable code

denom = (mu - alpha*r0)/6.0_rk
a2 = 0.5_rk*u/denom
a1 = r0/denom
a0 = -dt/denom

q = (a1 - (a2**2)/3.0_rk)/3.0_rk
r = (a1*a2 - 3.0_rk*a0)/6.0_rk - (a2**3)/27.0_rk
sq2 = q**3 + r**2

if(sq2 >= 0.0_rk) then

  sq = sqrt(sq2)

  if((r + sq) <= 0.0_rk) then

    p1 = -(-(r + sq))**(1.0_rk/3.0_rk)

  else

    p1 = (r + sq)**(1.0_rk/3.0_rk)

  end if

  if((r - sq) <= 0.0_rk) then

    p2 =  -(-(r - sq))**(1.0_rk/3.0_rk)

  else

    p2 = (r - sq)**(1.0_rk/3.0_rk)

  end if

  iflg = 0
  s = p1 + p2 - a2/3.0_rk

else

  iflg = 1
  s = 0.0_rk

end if

return
end subroutine drift_kepu_p3solve