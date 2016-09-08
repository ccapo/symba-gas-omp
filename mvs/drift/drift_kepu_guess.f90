subroutine drift_kepu_guess(dt, r0, mu, alpha, u, s)
!*************************************************************************
!                        DRIFT_KEPU_GUESS.F
!*************************************************************************
! Initial guess for solving kepler's equation using universal variables.
!
!             Input:
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 s             ==>  initial guess for the value of
!                                    universal variable
!
! Author:  Hal Levison & Martin Duncan
! Date:    3/12/93
! Last revision: April 6/93
use module_swift
use module_interfaces, except_this_one => drift_kepu_guess
implicit none

! Inputs:
real(rk) :: dt, r0, mu, alpha, u

! Inputs and Outputs:
real(rk) :: s

! Internals:
integer(ik) :: iflg
real(rk) :: y, sy, cy, sigma, es
real(rk) :: x, a
real(rk) :: en, ec, e

!----
! Executable code

if(alpha > 0.0_rk) then

  ! find initial guess for elliptic motion
  if(dt/r0 <= 0.4_rk)  then

    s = dt/r0 - (dt*dt*u)/(2.0_rk*r0*r0*r0)
    return

  else

    a = mu/alpha
    en = sqrt(mu/(a*a*a))
    ec = 1.0_rk - r0/a
    es = u/(en*a*a)
    e = sqrt(ec*ec + es*es)
    y = en*dt - es
    call orbel_scget(y, sy, cy)
    sigma = sign(1.0_rk, (es*cy + ec*sy))
    x = y + sigma*0.85_rk*e
    s = x/sqrt(alpha)

  end if

else

  ! find initial guess for hyperbolic motion
  call drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflg)
  if(iflg /= 0) s = dt/r0

end if

return
end subroutine drift_kepu_guess
