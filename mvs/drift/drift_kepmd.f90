subroutine drift_kepmd(dm, es, ec, x, s, c)
!********************************************************************
!                  DRIFT_KEPMD
!********************************************************************
!  Subroutine for solving kepler's equation in difference form for an
!  ellipse, given SMALL dm and SMALL eccentricity.  See DRIFT_DAN.F
!  for the criteria.
!  WARNING - BUILT FOR SPEED : DOES NOT CHECK HOW WELL THE ORIGINAL
!  EQUATION IS SOLVED! (CAN DO THAT IN THE CALLING ROUTINE BY
!  CHECKING HOW CLOSE (x - ec*s +es*(1.-c) - dm) IS TO ZERO.
!
!	Input:
!	    dm		==> increment in mean anomaly M (real(rk) :: scalar)
!	    es, ec       ==> ecc. times sin and cos of E_0 (real(rk) :: scalars)
!
!       Output:
!            x          ==> solution to Kepler's difference eqn (real(rk) :: scalar)
!            s, c        ==> sin and cosine of x (real(rk) :: scalars)
use module_swift
implicit none

! Inputs
real(rk) :: dm, es, ec

! Outputs
real(rk) :: x, s, c

! Internals
real(rk), parameter :: A0 = 39916800.0_rk, A1 = 6652800.0_rk, A2 = 332640.0_rk, A3 = 7920.0_rk,  A4 = 110.0_rk
real(rk) :: dx
real(rk) :: fac1, fac2, q, y
real(rk) :: f, fp, fpp, fppp

! calc initial guess for root
fac1 = 1.0_rk/(1.0_rk - ec)
q = fac1*dm
fac2 = es*es*fac1 - ec/3.0_rk
x = q*(1.0_rk - 0.5_rk*fac1*q*(es - q*fac2))

! excellent approx. to sin and cos of x for small x.
y = x*x
s = x*(A0 - y*(A1 - y*(A2 - y*(A3 - y*(A4 - y)))))/A0
c = sqrt(1.0_rk - s*s)

! Compute better value for the root using quartic Newton method
f = x - ec*s + es*(1.0_rk - c) - dm
fp = 1.0_rk - ec*c + es*s
fpp = ec*s + es*c
fppp = ec*c - es*s
dx = -f/fp
dx = -f/(fp + 0.5_rk*dx*fpp)
dx = -f/(fp + 0.5_rk*dx*fpp + dx*dx*fppp/6.0_rk)
x = x + dx

! excellent approx. to sin and cos of x for small x.
y = x*x
s = x*(A0 - y*(A1 - y*(A2 - y*(A3 - y*(A4 - y)))))/A0
c = sqrt(1.0_rk - s*s)

return
end subroutine drift_kepmd