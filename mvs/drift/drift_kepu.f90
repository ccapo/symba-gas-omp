subroutine drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg)
!*************************************************************************
!                        DRIFT_KEPU.F
!*************************************************************************
! subroutine for solving kepler's equation using universal variables.
!
!             Input:
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 fp            ==>  f' from p170
!                                       (real scalors)
!                 c1, c2, c3      ==>  c's from p171-172
!                                       (real scalors)
!                 iflg          ==>  =0 if converged; !=0 if not
!
! Author:  Hal Levison
! Date:    2/3/93
! Last revision: 2/3/93
use module_swift
use module_interfaces, except_this_one => drift_kepu
implicit none

! Inputs:
real(rk) :: dt, r0, mu, alpha, u

! Outputs:
integer(ik) :: iflg
real(rk) :: fp, c1, c2, c3

! Internals:
real(rk) :: s, st, fo, fn

!----
! Executable code

call drift_kepu_guess(dt, r0, mu, alpha, u, s)

! store initial guess for possible use later in
! laguerre's method,  in case newton's method fails.
st = s

call drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg)

if(iflg /= 0) then

  call drift_kepu_fchk(dt, r0, mu, alpha, u, st, fo)
  call drift_kepu_fchk(dt, r0, mu, alpha, u, s, fn)
  if(abs(fo) < abs(fn)) s = st

  call drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg)

endif

return
end subroutine drift_kepu
