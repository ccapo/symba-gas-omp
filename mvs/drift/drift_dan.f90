subroutine drift_dan(mu, x0, y0, z0, vx0, vy0, vz0, dt0, iflg)
!*************************************************************************
!                        DRIFT_DAN.F
!*************************************************************************
! This subroutine does the Danby and decides which vbles to use
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 x0, y0, z0         ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 vx0, vy0, vz0      ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 dt0            ==>  time step
!             Output:
!                 x0, y0, z0         ==>  final position in jacobi coord 
!                                       (real scalars)
!                 vx0, vy0, vz0      ==>  final position in jacobi coord 
!                                       (real scalars)
!                 iflg             ==>  integer flag (zero if satisfactory)
!					      (non-zero if nonconvergence)
!
! Authors:  Hal Levison & Martin Duncan  
! Date:    2/10/93
! Last revision: April 6/93 - MD adds dt and keeps dt0 unchanged
use module_swift
use module_interfaces, except_this_one => drift_dan
implicit none

! Inputs Only: 
real(rk) :: mu, dt0

! Inputs and Outputs:
real(rk) :: x0, y0, z0
real(rk) :: vx0, vy0, vz0

! Output
integer(ik) :: iflg

! Internals:
real(rk) :: x, y, z, vx, vy, vz, dt
real(rk) :: f, g, fdot, c1, c2
real(rk) :: c3, gdot
real(rk) :: u, alpha, fp, r0, v0s
real(rk) :: a, asq, en
real(rk) :: dm, ec, es, esq, xkep
real(rk) :: fchk, s, c

!----
! Executable code 

! Set dt = dt0 to be sure timestep is not altered while solving for new coords.
dt = dt0
iflg = 0
r0 = sqrt(x0**2 + y0**2 + z0**2)
v0s = vx0**2 + vy0**2 + vz0**2
u = x0*vx0 + y0*vy0 + z0*vz0
alpha = 2.0_rk*mu/r0 - v0s

if(alpha > 0.0_rk) then

  a = mu/alpha
  asq = a*a
  en = sqrt(mu/(a*asq))
  ec = 1.0_rk - r0/a
  es = u/(en*asq)
  esq = ec*ec + es*es
  dm = dt*en - TWOPI*int((dt*en)/TWOPI)
  dt = dm/en

  !if((esq <= 0.36_rk) .and. (dm*dm <= 0.16_rk) .and. (esq*dm*dm < 0.0016_rk)) then

  if((dm*dm > 0.16_rk) .or. (esq > 0.36_rk)) goto 100

  if(esq*dm*dm < 0.0016_rk) then

    call drift_kepmd(dm, es, ec, xkep, s, c)
    fchk = (xkep - ec*s + es*(1.0_rk - c) - dm)

    if(fchk*fchk > DANBYB*DANBYB) then

      iflg = 1
      return

    end if

    fp = 1.0_rk - ec*c + es*s
    f = (a/r0)*(c - 1.0_rk) + 1.0_rk
    g = dt + (s - xkep)/en
    fdot = -(a/(r0*fp))*en*s
    gdot = (c - 1.0_rk)/fp + 1.0_rk

    x = x0*f + vx0*g
    y = y0*f + vy0*g
    z = z0*f + vz0*g
    vx = x0*fdot + vx0*gdot
    vy = y0*fdot + vy0*gdot
    vz = z0*fdot + vz0*gdot

    x0 = x
    y0 = y
    z0 = z
    vx0 = vx
    vy0 = vy
    vz0 = vz

    iflg = 0
    return

  end if

end if

100   call drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg)

if(iflg == 0) then

  f = 1.0_rk - (mu/r0)*c2
  g = dt - mu*c3
  fdot = -(mu/(fp*r0))*c1
  gdot = 1.0_rk - (mu/fp)*c2

  x = x0*f + vx0*g
  y = y0*f + vy0*g
  z = z0*f + vz0*g
  vx = x0*fdot + vx0*gdot
  vy = y0*fdot + vy0*gdot
  vz = z0*fdot + vz0*gdot

  x0 = x
  y0 = y
  z0 = z
  vx0 = vx
  vy0 = vy
  vz0 = vz

end if

return
end subroutine drift_dan
