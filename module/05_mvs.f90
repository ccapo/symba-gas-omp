module mvs
! Module for routines in mvs directory
use swift
use orbel, only: orbel_scget
implicit none

contains
!!
!!
  subroutine drift_dan(mu, x0, y0, z0, vx0, vy0, vz0, dt0, iflg)
  !-------------------------------------------------------------------------
  !                        DRIFT_DAN.F
  !-------------------------------------------------------------------------
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
!!
  subroutine drift_kepmd(dm, es, ec, x, s, c)
  !--------------------------------------------------------------------
  !                  DRIFT_KEPMD
  !--------------------------------------------------------------------
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
!!
  subroutine drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg)
  !-------------------------------------------------------------------------
  !                        DRIFT_KEPU.F
  !-------------------------------------------------------------------------
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
!!
  subroutine drift_kepu_fchk(dt, r0, mu, alpha, u, s, f)
  !-------------------------------------------------------------------------
  !                        DRIFT_KEPU_FCHK.F
  !-------------------------------------------------------------------------
  ! Returns the value of the function f of which we are trying to find the root
  ! in universal variables.
  !
  !             Input:
  !                 dt            ==>  time step (real scalar)
  !                 r0            ==>  Distance between `Sun' and particle
  !                                     (real scalar)
  !                 mu            ==>  Reduced mass of system (real scalar)
  !                 alpha         ==>  Twice the binding energy (real scalar)
  !                 u             ==>  Vel. dot radial vector (real scalar)
  !                 s             ==>  Approx. root of f
  !             Output:
  !                 f             ==>  function value ( = 0 if O.K.) (integer)
  !
  ! Author:  Martin Duncan
  ! Date:    March 12/93
  ! Last revision: March 12/93
  implicit none

  ! Inputs:
  real(rk) :: dt, r0, mu, alpha, u, s

  ! Outputs:
  real(rk) :: f

  ! Internals:
  real(rk) :: x, c0, c1, c2, c3

  !----
  ! Executable code

  x = s*s*alpha
  call drift_kepu_stumpff(x, c0, c1, c2, c3)
  c1 = c1*s
  c2 = c2*s*s
  c3 = c3*s*s*s
  f = r0*c1 + u*c2 + mu*c3 - dt

  return
  end subroutine drift_kepu_fchk
!!
  subroutine drift_kepu_guess(dt, r0, mu, alpha, u, s)
  !-------------------------------------------------------------------------
  !                        DRIFT_KEPU_GUESS.F
  !-------------------------------------------------------------------------
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
!!
  subroutine drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg)
  !-------------------------------------------------------------------------
  !                        DRIFT_KEPU_LAG.F
  !-------------------------------------------------------------------------
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
!!
  subroutine drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg)
  !-------------------------------------------------------------------------
  !                        DRIFT_KEPU_NEW.F
  !-------------------------------------------------------------------------
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
!!
  subroutine drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflg)
  !-------------------------------------------------------------------------
  !                        DRIFT_KEPU_P3SOLVE.F
  !-------------------------------------------------------------------------
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
!!
  subroutine drift_kepu_stumpff(x, c0, c1, c2, c3)
  !-------------------------------------------------------------------------
  !                        DRIFT_KEPU_STUMPFF.F
  !-------------------------------------------------------------------------
  ! subroutine for the calculation of stumpff functions
  ! see Danby p.172  equations 6.9.15
  !
  !             Input:
  !                 x             ==>  argument
  !             Output:
  !                 c0, c1, c2, c3   ==>  c's from p171-172
  !                                       (real scalors)
  ! Author:  Hal Levison
  ! Date:    2/3/93
  ! Last revision: 2/3/93
  implicit none

  ! Inputs:
  real(rk) :: x

  ! Outputs:
  real(rk) :: c0, c1, c2, c3

  ! Internals:
  real(rk), parameter :: xm = 0.1_rk
  integer(ik) :: n, i

  !----
  ! Executable code

  n = 0
  do while(abs(x) >= xm)

    n = n + 1
    x = x/4.0_rk

  end do

  c2 = (1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x/182.0_rk)/132.0_rk)/90.0_rk)/56.0_rk)/30.0_rk) &
       /12.0_rk)/2.0_rk
  c3 = (1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x/210.0_rk)/156.0_rk)/110.0_rk)/72.0_rk)/42.0_rk) &
       /20.0_rk)/6.0_rk
  c1 = 1.0_rk - x*c3
  c0 = 1.0_rk - x*c2

  if(n /= 0) then

    do i = n, 1, -1

      c3 = (c2 + c0*c3)/4.0_rk
      c2 = (c1**2)/2.0_rk
      c1 = c0*c1
      c0 = 2.0_rk*c0**2 - 1.0_rk
      x = 4.0_rk*x

    end do

  end if

  return
  end subroutine drift_kepu_stumpff
!!
  subroutine drift_one(mu, x, y, z, vx, vy, vz, dt, iflg)
  !-------------------------------------------------------------------------
  !                        DRIFT_ONE.F
  !-------------------------------------------------------------------------
  ! This subroutine does the danby-type drift for one particle,  using
  ! appropriate vbles and redoing a drift if the accuracy is too poor
  ! (as flagged by the integer iflg).
  !
  !             Input:
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 mu            ==>  mass of central body (real scalar)
  !                 x, y, z         ==>  initial position in jacobi coord
  !                                    (real scalar)
  !                 vx, vy, vz      ==>  initial position in jacobi coord
  !                                    (real scalar)
  !                 dt            ==>  time step
  !             Output:
  !                 x, y, z         ==>  final position in jacobi coord
  !                                       (real scalars)
  !                 vx, vy, vz      ==>  final position in jacobi coord
  !                                       (real scalars)
  !                 iflg          ==>  integer (zero for successful step)
  !
  ! Authors:  Hal Levison & Martin Duncan
  ! Date:    2/10/93
  ! Last revision: 2/10/93
  implicit none

  ! Inputs Only:
  real(rk) :: mu, dt

  ! Inputs and Outputs:
  real(rk) :: x, y, z
  real(rk) :: vx, vy, vz

  ! Output
  integer(ik) :: iflg

  ! Internals:
  integer(ik) :: i
  real(rk) :: dttmp

  !----
  ! Executable code

  call drift_dan(mu, x, y, z, vx, vy, vz, dt, iflg)

  if(iflg /= 0) then

    dttmp = 0.1_rk*dt

    do i = 1, 10

      call drift_dan(mu, x, y, z, vx, vy, vz, dttmp, iflg)
      if(iflg /= 0) return ! Abandon all hope ye who enter here

    end do

  end if

  return
  end subroutine drift_one
!!
  subroutine getacch_ir3(nbod, istart, x, y, z, ir3, ir)
  !-------------------------------------------------------------------------
  !                        GETACCH_IR3.F
  !-------------------------------------------------------------------------
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
!!
  subroutine kickvh(nbod, vxh, vyh, vzh, axh, ayh, azh, dt)
  !-------------------------------------------------------------------------
  !                        KICKVH.F
  !-------------------------------------------------------------------------
  ! To kick the velocity components vxh(*) by axh(*)*dt
  !
  !             Input:
  !                 nbod          ==>  number of bodies (int scalar)
  !                 vxh, vyh, vzh   ==>  initial velocity in helio coord
  !                                    (real arrays)
  !                 axh, ayh, azh   ==>  acceleration in helio coord
  !                                    (real arrays)
  !                 dt            ==>  time step
  !             Output:
  !                 vxh, vyh, vzh   ==>  final velocity in helio coord
  !                                    (real arrays)
  !
  !     ALGORITHM: Obvious
  !     REMARKS:  Only alters particles 2 thru nbod since Sun is #1
  !
  !     AUTHOR:  M. Duncan.
  !     DATE WRITTEN:  Feb. 2,  1993.
  !     REVISIONS: 2/18/93   HFL
  implicit none

  ! Inputs Only:
  integer(ik) :: nbod
  real(rk) :: axh(nbod), ayh(nbod), azh(nbod)
  real(rk) :: dt

  ! Inputs and Output:
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

  ! Internals:
  integer(ik) :: i

  !----
  ! Executable code

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(dt) &
  !$OMP SHARED(nbod, vxh, vyh, vzh, axh, ayh, azh)
  do i = 2, nbod

    vxh(i) = vxh(i) + axh(i)*dt
    vyh(i) = vyh(i) + ayh(i)*dt
    vzh(i) = vzh(i) + azh(i)*dt

  end do
  !$OMP END PARALLEL DO

  return
  end subroutine kickvh
!!
!!
end module mvs