subroutine orbel_xv2el(x, y, z, vx, vy, vz, gmsum, ialpha, a, e, inc, capom, omega, capm)
!**********************************************************************
!	               ORBEL_XV2EL.F
!**********************************************************************
!     PURPOSE:  Given the cartesian position and velocity of an orbit, 
!       compute the osculating orbital elements.
!
!       input:
!            x, y, z    ==>  position of object (real scalars)
!            vx, vy, vz ==>  velocity of object (real scalars)
!            gmsum       ==> G*(M1+M2) (real scalar)
!
!       Output:
!	     ialpha   ==> conic section type ( see PURPOSE,  integer scalar)
!	     a        ==> semi-major axis or pericentric distance if a parabola
!                          (real scalar)
!            e        ==> eccentricity (real scalar)
!            inc      ==> inclination  (real scalar)
!            capom    ==> longitude of ascending node (real scalar)
!	     omega    ==> argument of perihelion (real scalar)
!	     capm     ==> mean anomoly(real scalar)
!
!     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech."
!     REMARKS:  If the inclination INC is less than TINY,  we
!       arbitrarily choose the longitude of the ascending node LGNODE
!       to be 0.0 (so the ascending node is then along the X axis).  If
!       the  eccentricity E is less than SQRT(TINY),  we arbitrarily
!       choose the argument of perihelion to be 0.
!     AUTHOR:  M. Duncan.
!     DATE WRITTEN:  May 8, 1992.
!     REVISIONS: 2/4/2K
use module_swift
implicit none

! Inputs Only:
real(rk) :: x, y, z, vx, vy, vz, gmsum

! Outputs
integer(ik) :: ialpha
real(rk) :: a, e, inc, capom, omega, capm

! Internals:
real(rk) :: hx, hy, hz, h2, h, r, v2, v, vdotr, energy, fac, face, cape, capf, tmpf
real(rk) :: cw, sw, w, u

!----
! Executable code

! Compute the angular momentum H, and thereby the inclination INC.
hx = y*vz - z*vy
hy = z*vx - x*vz
hz = x*vy - y*vx
h2 = hx*hx + hy*hy +hz*hz
h = sqrt(h2)
if(hz > h) hz = h ! Hal's fix
inc = acos(hz/h)

! Compute longitude of ascending node CAPOM and the argument of latitude u.
fac = sqrt(hx**2 + hy**2)/h

if(fac < TINY) then

  capom = 0.0_rk
  u = atan2(y, x)
  if(abs(inc - PI) < 10.0_rk*TINY) u = -u

else

  capom = atan2(hx, -hy)
  u = atan2(z/sin(inc), x*cos(capom) + y*sin(capom))

end if

if(capom < 0.0_rk) capom = capom + TWOPI
if(u < 0.0_rk) u = u + TWOPI

! Compute the radius R and velocity squared V2,  and the dot
! product RDOTV, the energy per unit mass ENERGY.
r = sqrt(x*x + y*y + z*z)
v2 = vx*vx + vy*vy + vz*vz
v = sqrt(v2)
vdotr = x*vx + y*vy + z*vz
energy = 0.5_rk*v2 - gmsum/r

! Determine type of conic section and label it via IALPHA
if(abs(energy*r/gmsum) < sqrt(TINY)) then

  ialpha = 0

else

  if(energy < 0.0_rk) ialpha = -1
  if(energy > 0.0_rk) ialpha = 1

end if

! Depending on the conic type,  determine the remaining elements

! ELLIPSE
if(ialpha == -1) then

  a = -0.5_rk*gmsum/energy
  fac = 1.0_rk - h2/(gmsum*a)

  if(fac > TINY) then

    e = sqrt(fac)
    face = (a - r)/(a*e)

    ! Apr. 16/93 : watch for case where face is slightly outside unity
    if(face > 1.0_rk) then

      cape = 0.0_rk

    else

      if(face > -1.0_rk) then

        cape = acos(face)

      else

        cape = PI

      end if

    end if

    if(vdotr < 0.0_rk) cape = TWOPI - cape
    cw = (cos(cape) - e)/(1.0_rk - e*cos(cape))
    sw = sqrt(1.0_rk - e*e)*sin(cape)/(1.0_rk - e*cos(cape))
    w = atan2(sw, cw)
    if(w < 0.0_rk) w = w + TWOPI

  else

    e = 0.0_rk
    w = u
    cape = u

  end if

  capm = cape - e*sin(cape)
  omega = u - w
  if(omega < 0.0_rk) omega = omega + TWOPI
  omega = omega - TWOPI*int(omega/TWOPI)

end if

! HYPERBOLA:
if(ialpha == 1) then

  a = 0.5_rk*gmsum/energy
  fac = h2/(gmsum*a)

  if(fac > TINY) then

    e = sqrt(1.0_rk + fac)
    tmpf = (a + r)/(a*e)
    if(tmpf < 1.0_rk) tmpf = 1.0_rk
    capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_rk))
    if(vdotr < 0.0_rk) capf = -capf
    cw = (e - cosh(capf))/(e*cosh(capf) - 1.0_rk)
    sw = sqrt(e*e - 1.0_rk)*sinh(capf)/(e*cosh(capf) - 1.0_rk)
    w = atan2(sw, cw)
    if(w < 0.0_rk) w = w + TWOPI

  else

    ! we only get here if a hyperbola is essentially a parabola
    ! so we calculate e and w accordingly to avoid singularities
    e = 1.0_rk
    tmpf = 0.5_rk*h2/gmsum
    w = acos(2.0_rk*tmpf/r - 1.0_rk)
    if(vdotr < 0.0_rk) w = TWOPI - w
    tmpf = (a + r)/(a*e)
    capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_rk))

  end if

  capm = e*sinh(capf) - capf
  omega = u - w
  if(omega < 0.0_rk) omega = omega + TWOPI
  omega = omega - TWOPI*int(omega/TWOPI)

end if

! PARABOLA: (NOTE - in this case we use "a" to mean pericentric distance)
if(ialpha == 0) then

  a = 0.5_rk*h2/gmsum
  e = 1.0_rk
  w = acos(2.0_rk*a/r - 1.0_rk)
  if(vdotr < 0.0_rk) w = TWOPI - w
  tmpf = tan(0.5_rk*w)
  capm = tmpf*(1.0_rk + tmpf*tmpf/3.0_rk)
  omega = u - w
  if(omega < 0.0_rk) omega = omega + TWOPI
  omega = omega - TWOPI*int(omega/TWOPI)

end if

return
end subroutine orbel_xv2el