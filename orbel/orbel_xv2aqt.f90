subroutine orbel_xv2aqt(x, y, z, vx, vy, vz, gmsum, ialpha, a, q, capm, tperi)
!***********************************************************************
!	               ORBEL_XV2AQT.F
!***********************************************************************
!     PURPOSE:  Given the cartesian position and velocity of an orbit, 
!       compute the osculating orbital elements a,  e,  and q only.
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
!            q        ==> perihelion distance (real scalar); q = a(1 - e)
!          capm       ==> mean anomoly(real scalar)
!          tperi      ==> time to next or last perihelion,  which ever is less
!                         (real scalar)
!
!     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech."
!     REMARKS: Based on M. Duncan's orbel_xv2el.f
!      This routine is generally applied to study (hyperbolic) close
!       encounters of test particles with planets.
!     AUTHOR:  Hal Levison
!     DATE WRITTEN:  8/7/01
!     REMARKS: The tperi may not be correct for parabolic orbits.
!              I Think it is OK but beware!
!     REVISIONS:
use module_swift
implicit none

! Inputs Only:
real(rk) :: x, y, z, vx, vy, vz, gmsum

! Outputs
integer(ik) :: ialpha
real(rk) :: a, q, capm, tperi

! Internals:
real(rk) :: hx, hy, hz, h2, r, v2, energy, fac, vdotr, cape, e
real(rk) :: capf, tmpf, meanmo, face, w

!----
! Executable code

! Compute the angular momentum H, and thereby the inclination INC
hx = y*vz - z*vy
hy = z*vx - x*vz
hz = x*vy - y*vx
h2 = hx*hx + hy*hy + hz*hz

! Compute the radius R and velocity squared V2, and the dot
! product RDOTV, the energy per unit mass ENERGY.
r = sqrt(x*x + y*y + z*z)
v2 = vx*vx + vy*vy + vz*vz
energy = 0.5_rk*v2 - gmsum/r
vdotr = x*vx + y*vy + z*vz

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

  else

    e = 0.0_rk
    cape = 0.0_rk

  end if

  capm = cape - e*sin(cape)
  q = a*(1.0_rk - e)

end if

! HYPERBOLA
if(ialpha == 1) then

  a = 0.5_rk*gmsum/energy
  fac = h2/(gmsum*a)

  if(fac > TINY) then

    e = sqrt(1.0_rk + fac)

    ! have to insert minus sign in expression for q because this code
    ! takes a > 0,  even for a hyperbola
    q = -a*(1.0_rk - e)

    tmpf = (a + r)/(a*e)
    if(tmpf < 1.0_rk) tmpf = 1.0_rk

    capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_rk))
    if(vdotr < 0.0_rk) capf = - capf

  else

    ! we only get here if a hyperbola is essentially a parabola
    ! so we calculate e accordingly to avoid singularities
    e = 1.0_rk
    q = 0.5_rk*h2/gmsum

    tmpf = (a + r)/(a*e)
    capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_rk))

  end if

  capm = e*sinh(capf) - capf

end if

! PARABOLA: (NOTE - in this case "a",  which is formally infinite,
!         is arbitrarily set equal to the pericentric distance q).
if(ialpha == 0) then

  a = 0.5_rk*h2/gmsum
  e = 1.0_rk
  q = a
  w = acos(2.0_rk*a/r - 1.0_rk)
  if(vdotr < 0.0_rk) w = TWOPI - w
  tmpf = tan(0.5_rk*w)
  capm = tmpf*(1.0_rk + tmpf*tmpf/3.0_rk)

end if

meanmo = sqrt(gmsum/a**3)
if((capm < PI) .or. (ialpha >= 0)) then

  tperi = -1.0_rk*capm/meanmo

else

  tperi = -1.0_rk*(capm - TWOPI)/meanmo

end if

return
end subroutine orbel_xv2aqt