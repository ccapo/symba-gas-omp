subroutine util_mass_peri(iflg, nbod, x, y, z, vx, vy, vz, mass, isperi, peri, lperi)
!-------------------------------------------------------------------------
!                            UTIL_MASS_PERI.F
!-------------------------------------------------------------------------
! This subroutine determines whether peri of a planet has taken place
!
!             Input:
!                 iflg           ==>  = 0 if first step; = 1 not (int scalar)
!                 nbod           ==>  number of bodies (int scalar)
!                 x, y, z          ==>  heliocentric position of planets
!                                       (real arrays)
!                 vx, vy, vz       ==>  heliocentric velcocities of planets
!                                       (real arrays)
!                 mass           ==>  mass of the bodies (real array)
!
!             Output:
!                 isperi         ==> = 0 if tp went through peri
!                                    =-1 if tp pre peri
!                                    = 1 if tp post peri
!                                         (integer array)
!                 peri           ==> set to pericenter dist. if isperi=0
!                                         (real array)
!                lperi           ==> set to .true. if isperi=0
!                                         (logical*2 array)
!
!
! Remarks: Based on util_peri.f
! Authors:  Hal Levison 
! Date:    12/30/96
! Last revision:
!-------------------------------------------------------------------------
use module_swift
use module_interfaces, except_this_one => util_mass_peri
implicit none

! Inputs Only: 
integer(ik) :: nbod, iflg
real(rk) :: x(nbod), y(nbod), z(nbod), mass(nbod)
real(rk) :: vx(nbod), vy(nbod), vz(nbod)

! Outputs:
real(rk) :: peri(nbod)
integer(ik) :: isperi(nbod)
logical(lk) :: lperi(nbod)

! Internals
integer(ik) :: i, ialpha
real(rk) :: mu, gm, a, e, vdotr

!-----------------!
! Executable code !
!-----------------!

if(iflg == 0) then    ! are we just setting thing up?

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, vdotr) SHARED(nbod, x, y, z, vx, vy, vz, isperi)
  do i = 2, nbod

    vdotr = x(i)*vx(i) + y(i)*vy(i) + z(i)*vz(i)

    if(vdotr > 0.0_rk) then

      isperi(i) = 1

    else

      isperi(i) = -1

    end if

  end do
  !$OMP END PARALLEL DO

else

  mu = mass(1)

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, vdotr, gm, ialpha, a, e) FIRSTPRIVATE(mu) &
  !$OMP SHARED(nbod, mass, x, y, z, vx, vy, vz, isperi, lperi, peri)
  do i = 2, nbod

    vdotr = x(i)*vx(i) + y(i)*vy(i) + z(i)*vz(i)

    if(isperi(i) == -1) then         ! was coming in

      if(vdotr < 0.0_rk) then        ! still coming in

        isperi(i) = -1

      else                           ! turned around

        isperi(i) = 0
        lperi(i) = .true.
        gm = mu + mass(i)
        call orbel_xv2aeq(x(i), y(i), z(i), vx(i), vy(i), vz(i), gm, ialpha, a, e, peri(i))

      end if

    else

      if(vdotr < 0.0_rk) then     ! coming in

        isperi(i) = -1

      else

        isperi(i) = 1             ! going out

      end if

    end if

  end do
  !$OMP END PARALLEL DO

end if

return
end subroutine util_mass_peri
