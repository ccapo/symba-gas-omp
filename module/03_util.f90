module util
! Module for routines contained in the util directory
use swift
use orbel

contains
!!
!!
  function util_disk_mass(sigma_0, a_in, a_out, p, a_snow, f_snow) result(m_disk)
  !--------------------------------------------------------------------------
  !			UTIL_DISK_MASS.F90
  !--------------------------------------------------------------------------
  ! Computes the disk mass in solids
  !
  ! Input:  sigma_0 ==> Surface density @ 1.0 AU [g/cm^2]
  !         a_in    ==> Inner radius of the disk [AU]
  !         a_out   ==> Outer radius of the disk [AU]
  !         p       ==> Power-law exponent
  !         a_snow  ==> Location of the snow line [AU]
  !         f_snow  ==> Material enhancement factor
  !
  ! Output: m_disk ==> Mass of disk in solids [MEarth]
  !
  ! By: Chris Capobianco
  ! Date: 05/27/07
  !--------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: sigma_0, a_in, a_out, p, a_snow, f_snow

  ! Output variables
  real(rk) :: m_disk

  ! Internal variables
  real(rk) :: sigma, r_0, p2, a_is, a_os

  ! Stop and print warning if p >= 2.0 and a_in = 0.0
  if((p >= 2.0_rk) .and. (a_in == 0.0_rk)) then

    write(*,'(a)')       ' Disk mass does not converge!'
    write(*,'(a,g14.6)') ' Power-law = ', p
    write(*,'(a,g14.6)') ' a_in (AU) = ', a_in
    stop

  end if

  ! Convert sigma_0 from [g/cm^2] to [M_Earth/AU^2]
  sigma = 3.749e-2_rk*sigma_0
  r_0 = 1.0_rk
  p2 = 2.0_rk - p
  a_is = a_in/a_snow
  a_os = a_out/a_snow

  if(p /= 2.0_rk) then

    if((a_in < a_snow) .and. (a_out < a_snow)) then

      m_disk = TWOPI*sigma*r_0**p*(a_out**p2 - a_in**p2)/p2

    else if((a_in < a_snow) .and. (a_out > a_snow)) then

      m_disk = TWOPI*sigma*r_0**p*a_snow**p2*(1.0_rk + f_snow*(a_os**p2 - 1.0_rk) - a_is**p2)/p2

    else ! a_in > a_snow and a_out > a_snow

      m_disk = TWOPI*sigma*r_0**p*f_snow*(a_out**p2 - a_in**p2)/p2

    end if

  else

    if((a_in < a_snow) .and. (a_out < a_snow)) then

      m_disk = TWOPI*sigma*r_0**p*log(a_out/a_in)

    else if((a_in < a_snow) .and. (a_out > a_snow)) then

      m_disk = TWOPI*sigma*r_0**p*(f_snow*log(a_os) - log(a_is))

    else ! a_in > a_snow and a_out > a_snow

      m_disk = TWOPI*sigma*r_0**p*f_snow*log(a_out/a_in)

    end if

  end if

  return
  end function util_disk_mass
!!
  subroutine util_exit(iflag)
  !-------------------------------------------------------------------------
  !    UTIL_EXIT.F90
  !-------------------------------------------------------------------------
  ! Exits program
  !
  ! Input: iflag ==> Status of exit
  !                   = 0 if normal exit
  !                   = 1 if exit because error
  !
  ! Authors: Hal Levison 
  ! Date: 08/06/93
  !-------------------------------------------------------------------------
  implicit none

  ! Inputs: 
  integer(ik) :: iflag

  !-----------------!
  ! Executable Code !
  !-----------------!

  if(iflag == 0) then

    write(*,'(/,a,f3.1,a)') 'Normal termination of SWIFT (Version: ', VER_NUM ,')'
    write(*,'(a)')          '------------------------------------------'

  else

    write(*,'(/,a,f3.1,a)') 'Terminating SWIFT (Version: ', VER_NUM ,') due to ERROR!!!'
    write(*,'(a)')          '-----------------------------------------------'

  end if

  stop
  end subroutine util_exit
!!
  subroutine util_hills1(mstar, mpl, xh, yh, zh, vxh, vyh, vzh, rhill)
  !-------------------------------------------------------------------------
  !    UTIL_HILLS1.F90
  !-------------------------------------------------------------------------
  ! This subroutine calculates the hill's sphere for the planets
  !
  !             Input:
  !                 mstar          ==>  mass of sun (real scalar)
  !                 mpl           ==>  mass of sun (real scalar)
  !                 xh, yh, zh      ==>  position of pl in helio coord 
  !                                    (real scalars)
  !                 vxh, vyh, vzh   ==>  velocity of pl in helio coord 
  !                                    (real scalars)
  !             Output:
  !                  rhill        ==>  the radius of planet's hill's sphere 
  !                                    (real scalar)
  !
  !
  ! Remarks: Based on util_hill
  ! Authors:  Hal Levison 
  ! Date:    1/8/97
  ! Last revision:
  !-------------------------------------------------------------------------
  implicit none

  ! Inputs:
  real(rk) :: mstar, mpl, xh, yh, zh
  real(rk) :: vxh, vyh, vzh

  ! Outputs
  real(rk) :: rhill

  ! Internals
  real(rk) :: mu, energy, ap, r, v2

  !-----------------!
  ! Executable code !
  !-----------------!

  mu = mstar*mpl/(mstar + mpl)
  r = sqrt(xh*xh + yh*yh + zh*zh)
  v2 = vxh*vxh + vyh*vyh + vzh*vzh
  energy = 0.5_rk*mu*v2 - mstar*mpl/r
  ap = -mstar*mpl/(2.0_rk*energy)
  rhill = ap*(mu/(3.0_rk*mstar))**(1.0_rk/3.0_rk)

  return
  end subroutine util_hills1
!!
  subroutine util_hills(nbod, mass, xh, yh, zh, vxh, vyh, vzh, r2hill)
  !-------------------------------------------------------------------------
  !    UTIL_HILLS.F90
  !-------------------------------------------------------------------------
  ! This subroutine calculates the hill's sphere for the planets
  !
  !             Input:
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 xh, yh, zh      ==>  initial position in helio coord 
  !                                    (real arrays)
  !                 vxh, vyh, vzh   ==>  initial velocity in helio coord 
  !                                    (real arrays)
  !             Output:
  !                  r2hill       ==>  the SQUARE of the planet's hill's sphere 
  !                                    (real array)
  !
  !
  ! Remarks: 
  ! Authors:  Hal Levison 
  ! Date:    2/19/93
  ! Last revision: 1/6/97
  !-------------------------------------------------------------------------
  implicit none

  ! Inputs: 
  integer(ik) :: nbod
  real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)

  ! Outputs
  real(rk) :: r2hill(nbod)

  ! Internals
  integer(ik) :: i
  real(rk) :: mass0, mu, energy, ap, rhill, r, v2

  !-----------------!
  ! Executable Code !
  !-----------------!

  r2hill(1) = 0.0_rk
  mass0 = mass(1)

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, mu, r, v2, energy, ap, rhill) &
  !$OMP FIRSTPRIVATE(mass0) SHARED(nbod, mass, xh, yh, zh, vxh, vyh, vzh, r2hill)
  do i = 2, nbod

    if(mass(i) /= 0.0_rk) then

      mu = mass0*mass(i)/(mass0 + mass(i))
      r = sqrt(xh(i)*xh(i) + yh(i)*yh(i) + zh(i)*zh(i))
      v2 = vxh(i)*vxh(i) + vyh(i)*vyh(i) + vzh(i)*vzh(i)
      energy = 0.5_rk*mu*v2 - mass0*mass(i)/r
      ap = -mass0*mass(i)/(2.0_rk*energy)
      rhill = ap*(mu/(3.0_rk*mass0))**(1.0_rk/3.0_rk)
      r2hill(i) = rhill*rhill

    else

      r2hill(i) = 0.0_rk

    end if

  end do
  !$OMP END PARALLEL DO

  return
  end subroutine util_hills
!!
  function util_kahan_sum(xsum_current, xi, xerror) result(xsum_new)
  !-------------------------------------------------------------------------------------------
  !         UTIL_KAHAN_SUM.F90
  !-------------------------------------------------------------------------------------------
  ! Sums two floating point scalars more accurately utilitizing the Kahan summation formula
  ! This function is designed to be used inside a do or while loop, where the initial value of
  ! of xsum_current is initialized appropriately and the initial value of xerror is 0.0_rk
  !
  ! N.B. Use this function if the summation is being performed for more than *three* terms
  !
  ! Input:  xsum_current - Current value of the sum
  !         xi           - i-th term to be added to the sum
  !         xerror       - Error term from the previous term of the sum
  !
  ! Output: xsum_new     - The updated value of the sum
  !         xerror       - The error term for this term of the sum
  !
  ! By: Chris Capobianco
  ! Date: 05/04/09
  !-------------------------------------------------------------------------------------------
  implicit none

  ! Input/Output variables
  real(rk), intent(in) :: xsum_current, xi
  real(rk), intent(inout) :: xerror
  real(rk) :: xsum_new

  ! Internal variables
  real(rk) :: low_bits

  low_bits = xi - xerror
  xsum_new = xsum_current + low_bits
  xerror = (xsum_new - xsum_current) - low_bits

  return
  end function util_kahan_sum
!!
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
!!
  function util_power_law(p, xmin, xmax) result(x)
  !------------------------------------------------------------------
  !    UTIL_POWER_LAW.F90
  !------------------------------------------------------------------
  ! Returns a random number drawn from a power-law distribution
  ! with exponent p, where xmin and xmax are the minimum and
  ! maximum range for the random number.
  !
  ! Input:  p    ==> Exponent of power-law distribution
  !         xmin ==> Minimum value for random number
  !         xmax ==> Maximum value for random number
  !
  ! Output: x    ==> Random number drawn from power-law distrubution
  !
  ! By: Chris Capobianco
  ! Date: 02/04/09
  !------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: p, xmin, xmax

  ! Output variable
  real(rk) :: x

  ! Internal variables
  real(rk) :: p1, ip1, xi

  p1 = 1.0_rk - p
  ip1 = 1.0_rk/p1
  xi = util_randomu()

  if(p /= 1.0_rk) then

    x = (xmin**p1 + xi*(xmax**p1 - xmin**p1))**ip1

  else

    x = xmin*(xmax/xmin)**xi

  end if

  return
  end function util_power_law
!!
  function util_randomu() result(ran)
  !--------------------------------------------------------------
  !    UTIL_RANDOMU.F90
  !--------------------------------------------------------------
  ! Generates a uniform random number between (0.0,1.0) using
  ! Intel Fortran 90 random number generator
  !--------------------------------------------------------------
  implicit none

  ! Output variable
  real(rk) :: ran

  call random_number(ran)

  return
  end function util_randomu
!!
  function util_rayleigh(rms, xmin, xmax) result(x)
  !--------------------------------------------------------------
  !    UTIL_RAYLEIGH.F90
  !--------------------------------------------------------------
  ! Returns a random number drawn from a rayleigh distribution
  ! with dispersion rms, where xmin and xmax are the minimum and
  ! maximum range for the random number.
  !--------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: rms, xmin, xmax

  ! Output variable
  real(rk) :: x

  ! Internal variables
  real(rk) :: fmin, fmax, xi

  fmin = -exp(-0.5_rk*(xmin/rms)**2)
  fmax = -exp(-0.5_rk*(xmax/rms)**2)
  xi = util_randomu()

  x = rms*sqrt(-2.0_rk*log((fmin - fmax)*xi - fmin))

  return
  end function util_rayleigh
!!
  subroutine util_version
  !-------------------------------------------------------------------------
  !    UTIL_VERSION.F90
  !-------------------------------------------------------------------------
  ! Prints version of SWIFT and contact information
  !
  ! Authors: Hal Levison
  ! Date: 02/21/94
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-------------------------------------------------------------------------
  implicit none

  !-----------------!
  ! Executable code !
  !-----------------!

  write(*,'(a)')          "!---------------------------------------------------------!"
  write(*,'(a)')          "!                                                         !"
  write(*,'(a,f3.1,a)')   "! SWIFT (Version: ", VER_NUM, ")                                    !"
  write(*,'(a)')          "!                                                         !"
  write(*,'(a)')          "!---------------------------------------------------------!"
  write(*,'(a)')          "!                                                         !"
  write(*,'(a)')          "! Authors:                                                !"
  write(*,'(a)')          "!  Martin Duncan: Queen's University                      !"
  write(*,'(a)')          "!  Hal Levison: Southwest Research Institute              !"
  write(*,'(a)')          "!                                                         !"
  write(*,'(a)')          "! Please address any comments or questions to:            !"
  write(*,'(a)')          "!  Hal Levison                                            !"
  write(*,'(a)')          "!  Geophysical, Astrophysical & Planetary Sciences        !"
  write(*,'(a)')          "!  Southwest Research Institute                           !"
  write(*,'(a)')          "!  1050 Walnut St.                                        !"
  write(*,'(a)')          "!  Suite 429                                              !"
  write(*,'(a)')          "!  Boulder, Co 80302                                      !"
  write(*,'(a)')          "!  (303) 546-0290                                         !"
  write(*,'(a)')          "!  Fax: (303) 546-9687                                    !"
  write(*,'(a)')          "!  (D)  swri::levison                                     !"
  write(*,'(a)')          "!  (I)  hal@gort.space.swri.edu                           !"
  write(*,'(a)')          "!                                                         !"
  write(*,'(a,/)')        "!---------------------------------------------------------!"

  return
  end subroutine util_version
!!
!!
end module util