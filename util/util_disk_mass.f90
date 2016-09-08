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
use module_swift
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
