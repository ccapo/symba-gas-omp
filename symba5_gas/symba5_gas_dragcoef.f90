subroutine symba5_gas_dragcoef(ma,kn,re,cdr)
!**************************************************************
!     dragcoef.f Compute the drag coefficient C_D as a function of
!     Knudsen number and Mach number via the Reynolds number
!     The routine returns M, K, R and C_D.
!
!  Version from Martin on 7/11/07
use module_swift
implicit none

real(rk) :: ma, kn, re, cdr

!-----
! Executable code

re = 4.448823857787100274_rk*ma/kn    ! Reynolds number

if(kn < 1.0_rk) then

  if(ma >= 1.0_rk) then

    cdr = 2.0_rk

  else

    if(re > 1.0e3_rk) then

      cdr = 0.44_rk + 1.56_rk*ma**2

    else

      cdr = 2.0_rk*ma**2 + 24.0_rk*(1.0_rk - ma**2)*(1.0_rk + 0.15_rk*re**0.687_rk)/re

    end if

  end if

else

  if(ma < 1.8_rk) then

    cdr = 3.6_rk/ma

  else

    cdr = 2.0_rk

  end if

end if

return
end subroutine symba5_gas_dragcoef
