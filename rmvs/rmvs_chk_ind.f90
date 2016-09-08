subroutine rmvs_chk_ind(xr, yr, zr, vxr, vyr, vzr, dt, r2crit, r2critp, iflag)
!*************************************************************************
!                            RMVS_CHK_IND.F
!*************************************************************************
!  Subroutine to check if a test particle and planet
!  are having or **will** have an encounter
!  in the next timestep.
!
!             Input:
!                 xr, yr, zr     ==>  relative position of tp wrt planet
!                                   (real scalar)
!                 vxr, vyr, vzr  ==>  relative velocity of tp wrt planet
!                                   (real scalor)
!                 dt           ==>  time step (real scalor)
!                 r2crit       ==> boundary of outer enc region
!                                   (real scalor)
!                 r2critp      ==> boundary of inner (planocentric) enc region
!                                   (real scalor)
!             Output:
!                 iflag        ==> encounter?  =  0 no
!                                              =  1 yes,  in outer region
!                                              = -1 yes,  in inner region
!
!
! Remarks: Based on Hal's wiscl_fk.f' but originally written by Martin Duncan
! Authors:  Hal Levison
! Date:    2/19/93
! Last revision:
use module_swift
use module_rmvs
implicit none

! Inputs:
real(rk) :: xr, yr, zr, vxr, vyr, vzr, dt, r2crit, r2critp

! Outputs
integer(ik) :: iflag

! Internals
real(rk) :: r2, v2, vdotr, tmin, r2min

!-----
! Executable code

! First check if we're already in the encounter region. If so return with flag set to one.
r2 = xr**2 + yr**2 + zr**2
if(r2 <= r2critp) then

  iflag = -1
  return

end if

! If we're heading outward, use r2 to calc iflag
vdotr = xr*vxr + yr*vyr + zr*vzr
if(vdotr > 0.0_rk) then

  if(r2 >= r2crit) then

    iflag = 0

  else

    iflag = 1

  end if

  return

end if

! We are converging, so we need to calc. the minimum separation attained in time dt.
v2 = vxr**2 + vyr**2 + vzr**2
tmin = -vdotr/v2

if(tmin < dt) then

  r2min = r2 - (vdotr**2)/v2

else

  r2min = r2 + 2.0_rk*vdotr*dt + v2*dt**2

end if

r2min = min(r2min, r2) ! really make sure

if(r2min <= r2critp) then

  iflag = -1

else if(r2min <= r2crit) then

  iflag = 1

else

  iflag = 0

end if

return
end subroutine rmvs_chk_ind
