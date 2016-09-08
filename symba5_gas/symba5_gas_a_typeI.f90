subroutine symba5_gas_a_typeI(time,dt,nbodm,mass,rgap,wgap,xh,yh,zh,vxh,vyh,vzh,&
           sig0,spower,zscale,rgi,rgf,tgdecay,ca,ce)
!-------------------------------------------------------------------------------------
!
! Type - I migration
!
! By: Martin Duncan & Hal Levison
! Date: ?
!-------------------------------------------------------------------------------------
use module_swift
use module_interfaces, except_this_one => symba5_gas_a_typeI
implicit none

! Inputs Only:
integer(ik) :: nbodm
real(rk) :: mass(nbodm),time,dt
real(rk) :: xh(nbodm),yh(nbodm),zh(nbodm)
real(rk) :: spower,sig0,zscale
real(rk) :: rgi,rgf,tgdecay,drsq,ca,ce
real(rk) :: rgap(nbodm),wgap(nbodm)

! Inputs/Outputs:
real(rk) :: vxh(nbodm),vyh(nbodm),vzh(nbodm)

! Internals
integer(ik) :: i,j,ialpha
real(rk) :: po2,znaught,den,gm,a,e,qdum,sig1
real(rk) :: hovera,ommin1,mdisk1,ta,te,erat
real(rk) :: vdotr,r2,faca,ti,mdisk,mrat,cp,sp,po,aphi,fgap
real(rk) :: atypeIx,atypeIy,atypeIz

!-----------------
! Executable code
!-----------------

sig1 = sig0*exp(-time/tgdecay) ! Exponential decay

do i = 2, nbodm

  po2 = xh(i)**2 + yh(i)**2

  if((po2 >= rgi**2) .or. (po2 <= rgf**2)) then

    gm = mass(1) + mass(i)
    call orbel_xv2aeq(xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i),gm,ialpha,a,e,qdum)

    mrat = mass(i)/mass(1)
    hovera = zscale*a**0.25_rk
    ommin1 = sqrt(a**3/gm)
    mdisk = PI*sig1*(a**(2.0_rk - spower))/mass(1)
    ta = (ommin1*hovera**2)/(ca*mrat*mdisk)
    te = (ommin1*hovera**4)/(ce*mrat*mdisk)

    ! Now multiply by the eccentricity correction from Pap. and Larwood
    if(rgap(i) <= 0.0_rk) then

      erat = e/hovera
      ta = ta*(1.0_rk + (erat/1.3_rk)**5)/(1.0_rk - (erat/1.1_rk)**4)
      te = te*(1.0_rk + 0.25_rk*erat**3)
      fgap = 1.0_rk

    else

      fgap = (a - rgap(i))/(abs(a - rgap(i)) + a*wgap(i))

    end if

    ! For *now* assume timescale for vertical damping ti equals te
    ! even though Pap. and Larwood don't discuss inclination dependent
    ! term in te. What's done here could generate anomolous
    ! inclination damping timescale if inc takes body above scale height
    ! and ecc remains less than hovera. Can we legitimately use parameter
    ! like sin(inc)/hovera to modify ti as we do for te with e/hovera??

    ti = te

    ! Now get the acceleration components as in Pap. and Larwood
    vdotr = xh(i)*vxh(i) + yh(i)*vyh(i) + zh(i)*vzh(i)
    r2 = xh(i)**2 + yh(i)**2 + zh(i)**2
    faca = 2.0_rk*vdotr/(r2*te)

    atypeIx = -vxh(i)*fgap/ta - faca*xh(i)
    atypeIy = -vyh(i)*fgap/ta - faca*yh(i)
    atypeIz = -vzh(i)*fgap/ta - faca*zh(i) - 2.0_rk*vzh(i)/ti

    ! Apply the Type-I acceleration
    vxh(i) = vxh(i) + atypeIx*dt
    vyh(i) = vyh(i) + atypeIy*dt
    vzh(i) = vzh(i) + atypeIz*dt

  end if

end do

return
end subroutine symba5_gas_a_typeI
