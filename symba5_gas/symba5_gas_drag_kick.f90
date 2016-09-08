subroutine symba5_gas_drag_kick(time,nbod,nbodm,mass,rhill,xh,yh,zh,vxh,vyh,vzh,dt)
!
! Drag Kick
!
use module_swift
use module_interfaces, except_this_one => symba5_gas_drag_kick
implicit none

! Inputs Only:
integer(ik) :: nbod,nbodm
real(rk) :: mass(nbod),rhill(nbod),dt,time
real(rk) :: xh(nbod),yh(nbod),zh(nbod)

! Inputs and Outputs:
real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)

! Internals
integer(ik), save :: i1st = 0
integer(ik) :: i,j,ierr,iencounter(nbod),ngap,igap,ialpha
real(rk) :: fwgap,dr2,gm,apl,epl,qpl,xhi,yhi,zhi,rhilli2
real(rk), save :: deng0,gpower,deng0s,zscale,rgi,rgf,tgdecay,sig0,spower,ca,ce,rgap(NPLMAX),wgap(NPLMAX)

! rcomet and dragc are computed in gen_ray_size, and stored in the sizedist common block
!real(rk), save :: rcomet, dcomet, dragc
!save rcomet,dcomet,dragc

! rcomet and dragc are now variables, dcomet is no longer required
namelist / gas_params / deng0s, gpower, zscale, rgi, rgf, tgdecay, ca, ce, ngap, igap, fwgap

!-----------------
! Executable code
!-----------------

if(i1st == 0) then

  i1st = 1

  open(unit = 7, file = "param.in", status = 'old')
  read(unit = 7, nml = gas_params)
  close(unit = 7)

  deng0s = deng0s*1.68314195e6_rk ! convert to Solar mass/AU^3
  deng0 = deng0s*mass(1)          ! correct for units

  ! dragc is now computed in gen_ray_size *** 04/24/08 -- CCC ***
  !
  ! Drag coefficient 3*C_D/(8*rho_comet*r_comet) and proper unit conversion.
  ! dragc = 0.84423_rk/(rcomet*dcomet) !(0.375/(1.683d6*mstar/1.49587d8)

  sig0 = deng0*zscale*sqrt(PI)
  spower = gpower - 1.25_rk

  ! Now get the gap information
  !OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbodm, rgap, wgap)
  do i = 1, nbodm

    rgap(i) = -1.0_rk
    wgap(i) = -1.0_rk

  end do
  !OMP END PARALLEL DO

  if(ngap >= 1) then

    write(*,'(a)') ' Gap information:'
    write(*,'(a)') 'id, gap location (AU), fractional gap width'
    ialpha = -1

    do i = 1, ngap

      j = igap
      gm = mass(1) + mass(j)
      call orbel_xv2aeq(xh(j),yh(j),zh(j),vxh(j),vyh(j),vzh(j),gm,ialpha,apl,epl,qpl)
      write(*,*) j, apl, fwgap
      rgap(i) = apl
      wgap(i) = fwgap

    end do

  end if

end if

! Initialize the ecounter flag
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(j) SHARED(nbod, nbodm, iencounter)
do j = nbodm + 1, nbod

  iencounter(j) = 0

end do
!$OMP END PARALLEL DO

do i = 2, nbodm

  ! Extract a copy of the relevant information for body i
  xhi = xh(i); yhi = yh(i); zhi = zh(i)
  rhilli2 = rhill(i)**2

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(j, dr2) &
  !$OMP FIRSTPRIVATE(xhi, yhi, zhi, rhilli2) SHARED(nbod, nbodm, xh, yh, zh, iencounter)
  do j = nbodm + 1, nbod

    dr2 = (xh(j) - xhi)**2 + (yh(j) - yhi)**2 + (zh(j) - zhi)**2
    if(dr2 < rhilli2) iencounter(j) = 1

  end do
  !$OMP END PARALLEL DO

end do

! Apply aerodynamic gas drag on planetesimals
call symba5_gas_a_drag(time,dt,nbod,nbodm,iencounter,mass(1),xh,yh,zh,vxh,vyh,vzh,&
     deng0,gpower,deng0s,zscale,rgi,rgf,tgdecay)

! Apply Type - I drag on embryos
call symba5_gas_a_typeI(time,dt,nbodm,mass,rgap,wgap,xh,yh,zh,vxh,vyh,vzh,&
     sig0,spower,zscale,rgi,rgf,tgdecay,ca,ce)

return
end subroutine symba5_gas_drag_kick
