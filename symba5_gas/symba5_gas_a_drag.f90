!------------------------------------------------
!subroutine a_drag(time,nbodm,nbod,iencounter,mstar,x,y,z,vx
!     &     ,vy,vz,agasx,agasy,agasz,deng0,gpower,dragc,deng0s,rcomet
!     &     ,zscale,dcomet,rgi,rgf,tgdecay)
subroutine symba5_gas_a_drag(time,dt,nbod,nbodm,iencounter,mstar,x,y,z,vx,vy,vz,deng0,gpower,deng0s,zscale,rgi,rgf,tgdecay)
!-------------------------------------------------------------------------------------
!
! aero gas drag
!
!-------------------------------------------------------------------------------------
use module_swift
use module_interfaces, except_this_one => symba5_gas_a_drag
implicit none

! Inputs Only:
integer(ik) :: nbod,nbodm,iencounter(nbod)
real(rk) :: mstar,time,dt
real(rk) :: x(nbod),y(nbod),z(nbod)
real(rk) :: deng0,gpower,deng0s,zscale
real(rk) :: rgi,rgf,tgdecay
!real(rk) :: rcomet,dcomet,dragc
real(rk) :: rcomet(NTPMAX),dragc(NTPMAX)

! Inputs/Outputs:
real(rk) :: vx(nbod),vy(nbod),vz(nbod)

! Internals
logical(lk), save :: lfirst = .true.
integer(ik) :: i,j
real(rk) :: fac,po,znaught,den,eta,vkep,vfac,vgasx,vgasy,vgasz
real(rk) :: vrel,mach,knudsen,reynolds,c_d,gdrag
real(rk) :: agasx,agasy,agasz

! Size distribution common block
common / sizedist / rcomet, dragc

!-----------------
! Executable code
!-----------------

fac = exp(-time/tgdecay) ! Exponential decay

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(fac,mstar,rgi,rgf,zscale,deng0,gpower,dt) &
!$OMP PRIVATE(i,po,znaught,den,eta,vkep,vfac,vgasx,vgasy,vgasz,vrel,mach,knudsen,reynolds,c_d,gdrag,agasx,agasy,agasz) &
!$OMP SHARED(lfirst,nbod,nbodm,rcomet,dragc,x,y,z,vx,vy,vz,iencounter)
do i = nbodm + 1, nbod

  po = sqrt(x(i)**2 + y(i)**2)

  if((iencounter(i) == 0) .and. (po >= rgi) .and. (po <= rgf)) then

    ! Scale height of the gas (Hayashi and co.)
    znaught = zscale*po**1.25_rk

    ! Density follows power law model: po^(-gpower)*exp(-(z/z0)^2)
    den = deng0*(po**(-gpower))*exp(-(z(i)/znaught)**2)

    ! Pressure of gas parameter eta (e.g. Kokubo and Ida)
    eta = 6.0e-4_rk*(gpower + 0.5_rk)*sqrt(po)

    vkep = sqrt(mstar/sqrt(x(i)**2 + y(i)**2 + z(i)**2))

    ! Gas moves slower than Kepler velocity
    vfac = vkep*sqrt(1.0_rk - 2.0_rk*eta) ! Adachi et al., 1976
    vgasx = -vfac*y(i)/po                 ! brakes vx
    vgasy = vfac*x(i)/po                  ! brakes vy
    vgasz = 0.0_rk                        ! no braking in vz
    vrel = sqrt((vx(i) - vgasx)**2 + (vy(i) - vgasy)**2 + (vz(i) - vgasz)**2)

    mach = 3.32126045_rk*vrel*po**0.25_rk     ! Mach number
    knudsen = 1.1097308e-6_rk/(den*rcomet(i)) ! Knudsen number

    call symba5_gas_dragcoef(mach,knudsen,reynolds,c_d)

    ! Fix C_D = {0.5, 1.0, 2.0}
    !c_d = 0.5_rk
    !OMP MASTER
    !if(lfirst) then
    !
    !  write(*,'(a,f3.1)') "The fixed value of the drag coefficient is C_D = ", c_d
    !  lfirst = .false.
    !
    !end if
    !OMP END MASTER

    gdrag = dragc(i)*c_d*vrel*den*fac

    agasx = -gdrag*(vx(i) - vgasx)
    agasy = -gdrag*(vy(i) - vgasy)
    agasz = -gdrag*(vz(i) - vgasz)

    ! Apply the aerodynamic drag acceleration
    vx(i) = vx(i) + agasx*dt
    vy(i) = vy(i) + agasy*dt
    vz(i) = vz(i) + agasz*dt

  end if

end do
!$OMP END PARALLEL DO

return
end subroutine symba5_gas_a_drag
