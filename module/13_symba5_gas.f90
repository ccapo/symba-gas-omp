!!
!!
module symba5_gas
! Module for routines in symba5_gas directory
use swift
use orbel
use symba5
implicit none

contains
!!
!!
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
!!
  subroutine symba5_gas_a_typeI(time,dt,nbodm,mass,rgap,wgap,xh,yh,zh,vxh,vyh,vzh,&
             sig0,spower,zscale,rgi,rgf,tgdecay,ca,ce)
  !-------------------------------------------------------------------------------------
  !
  ! Type - I migration
  !
  ! By: Martin Duncan & Hal Levison
  ! Date: ?
  !-------------------------------------------------------------------------------------
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
!!
  subroutine symba5_gas_dragcoef(ma,kn,re,cdr)
  !**************************************************************
  !     dragcoef.f Compute the drag coefficient C_D as a function of
  !     Knudsen number and Mach number via the Reynolds number
  !     The routine returns M, K, R and C_D.
  !
  !  Version from Martin on 7/11/07
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
!!
  subroutine symba5_gas_drag_kick(time,nbod,nbodm,mass,rhill,xh,yh,zh,vxh,vyh,vzh,dt)
  !
  ! Drag Kick
  !
  implicit none

  ! Inputs Only:
  integer(ik) :: nbod,nbodm
  real(rk) :: mass(nbod),rhill(nbod),dt,time
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)

  ! Inputs and Outputs:
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)

  ! Internals
  logical(lk) :: lgas
  integer(ik), save :: i1st = 0
  integer(ik) :: i,j,ierr,iencounter(nbod),ngap,igap,ialpha
  real(rk) :: fwgap,dr2,gm,apl,epl,qpl,xhi,yhi,zhi,rhilli2
  real(rk), save :: deng0,gpower,deng0s,zscale,rgi,rgf,tgdecay,sig0,spower,ca,ce,rgap(NPLMAX),wgap(NPLMAX)

  ! rcomet and dragc are computed in gen_ray_size, and stored in the sizedist common block
  !real(rk), save :: rcomet, dcomet, dragc
  !save rcomet,dcomet,dragc

  ! rcomet and dragc are now variables, dcomet is no longer required
  namelist / gas_params / lgas, deng0s, gpower, zscale, rgi, rgf, tgdecay, ca, ce, ngap, igap, fwgap

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
!!
  subroutine symba5_gas_step_helio(i1st,nbod,nbodm,mass,j2rp2,j4rp4,lgas,xh,yh,zh,vxh,vyh,vzh,dt)
  !-------------------------------------------------------------------------
  !			SYMBA5_GAS_STEP_HELIO.F90
  !-------------------------------------------------------------------------
  ! This subroutine takes a step in helio coord.
  ! Does a KICK than a DRIFT than a KICK.
  ! ONLY DOES MASSIVE PARTICLES
  !
  !             Input:
  !                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
  !                                     (real scalars)
  !                 xh,yh,zh      ==>  initial position in helio coord
  !                                    (real arrays)
  !                 vxh,vyh,vzh   ==>  initial velocity in helio coord
  !                                    (real arrays)
  !                 dt            ==>  time step
  !             Output:
  !                 xh,yh,zh      ==>  final position in helio coord
  !                                       (real arrays)
  !                 vxh,vyh,vzh   ==>  final velocity in helio coord
  !                                       (real arrays)
  ! Remarks: Based on helio_step_pl.f but does not pass the intermediate
  !          positions and velocities back for the TP to use.
  ! Authors:  Hal Levison
  ! Date:    3/20/97
  ! Last revision: 12/13/00
  implicit none

  ! Inputs Only:
  logical(lk) :: lgas
  integer(ik) :: nbod,i1st,nbodm
  real(rk) :: mass(nbod),dt,j2rp2,j4rp4

  ! Inputs and Outputs:
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)

  ! Internals:
  integer(ik) :: i1stloc
  real(rk) :: dth,msys
  real(rk) :: axh(nbod),ayh(nbod),azh(nbod)
  real(rk), save :: vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX)
  real(rk) :: ptxb,ptyb,ptzb ! Not used here
  real(rk) :: ptxe,ptye,ptze

  !-----------------
  ! Executable code
  !-----------------

  dth = 0.5_rk*dt

  i1stloc = i1st
  if(lgas) i1st = 0 ! Forces converting heliocentric to barycentric velocities when gas drag is present

  if(i1st == 0) then

    ! Convert vel to bary to jacobi coords
    call coord_vh2b(nbod,mass,vxh,vyh,vzh,vxb,vyb,vzb,msys)
    i1st = 1 ! Turn this off
    !write(100,'(a,3(1x,1pe22.15))') 'vxh(1:3) = ', vxh(1:3)
    !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  end if

  ! Do the linear drift due to momentum of the Sun
  call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,xh,yh,zh,ptxb,ptyb,ptzb)
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  ! Get the accelerations in helio frame. if frist time step
  call symba5_helio_getacch(i1stloc,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,axh,ayh,azh)
  !if(i1stloc == 0) write(100,'(a,3(1x,1pe22.15))') 'axh(1:3) = ', axh(1:3)

  ! Apply a heliocentric kick for a half dt
  call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  ! Drift in helio coords for the full step
  call helio_drift(nbod,mass,xh,yh,zh,vxb,vyb,vzb,dt)
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  ! Get the accelerations in helio frame. if first time step
  call symba5_helio_getacch(0,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,axh,ayh,azh)
  !write(100,'(a,3(1x,1pe22.15))') 'axh(1:3) = ', axh(1:3)

  ! Apply a heliocentric kick for a half dt
  call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  ! Do the linear drift due to momentum of the Sun
  call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,xh,yh,zh,ptxe,ptye,ptze)
  !write(100,'(a,3(1x,1pe22.15))') 'vxb(1:3) = ', vxb(1:3)

  ! Convert back to helio velocities
  call coord_vb2h(nbod,mass,vxb,vyb,vzb,vxh,vyh,vzh)
  !write(100,'(a,3(1x,1pe22.15))') 'vxh(1:3) = ', vxh(1:3)
  !write(100,*) '--------'

  return
  end subroutine symba5_gas_step_helio
!!
  subroutine symba5_gas_step_pl(i1st,time,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,dt,lclose,lgas,rpl,isenc,&
             mergelst,mergecnt,iecnt,eoff,rhill,mtiny,ibound)
  !---------------------------------------------------------------------
  !			SYMBA5_GAS_STEP_PL.F90
  !---------------------------------------------------------------------
  !
  !             Input:
  !                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
  !                 time          ==>  current time (real scalar)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 nbodm  ==>  location of the last massie body
  !                                    (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
  !                                     (real scalars)
  !                 xh,yh,zh      ==>  initial position in helio coord
  !                                    (real arrays)
  !                 vxh,vyh,vzh   ==>  initial velocity in helio coord
  !                                    (real arrays)
  !                 dt            ==>  time step
  !                 lclose        ==> .true. --> check for close encounters
  !                                      (logical*2 scalar)
  !                 rpl           ==>  physical size of a planet.
  !                                    (real array)
  !                 eoff          ==>  Energy offset (real scalar)
  !                 rhill         ==>  size of planet's hills sphere
  !                                    (real array)
  !                 mtiny         ==>  Small mass  (real array)
  !             Output:
  !                 xh,yh,zh      ==>  final position in helio coord
  !                                       (real arrays)
  !                 vxh,vyh,vzh   ==>  final velocity in helio coord
  !                                       (real arrays)
  !                 rpl           ==>  Recalculated physical size of a planet.
  !                                    if merger happened (real array)
  !                 nbod          ==>  Recalculated number of massive bodies
  !                                    if merger happened (int scalar)
  !                 mass          ==>  Recalculated mass of bodies
  !                                    if merger happened (real array)
  !                 isenc         ==>  0 --> No encounter during last dt
  !                                    1 --> There was encounters
  !                                     (integer scalar)
  !                 mergelst      ==>  list of mergers (int array)
  !                 mergecnt      ==>  count of mergers (int array)
  !                 iecnt         ==>  Number of encounters (int*2 array)
  !                 eoff          ==>  Energy offset (real scalar)
  !                 rhill         ==>  size of planet's hills sphere
  !                                    (real array)
  !
  ! Remarks: Based on symba2_step_pl.f
  ! Authors:  Hal Levison
  ! Date:    11/27/97
  ! Last revision:
  !$ use omp_lib
  implicit none

  ! Inputs Only:
  logical(lk) :: lclose,lgas
  integer(ik) :: nbod,i1st,nbodm
  real(rk) :: mass(nbod),dt,time,j2rp2,j4rp4,mtiny

  ! Inputs and Outputs:
  integer(ik) :: ibound(nbod)
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: rpl(nbod),rhill(nbod),eoff

  ! Outputs only
  integer(ik) :: isenc
  integer(ik) :: iecnt(nbod),ielev(nbod)
  integer(ik) :: mergelst(NENMAX,2),mergecnt

  ! Internals
  logical(lk) :: lvdotr ! Not used in the routine
  integer(ik) :: i, j, k, ieflg
  integer(ik) :: ielst(NENMAX,2), ielc
  integer(ik) :: id, ielc_shared(nthreads), istart(nthreads), ielst_shared(nbod)
  real(rk) :: dth, rhillij(2), xhij(2), yhij(2), zhij(2), vxhij(2), vyhij(2), vzhij(2)

  !-----------------
  ! Executable code
  !-----------------

  ! Initialize relevant variables
  isenc = 0
  ielc = 0

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, iecnt, ielev)
  do i = 1, nbod

    iecnt(i) = 0
    ielev(i) = -1

  end do
  !$OMP END PARALLEL DO

  ! Check for encounters
  if(lclose) then

    do i = 2, nbodm

      ! Extract the information for body i, and make a copy available for each thread
      rhillij(1) = rhill(i)
      xhij(1) = xh(i); yhij(1) = yh(i); zhij(1) = zh(i)
      vxhij(1) = vxh(i); vyhij(1) = vyh(i); vzhij(1) = vzh(i)

      !$OMP PARALLEL DEFAULT(NONE) PRIVATE(j, id, ieflg) &
      !$OMP FIRSTPRIVATE(i, xhij, yhij, zhij, vxhij, vyhij, vzhij, rhillij, dt) &
      !$OMP SHARED(nbod, nthreads, rhill, xh, yh, zh, vxh, vyh, vzh, istart, ielc_shared, ielst_shared, iecnt, ielev)

      id = 1                              ! Set thread identifier for *serial* case
      !$ id = omp_get_thread_num() + 1    ! Set thread identifier for *parallel* case
      ielc_shared(id) = 0                 ! Initialize encounter counter for each thread
      istart(id) = (id - 1)*nbod/nthreads ! Set the starting position in ielst_shared for each thread

      !$OMP DO SCHEDULE(STATIC)
      do j = i + 1, nbod

        rhillij(2) = rhill(j)
        xhij(2) = xh(j); yhij(2) = yh(j); zhij(2) = zh(j)
        vxhij(2) = vxh(j); vyhij(2) = vyh(j); vzhij(2) = vzh(j)
        call symba5_chk2(rhillij, xhij, yhij, zhij, vxhij, vyhij, vzhij, dt, 0, ieflg)

        if(ieflg /= 0) then

          iecnt(j) = iecnt(j) + 1
          ielev(j) = 0
          ielc_shared(id) = ielc_shared(id) + 1 ! Update counter *before* storing index j since istart starts at 0
          ielst_shared(istart(id) + ielc_shared(id)) = j

        end if

      end do
      !$OMP END DO NOWAIT

      !$OMP END PARALLEL

      ! If there any encounters, append to the global encounter list for body i
      if(sum(ielc_shared) > 0) then

        ! Set global encounter flag
        isenc = 1

        ! Set the recursion level to 0 for body i
        ielev(i) = 0

        ! Search through each thread's local encounter list
        do k = 1, nthreads

          ! If an encounter was found for thread k, append ielc_shared(k) entries to ielst
          if(ielc_shared(k) > 0) then

            do j = 1, ielc_shared(k)

              iecnt(i) = iecnt(i) + 1
              ielst(ielc + j, 1) = i
              ielst(ielc + j, 2) = ielst_shared(istart(k) + j)

            end do

            ! Update the global encounter counter
            ielc = ielc + ielc_shared(k)

            if(ielc > NENMAX) then

              write(*,'(a)') 'ERROR: Encounter matrix is filled.'
              write(*,'(a)') 'STOPPING'
              call util_exit(FAILURE)

            end if

          end if

        end do

      end if

    end do

  end if

  !isenc = 0

  dth = 0.5_rk*dt
  if(lgas) call symba5_gas_drag_kick(time,nbod,nbodm,mass,rhill,xh,yh,zh,vxh,vyh,vzh,dth)

  ! Do a step
  if(isenc == 0) then

    !write(100,*) '--------'
    !write(100,*) time, ielc
    !write(100,*) 'symba5_gas_step_helio'
    call symba5_gas_step_helio(i1st,nbod,nbodm,mass,j2rp2,j4rp4,lgas,xh,yh,zh,vxh,vyh,vzh,dt)
    mergecnt = 0

  else

    !write(100,*) '--------'
    !write(100,*) time, ielc
    !write(100,*) 'symba5_step_interp'
    call symba5_step_interp(time,ielev,nbod,nbodm,mass,rhill,j2rp2,j4rp4,rpl,xh,yh,zh,vxh,vyh,vzh,dt,&
         mergelst,mergecnt,eoff,ielc,ielst,mtiny,ibound)
    i1st = 0

  end if

  if(lgas) call symba5_gas_drag_kick(time,nbod,nbodm,mass,rhill,xh,yh,zh,vxh,vyh,vzh,dth)

  ! Print number of encounters and mergers found in this time step
  !write(100,*) time, ielc, mergecnt
  !do i = 1, ielc
  !  write(100,'(2i9)') ielst(i,:)
  !end do

  return
  end subroutine symba5_gas_step_pl
!!
!!
end module symba5_gas