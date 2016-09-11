program swift_symba5
!**********************************************************************
!		      SWIFT_SYMBA5.F
!**********************************************************************
!
!                 To run, need 2 input files. The code prompts for
!                 the file names, but examples are :
!
!                   parameter file like       param.in
!		    planet file like          pl.in
!
!  NOTE:  No test particles in this code and the massive bodies
!         are dimensioned at NTPMAX
!  Created by Hal and Martin September 2006.
!
! Authors:  Hal Levison \& Martin Duncan
! Date:    11/21/96
! Last revision: 12/27/96
use swift
use util
use io
use discard
use symba5
implicit none

integer(ik) :: nbod,i1st,i,id,nbodm,nbodo
integer(ik) :: iflgchk,iub,iuj,iud,iue,ium

real(rk) :: mass(NTPMAX),j2rp2,j4rp4
real(rk) :: xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
real(rk) :: vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)

real(rk) :: xht(1),yht(1),zht(1)       ! Dummy for the io
real(rk) :: vxht(1),vyht(1),vzht(1)
integer(ik) :: ntp,istat(1),istat_r(1,NSTAT)

real(rk) :: t0,tstop,dt,dtout,dtdump
real(rk) :: t,tout,tdump,tfrac
real(rk) :: eoff = 0.0_rk
real(rk) :: rpl(NTPMAX),rhill(NTPMAX)

real(rk) :: rmin,rmax,rmaxu,qmin,mtiny
logical(lk) :: lclose, lgas
integer(ik) :: mergelst(NTPMAX,2),mergecnt
integer(ik) :: iecnt(NTPMAX), ibound(NTPMAX), isenc

character(len = 24) :: outfile, fopenstat
character(len = 24) :: inparfile = "param.in", inplfile = "pl.in"
character(len = 24) :: dparfile = "dump_param.dat", dplfile = "dump_pl.dat"

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
real(rk) :: ke,pot,energy,eltot(3)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----
! Executable code

ntp = 0

! print version number
call util_version

! Get data for the run and the test particles
call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,lgas,mtiny,outfile,fopenstat)

call io_init_pl_symba(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)

! Calculate the location of the last massive particle
call symba5_nbodm(nbod, mass, mtiny, nbodm)

! Initialize initial time and times for first output and first dump
t = t0
tout = t0 + dtout
tdump = t0 + dtdump

iub = 20
iuj = 30
iud = 40
iue = 60
ium = 21

! Do the initial io write
if(btest(iflgchk,0)) then ! bit 0 is set

  call io_write_frame(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
  call io_write_mass(t0,nbod,mass,outfile,ium,fopenstat)

end if

if(btest(iflgchk,1)) then ! bit 1 is set

  call io_write_frame_r(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat_r,outfile,iub,fopenstat)
  call io_write_mass_r(t0,nbod,mass,outfile,ium,fopenstat)

end if

if(btest(iflgchk,2)) then ! bit 2 is set

  eoff = 0.0_rk
  call anal_energy_discard5(1,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,ke,pot,energy,eltot)
  call anal_energy_write(t,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,iue,fopenstat,eoff)

else

  call anal_energy_discard5(-1,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,ke,pot,energy,eltot)

end if

! must initize discard io routine
if(btest(iflgchk,4)) then ! bit 4 is set

  call io_discard_mass(0,t,0,mass(1),rpl(1),xh(1),yh(1),zh(1),vxh(1),vyh(1),vzh(1),iud,-1,fopenstat)

end if

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, ibound)
do i = 1, nbod

  ibound(i) = 0

end do
!$OMP END PARALLEL DO

i1st = 0
!*************** The Big Loop **********************************
write(*,'(/,a,/)') ' ************** MAIN LOOP ******************'

do while((t <= tstop) .and. (nbod > 1))

  call symba5_step_pl(i1st,t,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,dt,lclose,rpl,isenc,&
       mergelst,mergecnt,iecnt,eoff,rhill,mtiny,ibound)

  t = t + dt

  if(btest(iflgchk,4)) then ! bit 4 is set

    nbodo = nbod
    call discard_massive5(t,dt,nbod,nbodm,mass,xh,yh,zh,vxh,vyh,vzh,rmin,rmax,rmaxu,qmin,lclose,rpl,rhill,isenc,mergelst,mergecnt,&
         iecnt,eoff,i1st,ibound)
    if(nbodo /= nbod) call symba5_nbodm(nbod,mass,mtiny,nbodm)

  end if


  ! if it is time, output orb. elements,
  if(t >= tout) then

    if(btest(iflgchk,0)) then ! bit 0 is set

      call io_write_frame(t,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
      call io_write_mass(t,nbod,mass,outfile,ium,fopenstat)

    end if

    if(btest(iflgchk,1)) then ! bit 1 is set
    
      call io_write_frame_r(t,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat_r,outfile,iub,fopenstat)
      call io_write_mass_r(t,nbod,mass,outfile,ium,fopenstat)
    
    end if

    tout = tout + dtout

  end if

  ! If it is time, do a dump
  if(t >= tdump) then

    tfrac = (t - t0)/(tstop - t0)
    write(*,'(a,1p1e12.5,a,0pf5.3,a,i6)') ' Time = ', t ,'; Fraction complete = ', tfrac , '; Number of bodies = ', nbod

    !call io_dump_pl_symba(dplfile,nbod,mass,xh,yh,zh,vxh,vyh,vzh,lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
    call io_dump_param(dparfile,t,tstop,dt,dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)

    tdump = tdump + dtdump

    if(btest(iflgchk,2)) then ! bit 2 is set

      call anal_energy_write(t,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,iue,fopenstat,eoff)

    end if

  end if

end do
!********** End Of The Big Loop **********************************

! Do a final dump for possible resumption later
t = tstop
call io_dump_pl_symba(dplfile,nbod,mass,xh,yh,zh,vxh,vyh,vzh,lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
call io_dump_param(dparfile,t,tstop,dt,dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)

call util_exit(SUCCESS)
end program swift_symba5
