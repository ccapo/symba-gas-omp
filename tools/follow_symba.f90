program follow_symba
!--------------------------------------------------------------------------
!			FOLLOW_SYMBA.F90
!--------------------------------------------------------------------------
!
! Extracts information on a particle in a binary file, and writes it to
! an ascii file
!
!--------------------------------------------------------------------------
use module_swift
use module_interfaces
implicit none

real(rk) :: mass(NTPMAX),dr
real(rk) :: xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
real(rk) :: vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)

integer(ik) :: nbodm,nbod,nbod0,ierr,ifol,istep,arg_num,iargc
integer(ik) :: iflgchk,iu,i,id,nleft,ium

real(rk) :: t0,tstop,dt,dtout,dtdump
real(rk) :: t,tmax,tm,mtiny

real(rk) :: rmin,rmax,rmaxu,qmin,rpl(NTPMAX),rhill(NTPMAX)
logical(lk) :: lclose,lgas
real(rk) :: a,e,inc,capom,omega,capm,j2rp2,j4rp4
real(rk) :: peri,apo,tg,obar

integer(ik) :: plist(NTPMAX),ifoln

character(len = 24) :: outfile, inparfile, inplfile, fopenstat, buf, ffile
character(len = 24) :: form = "UNFORMATTED", state = "OLD"

! Get data for the run and the test particles
inparfile = 'param.in'
call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,lgas,mtiny,outfile,fopenstat)

! Get planet data
inplfile = 'pl.in'
call io_init_pl_symba(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)

iu = 20
ium = 30

dr = 180.0_rk/PI

if(btest(iflgchk,0)) then

  write(*,'(a,/)') ' Reading an FXDR binary file'

else if(btest(iflgchk,1)) then

  write(*,'(a,/)') ' Reading an real*4 binary file'

else

  write(*,'(a,/)') ' ERROR: no binary file format specified in param file'
  stop

end if

! Get the number of command line arguments
arg_num = iargc()

! Read in the command line argument
if(arg_num == 1) then

  call getarg(1,buf)
  read(buf,*) ifol
  ifol = abs(ifol)
  write(*,'(a,i6)') ' Following particle number = ', ifol

else

  write(*,'(a)', advance = 'no') ' Input the particle number to follow: '
  read(*,*) ifol
  ifol = abs(ifol)

end if

if(btest(iflgchk,0))  then ! bit 0 is set

  call io_open_fxdr(outfile, 'r', .true., iu, ierr)
  call io_open_fxdr('mass.' // outfile, 'r', .true., ium, ierr)

else

  call io_open(iu, outfile, state, form, ierr)
  call io_open(ium, 'mass.' // outfile, state, form, ierr)

end if

if((0 <= ifol) .and. (ifol <= 9))         write(ffile, '(a11,i1,a4)') "follow_0000", ifol, ".dat"
if((10 <= ifol) .and. (ifol <= 99))       write(ffile, '(a10,i2,a4)') "follow_000", ifol, ".dat"
if((100 <= ifol) .and. (ifol <= 999))     write(ffile, '(a9, i3,a4)') "follow_00", ifol, ".dat"
if((1000 <= ifol) .and. (ifol <= 9999))   write(ffile, '(a8, i4,a4)') "follow_0", ifol, ".dat"
if((10000 <= ifol) .and. (ifol <= 99999)) write(ffile, '(a7, i5,a4)') "follow_", ifol, ".dat"
!write(ffile,'(a7,i1,a4)') "follow_", ifol, ".dat"
open(unit = 7, file = ffile, status = "UNKNOWN")

nbod0 = nbod

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, plist)
do i = 1, nbod

  plist(i) = i

end do
!$OMP END PARALLEL DO

call follow_plist(0,ifol,ifoln,plist,nbod0,tg,tstop)

ifoln = ifol

write(*,'(/,a)') ' Output format:'
write(*,'(a)')   ' 1  2   3  4  5    6      7      8     9     10  11  12  '
write(*,'(a,/)') ' t, id, a, e, inc, capom, omega, capm, peri, apo, M, obar'

tmax = t0
! 1    continue
do

  if(btest(iflgchk,0)) then ! bit 0 is set

    ierr = io_read_hdr(iu,t,nbod,nleft)
    if(ierr /= 0) exit !goto 2
    ierr = io_read_mass(tm,nbodm,mass,ium)
    if(ierr /= 0) exit !goto 2

  else

    ierr = io_read_hdr_r(iu,t,nbod,nleft)
    if(ierr /= 0) exit !goto 2
    ierr = io_read_mass_r(tm,nbodm,mass,ium)
    if(ierr /= 0) exit !goto 2

  end if


  if(nbodm /= nbod) then

    write(*,'(a,2i6)') ' Error 1: ', nbod, nbodm
    stop

  end if

  if(tm /= t) then

    write(*,'(a,2(1pe11.5))') ' Error 2: ', t, tm
    stop

  end if

  do while(t >= tg)

    call follow_plist(1,ifol,ifoln,plist,nbod0,tg,tstop)

  end do

  istep = 0

  do i = 2, nbod

    if(btest(iflgchk,0))  then ! bit 0 is set

      ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm)

    else

      ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm)

    end if

    if(ierr /= 0) exit !goto 2

    if(abs(id) == ifoln) then

      istep = 1
      inc = inc*dr
      obar = capom + omega
      if(obar >= TWOPI) obar = obar - TWOPI
      obar = obar*dr
      capom = capom*dr
      omega = omega*dr
      capm = capm*dr
      peri = a*(1.0_rk - e)
      apo = a*(1.0_rk + e)
      write(7,'(1x,e15.7,1x,i6,1x,f10.4,1x,f8.5,4(1x,f9.4),2(1x,f10.4),1e13.5,1x,f10.4)') t,ifoln,a,e,inc,capom,omega,capm,peri, &
                                                                                          apo,mass(abs(id)),obar
      tmax = t

    end if

  end do

  ! Did not find particle this times step
  if(istep == 0) exit !goto 2

!  goto 1
end do

! 2    continue

write(*,'(a,1pe11.5)') ' Tmax = ', tmax

end program follow_symba
!!
subroutine follow_plist(iflg,ifol,ifoln,plist,nbod,tg,tstop)
!--------------------------------------------------------------------
!
! Extracts information about why a particle was discarded from the
! discard_mass.dat file
!
!--------------------------------------------------------------------
use module_swift
implicit none

real(rk) :: tg,tstop
integer(ik) :: nbod
integer(ik) :: iflg,ifol,ifoln,plist(nbod)
integer(ik) :: ig,im,idum,i,ierr
integer(ik), save :: iwhy

if(iflg == 0) then

  open(unit = 2, file = 'discard_mass.dat', status = 'old', iostat = ierr)

  if(ierr /= 0) then

    write(*,'(a)') 'Could not open discard_mass.dat'
    tg = 5.0_rk*tstop
    return ! <==== NOTE!!!!

  end if

  read(2,*,iostat = ierr) tg, iwhy

  if(ierr /= 0) then

    write(*,*) 'Could not read discard_mass.dat'
    tg = 5.0_rk*tstop
    return ! <==== NOTE!!!!

  end if

  ifoln = ifol
  return ! <==== NOTE!!!!

endif

if(iwhy == 2) then

  read(2,*) idum, im
  read(2,fmt = '(1x)')
  read(2,fmt = '(1x)')
  read(2,*) idum,ig
  call left_reorder(ig,im,nbod,plist)

  do i = 1, 5

    read(2,fmt = '(1x)')

  end do

else

  read(2,*) idum, ig
  im = -1
  call left_reorder(ig,im,nbod,plist)
  read(2,fmt = '(1x)')
  read(2,fmt = '(1x)')

end if

read(2,*,iostat = ierr) tg, iwhy
if(ierr /= 0) tg = 5.0_rk*tstop

ifoln = plist(ifol)

return
end subroutine follow_plist
!!
subroutine left_reorder(ig,im,nbod,plist)
!-------------------------------------------------------------------
!
! Reorders plist, some how, I guess ...
!
!-------------------------------------------------------------------
use module_swift
implicit none

integer(ik) :: i, ig, im, nbod, plist(nbod)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(im, ig) SHARED(nbod, plist)
do i = 1, nbod

  if(plist(i) == ig) then

    if(im > 0) then

      plist(i) = im

    else

      plist(i) = -1

    end if

  end if

end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(ig) SHARED(nbod, plist)
do i = 1, nbod

  if(plist(i) > ig) plist(i) = plist(i) - 1

end do
!$OMP END PARALLEL DO

return
end subroutine left_reorder
