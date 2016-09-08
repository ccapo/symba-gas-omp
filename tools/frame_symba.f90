program frame_symba
!----------------------------------------------------------------------------
!				FRAME_SYMBA.F90
!----------------------------------------------------------------------------
!
! Extracts information on all particles in a binary file, for certain times,
! and writes to separate ascii files
!
!----------------------------------------------------------------------------
use module_swift
use module_interfaces
implicit none

real(rk) :: mass(NTPMAX), rdrag(NTPMAX)
real(rk) :: xh(NTPMAX), yh(NTPMAX), zh(NTPMAX)
real(rk) :: vxh(NTPMAX), vyh(NTPMAX), vzh(NTPMAX)

integer(ik) :: arg_num, iargc, system, status
integer(ik) :: nframe, nf
integer(ik) :: nbod, ierr, ialpha, nbodm
integer(ik) :: iflgchk, iu, nleft, i, j, id, ium, iur

real(rk) :: t0, tstop, dt, dtout, dtdump
real(rk) :: t, gm, mtiny, dtdelay = 0.1_rk

real(rk) :: rmin, rmax, rmaxu, qmin, rpl(NTPMAX), rhill(NTPMAX)
real(rk) :: a, e, inc, capom, omega, capm, j2rp2, j4rp4, temp
real(rk) :: mscale, tframe, dtframe, amin, amax, mpl, fg, fs, rplsml
character(len = 80) :: outfile, inparfile, inplfile, fopenstat
character(len = 80) :: buf, cmd, fn, tf, aminf, amaxf, mf, gf, sf, rf
character(len = 24) :: form = "UNFORMATTED", state = "OLD"
logical(lk) :: lclose, lgas

 1000   format(i5,1x,e13.5,1x,6(f10.5,1x))

! Get data for the run and the test particles
inparfile = 'param.in'
call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,lgas,mtiny,outfile,fopenstat)

! Get planet data
inplfile = 'pl.in'
call io_init_pl_symba(inplfile, lclose, iflgchk, nbod, mass, xh, yh, zh, vxh, vyh, vzh, rpl, rhill, j2rp2, j4rp4)

! Calculate the location of the last massive particle
call symba5_nbodm(nbod, mass, mtiny, nbodm)

! Get planetesimal radii
rdrag(1:nbodm) = 0.0_rk
open(unit = 10, file = 'drag_radius.dat', status = 'old')
do i = nbodm + 1, nbod
  read(10,'(i6,2(1x,1pe22.15))') id, rdrag(i), temp 
end do
close(unit = 10)

! Useful definitions
mscale = mass(1)*3.040432646e-06_rk
tframe = t0
nf = 0

! Get the number of command line arguments
arg_num = iargc()

! Read in the command line arguments
if(arg_num == 7) then

  call getarg(1, buf)
  read(buf,*) nframe

  call getarg(2, buf)
  read(buf,*) amin

  call getarg(3, buf)
  read(buf,*) amax

  call getarg(4, buf)
  read(buf,*) mpl

  call getarg(5, buf)
  read(buf,*) fg

  call getarg(6, buf)
  read(buf,*) fs

  call getarg(7, buf)
  read(buf,*) rplsml

else

  write(*, '(a)', advance = 'no') 'Input the number of frames: '
  read(*,*) nframe

  write(*, '(a)', advance = 'no') 'Input the Min. semi-major axis (AU): '
  read(*,*) amin

  write(*, '(a)', advance = 'no') 'Input the Max. semi-major axis (AU): '
  read(*,*) amax

  write(*, '(a)', advance = 'no') 'Input the mass for the embryo (M_Earth): '
  read(*,*) mpl

  write(*, '(a)', advance = 'no') 'Input the gas scale factor (MMSN): '
  read(*,*) fg

  write(*, '(a)', advance = 'no') 'Input the solid scale factor (MMSN): '
  read(*,*) fs

  write(*, '(a)', advance = 'no') 'Input the planetesimal radius (km): '
  read(*,*) rplsml

end if

if(nframe <= 1) nframe = 2
dtframe = (tstop - t0)/real(nframe - 1, rk)
if(modulo(dtframe,dtout) /= 0) then
  if( floor(dtframe/dtout) <= 0 ) then
    dtframe = dtout*(floor(dtframe/dtout) + 1)
  else
    dtframe = dtout*floor(dtframe/dtout)
  end if
end if
write(*,'(a,i5,a,1pe9.3,a,/)') 'Creating ', nframe, ' snapshots from t0 to tstop with step size of ', dtframe, ' years'

write(tf,'(1pe9.3)') t0
write(aminf,'(f4.1)') amin
write(amaxf,'(f4.1)') amax
write(mf,'(f5.2)') mpl
write(gf,'(f4.1)') fg
write(sf,'(f4.1)') fs
write(rf,'(1pe9.3)') rplsml

open(7, file = 'plsml.dat', status = 'unknown')
open(8, file = 'embryo.dat', status = 'unknown')

do i = 2, nbod

  gm = mass(1) + mass(i)
  call orbel_xv2el(xh(i), yh(i), zh(i), vxh(i), vyh(i), vzh(i), gm, ialpha, a, e, inc, capom, omega, capm)
  if(mass(i) < mtiny) then

    write(7, '(i5,8(1x,1pe22.15))') i, mass(i)/mscale, a, e, inc, omega, capom, capm, rdrag(i)

  else

    write(8, '(i5,8(1x,1pe22.15))') i, mass(i)/mscale, a, e, inc, omega, capom, capm, rdrag(i)

  end if

end do

close(unit = 7)

close(unit = 8)

if((0 <= nf) .and. (nf <= 9))       write(fn, '("000", i1)') nf
if((10 <= nf) .and. (nf <= 99))     write(fn, '("00",  i2)') nf
if((100 <= nf) .and. (nf <= 999))   write(fn, '("0",   i3)') nf
if((1000 <= nf) .and. (nf <= 9999)) write(fn, '(i4)') nf
write(tf, '(1pe9.3)') t0
write(cmd,*) './@plot_frame ' // trim(fn) // ' ' // trim(tf) // ' ' // trim(aminf) // ' ' // trim(amaxf) // &
             ' ' // trim(mf) // ' ' // trim(gf) // ' ' // trim(sf) // ' ' // trim(rf)
status = system(trim(cmd))

tframe = tframe + dtframe

iu = 20
ium = 30

if(btest(iflgchk,0)) then ! bit 0 is set

  !write(*,'(a)') 'Reading an FXDR binary file'
  call io_open_fxdr(outfile, 'r', .true., iu, ierr)
  call io_open_fxdr('mass.' // outfile, 'r', .true., ium, ierr)
  call io_open_fxdr('radius.' // outfile, 'r', .true., iur, ierr)

else if(btest(iflgchk,1)) then

  !write(*,'(a)') 'Reading an real*4 binary file'
  call io_open(iu, outfile, state, form, ierr)
  call io_open(ium, 'mass.' // outfile, state, form, ierr)
  call io_open(ium, 'radius.' // outfile, state, form, ierr)

else

  write(*,'(a)') 'ERROR: no binary file format specified in param file'
  stop

end if

outer: do

  open(7, file = 'plsml.dat', status = 'unknown')
  open(8, file = 'embryo.dat', status = 'unknown')

  inner: do

    if(btest(iflgchk,0)) then ! bit 0 is set

      ierr = io_read_hdr(iu, t, nbod, nleft)
      if(ierr /= 0) then

        write(*,'(a)') 'Did not find time!'
        stop

      end if

      ierr = io_read_mass(t, nbodm, mass, ium)
      if(ierr /= 0) then

        write(*,'(a)') 'Did not find time!'
        stop

      end if

      ierr = io_read_radius(t, nbodm, rdrag, iur)
      if(ierr /= 0) then

        write(*,'(a)') 'Did not find time!'
        stop

      end if

    else

      ierr = io_read_hdr_r(iu, t, nbod, nleft)
      if(ierr /= 0) then

        write(*,'(a)') 'Did not find time!'
        stop

      end if

      ierr = io_read_mass_r(t, nbodm, mass, ium)
      if(ierr /= 0) then

        write(*,'(a)') 'Did not find time!'
        stop

      end if

      ierr = io_read_mass_r(t, nbodm, rdrag, iur)
      if(ierr /= 0) then

        write(*,'(a)') 'Did not find time!'
        stop

      end if

    end if

    if(nbodm /= nbod) then

      write(*,'(a,2i7)') 'Error: ', nbod, nbodm
      stop

    end if

    do i = 2, nbod

      if(btest(iflgchk,0)) then ! bit 0 is set

        ierr = io_read_line(iu, id, a, e, inc, capom, omega, capm)

      else

        ierr = io_read_line_r(iu, id, a, e, inc, capom, omega, capm)

      end if

      if(ierr /= 0) then

        write(*,'(a)') 'Did not find time!'
        stop

      end if

      if(t == tframe) then

        if(mass(i) < mtiny) then

          write(7, '(i5,8(1x,1pe22.15))') i, mass(i)/mscale, a, e, inc, omega, capom, capm, rdrag(i)

        else

          write(8, '(i5,8(1x,1pe22.15))') i, mass(i)/mscale, a, e, inc, omega, capom, capm, rdrag(i)

        end if

      end if

    end do

    if(t == tframe) then

      close(unit = 7)
      close(unit = 8)

      nf = nf + 1
      if((1 <= nf) .and. (nf <= 9))       write(fn, '("000", i1)') nf
      if((10 <= nf) .and. (nf <= 99))     write(fn, '("00",  i2)') nf
      if((100 <= nf) .and. (nf <= 999))   write(fn, '("0",   i3)') nf
      if((1000 <= nf) .and. (nf <= 9999)) write(fn, '(i4)') nf

      write(tf, '(1pe9.3)') t
      write(cmd,*) './@plot_frame ' // trim(fn) // ' ' // trim(tf) // ' ' // trim(aminf) // ' ' // trim(amaxf) // &
                   ' ' // trim(mf) // ' ' // trim(gf) // ' ' // trim(sf) // ' ' // trim(rf)
      status = system(trim(cmd))

      tframe = tframe + dtframe

      if(tframe >= tstop) then

        status = system("rm plsml.dat embryo.dat")

        write(*,'(a)') '*** Make Movie ***'
        !write(cmd,'(a)') "./@make_movie"
        write(cmd,'(a15,f6.1,a)') "convert -delay ", 100.0_rk*dtdelay, " *.png movie.gif && rm *.png &> /dev/null"
        status = system(trim(cmd))

        !write(*,'(a)') '*** Remove PNG ***'
        !status = system("rm *.png")
        exit outer

      else

        cycle outer

      end if

    end if

  end do inner

end do outer

write(*,'(a)') '*** Done ***'

stop
end program frame_symba
