module io
! Module for io routines
!
! Author:  Hal Levison
! Date:    2/21/94
! Last revision: 2/24/94
use swift
use fxdr
use util
use orbel
implicit none

! Number of bytes in iflgchk
integer(ik), parameter :: IO_NBITS = 6

! Bit 0 set ==> write big binary data file
! Bit 1 set ==> write real*4 binary file rather than int*2: ignored if bit0 == F
! Bit 2 set ==> calc energy of system wrt time
! Bit 3 set ==> calc jacobi of the test particles
! Bit 4 set ==> check if particles are removed
! Bit 5 set ==> include J2 and J4 terms
! Bit 6 set ==> symba5

contains
!!
!!
  subroutine io_discard_mass(init,time,id,m1,r1,x1,y1,z1,vx1,vy1,vz1,iu,iwhy,fopenstat)
  !------------------------------------------------------------------------------------
  !				IO_DISCARD_MASS.F90
  !------------------------------------------------------------------------------------
  ! Write out information about a discarded massive body.
  !
  !             Input:
  !                 init          ==>  initiize flag if = 0 initialize and return
  !                                                     = 1 run through
  !                 id            ==> particle number (int scalar)
  !                 time          ==>  current time (real scalar)
  !                 m1            ==>  Mass of pl (real scalar)
  !                 r1            ==>  Radius of pl 2 (real scalar)
  !                 x1,y1,z1      ==>  current position of pl 1 in helio coord
  !                                    (real scalar)
  !                 vx1,vy1,vz1   ==>  current velocity of pl 1 in helio coord
  !                                    (real scalar)
  !                 iu            ==> IO unit (int scalar)
  !                 iwhy          ==> reason for discard (int scalar)
  !                 fopenstat     ==>  The status flag for the open
  !                                      statements of the output files.
  !                                          (character*80)
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    12/30/96
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: iwhy,iu,init,id
  real(rk) :: time
  real(rk) :: m1,r1
  real(rk) :: x1,y1,z1
  real(rk) :: vx1,vy1,vz1
  character(len = *) :: fopenstat

  ! Internals
  integer(ik) :: ierr

  !-----------------
  ! Executable code
  !-----------------

  if(init == 0) then

    ! Replaced fopenstat with 'unknown' for initialization call
    !call io_open(iu,'discard_mass.dat',fopenstat,'FORMATTED',ierr)
    call io_open(iu,'discard_mass.dat','unknown','FORMATTED',ierr)

    ! If there was an error and fopenstat='append' then
    ! try to open as **unknown** (CCC -- 09/05/07)
    if(ierr /= 0) then

      if((fopenstat == 'append') .or. (fopenstat == 'APPEND')) then

        call io_open(iu,'discard_mass.dat','unknown','FORMATTED',ierr)

      end if

    end if

    if(ierr /= 0) then

      write(*,'(a)') ' SWIFT ERROR: in io_discard_mass: '
      write(*,'(a)') '  Could not open discard output file'
      call util_exit(FAILURE)

    end if

    return ! <== NOTE!!!!

  else

    call io_open(iu,'discard_mass.dat','append','FORMATTED',ierr)

  end if

  write(iu,'(1x,1pe22.15,1x,i4)') time,iwhy
  write(iu,'("-1",1x,i6,1x,2(1x,1pe22.15))') id,m1,r1
  write(iu,'(3(1x,1pe22.15))') x1,y1,z1
  write(iu,'(3(1x,1pe22.15))') vx1,vy1,vz1

  close(unit = iu)

  return
  end subroutine io_discard_mass
!!
  subroutine io_discard_merge(time,ip1,ip2,m1,r1,x1,y1,z1,vx1,vy1,vz1, &
             m2,r2,x2,y2,z2,vx2,vy2,vz2,mn,rn,xn,yn,zn,vxn,vyn,vzn)
  !-------------------------------------------------------------------------------
  !				IO_DISCARD_MERGE.F90
  !-------------------------------------------------------------------------------
  ! Write out information about a merger.
  !
  !             Input:
  !                 time          ==>  current time (real scalar)
  !                 ip1,ip2       ==>  planets to merge (real scalar)
  !                 m1            ==>  Mass of pl 1 (real scalar)
  !                 r1            ==>  Radius of pl 1 (real scalar)
  !                 x1,y1,z1      ==>  current position of pl 1 in helio coord
  !                                    (real arrays)
  !                 vx1,vy1,vz1   ==>  current velocity of pl 1 in helio coord
  !                                    (real arrays)
  !                 m2            ==>  Mass of pl 2 (real scalar)
  !                 r2            ==>  Radius of pl 2 (real scalar)
  !                 x2,y2,z2      ==>  current position of pl 2 in helio coord
  !                                    (real arrays)
  !                 vx2,vy2,vz2   ==>  current velocity of pl 2 in helio coord
  !                                    (real arrays)
  !                 mn            ==>  Mass of new pl  (real scalar)
  !                 rn            ==>  Radius of new pl (real scalar)
  !                 xn,yn,zn      ==>  current position of new pl in helio coord
  !                                    (real arrays)
  !                 vxn,vyn,vzn   ==>  current velocity of new pl in helio coord
  !                                    (real arrays)
  !                 nleft           ==>  number of active test bodies(int scalar)
  !
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    12/30/96
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: ip1,ip2
  real(rk) :: time
  real(rk) :: m1,r1
  real(rk) :: x1,y1,z1
  real(rk) :: vx1,vy1,vz1
  real(rk) :: m2,r2
  real(rk) :: x2,y2,z2
  real(rk) :: vx2,vy2,vz2
  real(rk) :: mn,rn
  real(rk) :: xn,yn,zn
  real(rk) :: vxn,vyn,vzn

  ! Internals
  integer(ik) :: ierr,iu

  !-----------------
  ! Executable code
  !-----------------

  iu = 40

  call io_open(iu,'discard_mass.dat','append','FORMATTED',ierr)

  write(iu,'(1x,1pe22.15,"  2")') time

  write(iu,'("-1",1x,i6,1x,2(1x,1pe22.15))') ip1,m1,r1
  write(iu,'(3(1x,1pe22.15))') x1,y1,z1
  write(iu,'(3(1x,1pe22.15))') vx1,vy1,vz1

  write(iu,'("-1",1x,i6,1x,2(1x,1pe22.15))') ip2,m2,r2
  write(iu,'(3(1x,1pe22.15))') x2,y2,z2
  write(iu,'(3(1x,1pe22.15))') vx2,vy2,vz2

  write(iu,'("+1",1x,i6,1x,2(1x,1pe22.15))') ip1,mn,rn
  write(iu,'(3(1x,1pe22.15))') xn,yn,zn
  write(iu,'(3(1x,1pe22.15))') vxn,vyn,vzn

  close(unit = iu)

  return
  end subroutine io_discard_merge
!!
  subroutine io_dump_param(dparfile,t,tstop,dt,dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
  !-----------------------------------------------------------------------------------------------------
  !					IO_DUMP_PARAM.F90
  !-----------------------------------------------------------------------------------------------------
  ! IO_DUMP_PARAM dumps out the parameters for the integration.
  !
  !      Input:
  !       dparfile      ==>  Name of file to write to (character*80)
  !            t0       ==> Initial time (real scalar)
  !            tstop    ==> final time (real scalar)
  !            dt       ==> time step  (real scalar)
  !            dtout    ==> time between binary outputs (real scalar)
  !            dtdump   ==> time between dumps  (real scalar)
  !            iflgchk  ==>  =0 don't run diagnostic routines
  !                         !=0 run them
  !      rmin,rmax      ==>  maximum and min distance from Sun
  !                                if <0  then don't check
  !                                    (real scalar)
  !       rmaxu         ==>  maximum distance from Sun in not bound
  !                                 if <0  then don't check
  !                                      (real scalar)
  !       qmin          ==> Smallest perihelion distance
  !                                 if <0  then don't check
  !                                      (real scalar)
  !       lclose        ==> .true. --> discard particle if it gets
  !                                    too close to a planet. Read in that
  !                                    distance in io_init_pl
  !                                      (logical*2 scalar)
  !       outfile       ==>  Name of binary output file (character*80)
  !
  ! Remarks:
  ! Authors:  Martin Duncan
  ! Date:    3/2/93
  ! Last revision:  5/10/94 HFL
  implicit none

  ! Input variables
  logical(lk) :: lclose
  integer(ik) :: iflgchk
  real(rk) :: t, tstop, dt, dtout, dtdump, rmin, rmax, rmaxu, qmin
  character(len = 24) :: outfile, dparfile

  ! Internal variables
  integer(ik) :: ierr

  ! Simulation variables
  logical(lk) :: lenergy, ljacobi
  real(rk) :: t0
  character(len = 24) :: outtype, fopenstat

  ! Stellar variables
  logical(lk) :: loblate
  real(rk) :: j2, j4

  ! SyMBA variable
  real(rk) :: mtiny

  ! Gas disk variables
  logical(lk) :: lgas
  integer(ik) :: ngap, igap
  real(rk) :: deng0s, gpower, zscale, rgi, rgf, tgdecay, ca, ce, fwgap

  namelist / sim_params / t0, tstop, dt, dtout, dtdump, outfile, outtype, fopenstat, lenergy, ljacobi
  namelist / stellar_params / loblate, j2, j4, lclose, rmin, rmax, rmaxu, qmin
  namelist / symba_params / mtiny
  namelist / gas_params / lgas, deng0s, gpower, zscale, rgi, rgf, tgdecay, ca, ce, ngap, igap, fwgap

  !-----------------
  ! Executable code
  !-----------------

  ! Open and read in parameters
  open(unit = 7, file = "param.in", status = 'old')
  read(unit = 7, nml = sim_params)
  read(unit = 7, nml = stellar_params)
  read(unit = 7, nml = symba_params)
  read(unit = 7, nml = gas_params)
  close(unit = 7)

  ! Ensure integration time, along with output intervals are integer multiples of the timestep
  if(modulo(tstop,dt) /= 0) tstop = dt*(floor(tstop/dt) + 1)

  if(modulo(dtout,dt) /= 0) dtout = dt*(floor(dtout/dt) + 1)

  if(modulo(dtdump,dt) /= 0) dtdump = dt*(floor(dtdump/dt) + 1)

  ! Open parameter data file for the dump
  t0 = t
  fopenstat = 'APPEND'
  open(unit = 8, file = dparfile, status = 'unknown')
  write(unit = 8, nml = sim_params)
  write(unit = 8, nml = stellar_params)
  write(unit = 8, nml = symba_params)
  write(unit = 8, nml = gas_params)
  close(unit = 8)

  ! Put mtiny and mpuny back into simulation units
  mtiny = mtiny*MEARTH

  return
  end subroutine io_dump_param
!!
  subroutine io_dump_pl_symba(dplfile,nbod,mass,xh,yh,zh,vxh,vyh,vzh, &
                              lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
  !-----------------------------------------------------------------------------
  !			IO_DUMP_PL_SYMBA.F90
  !-----------------------------------------------------------------------------
  ! Dumps the data for the Sun and planets
  !
  !             Input:
  !                 dplfile       ==>  Name of file to write to (character*80)
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 xh,yh,zh      ==>  initial position in Helio coord
  !                                    (real arrays)
  !                 vxh,vyh,vzh   ==>  initial position in Helio coord
  !                                    (real arrays)
  !                 lclose        ==> .true. --> discard particle if it gets
  !                                    too close to a planet. Read in that
  !                                    distance in io_init_pl
  !                                      (logical*2 scalar)
  !                 iflgchk       ==>  bit 5 set ==>  include J2 and J4 terms
  !                 rpl           ==>  physical size of planet
  !                                    (real array)
  !                 rhill         ==>  size of planet's hills sphere
  !                                    (real array)
  !                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
  !                                     (real scalars)
  !
  !
  ! Remarks: Based on io_dump_pl.f
  ! Authors:  Hal Levison
  ! Date:    1/8/97
  ! Last revision:
  implicit none

  ! Input
  logical(lk) :: lclose
  integer(ik) :: nbod,iflgchk
  real(rk) :: mass(nbod),rpl(nbod),j2rp2,j4rp4
  real(rk) :: xh(nbod),yh(nbod),zh(nbod),rhill(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  character(len = *) :: dplfile

  ! Internal
  integer(ik) :: j,ierr

  !-----------------
  ! Executable code
  !-----------------

  call io_open(7,dplfile,'unknown','formatted',ierr)

  write(7,'(1x,i6)') nbod

  if(btest(iflgchk,5)) then ! bit 5 is set

    write(7,'(3(1x,1pe22.15))') mass(1),j2rp2,j4rp4
    write(7,'(3(1x,1pe22.15))') xh(1),yh(1),zh(1)
    write(7,'(3(1x,1pe22.15))') vxh(1),vyh(1),vzh(1)

  else

    write(7,'(1x,1pe22.15)') mass(1)
    write(7,'(3(1x,1pe22.15))') xh(1),yh(1),zh(1)
    write(7,'(3(1x,1pe22.15))') vxh(1),vyh(1),vzh(1)

  end if

  if(lclose) then

    do j = 2, nbod

      write(7,'(3(1x,1pe22.15))') mass(j),rhill(j),rpl(j)
      write(7,'(3(1x,1pe22.15))') xh(j),yh(j),zh(j)
      write(7,'(3(1x,1pe22.15))') vxh(j),vyh(j),vzh(j)

    end do

  else

    do j = 2, nbod

      write(7,'(2(1x,1pe22.15))') mass(j),rhill(j)
      write(7,'(3(1x,1pe22.15))') xh(j),yh(j),zh(j)
      write(7,'(3(1x,1pe22.15))') vxh(j),vyh(j),vzh(j)

    end do

  end if

  close(unit = 7)

  return
  end subroutine io_dump_pl_symba
!!
  subroutine io_energy_write(i1st,t,energy,eltot,iu,fopenstat)
  !--------------------------------------------------------------------------
  !			IO_ENERGY_WRITE.F90
  !--------------------------------------------------------------------------
  ! Does the write for anal_jacobi_write
  !
  !      Input:
  !            i1st           ==>  =0 if first write, =1 if not (int scalar)
  !            t              ==>  current time (real scalar)
  !            energy         ==>  Total energy
  !            eltot          ==>  components of total angular momentum
  !                               (real array)
  !            iu             ==>  unit to write to
  !            fopenstat      ==>  The status flag for the open
  !                                statements of the output files.
  !                                          (character*80)
  !
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    2/21/94
  ! Last revision: 3/4/94
  implicit none

  ! Inputs:
  integer(ik) :: iu,i1st
  real(rk) :: t,energy,eltot(NDIM)
  character(len = *) :: fopenstat

  ! Internals
  integer(ik) :: ierr
  character(len = 24) :: outfile = "energy.dat"

  !-----------------
  ! Executable code
  !-----------------

  if(i1st == 0) then

    call io_open(iu,outfile,fopenstat,'FORMATTED',ierr)

    if(ierr /= 0) then

      write(*,'(a)') ' SWIFT ERROR: in anal_energy_write'
      write(*,'(a)') '  Could not open energy.dat'
      call util_exit(1)

    end if

  else

    call io_open(iu,outfile,'append','FORMATTED',ierr)

  end if

  write(iu,'(1x,1pe13.6,4(1x,1pe22.15))') t,energy,eltot

  close(unit = iu)

  return
  end subroutine io_energy_write
!!
  subroutine io_init_param(infile,t0,tstop,dt,dtout,dtdump,iflgchk, &
             rmin,rmax,rmaxu,qmin,lclose,lgas,mtiny,outfile,fopenstat)
  !-----------------------------------------------------------------------------
  !				IO_INIT_PARAM.F90
  !-----------------------------------------------------------------------------
  ! INIT_PARAM reads in the parameters for the integration.
  !
  !      Input:
  !            infile   ==> File name to read from (character*80)
  !
  !      Output:
  !            t0       ==> Initial time (real scalar)
  !            tstop    ==> final time (real scalar)
  !            dt       ==> time step  (real scalar)
  !            dtout    ==> time between binary outputs (real scalar)
  !            dtdump   ==> time between dumps  (real scalar)
  !            iflgchk  ==>  =0 don't run diagnostic routines
  !                          bit 0 set ==>  write xdr*4 binary data file
  !                          bit 1 set ==>  write real*4 binary file
  !                          bit 2 set ==>  calc energy of system wrt time
  !                          bit 3 set ==>  calc jacobi of the test particles
  !                          bit 4 set ==>  check if particles are removed
  !                          bit 5 set ==>  include J2 and J4 terms
  !      rmin,rmax      ==>  maximum and min distance from Sun
  !                                if <0  then don't check
  !                                    (real scalar)
  !      rmaxu          ==>  maximum distance from Sun in not bound
  !                                 if <0  then don't check
  !                                      (real scalar)
  !       qmin          ==> Smallest perihelion distance
  !                                 if <0  then don't check
  !                                      (real scalar)
  !       lclose        ==> .true. --> discard particle if it gets
  !                                    too close to a planet. Read in that
  !                                    distance in io_init_pl
  !                                      (logical*2 scalar)
  !       outfile       ==>  Name of binary output file (character*80)
  !       fopenstat     ==>  The status flag for the open statements of the
  !                          output files.  Must be one of the following:
  !                                 new      (die if the file exists)
  !                                 append   (add to what is there)
  !                                 unknown  (just write over what is there)
  !                                 (character*80)
  !
  !
  ! Remarks:
  ! Authors:  Martin Duncan
  ! Date:    3/2/93
  ! Last revision:  5/10/94  HFL
  ! use omp_lib
  implicit none

  ! Input variables
  character(len = *) :: infile

  ! Output variables
  logical(lk) :: lclose, lgas
  integer(ik) :: iflgchk
  real(rk) :: t0, tstop, dt
  real(rk) :: dtout, dtdump
  real(rk) :: rmin, rmax, rmaxu, qmin, mtiny
  character(len = 24) :: outfile, fopenstat

  ! Internal variables
  logical(lk) :: lflg(0:IO_NBITS - 1), lenergy, ljacobi, loblate
  integer(ik) :: i, ierr, ngap, igap
  real(rk) :: j2, j4, deng0s, gpower, zscale, rgi, rgf, tgdecay, ca, ce, fwgap
  character(len = 24) :: outtype

  namelist / sim_params / t0, tstop, dt, dtout, dtdump, outfile, outtype, fopenstat, lenergy, ljacobi
  namelist / stellar_params / loblate, j2, j4, lclose, rmin, rmax, rmaxu, qmin
  namelist / symba_params / mtiny
  namelist / gas_params / lgas, deng0s, gpower, zscale, rgi, rgf, tgdecay, ca, ce, ngap, igap, fwgap

  !-----------------
  ! Executable code
  !-----------------

  ! Open and read in parameters
  open(unit = 7, file = infile, status = 'old')
  read(unit = 7, nml = sim_params)
  read(unit = 7, nml = stellar_params)
  read(unit = 7, nml = symba_params)
  read(unit = 7, nml = gas_params)
  close(unit = 7)

  ! Ensure integration time, along with output intervals are integer multiples of the timestep
  if(modulo(tstop,dt) /= 0) then

    tstop = dt*(floor(tstop/dt) + 1)
    write(*,'(a)') ' *** tstop modified to be an integer multiple of dt ***'

  end if

  if(modulo(dtout,dt) /= 0) then

    dtout = dt*(floor(dtout/dt) + 1)
    write(*,'(a)') ' *** dtout modified to be an integer multiple of dt ***'

  end if

  if(modulo(dtdump,dt) /= 0) then

    dtdump = dt*(floor(dtdump/dt) + 1)
    write(*,'(a)') ' *** dtdump modified to be an integer multiple of dt ***'

  end if

  write(*,'(/,a)') ' Simulation parameters:'
  write(*,'(a)') ' t0 (years), tstop (years), dt (years):'
  write(*,'(2(1pe12.5),e12.5,/)') t0, tstop, dt
  write(*,'(a)') ' dtout (years), dtdump (years):'
  write(*,'(2(1pe12.5),/)') dtout, dtdump

  ! Set logical flags
  lflg(0:IO_NBITS - 1) = .false.
  if((outtype == "XDR") .or. (outtype == "xdr")) then

    lflg(0) = .true.
    lflg(1) = .false.

  else ! Default outtype = "REAL"

    lflg(0) = .false.
    lflg(1) = .true.

  end if

  if(lenergy) lflg(2) = .true.
  if(ljacobi) lflg(3) = .true.
  if(lclose) lflg(4) = .true.
  if(loblate) lflg(5) = .true.

  iflgchk = 0
  do i = 0, IO_NBITS - 1

    if(lflg(i)) iflgchk = ibset(iflgchk,i)

  end do

  write(*,'(a)') ' Logical flags:'
  write(*,'(6g2.0,a3,i3,/)') (lflg(i), i = IO_NBITS - 1, 0, -1),' = ', iflgchk

  if(btest(iflgchk,0) .and. btest(iflgchk,1)) then

    write(*,'(a)') ' SWIFT ERROR: in io_init_param:'
    write(*,'(a)') '  Invalid logical flags'
    write(*,'(a)') '  You cannot request that both a real and an integer binary file be written'
    call util_exit(1)

  end if

  if(btest(iflgchk,4)) then ! bit 4 is set

    write(*,'(a)') ' rmin (AU), rmax (AU), rmaxu (AU), qmin (AU), lclose:'
    write(*,'(4(1pe13.5),g2.0,/)') rmin, rmax, rmaxu, qmin, lclose

  else

    rmin = -1.0_rk
    rmax = -1.0_rk
    rmaxu = -1.0_rk
    qmin = -1.0_rk
    lclose = .false.

  end if

  if(btest(iflgchk,0) .and. btest(iflgchk,1)) write(*,'(a,a,/)') ' Outfile file is: ', outfile

  if((fopenstat(1:3) /= 'new') .and. (fopenstat(1:3) /= 'NEW') .and.         &
     (fopenstat(1:7) /= 'unknown') .and. (fopenstat(1:7) /= 'UNKNOWN') .and. &
     (fopenstat(1:6) /= 'append') .and. (fopenstat(1:6) /= 'APPEND')) then

    write(*,'(a)')   ' SWIFT ERROR: in io_init_param:'
    write(*,'(a,a)') '  Invalid status flag: ', fopenstat(1:7), ':'
    call util_exit(1)

  end if

  ! Display gas parameters
  if(lgas) then

    write(*,'(a,/)')       ' ---------------------------------------------'
    write(*,'(a)')         ' Gas Disk parameters:              '
    write(*,'(a,1pe11.5)') ' Density of gas at 1 AU (g/cm^3) = ', deng0s
    write(*,'(a,1pe11.5)') ' Power law index                 = ', gpower
    write(*,'(a,1pe11.5)') ' Gas scale height at 1 AU (AU)   = ', zscale
    write(*,'(a,1pe11.5)') ' Inner edge of the gas disk (AU) = ', rgi
    write(*,'(a,1pe11.5)') ' Outer edge of the gas disk (AU) = ', rgf
    write(*,'(a,1pe11.5)') ' Gas decay time scale (year)     = ', tgdecay
    write(*,'(a,1pe11.5)') ' Type-I drag efficiency [c_a]    = ', ca
    write(*,'(a,1pe11.5)') ' Type-I drag efficiency [c_e]    = ', ce
    write(*,'(/,a,/)')     ' ---------------------------------------------'

  end if

  ! Put mtiny in simulation units
  mtiny = mtiny*mearth

  ! Define the maximum number of threads
  nthreads = 2                        ! In the *serial* case
  ! nthreads = omp_get_max_threads() ! In the *parallel* case

  return
  end subroutine io_init_param
!!
  subroutine io_init_pl_symba(infile,lclose,iflgchk,nbod,mass, &
             xh,yh,zh,vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)
  !-----------------------------------------------------------------------------
  !				IO_INIT_PL_SYMBA.F90
  !-----------------------------------------------------------------------------
  ! IO_INIT_PL_SYMBA reads in the data for the Sun and planets for
  ! symba routines
  !
  !             Input:
  !                 infile        ==> File name to read from (character*80)
  !                 lclose        ==> .true. --> discard particle if it gets
  !                                    too close to a planet. Read in that
  !                                    distance in io_init_pl_symba
  !                                      (logical*2 scalar)
  !                 iflgchk        ==>  bit 5 set ==>  include J2 and J4 terms
  !
  !             Output:
  !                 nbod          ==>  number of massive bodies (int scalar)
  !                 mass          ==>  mass of bodies (real array)
  !                 xh,yh,zh      ==>  initial position in Helio coord
  !                                    (real arrays)
  !                 vxh,vyh,vzh   ==>  initial position in Helio coord
  !                                    (real arrays)
  !                 rpl           ==>  physical size of planet
  !                                    (real array)
  !                 rhill         ==>  size of planet's hills sphere
  !                                    (real array)
  !                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
  !                                     (real scalars)
  !
  ! Remarks: Based on io_init_pl
  ! Authors:  Hal Levison
  ! Date:    11/21/96
  ! Last revision: 1/10/97
  implicit none

  ! Input
  logical(lk) :: lclose
  integer(ik) :: iflgchk
  character(len = 24) :: infile

  ! Output
  integer(ik) :: nbod
  real(rk) :: mass(NTPMAX),rpl(NTPMAX),j2rp2,j4rp4
  real(rk) :: xh(NTPMAX),yh(NTPMAX),zh(NTPMAX),rhill(NTPMAX)
  real(rk) :: vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)

  ! Internal
  integer(ik) :: j,ierr,ibad
  real(rk) :: r2hill(NTPMAX),rhrat

  !-----------------
  ! Executable code
  !-----------------

  write(*,'(a22,a24,/)') ' Planet data file is: ', infile
  call io_open(7,infile,'old','formatted',ierr)

  ! Read number of planets
  read(7,*) nbod

  if(nbod > NTPMAX) then

    write(*,'(a)') ' SWIFT ERROR: in io_init_pl_symba: '
    write(*,'(a)') '  The number of massive bodies,',nbod,','
    write(*,'(a)') '  is too large, it must be less than',NTPMAX
    call util_exit(1)

  end if

  write(*,'(a,i6,/)') ' Total number of bodies (including the Sun) = ', nbod

  ! For each planet read mass, heliocentric position and velocity
  if(btest(iflgchk,5)) then ! bit 5 is set

    read(7,'(3(1x,1pe22.15))') mass(1), j2rp2, j4rp4
    read(7,'(3(1x,1pe22.15))') xh(1), yh(1), zh(1)
    read(7,'(3(1x,1pe22.15))') vxh(1), vyh(1), vzh(1)

  else

    read(7,'(1x,1pe22.15)') mass(1)
    read(7,'(3(1x,1pe22.15))') xh(1), yh(1), zh(1)
    read(7,'(3(1x,1pe22.15))') vxh(1), vyh(1), vzh(1)
    j2rp2 = 0.0_rk
    j4rp4 = 0.0_rk

  end if

  rpl(1) = 0.0_rk
  rhill(1) = 0.0_rk

  if((xh(1) /= 0.0_rk) .or. (yh(1) /= 0.0_rk) .or. (zh(1) /= 0.0_rk) .or. &
     (vxh(1) /= 0.0_rk) .or. (vyh(1) /= 0.0_rk) .or. (vzh(1) /= 0.0_rk)) then

    write(*,'(a)') ' SWIFT ERROR in io_init_pl_symba:'
    write(*,'(a)') '  Input MUST be in heliocentric coordinates'
    write(*,'(a)') '  Position and Velocity of First Massive body MUST be zero'
    call util_exit(1)

  end if

  if(lclose) then

    do j = 2, nbod

      read(7,'(3(1x,1pe22.15))') mass(j), rhill(j), rpl(j)
      read(7,'(3(1x,1pe22.15))') xh(j), yh(j), zh(j)
      read(7,'(3(1x,1pe22.15))') vxh(j), vyh(j), vzh(j)

    end do

  else

    do j = 2, nbod

      read(7,'(3(1x,1pe22.15))') mass(j), rhill(j)
      read(7,'(3(1x,1pe22.15))') xh(j), yh(j), zh(j)
      read(7,'(3(1x,1pe22.15))') vxh(j), vyh(j), vzh(j)

    end do

  end if

  close(unit = 7)

  ! Check to see if the hills spheres are alright
  !call util_hills(nbod, mass, xh, yh, zh, vxh, vyh, vzh, r2hill)
  !ibad = 0
  !do j = 2, nbod
  !
  !  rhrat = rhill(j)/sqrt(r2hill(j))
  !  if((rhrat > 2.0_rk) .or. (rhrat < 0.5_rk) ) ibad = ibad + 1
  !
  !end do
  !
  !if(ibad /= 0) then
  !
  !  write(*,'(a)') " Warning in io_init_pl_symba:"
  !  write(*,'(a39,i8,a8)') "  Hill's spheres are not consistent on ", ibad, " objects"
  !
  !end if

  return
  end subroutine io_init_pl_symba
!!
  subroutine io_open(iu,fname,fopenstat,format,ierr)
  !-------------------------------------------------------------------------
  !				IO_OPEN.F90
  !-------------------------------------------------------------------------
  ! open files
  !
  !             Input:
  !                 iu              ==>  unit number (integer scalar)
  !                 fname           ==>  file name (character*80)
  !                 fopenstat       ==>  The status flag for the open
  !                                      statements of the output files.
  !                                          (character*80)
  !                 format          ==>  format string (character*80)
  !             Output:
  !                 ierr            ==>  output from iostat
  !
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    3/3/94
  ! Last revision: 1/30/98
  implicit none

  ! Inputs:
  integer(ik) :: iu
  character(len = *) :: fname, fopenstat, format

  ! Outputs:
  integer(ik) :: ierr

  !-----------------!
  ! Executable code !
  !-----------------!

  if((fopenstat == 'append') .or. (fopenstat == 'APPEND')) then

    open(unit = iu, file = fname, status = 'old', POSITION = 'append', FORM = format, IOSTAT = ierr)

    if(ierr /= 0) then

      write(*,*) 'Warning:  Could not open ', fname, ' with position = append.'
      write(*,*) '          Will open as status = new'

      open(unit = iu, file = fname, status = 'new', form = format, iostat = ierr)

    end if

  else

    open(unit = iu, file = fname, status = fopenstat, form = format, iostat = ierr)

  end if


  return
  end subroutine io_open
!!
  subroutine io_open_fxdr(fname,fopenstat,lflg,iu,ierr)
  !-------------------------------------------------------------------------
  !				IO_OPEN_FXDR.F90
  !-------------------------------------------------------------------------
  ! open files using the fxdr routines
  !
  ! Input:
  !   fname           ==>  file name (character*80)
  !   fopenstat       ==>  The status flag for the open
  !                        statements of the output files.
  !                        (character*1)
  !   lflg            ==>  if .TRUE., then FXDR routines return even if
  !                        there is an I/O error.  If .FALSE., then routines
  !                        halt on I/O error (logical scalar)
  !
  ! Output:
  !   iu              ==>  unit number (integer scalar)
  !   ierr            ==>  output from iostat
  !
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    11/3/99
  ! Last revision:
  implicit none

  ! Inputs:
  logical(lk) :: lflg
  integer(ik) :: iu
  character(len = *) :: fname
  character(len = 1) :: fopenstat

  ! Outputs:
  integer(ik) :: ierr

  ! Internals:

  !-----------------!
  ! Executable code !
  !-----------------!

  iu = initxdr(fname, fopenstat, lflg)

  if(iu > 0) then

    ierr = 0

  else

    ierr = iu

  end if

  return
  end subroutine io_open_fxdr
!!
  function io_read_hdr(iu,time,nbod,nleft) result(read_hdr)
  !-------------------------------------------------------------------------
  !				IO_READ_HDR.F90
  !-------------------------------------------------------------------------
  ! NEW VERSION OF THIS, USES FXDR
  !
  !             Input:
  !                 iu            ==> unit number to write to
  !             Output:
  !                 time          ==> current time (real scalar)
  !                 nbod          ==> number of massive bodies (int scalar)
  !                 nleft         ==> number of active tp (int scalar)
  !
  !             Returns:
  !               io_read_hdr     ==> = 0 read ok
  !                                   !=0 read failed, set to iostat variable
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    11/2/99
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: iu

  ! Output
  integer(ik) :: nbod,nleft,read_hdr
  real(rk) :: time

  ! Internals
  integer(ik) :: nn(2)
  integer(ik) :: ierr
  real(real4) :: ttmp

  !-----------------
  ! Executable code
  !-----------------

  ierr = ixdrreal( iu, ttmp )
  !ierr = ixdrdouble( iu, time )
  read_hdr = ierr

  if(ierr /= 0) return

  ierr = ixdrimat( iu, 2, nn )
  read_hdr = ierr

  if(ierr /= 0) return

  nbod = nn(1)
  nleft = nn(2)
  time = ttmp

  return
  end function io_read_hdr
!!
  function io_read_hdr_r(iu,time,nbod,nleft) result(read_hdr_r)
  !-------------------------------------------------------------------------
  !			IO_READ_HDR_R.F90
  !-------------------------------------------------------------------------
  ! Read in header part of the real*4 file
  !
  !             Input:
  !                 iu           ==> unit number to write to
  !             Output:
  !                 time         ==> current time (real scalar)
  !                 nbod         ==> number of massive bodies (int scalar)
  !                 nleft        ==> number of active tp (int scalar)
  !
  !             Returns:
  !               io_read_hdr_r  ==> = 0 read ok
  !                                   != 0 read failed is set to iostat variable
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    2/22/94
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: iu

  ! Output
  integer(ik) :: nbod,nleft,read_hdr_r
  real(rk) :: time

  ! Internals
  integer(ik) :: ierr
  integer(integer2) :: nleft2,nbod2
  real(real4) :: ttmp

  !-----------------
  ! Executable code
  !-----------------

  read(iu,iostat=ierr) ttmp,nbod2,nleft2
  read_hdr_r = ierr
  if(ierr /= 0) return

  nbod = nbod2
  nleft = nleft2
  time = ttmp

  return
  end function io_read_hdr_r
!!
  function io_read_line(iu,id,a,e,inc,capom,omega,capm) result(read_line)
  !-------------------------------------------------------------------------
  !				IO_READ_LINE.F90
  !-------------------------------------------------------------------------
  ! NEW VERSION OF THIS, USES FXDR
  !
  !      Input:
  !            iu       ==> unit number to write to
  !      Output:
  !            id       ==> Particle id number (int scalar)
  !            a        ==> semi-major axis or pericentric distance if a parabola
  !                          (real scalar)
  !            e        ==> eccentricity (real scalar)
  !            inc      ==> inclination  (real scalar)
  !            capom    ==> longitude of ascending node (real scalar)
  !            omega    ==> argument of perihelion (real scalar)
  !            capm     ==> mean anomoly(real scalar)
  !       Returns:
  !      io_read_line    ==>   =0 read ok
  !                           !=0 read failed is set to iostat variable
  !
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    11/2/99
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: iu

  ! Output:
  integer(ik) :: id, read_line
  real(rk) :: a,e,inc,capom,omega,capm

  ! Internals
  integer(ik) :: ierr
  real(real4) :: orbel(6)
  !real(rk) :: orbel(6)

  !-----------------
  ! Executable code
  !-----------------

  ierr = ixdrint(iu, id)
  read_line = ierr
  if(ierr /= 0) return

  ierr = ixdrrmat( iu, 6, orbel )
  !ierr = ixdrdmat( iu, 6, orbel )
  read_line = ierr
  if(ierr /= 0) return

  a = orbel(1)
  e = orbel(2)
  inc = orbel(3)
  capom = orbel(4)
  omega = orbel(5)
  capm = orbel(6)

  return
  end function io_read_line
!!
  function io_read_line_r(iu,id,a,e,inc,capom,omega,capm) result(read_line_r)
  !-------------------------------------------------------------------------
  !			IO_READ_LINE_R.F90
  !-------------------------------------------------------------------------
  ! Read one line from real*4 binary file.
  !
  ! Input:
  !   iu     ==> unit number to write to
  !
  ! Output:
  !   a      ==> semi-major axis or pericentric distance if a parabola
  !              (real scalar)
  !   e      ==> eccentricity (real scalar)
  !   inc    ==> inclination  (real scalar)
  !   capom  ==> longitude of ascending node (real scalar)
  !   omega  ==> argument of perihelion (real scalar)
  !   capm   ==> mean anomoly(real scalar)
  !
  ! Returns:
  !   read_line_r ==> = 0 read ok
  !                  != 0 read failed is set to iostat variable
  !
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    2/22/94
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: iu

  ! Output:
  integer(ik) :: id, read_line_r
  real(rk) :: a,e,inc,capom,omega,capm

  ! Internals
  integer(ik) :: ierr
  integer(integer2) :: id2
  real(real4) :: a4,e4,inc4,capom4,omega4,capm4

  !-----------------
  ! Executable code
  !-----------------

  read(iu, iostat = ierr) id2, a4, e4, inc4, capom4, omega4, capm4
  read_line_r = ierr
  if(ierr /= 0) return

  id = id2
  a = a4
  e = e4
  inc = inc4
  capom = capom4
  capm = capm4
  omega = omega4

  return
  end function io_read_line_r
!!
  function io_read_mass(time,nbod,mass,iu) result(read_mass)
  !-------------------------------------------------------------------------
  !				IO_READ_MASS.F90
  !-------------------------------------------------------------------------
  ! NEW VERSION OF THIS, USES FXDR
  ! Read in the mass file.
  !
  ! Output:
  !  time  ==>  current time (real scalar)
  !  nbod  ==>  number of massive bodies (int scalar)
  !  mass  ==>  mass of bodies (real array)
  !  iu    ==>  unit number to read from
  !
  ! Returns:
  !  io_read_mass ==> = 0 read ok
  !                  != 0 read failed is set to iostat variable
  !
  ! Remarks: Based on io_read_frame
  ! Authors:  Hal Levison
  ! Date:    11/2/99
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: iu

  ! Outputs
  integer(ik) :: nbod, read_mass
  real(rk) :: mass(nbod),time

  ! Internals
  integer(ik) :: i,ierr
  real(real4) :: mass4(NTPMAX), ttmp

  !-----------------
  ! Executable code
  !-----------------

  ierr = ixdrreal( iu, ttmp )
  !ierr = ixdrdouble( iu, time )
  read_mass = ierr
  if(ierr /= 0) return
  time = ttmp

  ierr = ixdrint( iu, nbod )
  read_mass = ierr
  if(ierr /= 0) return

  ierr = ixdrrmat( iu, nbod, mass4 )
  !ierr = ixdrdmat( iu, nbod, mass )
  read_mass = ierr
  if(ierr /= 0) return

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, mass, mass4)
  do i = 1, nbod

    mass(i) = mass4(i)

  end do
  !$OMP END PARALLEL DO

  return
  end function io_read_mass
!!
  function io_read_mass_r(time,nbod,mass,iu) result(read_mass_r)
  !-------------------------------------------------------------------------
  !			IO_READ_MASS_R.F90
  !-------------------------------------------------------------------------
  ! Read in the mass file.
  !
  ! Output:
  !  time  ==>  current time (real scalar)
  !  nbod  ==>  number of massive bodies (int scalar)
  !  mass  ==>  mass of bodies (real array)
  !  iu    ==>  unit number to read from
  !
  ! Returns:
  !  read_mass  ==>  = 0 read ok
  !                 != 0 read failed is set to iostat variable
  !
  ! Remarks: Based on io_read_frame
  ! Authors:  Hal Levison
  ! Date:    1/9/97
  ! Last revision: 11/2/99
  implicit none

  ! Inputs:
  integer(ik) :: iu

  ! Outputs
  integer(ik) :: nbod, read_mass_r
  real(rk) :: mass(nbod),time

  ! Internals
  integer(ik) :: i,ierr
  integer(integer2) :: nbod2
  real(real4) :: mass4(NTPMAX),ttmp

  !-----------------
  ! Executable code
  !-----------------

  read(iu, iostat = ierr) ttmp,nbod2
  read_mass_r = ierr
  if(ierr /= 0) return

  read(iu, iostat = ierr) (mass4(i), i = 1, nbod2)
  read_mass_r = ierr
  if(ierr /= 0) return

  do i = 1, nbod2

    mass(i) = mass4(i)

  end do
  nbod = nbod2
  time = ttmp

  return
  end function io_read_mass_r
!!
  function io_read_radius(time,nbod,radius,iu) result(read_radius)
  !-------------------------------------------------------------------------
  !				IO_READ_RADIUS.F90
  !-------------------------------------------------------------------------
  ! NEW VERSION OF THIS, USES FXDR
  ! Read in the radius file.
  !
  ! Output:
  !  time  ==>  current time (real scalar)
  !  nbod  ==>  number of radiusive bodies (int scalar)
  !  radius  ==>  radius of bodies (real array)
  !  iu    ==>  unit number to read from
  !
  ! Returns:
  !  io_read_radius ==> = 0 read ok
  !                  != 0 read failed is set to iostat variable
  !
  ! Remarks: Based on io_read_frame
  ! Authors:  Hal Levison
  ! Date:    11/2/99
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: iu

  ! Outputs
  integer(ik) :: nbod, read_radius
  real(rk) :: radius(nbod),time

  ! Internals
  integer(ik) :: i,ierr
  real(real4) :: radius4(NTPMAX), ttmp

  !-----------------
  ! Executable code
  !-----------------

  ierr = ixdrreal( iu, ttmp )
  !ierr = ixdrdouble( iu, time )
  read_radius = ierr
  if(ierr /= 0) return
  time = ttmp

  ierr = ixdrint( iu, nbod )
  read_radius = ierr
  if(ierr /= 0) return

  ierr = ixdrrmat( iu, nbod, radius4 )
  !ierr = ixdrdmat( iu, nbod, radius )
  read_radius = ierr
  if(ierr /= 0) return

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, radius, radius4)
  do i = 1, nbod

    radius(i) = radius4(i)

  end do
  !$OMP END PARALLEL DO

  return
  end function io_read_radius
!!
  function io_read_radius_r(time,nbod,radius,iu) result(read_radius_r)
  !-------------------------------------------------------------------------
  !			IO_READ_RADIUS_R.F90
  !-------------------------------------------------------------------------
  ! Read in the radius file.
  !
  ! Output:
  !  time  ==>  current time (real scalar)
  !  nbod  ==>  number of radiusive bodies (int scalar)
  !  radius  ==>  radius of bodies (real array)
  !  iu    ==>  unit number to read from
  !
  ! Returns:
  !  read_radius  ==>  = 0 read ok
  !                 != 0 read failed is set to iostat variable
  !
  ! Remarks: Based on io_read_frame
  ! Authors:  Hal Levison
  ! Date:    1/9/97
  ! Last revision: 11/2/99
  implicit none

  ! Inputs:
  integer(ik) :: iu

  ! Outputs
  integer(ik) :: nbod, read_radius_r
  real(rk) :: radius(nbod),time

  ! Internals
  integer(ik) :: i,ierr
  integer(integer2) :: nbod2
  real(real4) :: radius4(NTPMAX),ttmp

  !-----------------
  ! Executable code
  !-----------------

  read(iu, iostat = ierr) ttmp,nbod2
  read_radius_r = ierr
  if(ierr /= 0) return

  read(iu, iostat = ierr) (radius4(i), i = 1, nbod2)
  read_radius_r = ierr
  if(ierr /= 0) return

  do i = 1, nbod2

    radius(i) = radius4(i)

  end do
  nbod = nbod2
  time = ttmp

  return
  end function io_read_radius_r
!!
  subroutine io_write_frame(time,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh, &
             xht,yht,zht,vxht,vyht,vzht,istat,oname,iu,fopenstat)
  !----------------------------------------------------------------------------
  !				IO_WRITE_FRAME.F90
  !----------------------------------------------------------------------------
  ! NEW VERSION OF THIS, USES FXDR
  !
  ! Input:
  !  time           ==>  current time (real scalar)
  !  nbod           ==>  number of massive bodies (int scalar)
  !  ntp            ==>  number of massive bodies (int scalar)
  !  mass           ==>  mass of bodies (real array)
  !  xh,yh,zh       ==>  current position in helio coord
  !                      (real arrays)
  !  vxh,vyh,vzh    ==>  current velocity in helio coord
  !                      (real arrays)
  !  xht,yht,zht    ==>  current part position in helio coord
  !                      (real arrays)
  !  vxht,vyht,vzht ==>  current velocity in helio coord
  !                      (real arrays)
  !  istat          ==>  status of the test paricles
  !  oname          ==>  output file name (character string)
  !  iu             ==>  unit number to write to
  !  fopenstat      ==>  The status flag for the open
  !                      statements of the output files.
  !                      (character*80)
  !
  ! Output:
  !  iu             ==>  unit number to write to.
  !                      FXDR changes this
  !
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    11/2/99
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: nbod,ntp,iu
  integer(ik) :: istat(ntp)
  real(rk) :: mass(nbod),time
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: xht(ntp),yht(ntp),zht(ntp)
  real(rk) :: vxht(ntp),vyht(ntp),vzht(ntp)
  character(len = *) :: oname,fopenstat

  ! Internals
  integer(ik) :: i,id
  integer(ik) :: ialpha, ierr
  integer(ik), save :: i1st = 0 ! = 0 first time through; = 1 after
  real(rk) :: gm,a,e,inc,capom,omega,capm

  !-----------------
  ! Executable code
  !-----------------

  ! If first time through open file
  if(i1st == 0) then

    if((fopenstat(1:6) == 'append') .or. (fopenstat(1:6) == 'APPEND')) then

      call io_open_fxdr(oname, 'a', .true., iu, ierr)

    else

      if((fopenstat(1:3) == 'new') .or. (fopenstat(1:3) == 'NEW')) then

        call io_open_fxdr(oname, 'r', .true., iu, ierr)

        if(iu >= 0) then

          write(*,'(a)') ' SWIFT ERROR: in io_write_frame: '
          write(*,'(a)') '  Binary output file exists'
          call util_exit(FAILURE)

        end if

      end if

      call io_open_fxdr(oname, 'w', .true., iu, ierr)

    end if

    if(ierr < 0) then

      write(*,'(a)') ' SWIFT ERROR: in io_write_frame: '
      write(*,'(a)') '  Could not open binary output file'
      call util_exit(FAILURE)

    end if

    i1st = 1

  else

    call io_open_fxdr(oname, 'a', .true., iu, ierr)

    if(ierr < 0) then

      write(*,'(a)') ' SWIFT ERROR: in io_write_frame: '
      write(*,'(a)') '  Could not open binary output file with append'
      call util_exit(FAILURE)

    end if

  end if

  call io_write_hdr(iu,time,nbod,ntp,istat)

  ! Write out planets
  do i = 2, nbod

    gm = mass(1) + mass(i)
    id = -i
    call orbel_xv2el(xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i),gm, &
         ialpha,a,e,inc,capom,omega,capm)
    call io_write_line(iu,id,a,e,inc,capom,omega,capm)
    !call io_write_line(iu,id,xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i))

  end do

  ! Write out test particles
  gm = mass(1)
  do i = 1, ntp

    if(istat(i) == 0) then

      call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),vzht(i),gm, &
           ialpha,a,e,inc,capom,omega,capm)
      call io_write_line(iu,i,a,e,inc,capom,omega,capm)

    end if

  end do

  ierr = ixdrclose(iu)
  if(ierr < 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_frame: '
    write(*,'(a)') '  Could not close binary output file'
    call util_exit(FAILURE)

  end if

  return
  end subroutine io_write_frame
!!
  subroutine io_write_frame_r(time,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh, &
             xht,yht,zht,vxht,vyht,vzht,istat,oname,iu,fopenstat)
  !-------------------------------------------------------------------------
  !			IO_WRITE_FRAME_R.F90
  !-------------------------------------------------------------------------
  ! Write out a whole frame to an real*4 binary file.
  ! both massive and test particles
  !
  ! Input:
  !  time           ==>  current time (real scalar)
  !  nbod           ==>  number of massive bodies (int scalar)
  !  ntp            ==>  number of massive bodies (int scalar)
  !  mass           ==>  mass of bodies (real array)
  !  xh,yh,zh       ==>  current position in helio coord
  !                      (real arrays)
  !  vxh,vyh,vzh    ==>  current velocity in helio coord
  !                      (real arrays)
  !  xht,yht,zht    ==>  current part position in helio coord
  !                      (real arrays)
  !  vxht,vyht,vzht ==>  current velocity in helio coord
  !                      (real arrays)
  !  istat          ==>  status of the test paricles
  !                      2d integer array)
  !                      istat(i,1) = 0 ==> active:  = 1 not
  !                      istat(i,2) = -1 ==> Danby did not work
  !  oname          ==>  output file name (character string)
  !  iu             ==>  unit number to write to
  !  fopenstat      ==>  The status flag for the open
  !                      statements of the output files.
  !                      (character*80)
  !
  ! Remarks: Based on io_write_frame
  ! Authors:  Hal Levison
  ! Date:    2/22/94
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: nbod,ntp,iu
  integer(ik) :: istat(ntp+1,NSTAT)
  real(rk) :: mass(nbod),time
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: xht(ntp),yht(ntp),zht(ntp)
  real(rk) :: vxht(ntp),vyht(ntp),vzht(ntp)
  character(len = *) :: oname,fopenstat

  ! Internals
  integer(ik) :: i,id
  integer(ik) :: ialpha,ierr
  integer(ik), save :: i1st = 0 ! = 0 first time through; = 1 after
  real(rk) :: gm,a,e,inc,capom,omega,capm

  !-----------------
  ! Executable code
  !-----------------

  ! If first time through open file
  if(i1st == 0) then

    call io_open(iu,oname,fopenstat,'UNFORMATTED',ierr)
    if(ierr /= 0) then

      write(*,*) ' SWIFT ERROR: in io_write_frame_r: '
      write(*,*) '  Could not open binary output file'
      call util_exit(FAILURE)

    end if

    i1st = 1

  else

    call io_open(iu,oname,'append','UNFORMATTED',ierr)

  end if

  call io_write_hdr_r(iu,time,nbod,ntp,istat)

  ! Write out planets
  do i = 2, nbod

    gm = mass(1) + mass(i)
    id = -i
    call orbel_xv2el(xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i),gm, &
         ialpha,a,e,inc,capom,omega,capm)
    call io_write_line_r(iu,id,a,e,inc,capom,omega,capm)

  end do

  ! Write out test particles
  gm = mass(1)
  do i = 1, ntp

    if(istat(i,1) == 0) then

      call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),vzht(i),gm, &
           ialpha,a,e,inc,capom,omega,capm)
      call io_write_line_r(iu,i,a,e,inc,capom,omega,capm)

    end if

  end do

  close(unit = iu)

  return
  end subroutine io_write_frame_r
!!
  subroutine io_write_hdr(iu,time,nbod,ntp,istat)
  !-------------------------------------------------------------------------
  !			IO_WRITE_HDR.F90
  !-------------------------------------------------------------------------
  ! NEW VERSION OF THIS, USES FXDR
  !
  ! Input:
  !  iu     ==>  unit number to write to
  !  time   ==>  current time (real scalar)
  !  nbod   ==>  number of massive bodies (int scalar)
  !  ntp    ==>  number of massive bodies (int scalar)
  !  istat  ==>  status of the test paricles
  !
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    11/2/99
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: nbod,ntp,istat(ntp+1,NSTAT),iu
  real(rk) :: time

  ! Internals
  integer(ik) :: i,ierr
  integer(ik) :: nleft,nn(2)
  real(real4) :: ttmp
  !real(rk) :: ttmp

  !-----------------
  ! Executable code
  !-----------------

  ! Calculate number of remaining test particles
  nleft = 0
  do i = 1, ntp

    if(istat(i,1) == 0) nleft = nleft + 1

  end do

  ttmp = time

  ierr = ixdrreal( iu, ttmp )
  !ierr = ixdrdouble( iu, ttmp )
  if(ierr > 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_hdr: '
    write(*,'(a)') '  Could not write time'
    call util_exit(1)

  end if

  nn(1) = nbod
  nn(2) = nleft
  ierr = ixdrimat( iu, 2, nn )
  if(ierr > 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_hdr: '
    write(*,'(a)') '  Could not write nbod and nleft'
    call util_exit(1)

  end if

  return
  end subroutine io_write_hdr
!!
  subroutine io_write_hdr_r(iu,time,nbod,ntp,istat)
  !-------------------------------------------------------------------------
  !				IO_WRITE_HDR_R.F90
  !-------------------------------------------------------------------------
  ! Write out header part of the real*4 binary file
  !
  ! Input:
  !  iu     ==>  unit number to write to
  !  time   ==>  current time (real scalar)
  !  nbod   ==>  number of massive bodies (int scalar)
  !  ntp    ==>  number of massive bodies (int scalar)
  !  istat  ==>  status of the test paricles
  !
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    2/22/94
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: nbod,ntp,istat(ntp+1,NSTAT),iu
  real(rk) :: time

  ! Internals
  integer(ik) :: i
  integer(integer2) :: nleft,nbod2
  real(real4) :: ttmp

  !-----------------
  ! Executable code
  !-----------------

  ! Calculate number of remaining test particles
  nleft = 0
  do i = 1, ntp

    if(istat(i,1) == 0) nleft = nleft + 1

  end do

  nbod2 = nbod
  ttmp = time

  write(iu) ttmp,nbod2,nleft

  return
  end subroutine io_write_hdr_r
!!
  subroutine io_write_line(iu,id,a,e,inc,capom,omega,capm)
  !-------------------------------------------------------------------------
  !			IO_WRITE_LINE.F90
  !-------------------------------------------------------------------------
  ! NEW VERSION OF THIS, USES FXDR
  !
  ! Input:
  !  iu       ==> unit number to write to
  !  a        ==> semi-major axis or pericentric distance if a parabola
  !               (real scalar)
  !  e        ==> eccentricity (real scalar)
  !  inc      ==> inclination  (real scalar)
  !  capom    ==> longitude of ascending node (real scalar)
  !  omega    ==> argument of perihelion (real scalar)
  !  capm     ==> mean anomoly(real scalar)
  !
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    11/2/99
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: iu,id
  real(rk) :: a,e,inc,capom,omega,capm

  ! Internals
  integer(ik) :: ierr
  real(real4) :: orbel(6)
  !real(rk) :: orbel(6)

  !-----------------
  ! Executable code
  !-----------------

  ierr = ixdrint(iu, id)
  if(ierr > 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_line: '
    write(*,'(a)') '  Could not write id'
    call util_exit(1)

  end if

  orbel(1) = a
  orbel(2) = e
  orbel(3) = inc
  orbel(4) = capom
  orbel(5) = omega
  orbel(6) = capm

  ierr = ixdrrmat( iu, 6, orbel )
  !ierr = ixdrdmat( iu, 6, orbel )
  if(ierr.gt.0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_line: '
    write(*,'(a)') '  Could not write orbit elements'
    call util_exit(1)

  end if

  return
  end subroutine io_write_line
!!
  subroutine io_write_line_r(iu,id,a,e,inc,capom,omega,capm)
  !-------------------------------------------------------------------------
  !			IO_WRITE_LINE_R.F90
  !-------------------------------------------------------------------------
  ! Write out one line to real*4 binary file.
  !
  ! Input:
  !  iu       ==> unit number to write to
  !  a        ==> semi-major axis or pericentric distance if a parabola
  !               real scalar)
  !  e        ==> eccentricity (real scalar)
  !  inc      ==> inclination  (real scalar)
  !  capom    ==> longitude of ascending node (real scalar)
  !  omega    ==> argument of perihelion (real scalar)
  !  capm     ==> mean anomoly(real scalar)
  !
  ! Remarks:
  ! Authors:  Hal Levison
  ! Date:    2/22/94
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: iu,id
  real(rk) :: a,e,inc,capom,omega,capm

  ! Internals
  integer(integer2) :: id2
  real(real4) :: a4,e4,inc4,capom4,omega4,capm4

  !-----------------
  ! Executable code
  !-----------------

  id2 = id

  a4 = a
  e4 = e
  inc4 = inc
  capom4 = capom
  capm4 = capm
  omega4 = omega

  write(iu) id2,a4,e4,inc4,capom4,omega4,capm4

  return
  end subroutine io_write_line_r
!!
  subroutine io_write_mass(time,nbod,mass,oname,iu,fopenstat)
  !-------------------------------------------------------------------------
  !				IO_WRITE_MASS.F90
  !-------------------------------------------------------------------------
  ! NEW VERSION OF THIS, USES FXDR
  ! write out masses
  !
  ! Input:
  !  time       ==>  current time (real scalar)
  !  nbod       ==>  number of massive bodies (int scalar)
  !  mass       ==>  mass of bodies (real array)
  !  oname      ==>  output file name (character string)
  !  iu         ==>  unit number to write to
  !  fopenstat  ==>  The status flag for the open
  !                  statements of the output files.
  !                  (character*80)
  !
  ! Remarks: Based on io_write_frame
  ! Authors:  Hal Levison
  ! Date:    11/2/99
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: nbod,iu
  real(rk) :: mass(nbod),time
  character(len = *) :: oname,fopenstat

  ! Internals
  integer(ik) :: ierr,i
  integer(ik), save :: i1st = 0 ! i1st = 0 first time through; i1st = 1 after
  real(real4) :: mass4(nbod), ttmp

  !-----------------
  ! Executable code
  !-----------------

  ! If first time through open file
  if(i1st == 0) then

    if((fopenstat == 'append') .or. (fopenstat == 'APPEND')) then

      call io_open_fxdr('mass.' // oname, 'a', .true., iu, ierr)

    else

      if((fopenstat == 'new') .or. (fopenstat == 'NEW')) then

        call io_open_fxdr('mass.' // oname, 'w', .true., iu, ierr)

        if(iu >= 0) then

          write(*,'(a)') ' SWIFT ERROR: in io_write_mass:'
          write(*,'(a)') '  Binary output file exists'
          call util_exit(FAILURE)

        end if

      end if

      call io_open_fxdr('mass.' // oname, 'w', .true., iu, ierr)

    end if

    if(ierr < 0) then

      write(*,'(a)') ' SWIFT ERROR: in io_write_mass:'
      write(*,'(a)') '  Could not open binary output file'
      call util_exit(FAILURE)

    end if

    i1st = 1

  else

    call io_open_fxdr('mass.' // oname, 'a', .true., iu, ierr)

    if(ierr < 0) then

      write(*,'(a)') ' SWIFT ERROR: in io_write_mass:'
      write(*,'(a)') '  Could not open binary output file with append'
      call util_exit(FAILURE)

    end if

  end if

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, mass, mass4)
  do i = 1, nbod

    mass4(i) = mass(i)

  end do
  !$OMP END PARALLEL DO

  ttmp = time

  ierr = ixdrreal( iu, ttmp )
  !ierr = ixdrdouble( iu, time )
  if(ierr < 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_mass:'
    write(*,'(a)') '  Could not write time'
    call util_exit(FAILURE)

  end if

  ierr = ixdrint( iu, nbod )
  if(ierr < 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_mass:'
    write(*,'(a)') '  Could not write nbod'
    call util_exit(FAILURE)

  end if

  ierr = ixdrrmat( iu, nbod, mass4 )
  !ierr = ixdrdmat( iu, nbod, mass )
  if(ierr < 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_mass:'
    write(*,'(a)') '  Could not write masses'
    call util_exit(FAILURE)

  end if

  ierr = ixdrclose(iu)
  if(ierr < 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_frame: '
    write(*,'(a)') '  Could not close mass output file'
    call util_exit(FAILURE)

  end if

  return
  end subroutine io_write_mass
!!
  subroutine io_write_mass_r(time,nbod,mass,oname,iu,fopenstat)
  !-------------------------------------------------------------------------
  !			IO_WRITE_MASS_R.F90
  !-------------------------------------------------------------------------
  ! Write out masses
  !
  ! Input:
  !  time       ==>  current time (real scalar)
  !  nbod       ==>  number of massive bodies (int scalar)
  !  mass       ==>  mass of bodies (real array)
  !  oname      ==>  output file name (character string)
  !  iu         ==>  unit number to write to
  !  fopenstat  ==>  The status flag for the open
  !                  statements of the output files.
  !                  (character*80)
  !
  ! Remarks: Based on io_write_frame
  ! Authors:  Hal Levison
  ! Date:    1/9/97
  ! Last revision: 11/2/99
  implicit none

  ! Inputs:
  integer(ik) :: nbod, iu
  real(rk) :: mass(nbod), time
  character(len = *) :: oname, fopenstat

  ! Internals
  integer(ik) :: i, ierr
  integer(ik), save :: i1st = 0 ! i1st = 0 first time through; i1st = 1 after
  integer(integer2) :: nbod2
  real(real4) :: mass4(nbod), ttmp

  !-----------------!
  ! Executable code !
  !-----------------!

  ! If first time through open file
  if(i1st == 0) then

    call io_open(iu,'mass.' // oname,fopenstat,'UNFORMATTED',ierr)

    if(ierr /= 0) then

      write(*,'(a)') ' SWIFT ERROR: in io_write_mass_r:'
      write(*,'(a)') '  Could not open binary output file'
      call util_exit(FAILURE)

    end if

    i1st = 1

  else

    call io_open(iu,'mass.' // oname,'append','UNFORMATTED',ierr)

  end if

  do i = 1, nbod

    mass4(i) = mass(i)

  end do
  nbod2 = nbod
  ttmp = time

  write(iu) ttmp,nbod2
  write(iu) (mass4(i), i = 1, nbod2)

  close(unit = iu)

  return
  end subroutine io_write_mass_r
!!
  subroutine io_write_radius(time,nbod,radius,oname,iu,fopenstat)
  !-------------------------------------------------------------------------
  !				IO_WRITE_RADIUS.F90
  !-------------------------------------------------------------------------
  ! NEW VERSION OF THIS, USES FXDR
  ! Write out radiii
  !
  ! Input:
  !  time       ==>  current time (real scalar)
  !  nbod       ==>  number of massive bodies (int scalar)
  !  radius     ==>  radii of bodies (real array)
  !  oname      ==>  output file name (character string)
  !  iu         ==>  unit number to write to
  !  fopenstat  ==>  The status flag for the open
  !                  statements of the output files.
  !                  (character*80)
  !
  ! Remarks: Based on io_write_frame
  ! Authors:  Hal Levison
  ! Date:    11/2/99
  ! Last revision:
  implicit none

  ! Inputs:
  integer(ik) :: nbod,iu
  real(rk) :: radius(nbod),time
  character(len = *) :: oname,fopenstat

  ! Internals
  integer(ik) :: ierr,i
  integer(ik), save :: i1st = 0 ! i1st = 0 first time through; i1st = 1 after
  real(real4) :: radius4(NTPMAX), ttmp

  !-----------------
  ! Executable code
  !-----------------

  ! If first time through open file
  if(i1st == 0) then

    if((fopenstat == 'append') .or. (fopenstat == 'APPEND')) then

      call io_open_fxdr('radius.' // oname, 'a', .true., iu, ierr)

    else

      if((fopenstat == 'new') .or. (fopenstat == 'NEW')) then

        call io_open_fxdr('radius.' // oname, 'w', .true., iu, ierr)

        if(iu >= 0) then

          write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
          write(*,'(a)') '  Binary output file exists'
          call util_exit(FAILURE)

        end if

      end if

      call io_open_fxdr('radius.' // oname, 'w', .true., iu, ierr)

    end if

    if(ierr < 0) then

      write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
      write(*,'(a)') '  Could not open binary output file'
      call util_exit(FAILURE)

    end if

    i1st = 1

  else

    call io_open_fxdr('radius.' // oname, 'a', .true., iu, ierr)

    if(ierr < 0) then

      write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
      write(*,'(a)') '  Could not open binary output file with append'
      call util_exit(FAILURE)

    end if

  end if

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, radius, radius4)
  do i = 1, nbod

    radius4(i) = radius(i)

  end do
  !$OMP END PARALLEL DO

  ttmp = time

  ierr = ixdrreal( iu, ttmp )
  !ierr = ixdrdouble( iu, time )
  if(ierr < 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
    write(*,'(a)') '  Could not write time'
    call util_exit(FAILURE)

  end if

  ierr = ixdrint( iu, nbod )
  if(ierr < 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
    write(*,'(a)') '  Could not write nbod'
    call util_exit(FAILURE)

  end if

  ierr = ixdrrmat( iu, nbod, radius4 )
  !ierr = ixdrdmat( iu, nbod, radius )
  if(ierr < 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_radius:'
    write(*,'(a)') '  Could not write radii'
    call util_exit(FAILURE)

  end if

  ierr = ixdrclose(iu)
  if(ierr < 0) then

    write(*,'(a)') ' SWIFT ERROR: in io_write_frame: '
    write(*,'(a)') '  Could not close radius output file'
    call util_exit(FAILURE)

  end if

  return
  end subroutine io_write_radius
!!
  subroutine io_write_radius_r(time,nbod,radius,oname,iu,fopenstat)
  !-------------------------------------------------------------------------
  !			IO_WRITE_RADIUS_R.F90
  !-------------------------------------------------------------------------
  ! Write out radiuses
  !
  ! Input:
  !  time       ==>  current time (real scalar)
  !  nbod       ==>  number of massive bodies (int scalar)
  !  radius     ==>  radii of bodies (real array)
  !  oname      ==>  output file name (character string)
  !  iu         ==>  unit number to write to
  !  fopenstat  ==>  The status flag for the open
  !                  statements of the output files.
  !                  (character*80)
  !
  ! Remarks: Based on io_write_frame
  ! Authors:  Hal Levison
  ! Date:    1/9/97
  ! Last revision: 11/2/99
  implicit none

  ! Inputs:
  integer(ik) :: nbod,iu
  real(rk) :: radius(nbod),time
  character(len = *) :: oname,fopenstat

  ! Internals
  integer(ik) :: i,ierr
  integer(ik), save :: i1st = 0 ! i1st = 0 first time through; i1st = 1 after
  integer(integer2) :: nbod2
  real(real4) :: radius4(NTPMAX), ttmp

  !-----------------
  ! Executable code
  !-----------------

  ! If first time through open file
  if(i1st == 0) then

    call io_open(iu,'radius.' // oname,fopenstat,'UNFORMATTED',ierr)

    if(ierr /= 0) then

      write(*,'(a)') ' SWIFT ERROR: in io_write_radius_r:'
      write(*,'(a)') '  Could not open binary output file'
      call util_exit(FAILURE)

    end if

    i1st = 1

  else

    call io_open(iu,'radius.' // oname,'append','UNFORMATTED',ierr)

  end if

  do i = 1, nbod

    radius4(i) = radius(i)

  end do
  nbod2 = nbod
  ttmp = time

  write(iu) ttmp,nbod2
  write(iu) (radius4(i), i = 1, nbod2)

  close(unit = iu)

  return
  end subroutine io_write_radius_r
!!
!!
end module io