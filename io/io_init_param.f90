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
use module_swift
use module_io
use module_interfaces, except_this_one => io_init_param
!$ use omp_lib
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
nthreads = 1                        ! In the *serial* case
!$ nthreads = omp_get_max_threads() ! In the *parallel* case

return
end subroutine io_init_param
