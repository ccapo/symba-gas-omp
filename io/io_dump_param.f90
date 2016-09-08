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
use module_swift
use module_io
use module_interfaces, except_this_one => io_dump_param
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
