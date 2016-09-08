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
use module_swift
use module_io
use module_interfaces, only: io_open, util_exit, util_hills
implicit none

! Input
logical(lk) :: lclose
integer(ik) :: iflgchk
character(len = *) :: infile

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
