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
use module_swift
use module_io
use module_interfaces, except_this_one => io_dump_pl_symba
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
