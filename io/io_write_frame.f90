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
use module_swift
use module_io
use module_fxdr
use module_interfaces, except_this_one => io_write_frame
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
