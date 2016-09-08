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
use module_swift
use module_io
use module_interfaces, except_this_one => io_energy_write
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
