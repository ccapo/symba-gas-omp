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
use module_swift
use module_io
use module_interfaces, only: io_open
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