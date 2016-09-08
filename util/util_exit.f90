subroutine util_exit(iflag)
!-------------------------------------------------------------------------
!    UTIL_EXIT.F90
!-------------------------------------------------------------------------
! Exits program
!
! Input: iflag ==> Status of exit
!                   = 0 if normal exit
!                   = 1 if exit because error
!
! Authors: Hal Levison 
! Date: 08/06/93
!-------------------------------------------------------------------------
use module_swift
implicit none

! Inputs: 
integer(ik) :: iflag

!-----------------!
! Executable Code !
!-----------------!

if(iflag == 0) then

  write(*,'(/,a,f3.1,a)') 'Normal termination of SWIFT (Version: ', VER_NUM ,')'
  write(*,'(a)')          '------------------------------------------'

else

  write(*,'(/,a,f3.1,a)') 'Terminating SWIFT (Version: ', VER_NUM ,') due to ERROR!!!'
  write(*,'(a)')          '-----------------------------------------------'

end if

stop
end subroutine util_exit
