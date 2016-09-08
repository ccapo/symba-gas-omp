subroutine util_version
!-------------------------------------------------------------------------
!    UTIL_VERSION.F90
!-------------------------------------------------------------------------
! Prints version of SWIFT and contact information
!
! Authors: Hal Levison
! Date: 02/21/94
! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
!-------------------------------------------------------------------------
use module_swift
implicit none

!-----------------!
! Executable code !
!-----------------!

write(*,'(a)')          "!---------------------------------------------------------!"
write(*,'(a)')          "!                                                         !"
write(*,'(a,f3.1,a)')   "! SWIFT (Version: ", VER_NUM, ")                                    !"
write(*,'(a)')          "!                                                         !"
write(*,'(a)')          "!---------------------------------------------------------!"
write(*,'(a)')          "!                                                         !"
write(*,'(a)')          "! Authors:                                                !"
write(*,'(a)')          "!  Martin Duncan: Queen's University                      !"
write(*,'(a)')          "!  Hal Levison: Southwest Research Institute              !"
write(*,'(a)')          "!                                                         !"
write(*,'(a)')          "! Please address any comments or questions to:            !"
write(*,'(a)')          "!  Hal Levison                                            !"
write(*,'(a)')          "!  Geophysical, Astrophysical & Planetary Sciences        !"
write(*,'(a)')          "!  Southwest Research Institute                           !"
write(*,'(a)')          "!  1050 Walnut St.                                        !"
write(*,'(a)')          "!  Suite 429                                              !"
write(*,'(a)')          "!  Boulder, Co 80302                                      !"
write(*,'(a)')          "!  (303) 546-0290                                         !"
write(*,'(a)')          "!  Fax: (303) 546-9687                                    !"
write(*,'(a)')          "!  (D)  swri::levison                                     !"
write(*,'(a)')          "!  (I)  hal@gort.space.swri.edu                           !"
write(*,'(a)')          "!                                                         !"
write(*,'(a,/)')        "!---------------------------------------------------------!"

return
end subroutine util_version
