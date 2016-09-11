module swift
! Module for SWIFT
!
! Author: Hal Levison
! Date: 2/2/93
! Last revision: 3/7/93
!              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
!                             - Uses new data types to simplify code, and
!                               make it easier to add new particle quantities
!                               or global parameters/flags
implicit none

! Type definitions

! Defines precision for integers and reals
integer, parameter :: integer2 = selected_int_kind(4)
integer, parameter :: integer4 = selected_int_kind(9)
integer, parameter :: real4 = selected_real_kind(6)
integer, parameter :: real8 = selected_real_kind(15)

! The user's choice between single and double precision for reals, and range for integers
integer, parameter :: ik = integer4
integer, parameter :: rk = real8
integer, parameter :: lk = kind(.true.)

! Version of SWIFT
real(rk), parameter :: VER_NUM = 2.5_rk

! Return flags
integer(ik), parameter :: SUCCESS = 0
integer(ik), parameter :: FAILURE = 1

! Maximum array size
integer(ik), parameter :: NDIM = 3        ! Number of spatial dimensions (update if string theory proves true)
integer(ik), parameter :: NPLMAX = 1023   ! max. number of planets, including the sun (2**7 - 1)
integer(ik), parameter :: NTPMAX = 524287 ! max. number of test particles (2**19 - 1)

! Size of the test particle integer status flag
integer(ik), parameter :: NSTATP = 3
integer(ik), parameter :: NSTAT = NSTATP + NPLMAX - 1 ! include one @ planet

! Size of the test particle integer status flag
integer(ik), parameter :: NSTATR = NSTAT              ! io_init_tp assumes nstat == nstatr

! Convergence criteria for danby
real(rk), parameter :: DANBYAC = 1.0e-14_rk, DANBYB = 1.0e-13_rk

! Loop limits in the Laguerre attempts
integer(ik), parameter :: NLAG1 = 50, NLAG2 = 400

! A small number
real(rk), parameter :: TINY_NUMBER = 4.0e-15_rk

! Trignometric stuff
real(rk), parameter :: PI = 3.14159265358979324_rk
real(rk), parameter :: TWOPI = 2.0_rk*PI
real(rk), parameter :: PIBY2 = 0.5_rk*PI
real(rk), parameter :: PI3BY2 = 1.5_rk*PI
real(rk), parameter :: DEGRAD = 180.0_rk/PI

! Simulation Parameters
real(rk), parameter :: MSUN = TWOPI**2
real(rk), parameter :: MEARTH = 3.0e-6_rk*MSUN

! Symbolic names for binary output file contents
integer(ik), parameter :: EL = 1
integer(ik), parameter :: XV = 2

! OpenMP Parameters
integer(ik), save :: nthreads = 1
integer(ik), parameter :: NTHRESHOLD = 1000

end module swift