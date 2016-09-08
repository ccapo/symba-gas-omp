function orbel_flon(e, capn) result(ea)
!**********************************************************************
!                    ORBEL_FLON.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                        capn ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!                  orbel_flon ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Uses power series for N in terms of F and Newton, s method
!     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 26, 1992.
!     REVISIONS:
use module_swift
implicit none

! Inputs Only:
real(rk) :: e, capn

! Output Only:
real(rk) :: ea

! Internals:
integer(ik), parameter :: IMAX = 10
real(rk), parameter :: a11 = 156.0_rk, a9 = 17160.0_rk, a7 = 1235520.0_rk, a5 = 51891840.0_rk, a3 = 1037836800.0_rk
real(rk), parameter :: b11 = 11.0_rk*a11, b9 = 9.0_rk*a9, b7 = 7.0_rk*a7, b5 = 5.0_rk*a5, b3 = 3.0_rk*a3
integer(ik) :: i, iflag
real(rk) :: a, b, sq, biga, bigb
real(rk) :: x, x2
real(rk) :: f, fp, dx
real(rk) :: diff
real(rk) :: a0, a1, b1

!----
! Executable code


! Function to solve "Kepler's eqn" for F (here called x) for given e and CAPN. Only good for smallish CAPN

iflag = 0
if(capn < 0.0_rk) then

  iflag = 1
  capn = -capn

end if

a1 = 6227020800.0_rk*(1.0_rk - 1.0_rk/e)
a0 = -6227020800.0_rk*capn/e
b1 = a1

! Set iflag nonzero if capn < 0.,  in which case solve for -capn
! and change the sign of the final answer for F.
! Begin with a reasonable guess based on solving the cubic for small F
a = 6.0_rk*(e - 1.0_rk)/e
b = -6.0_rk*capn/e
sq = sqrt(0.25_rk*b*b + a*a*a/27.0_rk)
biga = (sq - 0.5_rk*b)**(1.0_rk/3.0_rk)
bigb = -(sq + 0.5_rk*b)**(1.0_rk/3.0_rk)
x = biga + bigb
!write(*,'(a,1pe14.6)') 'cubic = ', x**3 + a*x + b
ea = x

! If capn is tiny (or zero) no need to go further than cubic even for e = 1.0
if(capn >= TINY) then

  do i = 1, IMAX

    x2 = x*x
    f = a0 + x*(a1 + x2*(a3 + x2*(a5 + x2*(a7 + x2*(a9 + x2*(a11 + x2))))))
    fp = b1 + x2*(b3 + x2*(b5 + x2*(b7 + x2*(b9 + x2*(b11 + 13.0_rk*x2)))))
    dx = -f/fp
    !write(*,'(a)')                    'i, dx, x, f:'
    !write(*,'(1x,i3,3(2x,1p1e22.15))') i, dx, x, f
    ea = x + dx

    ! If we have converged here there's no point in going on
    if(abs(dx) <= TINY) exit
    x = ea

  end do

end if

! Check if capn was originally negative
if(iflag == 1) then

  ea = -ea
  capn = -capn

end if

! Abnormal return here - we've gone thru the loop IMAX times without convergence
if(i >= imax) then

  write(*,'(a)') 'FLON: RETURNING WITHOUT COMPLETE CONVERGENCE'
  diff = e*sinh(ea) - ea - capn
  write(*,'(a)') 'N, F, ecc*sinh(F) - F - N: '
  write(*,'(3(1pe14.6))') capn, ea, diff

end if

return
end function orbel_flon