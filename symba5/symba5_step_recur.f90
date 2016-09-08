recursive subroutine symba5_step_recur(t, nbod, nbodm, mass, ireci, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, &
                     rpl, mergelst, mergecnt, dt0, eoff, lvdotr, ibound, ielc, ielst)
!------------------------------------------------------------------------------
!				SYMBA5_STEP_RECUR.F90
!------------------------------------------------------------------------------
!
!             Input:
!                 t             ==>  time (real Scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 nbodm         ==>  Location of last massive body(int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 ireci         ==>  Input recursion level  (integer scalar)
!                 ilevl         ==>  largest recursion level used
!                                    (integer array)
!                 ielev         ==>  The level that this particle should go
!                                             (int*2 array)
!                 j2rp2, j4rp4   ==>  J2*radii_pl^2 and J4*radii_pl^4
!                                     (real scalars)
!                 rhill         ==>  Hill sphere of planet (real Scalar)
!                 xh, yh, zh      ==>  initial position in helio coord
!                                    (real arrays)
!                 vxb, vyb, vzb   ==>  initial velocity in bari coord
!                                    (real arrays)
!                dt0            ==>  Global timestep  (real scalar)
!                rpl            ==>  physical size of a planet.
!                                    (real array)
!                eoff           ==>  Energy offset (real scalar)
!                ielc           ==>  number of encounters (integer*2 scalar)
!                ielst          ==>  list of ecnounters (2D integer*2 array)
!             Output:
!                 xh, yh, zh      ==>  final position in helio coord
!                                       (real arrays)
!                 vxb, vyb, vzb   ==>  final velocity in bari coord
!                                       (real arrays)
!             mergelst          ==>  list of mergers (int array)
!             mergecnt          ==>  count of mergers (int array)
!                 rpl           ==>  Recalculated physical size of a planet.
!                                    if merger happened (real array)
!                 mass          ==>  Recalculated mass of bodies
!                                    if merger happened (real array)
!                eoff           ==>  Energy offset (real scalar)
!                lvdotr         ==> vdotr relative flag
!                                   = .true. if i, j are receding
!                                   = .false is approaching
!                                     (2D logical*1 array)
!                                   Used internally,  but only need 1 copy.
!
! Remarks: If a merger occurs,  does not change nbod and puts the mass
!          of one of the particles to zero.
! Authors:  Hal Levison
! Date:    3/20/97
! Last revision: 5/13/99
use module_swift
use module_symba5
use module_interfaces, except_this_one =>  symba5_step_recur
implicit none

! Inputs Only:
integer(ik) :: nbod, ireci, nbodm
integer(ik) :: ielev(nbod)
integer(ik) :: ielst(NENMAX,2), ielc
real(rk) :: mass(nbod), dt0, rhill(nbod), t

! Inputs and Outputs:
logical(lk) :: lvdotr(NTPMAX)
integer(ik) :: ibound(nbod)
integer(ik) :: mergelst(NENMAX,2), mergecnt
real(rk) :: xh(nbod), yh(nbod), zh(nbod), eoff
real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod), rpl(nbod)

! Internals:
integer(ik) :: i, j, ie
integer(ik) :: icflg, it, irecp, ieflg
real(rk) :: dtl, dth, sgn

!----
! Executable code

dtl = dt0/real(NTENC**ireci, rk)
dth = 0.5_rk*dtl

if(dtl/dt0 <= TINY) then

  write(*,'(a)') 'Warning in SYMBA_STEP_RECUR:'
  write(*,'(a)') ' Local timestep too small '
  write(*,'(a)') ' Roundoff will be important!!!! '
  call util_exit(FAILURE)

end if

irecp = ireci + 1

if(ireci == 0) then

  ! Do we need to go deeper?
  icflg = 0
  do ie = 1, ielc

    i = ielst(ie,1)
    j = ielst(ie,2)

    if((ielev(i) >= ireci) .and. (ielev(j) >= ireci)) then

      call symba5_chk(rhill, nbod, i, j, mass, xh, yh, zh, vxb, vyb, vzb, dtl, irecp, ieflg, lvdotr(ie))

      if(ieflg /= 0) then

        icflg = 1
        ielev(i) = irecp
        ielev(j) = irecp

      end if

    end if

  end do

  sgn = 1.0_rk
  call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

  call symba5_enc_drift(nbod, ielc, ielst, ielev, ireci, mass, xh, yh, zh, vxb, vyb, vzb, dtl)

  if(icflg /= 0) call symba5_step_recur(t, nbod, nbodm, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, &
                      rpl, mergelst, mergecnt, dt0, eoff, lvdotr, ibound, ielc, ielst)

  sgn = 1.0_rk
  call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

  ! look for mergers
  do ie = 1, ielc

    i = ielst(ie,1)
    j = ielst(ie,2)

    if((ielev(i) >= ireci) .and. (ielev(j) >= ireci)) then

      call symba5_merge(t, dtl, nbod, nbodm, i, j, mass, xh, yh, zh, vxb, vyb, vzb, &
           ireci, lvdotr(ie), ibound, rpl, mergelst, mergecnt, rhill, eoff, ielc, ielst)

    end if

    if(ielev(i) == irecp) ielev(i) = ireci
    if(ielev(j) == irecp) ielev(j) = ireci

  end do

else

  do it = 1, NTENC

    ! Do we need to go deeper?
    icflg = 0
    do ie = 1, ielc

      i = ielst(ie,1)
      j = ielst(ie,2)

      if((ielev(i) >= ireci) .and. (ielev(j) >= ireci)) then

        call symba5_chk(rhill, nbod, i, j, mass, xh, yh, zh, vxb, vyb, vzb, dtl, irecp, ieflg, lvdotr(ie))

        if(ieflg /= 0) then

          icflg = 1
          ielev(i) = irecp
          ielev(j) = irecp

        end if

      end if

    end do

    sgn = 1.0_rk
    call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

    sgn = -1.0_rk
    call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

    call symba5_enc_drift(nbod, ielc, ielst, ielev, ireci, mass, xh, yh, zh, vxb, vyb, vzb, dtl)

    if(icflg /= 0) call symba5_step_recur(t, nbod, nbodm, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, &
                        rpl, mergelst, mergecnt, dt0, eoff, lvdotr, ibound, ielc, ielst)

    sgn = 1.0_rk
    call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

    sgn = -1.0_rk
    call symba5_kick(nbod, mass, irecp, ielev, rhill, xh, yh, zh, vxb, vyb, vzb, dth, sgn, ielc, ielst)

    ! look for mergers
    do ie = 1, ielc

      i = ielst(ie,1)
      j = ielst(ie,2)

      if((ielev(i) >= ireci) .and. (ielev(j) >= ireci)) then

        call symba5_merge(t, dtl, nbod, nbodm, i, j, mass, xh, yh, zh, vxb, vyb, vzb, &
             ireci, lvdotr(ie), ibound, rpl, mergelst, mergecnt, rhill, eoff, ielc, ielst)

      end if

      if(ielev(i) == irecp) ielev(i) = ireci
      if(ielev(j) == irecp) ielev(j) = ireci

    end do

  end do

end if

return
end subroutine symba5_step_recur
