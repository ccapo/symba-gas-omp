subroutine symba5_step_pl(i1st, time, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, dt, &
           lclose, rpl, isenc, mergelst, mergecnt, iecnt, eoff, rhill, mtiny, ibound)
!----------------------------------------------------------------------------
!				SYMBA5_STEP_PL.F90
!----------------------------------------------------------------------------
!
!             Input:
!                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
!                 time          ==>  current time (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 nbodm         ==>  location of the last massie body
!                                    (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 xh,yh,zh      ==>  initial position in helio coord
!                                    (real arrays)
!                 vxh,vyh,vzh   ==>  initial velocity in helio coord
!                                    (real arrays)
!                 dt            ==>  time step
!                 lclose        ==> .true. --> check for close encounters
!                                      (logical*2 scalar)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!                 eoff          ==>  Energy offset (real scalar)
!                 rhill         ==>  size of planet's hills sphere
!                                    (real array)
!                 mtiny         ==>  Small mass  (real array)
!             Output:
!                 xh,yh,zh      ==>  final position in helio coord
!                                       (real arrays)
!                 vxh,vyh,vzh   ==>  final velocity in helio coord
!                                       (real arrays)
!                 rpl           ==>  Recalculated physical size of a planet.
!                                    if merger happened (real array)
!                 nbod          ==>  Recalculated number of massive bodies
!                                    if merger happened (int scalar)
!                 mass          ==>  Recalculated mass of bodies
!                                    if merger happened (real array)
!                 isenc         ==>  0 --> No encounter during last dt
!                                    1 --> There was encounters
!                                     (integer scalar)
!                 mergelst      ==>  list of mergers (int array)
!                 mergecnt      ==>  count of mergers (int array)
!                 iecnt         ==>  Number of encounters (int*2 array)
!                 eoff          ==>  Energy offset (real scalar)
!                 rhill         ==>  size of planet's hills sphere
!                                    (real array)
!
! Remarks: Based on symba2_step_pl.f
! Authors:  Hal Levison
! Date:    11/27/97
! Last revision:
use module_swift
use module_symba5
use module_interfaces, except_this_one => symba5_step_pl
!$ use omp_lib
implicit none

! Inputs Only: 
logical(lk) :: lclose
integer(ik) :: nbod, i1st, nbodm
real(rk) :: mass(nbod), dt, time, j2rp2, j4rp4, mtiny

! Inputs and Outputs:
integer(ik) :: ibound(nbod)
real(rk) :: xh(nbod), yh(nbod), zh(nbod)
real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)
real(rk) :: rpl(nbod), rhill(nbod), eoff

! Outputs only
integer(ik) :: isenc
integer(ik) :: iecnt(nbod), ielev(nbod)
integer(ik) :: mergelst(NENMAX,2), mergecnt

! Internals
logical(lk) :: lvdotr            ! Not used in the routine
integer(ik) :: i, j, k, ieflg, irec
integer(ik) :: ielst(NENMAX,2), ielc
integer(ik) :: id, ielc_shared(nthreads), istart(nthreads), ielst_shared(nbod)
real(rk) :: rhillij(2), xhij(2), yhij(2), zhij(2), vxhij(2), vyhij(2), vzhij(2)

!----
! Executable code
!-----------------

! Initialize relevant variables
isenc = 0
ielc = 0
irec = 0

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, iecnt, ielev)
do i = 1, nbod

  iecnt(i) = 0
  ielev(i) = -1

end do
!$OMP END PARALLEL DO

! Check for encounters
if(lclose) then

  do i = 2, nbodm

    ! Extract the information for body i, and make a copy available for each thread
    rhillij(1) = rhill(i)
    xhij(1) = xh(i); yhij(1) = yh(i); zhij(1) = zh(i)
    vxhij(1) = vxh(i); vyhij(1) = vyh(i); vzhij(1) = vzh(i)

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(j, id, ieflg) &
    !$OMP FIRSTPRIVATE(i, nthreads, irec, xhij, yhij, zhij, vxhij, vyhij, vzhij, rhillij, dt) &
    !$OMP SHARED(nbod, rhill, xh, yh, zh, vxh, vyh, vzh, istart, ielc_shared, ielst_shared, iecnt, ielev)

    id = 1                              ! Set thread identifier for *serial* case
    !$ id = omp_get_thread_num() + 1    ! Set thread identifier for *parallel* case
    ielc_shared(id) = 0                 ! Initialize encounter counter for each thread
    istart(id) = (id - 1)*nbod/nthreads ! Set the starting position in ielst_shared for each thread

    !$OMP DO SCHEDULE(STATIC)
    do j = i + 1, nbod

      rhillij(2) = rhill(j)
      xhij(2) = xh(j); yhij(2) = yh(j); zhij(2) = zh(j)
      vxhij(2) = vxh(j); vyhij(2) = vyh(j); vzhij(2) = vzh(j)
      call symba5_chk2(rhillij, xhij, yhij, zhij, vxhij, vyhij, vzhij, dt, irec, ieflg)

      if(ieflg /= 0) then

        iecnt(j) = iecnt(j) + 1
        ielev(j) = 0
	ielc_shared(id) = ielc_shared(id) + 1 ! Update counter *before* storing index j since istart starts at 0
        ielst_shared(istart(id) + ielc_shared(id)) = j

      end if

    end do
    !$OMP END DO NOWAIT

    !$OMP END PARALLEL

    ! If there any encounters, append to the global encounter list for body i
    if(sum(ielc_shared) > 0) then

      ! Set global encounter flag
      isenc = 1

      ! Set the recursion level to 0 for body i
      ielev(i) = 0

      ! Search through each thread's local encounter list
      do k = 1, nthreads

        ! If an encounter was found for thread k, append ielc_shared(k) entries to ielst
        if(ielc_shared(k) > 0) then

          do j = 1, ielc_shared(k)

            iecnt(i) = iecnt(i) + 1
            ielst(ielc + j, 1) = i
            ielst(ielc + j, 2) = ielst_shared(istart(k) + j)

          end do

          ! Update the global encounter counter
          ielc = ielc + ielc_shared(k)

          if(ielc > NENMAX) then

            write(*,'(a)') 'ERROR: Encounter matrix is filled.'
            write(*,'(a)') 'STOPPING'
            call util_exit(FAILURE)

          end if

        end if

      end do

    end if

  end do

end if

! do a step
if(isenc == 0) then

  call symba5_step_helio(i1st, nbod, nbodm, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, dt)
  mergecnt = 0

else

  call symba5_step_interp(time, ielev, nbod, nbodm, mass, rhill, j2rp2, j4rp4, rpl, xh, yh, zh, vxh, vyh, vzh, &
       dt, mergelst, mergecnt, eoff, ielc, ielst, mtiny, ibound)
  i1st = 0

end if

! Print number of encounters and mergers found in this time step
!write(100,*) time, ielc, mergecnt
!do i = 1, ielc
!  write(100,'(2i9)') ielst(i,:)
!end do

return
end subroutine symba5_step_pl
