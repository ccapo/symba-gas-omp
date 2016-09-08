module module_interfaces
! Definition of interfaces of subroutines and functions used in swift2 package
implicit none
!!
!!
interface
  subroutine anal_energy(nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,ke,pot,energy,l)
  use module_swift
  implicit none
  integer(ik) :: nbod,nbodm
  real(rk) :: mass(nbod),j2rp2,j4rp4
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: energy,l(NDIM),ke,pot
  end subroutine anal_energy
end interface
!!
interface
  subroutine anal_energy_write(t,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,iu,fopenstat,eoff)
  use module_swift
  implicit none
  integer(ik) :: nbod,nbodm,iu
  real(rk) :: mass(nbod),t,j2rp2,j4rp4,eoff
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  character(len = 24) :: fopenstat
  end subroutine anal_energy_write
end interface
!!
interface
  subroutine anal_energy_discard5(iflg,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,ke,pot,energy,eltot) 
  use module_swift
  implicit none
  integer(ik) :: iflg, nbod, nbodm
  real(rk) :: mass(nbod), j2rp2, j4rp4
  real(rk) :: xh(nbod), yh(nbod), zh(nbod)
  real(rk) :: vxh(nbod), vyh(nbod), vzh(nbod)
  real(rk) :: energy, eltot(NDIM), ke, pot
  end subroutine anal_energy_discard5
end interface
!!
!!
interface
  subroutine coord_b2h(nbod,mass,xb,yb,zb,vxb,vyb,vzb,xh,yh,zh,vxh,vyh,vzh)
  use module_swift
  implicit none
  integer(ik) :: nbod
  real(rk) :: mass(nbod)
  real(rk) :: xb(nbod),yb(nbod),zb(nbod)
  real(rk) :: vxb(nbod),vyb(nbod),vzb(nbod)
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  end subroutine coord_b2h
end interface
!!
interface
  subroutine coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,xb,yb,zb,vxb,vyb,vzb,msys)
  use module_swift
  implicit none
  integer(ik) :: nbod
  real(rk) :: mass(nbod)
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: xb(nbod),yb(nbod),zb(nbod),msys
  real(rk) :: vxb(nbod),vyb(nbod),vzb(nbod)
  end subroutine coord_h2b
end interface
!!
interface
  subroutine coord_vb2h(nbod,mass,vxb,vyb,vzb,vxh,vyh,vzh)
  use module_swift
  implicit none
  integer(ik) :: nbod
  real(rk) :: mass(nbod)
  real(rk) :: vxb(nbod),vyb(nbod),vzb(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  end subroutine coord_vb2h
end interface
!!
interface
  subroutine coord_vh2b(nbod,mass,vxh,vyh,vzh,vxb,vyb,vzb,msys)
  use module_swift
  implicit none
  integer(ik) :: nbod
  real(rk) :: mass(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: vxb(nbod),vyb(nbod),vzb(nbod),msys
  end subroutine coord_vh2b
end interface
!!
!!
interface
  subroutine discard_mass_merge5(time,nbod,nbodm,ip1,ip2,mass,xh,yh,zh, &
             vxh,vyh,vzh,rpl,eoff,ielc,ielst)
  use module_swift
  use module_symba5
  implicit none
  integer(ik) :: ip1,ip2,nbod,nbodm
  integer(ik) :: ielst(NENMAX,2),ielc
  real(rk) :: mass(nbod),xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod),rpl(nbod)
  real(rk) :: eoff, time
  end subroutine discard_mass_merge5
end interface
!!
interface
  subroutine discard_mass_peri(time,nbod,iecnt,mass,xh,yh,zh,vxh,vyh,vzh,qmin,iwhy,lwhy,isperi)
  use module_swift
  implicit none
  logical(lk) :: lwhy
  integer(ik) :: nbod
  integer(ik) :: iecnt(nbod)
  real(rk) :: mass(nbod),time,qmin
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  integer(ik) :: iwhy(nbod)
  integer(ik) :: isperi(nbod)
  end subroutine discard_mass_peri
end interface
!!
interface
  subroutine discard_mass_reorder5(ip,nbod,mass,xh,yh,zh,vxh,vyh,vzh,rpl,rhill,isperih,ibound)
  use module_swift
  implicit none
  integer(ik) :: ip, nbod
  integer(ik) :: isperih(nbod), ibound(nbod)
  real(rk) :: mass(nbod),xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod),rpl(nbod)
  real(rk) :: rhill(nbod)
  end subroutine discard_mass_reorder5
end interface
!!
interface
  subroutine discard_massive5(time,dt,nbod,nbodm,mass,xh,yh,zh,vxh,vyh,vzh, &
             rmin,rmax,rmaxu,qmin,lclose,rpl,rhill,isenc,mergelst,    &
             mergecnt,iecnt,eoff,i1st,ibound)
  use module_swift
  use module_symba5
  implicit none
  logical(lk) :: lclose
  integer(ik) :: isenc
  integer(ik) :: mergelst(NENMAX,2),mergecnt
  integer(ik) :: iecnt(nbod)
  real(rk) :: time,dt
  real(rk) :: rmin,rmax,rmaxu,qmin
  integer(ik) :: nbod,nbodm,i1st,ibound(nbod)
  real(rk) :: mass(nbod),xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: eoff,rpl(nbod),rhill(nbod)
  end subroutine discard_massive5
end interface
!!
!!
interface
  subroutine helio_drift(nbod,mass,xh,yh,zh,vxb,vyb,vzb,dt)
  use module_swift
  implicit none
  integer(ik) :: nbod
  real(rk) :: mass(nbod),dt
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxb(nbod),vyb(nbod),vzb(nbod)
  end subroutine helio_drift
end interface
!!
interface
  subroutine helio_lindrift(nbod,mass,vxb,vyb,vzb,dt,xh,yh,zh,ptx,pty,ptz)
  use module_swift
  implicit none
  integer(ik) :: nbod
  real(rk) :: mass(nbod),dt
  real(rk) :: vxb(nbod),vyb(nbod),vzb(nbod)
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: ptx,pty,ptz
  end subroutine helio_lindrift
end interface
!!
!!
interface
  subroutine drift_dan(mu,x0,y0,z0,vx0,vy0,vz0,dt0,iflg)
  use module_swift
  implicit none
  real(rk) :: mu,dt0
  real(rk) :: x0,y0,z0
  real(rk) :: vx0,vy0,vz0
  integer(ik) :: iflg
  end subroutine drift_dan
end interface
!!
interface
  subroutine drift_kepmd(dm,es,ec,x,s,c)
  use module_swift
  implicit none
  real(rk) :: dm,es,ec
  real(rk) :: x,s,c
  end subroutine drift_kepmd
end interface
!!
interface
  subroutine drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
  use module_swift
  implicit none
  real(rk) :: dt,r0,mu,alpha,u
  integer(ik) :: iflg
  real(rk) :: fp,c1,c2,c3
  end subroutine drift_kepu
end interface
!!
interface
  subroutine drift_kepu_fchk(dt,r0,mu,alpha,u,s,f)
  use module_swift
  implicit none
  real(rk) :: dt,r0,mu,alpha,u,s
  real(rk) :: f
  end subroutine drift_kepu_fchk
end interface
!!
interface
  subroutine drift_kepu_guess(dt,r0,mu,alpha,u,s)
  use module_swift
  implicit none
  real(rk) :: dt,r0,mu,alpha,u
  real(rk) :: s
  end subroutine drift_kepu_guess
end interface
!!
interface
  subroutine drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
  use module_swift
  implicit none
  real(rk) :: s,dt,r0,mu,alpha,u
  integer(ik) :: iflg
  real(rk) :: fp,c1,c2,c3
  end subroutine drift_kepu_lag
end interface
!!
interface
  subroutine drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflgn)
  use module_swift
  implicit none
  real(rk) :: s,dt,r0,mu,alpha,u
  integer(ik) :: iflgn
  real(rk) :: fp,c1,c2,c3
  end subroutine drift_kepu_new
end interface
!!
interface
  subroutine drift_kepu_p3solve(dt,r0,mu,alpha,u,s,iflg)
  use module_swift
  implicit none
  real(rk) :: dt,r0,mu,alpha,u
  integer(ik) :: iflg
  real(rk) :: s
  end subroutine drift_kepu_p3solve
end interface
!!
interface
  subroutine drift_kepu_stumpff(x,c0,c1,c2,c3)
  use module_swift
  implicit none
  real(rk) :: x
  real(rk) :: c0,c1,c2,c3
  end subroutine drift_kepu_stumpff
end interface
!!
interface
  subroutine drift_one(mu,x,y,z,vx,vy,vz,dt,iflg)
  use module_swift
  implicit none
  real(rk) :: mu,dt
  real(rk) :: x,y,z
  real(rk) :: vx,vy,vz
  integer(ik) :: iflg
  end subroutine drift_one
end interface
!!
!!
interface
  subroutine getacch_ir3(nbod,istart,x,y,z,ir3,ir)
  use module_swift
  implicit none
  integer(ik) :: nbod,istart
  real(rk) :: x(nbod),y(nbod),z(nbod)
  real(rk) :: ir3(nbod), ir(nbod)
  end subroutine getacch_ir3
end interface
!!
!!
interface
  subroutine kickvh(nbod,vxh,vyh,vzh,axh,ayh,azh,dt)
  use module_swift
  implicit none
  integer(ik) :: nbod
  real(rk) :: axh(nbod),ayh(nbod),azh(nbod)
  real(rk) :: dt
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  end subroutine kickvh
end interface
!!
!!
interface
  subroutine io_discard_mass(init,time,id,m1,r1,x1,y1,z1,vx1,vy1,vz1,iu,iwhy,fopenstat)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iwhy,iu,init,id
  real(rk) :: time
  real(rk) :: m1,r1
  real(rk) :: x1,y1,z1
  real(rk) :: vx1,vy1,vz1
  character(len = 24) :: fopenstat
  end subroutine io_discard_mass
end interface
!!
interface
  subroutine io_discard_merge(time,ip1,ip2,m1,r1,x1,y1,z1,vx1,vy1,vz1, &
             m2,r2,x2,y2,z2,vx2,vy2,vz2,mn,rn,xn,yn,zn,vxn,vyn,vzn)
  use module_swift
  use module_io
  implicit none
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
  end subroutine io_discard_merge
end interface
!!
interface
  subroutine io_dump_param(dparfile,t,tstop,dt,dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
  use module_swift
  use module_io
  implicit none
  logical(lk) :: lclose
  integer(ik) :: iflgchk
  real(rk) :: t, tstop, dt, dtout, dtdump, rmin, rmax, rmaxu, qmin
  character(len = 24) :: outfile, dparfile
  end subroutine io_dump_param
end interface
!!
interface
  subroutine io_dump_pl_symba(dplfile,nbod,mass,xh,yh,zh,vxh,vyh,vzh, &
             lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
  use module_swift
  use module_io
  implicit none
  logical(lk) :: lclose
  integer(ik) :: nbod,iflgchk
  real(rk) :: mass(nbod),rpl(nbod),j2rp2,j4rp4
  real(rk) :: xh(nbod),yh(nbod),zh(nbod),rhill(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  character(len = 24) :: dplfile
  end subroutine io_dump_pl_symba
end interface
!!
interface
  subroutine io_energy_write(i1st,t,energy,eltot,iu,fopenstat)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu,i1st
  real(rk) :: t,energy,eltot(NDIM)
  character(len = *) :: fopenstat
  end subroutine io_energy_write
end interface
!!
interface
  subroutine io_init_param(infile,t0,tstop,dt,dtout,dtdump,iflgchk, &
             rmin,rmax,rmaxu,qmin,lclose,lgas,mtiny,outfile,fopenstat)
  use module_swift
  use module_io
  implicit none
  character(len = *) :: infile
  logical(lk) :: lclose,lgas
  integer(ik) :: iflgchk
  real(rk) :: t0, tstop, dt
  real(rk) :: dtout, dtdump
  real(rk) :: rmin, rmax, rmaxu, qmin, mtiny
  character(len = 24) :: outfile, fopenstat
  end subroutine io_init_param
end interface
!!
interface
  subroutine io_init_pl_symba(infile,lclose,iflgchk,nbod,mass, &
             xh,yh,zh,vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)
  use module_swift
  use module_io
  implicit none
  logical(lk) :: lclose
  integer(ik) :: iflgchk
  character(len = *) :: infile
  integer(ik) :: nbod
  real(rk) :: mass(NTPMAX),rpl(NTPMAX),j2rp2,j4rp4
  real(rk) :: xh(NTPMAX),yh(NTPMAX),zh(NTPMAX),rhill(NTPMAX)
  real(rk) :: vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)
  end subroutine io_init_pl_symba
end interface
!!
interface
  subroutine io_open(iu,fname,fopenstat,format,ierr)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu
  character(len = *) :: fname,fopenstat,format
  integer(ik) :: ierr
  end subroutine io_open
end interface
!!
interface
  subroutine io_open_fxdr(fname,fopenstat,lflg,iu,ierr)
  use module_swift
  use module_io
  implicit none
  logical(lk) :: lflg
  integer(ik) :: iu
  character(len = *) :: fname
  character(len = 1) :: fopenstat
  integer(ik) :: ierr
  end subroutine io_open_fxdr
end interface
!!
interface
  function io_read_hdr(iu,time,nbod,nleft) result(read_hdr)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu
  integer(ik) :: nbod,nleft,read_hdr
  real(rk) :: time
  end function io_read_hdr
end interface
!!
interface
  function io_read_hdr_r(iu,time,nbod,nleft) result(read_hdr_r)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu
  integer(ik) :: nbod,nleft,read_hdr_r
  real(rk) :: time
  end function io_read_hdr_r
end interface
!!
interface
  function io_read_line(iu,id,a,e,inc,capom,omega,capm) result(read_line)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu
  integer(ik) :: id, read_line
  real(rk) :: a,e,inc,capom,omega,capm
  end function io_read_line
end interface
!!
interface
  function io_read_line_r(iu,id,a,e,inc,capom,omega,capm) result(read_line_r)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu
  integer(ik) :: id, read_line_r
  real(rk) :: a,e,inc,capom,omega,capm
  end function io_read_line_r
end interface
!!
interface
  function io_read_mass(time,nbod,mass,iu) result(read_mass)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu
  integer(ik) :: nbod, read_mass
  real(rk) :: mass(nbod),time
  end function io_read_mass
end interface
!!
interface
  function io_read_mass_r(time,nbod,mass,iu) result(read_mass_r)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu
  integer(ik) :: nbod, read_mass_r
  real(rk) :: mass(nbod),time
  end function io_read_mass_r
end interface
!!
interface
  function io_read_radius(time,nbod,radius,iu) result(read_radius)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu
  integer(ik) :: nbod, read_radius
  real(rk) :: radius(nbod),time
  end function io_read_radius
end interface
!!
interface
  function io_read_radius_r(time,nbod,radius,iu) result(read_radius_r)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu
  integer(ik) :: nbod, read_radius_r
  real(rk) :: radius(nbod),time
  end function io_read_radius_r
end interface
!!
interface
  subroutine io_write_frame(time,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh, &
             xht,yht,zht,vxht,vyht,vzht,istat,oname,iu,fopenstat)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: nbod,ntp,iu
  integer(ik) :: istat(ntp)
  real(rk) :: mass(nbod),time
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: xht(ntp),yht(ntp),zht(ntp)
  real(rk) :: vxht(ntp),vyht(ntp),vzht(ntp)
  character(len = 24) :: oname,fopenstat
  end subroutine io_write_frame
end interface
!!
interface
  subroutine io_write_frame_r(time,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh, &
             xht,yht,zht,vxht,vyht,vzht,istat,oname,iu,fopenstat)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: nbod,ntp,iu
  integer(ik) :: istat(ntp+1,NSTAT)
  real(rk) :: mass(nbod),time
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: xht(ntp),yht(ntp),zht(ntp)
  real(rk) :: vxht(ntp),vyht(ntp),vzht(ntp)
  character(len = 24) :: oname,fopenstat
  end subroutine io_write_frame_r
end interface
!!
interface
  subroutine io_write_hdr(iu,time,nbod,ntp,istat)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: nbod,ntp,istat(ntp+1,NSTAT),iu
  real(rk) :: time
  end subroutine io_write_hdr
end interface
!!
interface
  subroutine io_write_hdr_r(iu,time,nbod,ntp,istat)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: nbod,ntp,istat(ntp+1,NSTAT),iu
  real(rk) :: time
  end subroutine io_write_hdr_r
end interface
!!
interface
  subroutine io_write_line(iu,id,a,e,inc,capom,omega,capm)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu,id
  real(rk) :: a,e,inc,capom,omega,capm
  end subroutine io_write_line
end interface
!!
interface
  subroutine io_write_line_r(iu,id,a,e,inc,capom,omega,capm)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: iu,id
  real(rk) :: a,e,inc,capom,omega,capm
  end subroutine io_write_line_r
end interface
!!
interface
  subroutine io_write_mass(time,nbod,mass,oname,iu,fopenstat)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: nbod,iu
  real(rk) :: mass(nbod),time
  character(len = 24) :: oname,fopenstat
  end subroutine io_write_mass
end interface
!!
interface
  subroutine io_write_mass_r(time,nbod,mass,oname,iu,fopenstat)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: nbod,iu
  real(rk) :: mass(nbod),time
  character(len = 24) :: oname,fopenstat
  end subroutine io_write_mass_r
end interface
!!
interface
  subroutine io_write_radius(time,nbod,radius,oname,iu,fopenstat)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: nbod,iu
  real(rk) :: radius(nbod),time
  character(len = 24) :: oname,fopenstat
  end subroutine io_write_radius
end interface
!!
interface
  subroutine io_write_radius_r(time,nbod,radius,oname,iu,fopenstat)
  use module_swift
  use module_io
  implicit none
  integer(ik) :: nbod,iu
  real(rk) :: radius(nbod),time
  character(len = 24) :: oname,fopenstat
  end subroutine io_write_radius_r
end interface
!!
!!
interface
  subroutine obl_acc(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,aoblx,aobly,aoblz)
  use module_swift
  implicit none
  integer(ik) :: nbod
  real(rk) :: j2rp2,j4rp4
  real(rk) :: mass(nbod)
  real(rk) :: xh(nbod),yh(nbod),zh(nbod),irh(nbod)
  real(rk) :: aoblx(nbod),aobly(nbod),aoblz(nbod)
  end subroutine obl_acc
end interface
!!
interface
  subroutine obl_pot(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,oblpot)
  use module_swift
  implicit none
  integer(ik) :: nbod
  real(rk) :: mass(nbod)
  real(rk) :: j2rp2,j4rp4
  real(rk) :: xh(nbod),yh(nbod),zh(nbod),irh(nbod)
  real(rk) :: oblpot
  end subroutine obl_pot
end interface
!!
!!
interface
  function orbel_eget(e,m) result(ea)
  use module_swift
  implicit none
  real(rk) :: e,m
  real(rk) :: ea
  end function orbel_eget
end interface
!!
interface
  function orbel_ehie(e,m) result(ea)
  use module_swift
  implicit none
  real(rk) :: e, m
  real(rk) :: ea
  end function orbel_ehie
end interface
!!
interface
  function orbel_ehybrid(e,m) result(ea)
  use module_swift
  implicit none
  real(rk) :: e,m
  real(rk) :: ea
  end function orbel_ehybrid
end interface
!!
interface
  subroutine orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,x,y,z,vx,vy,vz)
  use module_swift
  implicit none
  integer(ik) :: ialpha
  real(rk) :: gm,a,e,inc,capom,omega,capm
  real(rk) :: x,y,z,vx,vy,vz
  end subroutine orbel_el2xv
end interface
!!
interface
  function orbel_esolmd(e,m) result(ea)
  use module_swift
  implicit none
  real(rk) :: e,m
  real(rk) :: ea
  end function orbel_esolmd
end interface
!!
interface
  function orbel_fget(e,capn) result(ea)
  use module_swift
  implicit none
  real(rk) :: e,capn
  real(rk) :: ea
  end function orbel_fget
end interface
!!
interface
  function orbel_fhybrid(e,capn) result(ea)
  use module_swift
  implicit none
  real(rk) :: e,capn
  real(rk) :: ea
  end function orbel_fhybrid
end interface
!!
interface
  function orbel_flon(e,capn) result(ea)
  use module_swift
  implicit none
  real(rk) :: e,capn
  real(rk) :: ea
  end function orbel_flon
end interface
!!
interface
  subroutine orbel_scget(angle,sx,cx)
  use module_swift
  implicit none
  real(rk) :: angle
  real(rk) :: sx,cx
  end subroutine orbel_scget
end interface
!!
interface
  subroutine orbel_schget(angle,shx,chx)
  use module_swift
  implicit none
  real(rk) :: angle
  real(rk) :: shx,chx
  end subroutine orbel_schget
end interface
!!
interface
  subroutine orbel_xv2aeq(x,y,z,vx,vy,vz,gmsum,ialpha,a,e,q)
  use module_swift
  implicit none
  real(rk) :: x,y,z,vx,vy,vz,gmsum
  integer(ik) :: ialpha
  real(rk) :: a,e,q
  end subroutine orbel_xv2aeq
end interface
!!
interface
  subroutine orbel_xv2aqt(x,y,z,vx,vy,vz,gmsum,ialpha,a,q,capm,tperi)
  use module_swift
  implicit none
  real(rk) :: x,y,z,vx,vy,vz,gmsum
  integer(ik) :: ialpha
  real(rk) :: a,q,capm,tperi
  end subroutine orbel_xv2aqt
end interface
!!
interface
  subroutine orbel_xv2el(x,y,z,vx,vy,vz,gmsum,ialpha,a,e,inc,capom,omega,capm)
  use module_swift
  implicit none
  real(rk) :: x,y,z,vx,vy,vz,gmsum
  integer(ik) :: ialpha
  real(rk) :: a,e,inc,capom,omega,capm
  end subroutine orbel_xv2el
end interface
!!
interface
  function orbel_zget(q) result(ea)
  use module_swift
  implicit none
  real(rk) :: q
  real(rk) :: ea
  end function orbel_zget
end interface
!!
!!
interface
  subroutine rmvs_chk_ind(xr,yr,zr,vxr,vyr,vzr,dt,r2crit,r2critp,iflag)
  use module_swift
  use module_rmvs
  implicit none
  real(rk) :: xr,yr,zr,vxr,vyr,vzr,dt,r2crit,r2critp
  integer(ik) :: iflag
  end subroutine rmvs_chk_ind
end interface
!!
!!
interface
  subroutine symba5_chk(rhill,nbod,ip1,ip2,mass,xh,yh,zh,vxh,vyh,vzh,dt,irec,icflg,lvdotr)
  use module_swift
  use module_symba5
  implicit none
  integer(ik) :: nbod,irec,ip1,ip2
  real(rk) :: mass(nbod),xh(nbod),yh(nbod),zh(nbod),dt
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod),rhill(nbod)
  integer(ik) :: icflg
  logical(lk) :: lvdotr
  end subroutine symba5_chk
end interface
!!
interface
  subroutine symba5_chk2(rhill,xh,yh,zh,vxh,vyh,vzh,dt,irec,icflg)
  use module_swift
  use module_symba5
  implicit none
  integer(ik) :: irec
  real(rk) :: xh(2),yh(2),zh(2),dt
  real(rk) :: vxh(2),vyh(2),vzh(2),rhill(2)
  integer(ik) :: icflg
  end subroutine symba5_chk2
end interface
!!
interface
  subroutine symba5_getacch(nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,axh,ayh,azh,mtiny,ielc,ielst)
  use module_swift
  use module_symba5
  implicit none
  integer(ik) :: nbod,nbodm
  integer(ik) :: ielst(NENMAX,2),ielc
  real(rk) :: mass(nbod),j2rp2,j4rp4,mtiny
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: axh(nbod),ayh(nbod),azh(nbod)
  end subroutine symba5_getacch
end interface
!!
interface
  subroutine symba5_helio_drift(nbod,ielev,irec,mass,xh,yh,zh,vxb,vyb,vzb,dt)
  use module_swift
  use module_symba5
  implicit none
  integer(ik) :: nbod,irec
  integer(ik) :: ielev(nbod)
  real(rk) :: mass(nbod),dt
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxb(nbod),vyb(nbod),vzb(nbod)
  end subroutine symba5_helio_drift
end interface
!!
interface
  subroutine symba5_enc_drift(nbod,ielc,ielst,ielev,irec,mass,xh,yh,zh,vxb,vyb,vzb,dt)
  use module_swift
  use module_symba5
  implicit none
  integer(ik) :: nbod,irec,ielc
  integer(ik) :: ielst(NENMAX,2),ielev(nbod)
  real(rk) :: mass(nbod),dt
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxb(nbod),vyb(nbod),vzb(nbod)
  end subroutine symba5_enc_drift
end interface
!!
interface
  subroutine symba5_helio_getacch(iflg,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,axh,ayh,azh)
  use module_swift
  implicit none
  integer(ik) :: nbod,nbodm,iflg
  real(rk) :: mass(nbod),j2rp2,j4rp4
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: axh(nbod),ayh(nbod),azh(nbod)
  end subroutine symba5_helio_getacch
end interface
!!
interface
  subroutine symba5_kick(nbod,mass,irec,ielev,rhill,xh,yh,zh,vxb,vyb,vzb,dt,sgn,ielc,ielst)
  use module_swift
  use module_symba5
  implicit none
  integer(ik) :: nbod,irec
  integer(ik) :: ielev(nbod),ielst(NENMAX,2),ielc
  real(rk) :: mass(nbod),dt,rhill(nbod),sgn
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxb(nbod),vyb(nbod),vzb(nbod)
  end subroutine symba5_kick
end interface
!!
interface
  subroutine symba5_merge(t, dt, nbod, nbodm, ip1, ip2, mass, xh, yh, zh, vxb, vyb, vzb, &
             ireci, lvdotrold, ibound, rpl, mergelst, mergecnt, rhill, eoff, ielc, ielst)
  use module_swift
  use module_symba5
  implicit none
  integer(ik) :: nbod, nbodm, ireci, ip1, ip2
  real(rk) :: t, dt
  logical(lk) :: lvdotrold
  real(rk) :: mass(nbod), xh(nbod), yh(nbod), zh(nbod), eoff
  real(rk) :: vxb(nbod), vyb(nbod), vzb(nbod), rpl(nbod), rhill(nbod)
  integer(ik) :: mergelst(NENMAX,2), mergecnt
  integer(ik) :: ielst(NENMAX,2), ielc, ibound(nbod)
  end subroutine symba5_merge
end interface
!!
interface
  subroutine symba5_nbodm(nbod,mass,mtiny,nbodm)
  use module_swift
  use module_symba5
  integer(ik) :: nbod
  real(rk) :: mass(nbod), mtiny
  integer(ik) :: nbodm
  end subroutine symba5_nbodm
end interface
!!
interface
  subroutine symba5_step_helio(i1st,nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,dt)
  use module_swift
  implicit none
  integer(ik) :: nbod,i1st,nbodm
  real(rk) :: mass(nbod),dt,j2rp2,j4rp4
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  end subroutine symba5_step_helio
end interface
!!
interface
  subroutine symba5_step_interp(time,ielev,nbod,nbodm,mass,rhill,j2rp2,j4rp4, &
             rpl,xh,yh,zh,vxh,vyh,vzh,dt,mergelst,mergecnt,eoff,ielc,ielst,mtiny,ibound)
  use module_swift
  use module_symba5
  implicit none
  integer(ik) :: ielev(nbod),ielst(NENMAX,2),ielc
  real(rk) :: mass(nbod),dt,j2rp2,j4rp4,time,mtiny
  integer(ik) :: nbod,nbodm
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: rpl(nbod),rhill(nbod),eoff
  integer(ik) :: mergelst(NENMAX,2),mergecnt,ibound(nbod)
  end subroutine symba5_step_interp
end interface
!!
interface
  subroutine symba5_step_pl(i1st,time,nbod,nbodm,mass,j2rp2,j4rp4, &
             xh,yh,zh,vxh,vyh,vzh,dt,lclose,rpl,isenc,mergelst,    &
             mergecnt,iecnt,eoff,rhill,mtiny,ibound)
  use module_swift
  use module_symba5
  implicit none
  logical(lk) :: lclose
  integer(ik) :: nbod,i1st,nbodm
  real(rk) :: mass(nbod),dt,time,j2rp2,j4rp4,mtiny
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: rpl(nbod),rhill(nbod),eoff
  integer(ik) :: isenc
  integer(ik) :: iecnt(nbod),ielev(nbod)
  integer(ik) :: mergelst(NENMAX,2),mergecnt,ibound(nbod)
  end subroutine symba5_step_pl
end interface
!!
interface
  recursive subroutine symba5_step_recur(t,nbod,nbodm,mass,ireci,ielev, &
                       rhill,xh,yh,zh,vxb,vyb,vzb,rpl,mergelst,mergecnt,dt0, &
                       eoff,lvdotr,ibound,ielc,ielst)
  use module_swift
  use module_symba5
  implicit none
  integer(ik) :: nbod,ireci,nbodm
  integer(ik) :: ielev(nbod),ielst(NENMAX,2),ielc
  real(rk) :: mass(nbod),rhill(nbod),t,dt0
  logical(lk) :: lvdotr(NTPMAX)
  integer(ik) :: mergelst(NENMAX,2),mergecnt,ibound(nbod)
  real(rk) :: xh(nbod),yh(nbod),zh(nbod),eoff
  real(rk) :: vxb(nbod),vyb(nbod),vzb(nbod),rpl(nbod)
  end subroutine symba5_step_recur
end interface
!!
!!
interface
  subroutine symba5_gas_a_drag(time,dt,nbod,nbodm,iencounter,mstar,x,y,z,vx,vy,vz, &
             deng0,gpower,deng0s,zscale,rgi,rgf,tgdecay)
  use module_swift
  implicit none
  integer(ik) :: nbod,nbodm,iencounter(nbod)
  real(rk) :: mstar,time,dt
  real(rk) :: x(nbod),y(nbod),z(nbod)
  real(rk) :: vx(nbod),vy(nbod),vz(nbod)
  real(rk) :: deng0,gpower,deng0s,zscale
  real(rk) :: rgi,rgf,tgdecay,drsq
  end subroutine symba5_gas_a_drag
end interface
!!
interface
  subroutine symba5_gas_a_typeI(time,dt,nbodm,mass,rgap,wgap,xh,yh,zh,vxh,vyh,vzh, &
             sig0,spower,zscale,rgi,rgf,tgdecay,ca,ce)
  use module_swift
  implicit none
  integer(ik) :: nbodm
  real(rk) :: mass(nbodm),time,dt
  real(rk) :: xh(nbodm),yh(nbodm),zh(nbodm)
  real(rk) :: vxh(nbodm),vyh(nbodm),vzh(nbodm)
  real(rk) :: spower,sig0,zscale
  real(rk) :: rgi,rgf,tgdecay,drsq,ca,ce
  real(rk) :: rgap(nbodm),wgap(nbodm)
  end subroutine symba5_gas_a_typeI
end interface
!!
interface
  subroutine symba5_gas_dragcoef(ma,kn,re,cdr)
  use module_swift
  implicit none
  real(rk) :: ma, kn, re
  real(rk) :: cdr
  end subroutine symba5_gas_dragcoef
end interface
!!
interface
  subroutine symba5_gas_drag_kick(time,nbod,nbodm,mass,rhill,xh,yh,zh,vxh,vyh,vzh,dt)
  use module_swift
  implicit none
  integer(ik) :: nbod,nbodm
  real(rk) :: mass(nbod),rhill(nbod),dt,time
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  end subroutine symba5_gas_drag_kick
end interface
!!
interface
  subroutine symba5_gas_step_helio(i1st,nbod,nbodm,mass,j2rp2,j4rp4,lgas,xh,yh,zh,vxh,vyh,vzh,dt)
  use module_swift
  implicit none
  logical(lk) :: lgas
  integer(ik) :: nbod,i1st,nbodm
  real(rk) :: mass(nbod),dt,j2rp2,j4rp4
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  end subroutine symba5_gas_step_helio
end interface
!!
interface
  subroutine symba5_gas_step_pl(i1st,time,nbod,nbodm,mass,j2rp2,j4rp4, &
             xh,yh,zh,vxh,vyh,vzh,dt,lclose,lgas,rpl,isenc, &
             mergelst,mergecnt,iecnt,eoff,rhill,mtiny,ibound)
  use module_swift
  use module_symba5
  implicit none
  logical(lk) lclose,lgas
  integer(ik) :: nbod,i1st,nbodm
  real(rk) :: mass(nbod),dt,time,j2rp2,j4rp4,mtiny
  real(rk) :: xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: rpl(nbod),rhill(nbod),eoff
  integer(ik) :: isenc
  integer(ik) :: iecnt(nbod),ielev(nbod)
  integer(ik) :: mergelst(NENMAX,2),mergecnt,ibound(nbod)
  end subroutine symba5_gas_step_pl
end interface
!!
!!
interface
  function util_disk_mass(s_0, a_in, a_out, p, a_snow, f_snow) result(m_disk)
  use module_swift
  implicit none
  real(rk) :: s_0, a_in, a_out, p, a_snow, f_snow
  real(rk) :: m_disk
  end function util_disk_mass
end interface
!!
interface
  subroutine util_exit(iflag)
  use module_swift
  implicit none
  integer(ik) :: iflag
  end subroutine util_exit
end interface
!!
interface
  subroutine util_hills(nbod,mass,xh,yh,zh,vxh,vyh,vzh,r2hill)
  use module_swift
  implicit none
  integer(ik) :: nbod
  real(rk) :: mass(nbod),xh(nbod),yh(nbod),zh(nbod)
  real(rk) :: vxh(nbod),vyh(nbod),vzh(nbod)
  real(rk) :: r2hill(nbod)
  end subroutine util_hills
end interface
!!
interface
  subroutine util_hills1(mstar,mpl,xh,yh,zh,vxh,vyh,vzh,rhill)
  use module_swift
  implicit none
  real(rk) :: mstar,mpl,xh,yh,zh
  real(rk) :: vxh,vyh,vzh
  real(rk) :: rhill
  end subroutine util_hills1
end interface
!!
interface
  subroutine util_mass_peri(iflg,nbod,x,y,z,vx,vy,vz,mass,isperi,peri,lperi)
  use module_swift
  implicit none
  integer(ik) :: nbod,iflg
  real(rk) :: x(nbod),y(nbod),z(nbod),mass(nbod)
  real(rk) :: vx(nbod),vy(nbod),vz(nbod),gm
  logical(lk) :: lperi(nbod)
  integer(ik) :: isperi(nbod)
  real(rk) :: peri(nbod)
  end subroutine util_mass_peri
end interface
!!
interface
  function util_power_law(p, xmin, xmax) result(x)
  use module_swift
  implicit none
  real(rk) :: p, xmin, xmax
  real(rk) :: x
  end function util_power_law
end interface
!!
interface
  function util_randomu() result(ran)
  use module_swift
  implicit none
  real(rk) :: ran
  end function util_randomu
end interface
!!
interface
  function util_rayleigh(rms, xmin, xmax) result(x)
  use module_swift
  implicit none
  real(rk) :: rms, xmin, xmax
  real(rk) :: x
  end function util_rayleigh
end interface
!!
interface
  subroutine util_version
  use module_swift
  implicit none
  end subroutine util_version
end interface
!!
!!
end module module_interfaces
