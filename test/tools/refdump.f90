! ------------------------------------------------------------------------------
! This file is part of gfnff.
!
! Copyright (C) 2023-2026 Philipp Pracht
!
! gfnff is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! gfnff is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with gfnff.  If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------
!> Full precision reference dump: prints energies, gradients, stress tensors
!> and internal arrays at full precision so a refactor can be diffed against
!> the untouched code for bit-identical output; unit tests only check
!> tolerances, so a change that moves the last digit shows up here only.
!> Sections: MOL (molecular path), PBC (periodic path), SOLV (ALPB internals).
!> Usage: gfnff-refdump [jitter-amplitude]
program refdump
  use iso_fortran_env,only:wp => real64
  use gfnff_interface
  use coffeine,only:cnat => testnat,cat => testat,cxyz => testxyz
  use supermol,only:snat => testnat,sat => testat,sxyz => testxyz
  use sio2
  implicit none

  real(wp) :: jitamp
  character(len=32) :: arg

  jitamp = 0.002_wp
  if (command_argument_count() >= 1) then
    call get_command_argument(1,arg)
    read (arg,*) jitamp
  end if

  write (*,'(a)') '## MOL'
  call mol_run('caffeine  ',cnat,cat,cxyz,0,'')
  call mol_run('caffeine+1',cnat,cat,cxyz,1,'')
  call mol_run('caffeine-1',cnat,cat,cxyz,-1,'')
  call mol_run('caff-h2o  ',cnat,cat,cxyz,0,'h2o')
  call mol_run('supermol  ',snat,sat,sxyz,0,'')
  call mol_run('supermol-s',snat,sat,sxyz,0,'h2o')
  call mol_water_cluster()
  call mol_halogen()
  call mol_quartz()

  write (*,'(a)') '## PBC'
  call pbc_sio2(1)
  call pbc_sio2(2)
  call pbc_water(3)
  call pbc_water(4)
  call pbc_caffeine()

  write (*,'(a)') '## SOLV'
  call solv_run('solv-caff ',cnat,cat,cxyz)
  call solv_run('solv-super',snat,sat,sxyz)
  call solv_water()

contains

  subroutine mol_report(tag,nat,e,g,calc)
    !***********************************************************************
    !* Molecular dump: energy/gradient plus H/X-bond counts, bond-pair range
    !* and non-bonded exponent sum, the topology terms a setup bug moves.
    !***********************************************************************
    character(len=*),intent(in) :: tag
    integer,intent(in) :: nat
    real(wp),intent(in) :: e,g(3,nat)
    type(gfnff_data),intent(in) :: calc
    integer :: i
    write (*,'(a,a,1x,es24.16)') tag,' E    ',e
    write (*,'(a,a,1x,es24.16)') tag,' |g|  ',sqrt(sum(g*g))
    write (*,'(a,a,1x,es24.16)') tag,' netf ',sum(abs(sum(g,dim=2)))
    write (*,'(a,a,3(1x,i0))') tag,' hb   ',calc%nlist%nhb1,calc%nlist%nhb2,calc%nlist%nxb
    write (*,'(a,a,2(1x,i0))') tag,' bp   ',minval(calc%neigh%bpair),maxval(calc%neigh%bpair)
    write (*,'(a,a,1x,es24.16)') tag,' alp  ',sum(abs(calc%topo%alphanb))
    write (*,'(a,a,1x,es24.16)') tag,' chkg ',csum(reshape(g,[3*nat]))
    do i = 1,min(nat,4)
      write (*,'(a,a,i0,3(1x,es24.16))') tag,' g',i,g(1,i),g(2,i),g(3,i)
    end do
  end subroutine mol_report

  subroutine mol_run(tag,nat,at,xyz0,ichrg,solv)
    !***********************************************************************
    !* Runs 3 jittered singlepoints for one system/charge/solvent combination.
    !***********************************************************************
    character(len=*),intent(in) :: tag
    integer,intent(in) :: nat,at(nat),ichrg
    real(wp),intent(in) :: xyz0(3,nat)
    character(len=*),intent(in) :: solv
    integer :: io,istep
    real(wp) :: xyz(3,nat),g(3,nat),e
    type(gfnff_data) :: calc
    xyz = xyz0
    if (len_trim(solv) > 0) calc%solvent = solv
    call gfnff_initialize(nat,at,xyz,calc,ichrg=ichrg,iostat=io,printlevel=0)
    if (io /= 0) then
      write (*,'(a,a)') tag,' INIT FAILED'
      return
    end if
    do istep = 1,3
      call jitter(nat,xyz,3.0_wp*jitamp)
      call calc%singlepoint(nat,at,xyz,e,g,printlevel=0)
      call mol_report(tag,nat,e,g,calc)
    end do
    call calc%deallocate()
  end subroutine mol_run

  subroutine mol_water_cluster()
    !***********************************************************************
    !* Dense water cluster; H-bond lists are what a setup change disturbs most.
    !***********************************************************************
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    call water_cluster(3,5.4_wp,nat,at,xyz)
    call mol_run('water81   ',nat,at,xyz,0,'')
    call mol_run('water81-s ',nat,at,xyz,0,'h2o')
  end subroutine mol_water_cluster

  subroutine mol_halogen()
    !***********************************************************************
    !* Halogen-bond geometry: the other bond list a setup change can disturb.
    !***********************************************************************
    integer,parameter :: nat = 10
    integer :: at(nat)
    real(wp) :: xyz(3,nat)
    at = [53,53,8,6,8,1,1,7,1,1]
    xyz(:,1) = [0.0_wp,0.0_wp,0.0_wp]
    xyz(:,2) = [0.0_wp,0.0_wp,5.04_wp]
    xyz(:,3) = [0.0_wp,0.0_wp,10.4_wp]
    xyz(:,4) = [0.0_wp,0.0_wp,12.6_wp]
    xyz(:,5) = [0.0_wp,0.0_wp,14.8_wp]
    xyz(:,6) = [1.8_wp,0.0_wp,16.2_wp]
    xyz(:,7) = [-1.8_wp,0.4_wp,16.4_wp]
    xyz(:,8) = [0.3_wp,3.6_wp,13.1_wp]
    xyz(:,9) = [1.6_wp,4.4_wp,14.0_wp]
    xyz(:,10) = [-1.4_wp,4.6_wp,13.4_wp]
    call mol_run('halogen   ',nat,at,xyz,0,'')
  end subroutine mol_halogen

  subroutine mol_quartz()
    !***********************************************************************
    !* Many-element system without hydrogen, run through the molecular path.
    !***********************************************************************
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: lat(3,3)
    call super_sio2(2,nat,at,xyz,lat)
    call mol_run('sio2-mol  ',nat,at,xyz,0,'')
  end subroutine mol_quartz

  subroutine pbc_report(tag,nat,e,g,sig)
    !***********************************************************************
    !* Periodic dump: energy, gradient, stress, plus a few gradient entries
    !* that catch sign or index slips a norm would average away.
    !***********************************************************************
    character(len=*),intent(in) :: tag
    integer,intent(in) :: nat
    real(wp),intent(in) :: e,g(3,nat),sig(3,3)
    integer :: i
    write (*,'(a,a,1x,es24.16)') tag,' E     ',e
    write (*,'(a,a,1x,es24.16)') tag,' |g|   ',sqrt(sum(g*g))
    write (*,'(a,a,1x,es24.16)') tag,' netf  ',sum(abs(sum(g,dim=2)))
    write (*,'(a,a,1x,es24.16)') tag,' chkg  ',csum(reshape(g,[3*nat]))
    do i = 1,3
      write (*,'(a,a,i0,3(1x,es24.16))') tag,' sig',i,sig(i,1),sig(i,2),sig(i,3)
    end do
    do i = 1,min(nat,4)
      write (*,'(a,a,i0,3(1x,es24.16))') tag,' g',i,g(1,i),g(2,i),g(3,i)
    end do
  end subroutine pbc_report

  subroutine pbc_sio2(nr)
    !***********************************************************************
    !* Periodic SiO2 supercell (nr^3 cells); 3 geometry steps exercise the
    !* neighbour-list refresh path.
    !***********************************************************************
    integer,intent(in) :: nr
    integer :: nat,io,i
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),g(:,:)
    real(wp) :: lat(3,3),e,sig(3,3)
    character(len=16) :: tag
    type(gfnff_data) :: calc
    call super_sio2(nr,nat,at,xyz,lat)
    allocate (g(3,nat))
    write (tag,'(a,i0,a)') 'sio2-',nr,'x'
    call gfnff_initialize(nat,at,xyz,calc,lattice=lat,npbc=3,iostat=io,printlevel=0)
    if (io /= 0) then
      write (*,'(a,a)') trim(tag),' INIT FAILED'
      return
    end if
    do i = 1,3
      call jitter(nat,xyz,jitamp)
      call calc%singlepoint(nat,at,xyz,e,g,lattice=lat,sigma=sig,printlevel=0)
      call pbc_report(trim(tag),nat,e,g,sig)
    end do
    call calc%deallocate()
  end subroutine pbc_sio2

  subroutine pbc_water(nr)
    !***********************************************************************
    !* Periodic water box, side length nr*5.858: singlepoint plus stress.
    !***********************************************************************
    integer,intent(in) :: nr
    integer :: nat,io,i
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),g(:,:)
    real(wp) :: lat(3,3),e,sig(3,3)
    character(len=16) :: tag
    type(gfnff_data) :: calc
    call water_cluster(nr,5.858_wp,nat,at,xyz)
    lat = 0.0_wp
    lat(1,1) = nr*5.858_wp
    lat(2,2) = nr*5.858_wp
    lat(3,3) = nr*5.858_wp
    allocate (g(3,nat))
    write (tag,'(a,i0,a)') 'water-',nr,'^3'
    call gfnff_initialize(nat,at,xyz,calc,lattice=lat,npbc=3,iostat=io,printlevel=0)
    if (io /= 0) then
      write (*,'(a,a)') trim(tag),' INIT FAILED'
      return
    end if
    do i = 1,3
      call jitter(nat,xyz,jitamp)
      call calc%singlepoint(nat,at,xyz,e,g,lattice=lat,sigma=sig,printlevel=0)
      call pbc_report(trim(tag),nat,e,g,sig)
    end do
    call calc%deallocate()
  end subroutine pbc_water

  subroutine pbc_caffeine()
    !***********************************************************************
    !* Isolated molecule in a large box; result should track the molecular
    !* path closely.
    !***********************************************************************
    integer :: io
    real(wp) :: xyz(3,cnat),g(3,cnat),lat(3,3),e,sig(3,3)
    type(gfnff_data) :: calc
    xyz = cxyz
    lat = 0.0_wp
    lat(1,1) = 20.0_wp
    lat(2,2) = 20.0_wp
    lat(3,3) = 20.0_wp
    call gfnff_initialize(cnat,cat,xyz,calc,lattice=lat,npbc=3,iostat=io,printlevel=0)
    if (io /= 0) then
      write (*,'(a)') 'caffeine-box INIT FAILED'
      return
    end if
    call calc%singlepoint(cnat,cat,xyz,e,g,lattice=lat,sigma=sig,printlevel=0)
    call pbc_report('caffeine-box',cnat,e,g,sig)
    call calc%deallocate()
  end subroutine pbc_caffeine

  subroutine solv_run(tag,nat,at,xyz0)
    !***********************************************************************
    !* Reaches into ALPB internals (Born radii, SASA, derivative matrices):
    !* the energy alone would hide a change that cancels in the sum.
    !***********************************************************************
    character(len=*),intent(in) :: tag
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz0(3,nat)
    integer :: io,i
    real(wp) :: xyz(3,nat),g(3,nat),e
    type(gfnff_data) :: calc
    xyz = xyz0
    calc%solvent = 'h2o'
    call gfnff_initialize(nat,at,xyz,calc,ichrg=0,iostat=io,printlevel=0)
    if (io /= 0) then
      write (*,'(a,a)') tag,' INIT FAILED'
      return
    end if
    call calc%singlepoint(nat,at,xyz,e,g,printlevel=0)
    write (*,'(a,a,1x,es24.16)') tag,' E     ',e
    write (*,'(a,a,1x,es24.16)') tag,' e_es  ',calc%res%e_es
    write (*,'(a,a,1x,es24.16)') tag,' gborn ',calc%res%g_born
    write (*,'(a,a,1x,es24.16)') tag,' gsasa ',calc%res%g_sasa
    write (*,'(a,a,1x,es24.16)') tag,' ghb   ',calc%res%g_hb
    write (*,'(a,a,1x,es24.16)') tag,' gsolv ',calc%res%g_solv
    write (*,'(a,a,1x,es24.16)') tag,' gnorm ',calc%res%gnorm
    write (*,'(a,a,1x,es24.16)') tag,' chkq  ',csum(calc%nlist%q)
    associate (s => calc%solvation)
      if (allocated(s%brad)) write (*,'(a,a,1x,es24.16)') tag,' brad  ',csum(s%brad)
      if (allocated(s%sasa)) write (*,'(a,a,1x,es24.16)') tag,' sasa  ',csum(s%sasa)
      if (allocated(s%brdr)) write (*,'(a,a,1x,es24.16)') tag,' brdr  ', &
         & csum(reshape(s%brdr,[size(s%brdr)]))
      if (allocated(s%dsdr)) write (*,'(a,a,1x,es24.16)') tag,' dsdr  ', &
         & csum(reshape(s%dsdr,[size(s%dsdr)]))
      if (allocated(s%bornMat)) write (*,'(a,a,1x,es24.16)') tag,' bmat  ', &
         & csum(reshape(s%bornMat,[size(s%bornMat)]))
    end associate
    do i = 1,min(nat,4)
      write (*,'(a,a,i0,3(1x,es24.16))') tag,' g',i,g(1,i),g(2,i),g(3,i)
    end do
    call calc%deallocate()
  end subroutine solv_run

  subroutine solv_water()
    !***********************************************************************
    !* Solvated 16-water cluster (2x2x2), same dump as solv_run.
    !***********************************************************************
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    call water_cluster(2,5.4_wp,nat,at,xyz)
    call solv_run('solv-w24  ',nat,at,xyz)
  end subroutine solv_water

  pure function csum(a) result(s)
    !***********************************************************************
    !* Order-sensitive checksum; catches a changed summation order or a
    !* permuted list that a plain sum would hide.
    !***********************************************************************
    real(wp),intent(in) :: a(:)
    real(wp) :: s
    integer :: i
    s = 0.0_wp
    do i = 1,size(a)
      s = s+a(i)*real(i,wp)
    end do
  end function csum

  subroutine jitter(nat,xyz,d)
    !***********************************************************************
    !* Deterministic displacement (no RNG): identical geometries across builds.
    !***********************************************************************
    integer,intent(in) :: nat
    real(wp),intent(inout) :: xyz(3,nat)
    real(wp),intent(in) :: d
    integer :: k
    do k = 1,nat
      xyz(1,k) = xyz(1,k)+d*sin(0.37_wp*k)
      xyz(2,k) = xyz(2,k)+d*sin(0.71_wp*k)
      xyz(3,k) = xyz(3,k)+d*sin(1.13_wp*k)
    end do
  end subroutine jitter

  subroutine super_sio2(nr,nat,at,xyz,lat)
    !***********************************************************************
    !* Builds an nr x nr x nr SiO2 supercell from the unit cell.
    !***********************************************************************
    integer,intent(in) :: nr
    integer,intent(out) :: nat
    integer,allocatable,intent(out) :: at(:)
    real(wp),allocatable,intent(out) :: xyz(:,:)
    real(wp),intent(out) :: lat(3,3)
    integer :: ix,iy,iz,k,j
    real(wp) :: sh(3)
    nat = sio2nat*nr**3
    if (allocated(at)) deallocate (at)
    if (allocated(xyz)) deallocate (xyz)
    allocate (at(nat),xyz(3,nat))
    k = 0
    do ix = 0,nr-1
      do iy = 0,nr-1
        do iz = 0,nr-1
          sh = ix*sio2lattice(:,1)+iy*sio2lattice(:,2)+iz*sio2lattice(:,3)
          do j = 1,sio2nat
            k = k+1
            at(k) = sio2at(j)
            xyz(:,k) = sio2xyz(:,j)+sh
          end do
        end do
      end do
    end do
    lat(:,1) = nr*sio2lattice(:,1)
    lat(:,2) = nr*sio2lattice(:,2)
    lat(:,3) = nr*sio2lattice(:,3)
  end subroutine super_sio2

  subroutine water_cluster(nr,d,nat,at,xyz)
    !***********************************************************************
    !* Builds an nr x nr x nr grid of water molecules, spacing d.
    !***********************************************************************
    integer,intent(in) :: nr
    real(wp),intent(in) :: d
    integer,intent(out) :: nat
    integer,allocatable,intent(out) :: at(:)
    real(wp),allocatable,intent(out) :: xyz(:,:)
    real(wp),parameter :: rOH = 1.81_wp,hoh = 1.824_wp
    integer :: ix,iy,iz,k
    real(wp) :: o(3),ph
    nat = 3*nr**3
    if (allocated(at)) deallocate (at)
    if (allocated(xyz)) deallocate (xyz)
    allocate (at(nat),xyz(3,nat))
    k = 0
    do ix = 0,nr-1
      do iy = 0,nr-1
        do iz = 0,nr-1
          o = [ix*d,iy*d,iz*d]
          ph = 0.7_wp*(ix+2*iy+3*iz)
          k = k+1; at(k) = 8; xyz(:,k) = o
          k = k+1; at(k) = 1; xyz(:,k) = o+[rOH*cos(ph),rOH*sin(ph),0.2_wp]
          k = k+1; at(k) = 1; xyz(:,k) = o+[rOH*cos(ph+hoh),rOH*sin(ph+hoh),-0.2_wp]
        end do
      end do
    end do
  end subroutine water_cluster

end program refdump
