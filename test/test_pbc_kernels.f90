module test_pbc_kernels
  !> Guards for the periodic neighbour-list and Ewald kernels not covered by a reference value.
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use iso_fortran_env,only:wp => real64
  use gfnff_interface
  use sio2
!$ use omp_lib,only:omp_get_max_threads,omp_set_num_threads
  implicit none
  private

  public :: collect_pbc_kernels

contains

  subroutine collect_pbc_kernels(testsuite)
    !***********************************
    !* Registers the PBC kernel test cases with testdrive.
    !***********************************
    type(unittest_type),allocatable,intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("PBC HB list capacity is reused    ",test_list_capacity),  &
      new_unittest("PBC lists survive geometry steps  ",test_list_reuse),     &
      new_unittest("PBC halogen bonds without HB H    ",test_xb_without_h),   &
      new_unittest("PBC results independent of threads",test_thread_invariance), &
      new_unittest("PBC gradient vs finite differences",test_grad_fd),        &
      new_unittest("PBC gradient on a symmetric cell  ",test_grad_fd_quartz), &
      new_unittest("PBC stress vs finite differences  ",test_stress_fd),      &
      new_unittest("PBC stress on a reused calculator ",test_stress_reused),  &
      new_unittest("PBC bond pair matrix values       ",test_bpair_values),   &
      new_unittest("PBC charge on the second fragment ",test_frag_charge),    &
      new_unittest("HB to a linear two-neighbour N    ",test_hb_linear_n),    &
      new_unittest("HB list picks up an incoming donor",test_hb_list_growth)   &
    ]
  end subroutine collect_pbc_kernels

  subroutine water_box(nr,nat,at,xyz,lat)
    !***********************************
    !* Cubic box of nr**3 rotated water molecules; no symmetry the kernels could exploit.
    !***********************************
    integer,intent(in) :: nr
    integer,intent(out) :: nat
    integer,allocatable,intent(out) :: at(:)
    real(wp),allocatable,intent(out) :: xyz(:,:)
    real(wp),intent(out) :: lat(3,3)

    real(wp),parameter :: d = 5.858_wp    !> grid spacing, ~1 g/cm3
    real(wp),parameter :: rOH = 1.81_wp
    real(wp),parameter :: hoh = 1.824_wp
    integer :: ix,iy,iz,k
    real(wp) :: o(3),ph

    nat = 3*nr**3
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
    lat = 0.0_wp
    lat(1,1) = nr*d
    lat(2,2) = nr*d
    lat(3,3) = nr*d
  end subroutine water_box

  subroutine quartz_super(nr,nat,at,xyz,lat)
    !***********************************
    !* nr**3 alpha-quartz supercell; a symmetric crystal with many equidistant images.
    !***********************************
    integer,intent(in) :: nr
    integer,intent(out) :: nat
    integer,allocatable,intent(out) :: at(:)
    real(wp),allocatable,intent(out) :: xyz(:,:)
    real(wp),intent(out) :: lat(3,3)
    integer :: ix,iy,iz,k,j
    real(wp) :: sh(3)

    nat = sio2nat*nr**3
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
  end subroutine quartz_super

  subroutine step_geometry(nat,xyz,d)
    !***********************************
    !* Displace every atom a little, as an optimisation step would.
    !***********************************
    integer,intent(in) :: nat
    real(wp),intent(inout) :: xyz(3,nat)
    real(wp),intent(in) :: d
    integer :: k
    do k = 1,nat
      xyz(1,k) = xyz(1,k)+d*sin(0.37_wp*k)
      xyz(2,k) = xyz(2,k)+d*sin(0.71_wp*k)
      xyz(3,k) = xyz(3,k)+d*sin(1.13_wp*k)
    end do
  end subroutine step_geometry

  subroutine test_list_capacity(error)
    !***********************************
    !* HB list capacity must stay fixed while the candidate count drifts under motion.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io,istep,cap1,cap2,cap3
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    real(wp) :: lat(3,3),energy

    call water_box(5,nat,at,xyz,lat)
    allocate (grad(3,nat),source=0.0_wp)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lat,npbc=3,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    cap1 = 0; cap2 = 0; cap3 = 0
    do istep = 1,6
      call step_geometry(nat,xyz,0.002_wp)
      call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,lattice=lat, &
        &                    iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return
      if (istep == 1) then
        cap1 = size(calc%nlist%hblist1,dim=2)
        cap2 = size(calc%nlist%hblist2,dim=2)
        cap3 = size(calc%nlist%hblist3,dim=2)
      else
        if (size(calc%nlist%hblist1,dim=2) /= cap1 .or. &
          & size(calc%nlist%hblist2,dim=2) /= cap2 .or. &
          & size(calc%nlist%hblist3,dim=2) /= cap3) then
          call test_failed(error,"PBC HB lists were reallocated although the headroom sufficed")
          return
        end if
      end if
      if (calc%nlist%nhb1 > cap1 .or. calc%nlist%nhb2 > cap2 &
        & .or. calc%nlist%nxb > cap3) then
        call test_failed(error,"PBC HB list count exceeds its allocated size")
        return
      end if
    end do

    call calc%deallocate()
  end subroutine test_list_capacity

  subroutine test_list_reuse(error)
    !***********************************
    !* A calculator walked across small steps must agree with a fresh one at each geometry.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: walk,fresh
    integer :: nat,io,istep
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),gw(:,:),gf(:,:)
    real(wp) :: lat(3,3),ew,ef
    real(wp),parameter :: thr_e = 1.0e-5_wp
    real(wp),parameter :: thr_g = 1.0e-5_wp

    call water_box(3,nat,at,xyz,lat)
    allocate (gw(3,nat),gf(3,nat),source=0.0_wp)

    call gfnff_initialize(nat,at,xyz,walk,lattice=lat,npbc=3,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    do istep = 1,4
      call step_geometry(nat,xyz,0.01_wp)

      call gfnff_singlepoint(nat,at,xyz,walk,ew,gw,lattice=lat, &
        &                    iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return

      call gfnff_initialize(nat,at,xyz,fresh,lattice=lat,npbc=3,iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return
      call gfnff_singlepoint(nat,at,xyz,fresh,ef,gf,lattice=lat, &
        &                    iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return

      call check(error,ew,ef,thr=thr_e)
      if (allocated(error)) return
      if (any(abs(gw-gf) > thr_g)) then
        call test_failed(error,"PBC reused lists give a different gradient than a fresh setup")
        return
      end if

      call fresh%deallocate()
    end do

    call walk%deallocate()
  end subroutine test_list_reuse

  subroutine test_xb_without_h(error)
    !***********************************
    !* Halogen bonds must be found even when no HB-relevant hydrogen skips the H-bond scan.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer,parameter :: nat = 5
    integer :: at(nat),io
    real(wp) :: xyz(3,nat),lat(3,3),energy,grad(3,nat)

    at = [53,53,8,6,8]
    xyz(:,1) = [0.0_wp,0.0_wp,0.0_wp]     !> I
    xyz(:,2) = [0.0_wp,0.0_wp,5.04_wp]    !> I
    xyz(:,3) = [0.0_wp,0.0_wp,10.4_wp]    !> O, the halogen bond acceptor
    xyz(:,4) = [0.0_wp,0.0_wp,12.6_wp]    !> C
    xyz(:,5) = [0.0_wp,0.0_wp,14.8_wp]    !> O
    lat = 0.0_wp
    lat(1,1) = 16.0_wp
    lat(2,2) = 16.0_wp
    lat(3,3) = 22.0_wp

    call gfnff_initialize(nat,at,xyz,calc,lattice=lat,npbc=3,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,lattice=lat, &
      &                    iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,calc%topo%nathbH,0)
    if (allocated(error)) return

    if (calc%nlist%nxb <= 0) then
      call test_failed(error,"PBC halogen bond list is empty without HB hydrogen")
      return
    end if

    call calc%deallocate()
  end subroutine test_xb_without_h

  subroutine test_thread_invariance(error)
    !***********************************
    !* Energy, gradient, and stress must not depend on the OpenMP thread count.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io,nthreads
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),g1(:,:),gn(:,:)
    real(wp) :: lat(3,3),e1,en,s1(3,3),sn(3,3)
    real(wp),parameter :: thr = 1.0e-9_wp

    nthreads = 1
!$  nthreads = omp_get_max_threads()
    if (nthreads < 2) return   !> nothing to compare against

    call water_box(3,nat,at,xyz,lat)
    allocate (g1(3,nat),gn(3,nat),source=0.0_wp)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lat,npbc=3,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

!$  call omp_set_num_threads(1)
    call gfnff_singlepoint(nat,at,xyz,calc,e1,g1,lattice=lat,sigma=s1, &
      &                    iostat=io,printlevel=0)
!$  call omp_set_num_threads(nthreads)
    call check(error,io,0)
    if (allocated(error)) return

    call gfnff_singlepoint(nat,at,xyz,calc,en,gn,lattice=lat,sigma=sn, &
      &                    iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,e1,en,thr=thr)
    if (allocated(error)) return
    if (any(abs(g1-gn) > thr)) then
      call test_failed(error,"PBC gradient depends on the number of threads")
      return
    end if
    if (any(abs(s1-sn) > thr)) then
      call test_failed(error,"PBC stress depends on the number of threads")
      return
    end if

    call calc%deallocate()
  end subroutine test_thread_invariance

  subroutine test_grad_fd(error)
    !***********************************
    !* Central differences on a periodic water box, covering the CN and Ewald matrix gradients.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io,idof,iat,ic
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),grad(:,:),gdum(:,:)
    real(wp) :: lat(3,3),energy,ep,em,gnum
    real(wp),parameter :: h = 1.0e-4_wp
    real(wp),parameter :: thr = 1.0e-6_wp

    call water_box(2,nat,at,xyz,lat)
    allocate (grad(3,nat),gdum(3,nat),source=0.0_wp)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lat,npbc=3,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,lattice=lat, &
      &                    iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    do idof = 1,9
      iat = 1+mod(3*idof,nat)
      ic = 1+mod(idof,3)

      xyz(ic,iat) = xyz(ic,iat)+h
      call gfnff_singlepoint(nat,at,xyz,calc,ep,gdum,lattice=lat, &
        &                    iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return

      xyz(ic,iat) = xyz(ic,iat)-2.0_wp*h
      call gfnff_singlepoint(nat,at,xyz,calc,em,gdum,lattice=lat, &
        &                    iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return

      xyz(ic,iat) = xyz(ic,iat)+h

      gnum = (ep-em)/(2.0_wp*h)
      if (abs(gnum-grad(ic,iat)) > thr) then
        call test_failed(error,"PBC analytical gradient does not match numerical")
        return
      end if
    end do

    call calc%deallocate()
  end subroutine test_grad_fd

  subroutine test_grad_fd_quartz(error)
    !***********************************
    !* Gradient of an alpha-quartz supercell against central differences, at a
    !* much tighter tolerance than the water box test: the nearest Wigner-Seitz
    !* image is pinned to the setup geometry rather than tracked live, since
    !* doing the latter degrades consistency here from 4e-10 to 1.4e-7 (the
    !* threshold below), most visibly for this symmetric crystal.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer,parameter :: nrep = 2
    integer :: nat,io,idof,iat,ic
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),grad(:,:),gdum(:,:)
    real(wp) :: lat(3,3),energy,ep,em,gnum
    real(wp),parameter :: h = 1.0e-4_wp
    real(wp),parameter :: thr = 1.0e-8_wp

    call quartz_super(nrep,nat,at,xyz,lat)
    allocate (grad(3,nat),gdum(3,nat),source=0.0_wp)

    !> off the ideal crystal, so the image ties are generic
    call step_geometry(nat,xyz,0.05_wp)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lat,npbc=3,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,lattice=lat, &
      &                    iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    do idof = 1,12
      iat = 1+mod(5*idof,nat)
      ic = 1+mod(idof,3)

      xyz(ic,iat) = xyz(ic,iat)+h
      call gfnff_singlepoint(nat,at,xyz,calc,ep,gdum,lattice=lat, &
        &                    iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return

      xyz(ic,iat) = xyz(ic,iat)-2.0_wp*h
      call gfnff_singlepoint(nat,at,xyz,calc,em,gdum,lattice=lat, &
        &                    iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return

      xyz(ic,iat) = xyz(ic,iat)+h

      gnum = (ep-em)/(2.0_wp*h)
      if (abs(gnum-grad(ic,iat)) > thr) then
        call test_failed(error, &
          & "PBC gradient on a symmetric cell is not consistent with the energy")
        return
      end if
    end do

    call calc%deallocate()
  end subroutine test_grad_fd_quartz

  subroutine test_stress_fd(error)
    !***********************************
    !* Diagonal strain derivatives vs. central differences, covering CN/Ewald strain with hydrogen bonds.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc,calc_p,calc_m
    integer :: nat,io,ii
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),xyz_p(:,:),xyz_m(:,:),grad(:,:)
    real(wp) :: lat(3,3),lat_p(3,3),lat_m(3,3)
    real(wp) :: energy,ep,em,sigma(3,3),snum
    real(wp),parameter :: h = 1.0e-4_wp
    real(wp),parameter :: thr = 5.0e-5_wp

    call water_box(2,nat,at,xyz,lat)
    allocate (grad(3,nat),xyz_p(3,nat),xyz_m(3,nat),source=0.0_wp)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lat,npbc=3,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,lattice=lat,sigma=sigma, &
      &                    iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    do ii = 1,3
      lat_p = lat
      lat_p(ii,:) = lat(ii,:)+h*lat(ii,:)
      xyz_p = xyz
      xyz_p(ii,:) = xyz(ii,:)+h*xyz(ii,:)

      call gfnff_initialize(nat,at,xyz_p,calc_p,lattice=lat_p,npbc=3, &
        &                   iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return
      call gfnff_singlepoint(nat,at,xyz_p,calc_p,ep,grad,lattice=lat_p, &
        &                    iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return
      call calc_p%deallocate()

      lat_m = lat
      lat_m(ii,:) = lat(ii,:)-h*lat(ii,:)
      xyz_m = xyz
      xyz_m(ii,:) = xyz(ii,:)-h*xyz(ii,:)

      call gfnff_initialize(nat,at,xyz_m,calc_m,lattice=lat_m,npbc=3, &
        &                   iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return
      call gfnff_singlepoint(nat,at,xyz_m,calc_m,em,grad,lattice=lat_m, &
        &                    iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return
      call calc_m%deallocate()

      snum = (ep-em)/(2.0_wp*h)
      if (abs(sigma(ii,ii)-snum) > thr) then
        call test_failed(error,"PBC analytical stress does not match numerical")
        return
      end if
    end do

    call calc%deallocate()
  end subroutine test_stress_fd

  subroutine test_stress_reused(error)
    !***********************************
    !* Same strain check, but through one calculator handed a new lattice
    !* each step. That is the variable cell path a host like ASE takes, and
    !* the only test reaching the lattice-changed branch of
    !* gfnff_singlepoint, where the Wigner-Seitz images are refreshed.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io,ii
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),xyz_p(:,:),xyz_m(:,:),grad(:,:)
    real(wp) :: lat(3,3),lat_p(3,3),lat_m(3,3)
    real(wp) :: energy,ep,em,sigma(3,3),snum
    real(wp),parameter :: h = 1.0e-4_wp
    real(wp),parameter :: thr = 5.0e-5_wp

    call water_box(2,nat,at,xyz,lat)
    allocate (grad(3,nat),xyz_p(3,nat),xyz_m(3,nat),source=0.0_wp)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lat,npbc=3,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,lattice=lat,sigma=sigma, &
      &                    iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    do ii = 1,3
      lat_p = lat
      lat_p(ii,:) = lat(ii,:)+h*lat(ii,:)
      xyz_p = xyz
      xyz_p(ii,:) = xyz(ii,:)+h*xyz(ii,:)
      call gfnff_singlepoint(nat,at,xyz_p,calc,ep,grad,lattice=lat_p, &
        &                    iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return

      lat_m = lat
      lat_m(ii,:) = lat(ii,:)-h*lat(ii,:)
      xyz_m = xyz
      xyz_m(ii,:) = xyz(ii,:)-h*xyz(ii,:)
      call gfnff_singlepoint(nat,at,xyz_m,calc,em,grad,lattice=lat_m, &
        &                    iostat=io,printlevel=0)
      call check(error,io,0)
      if (allocated(error)) return

      snum = (ep-em)/(2.0_wp*h)
      if (abs(sigma(ii,ii)-snum) > thr) then
        call test_failed(error, &
          & "PBC stress on a reused calculator does not match numerical")
        return
      end if
    end do

    !> and the original lattice must reproduce the original energy
    call gfnff_singlepoint(nat,at,xyz,calc,ep,grad,lattice=lat, &
      &                    iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call check(error,ep,energy,thr=1.0e-10_wp)

    call calc%deallocate()
  end subroutine test_stress_reused

  subroutine test_bpair_values(error)
    !***********************************
    !* The bond pair matrix (shortest-path bond count, saturated at 5) is
    !* stored one byte wide. A too narrow type would wrap silently, and every
    !* consumer only compares against 1, 2 or 3, so the values are pinned
    !* here: the range, and three entries of a water box.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    real(wp) :: lat(3,3),energy

    call water_box(2,nat,at,xyz,lat)
    allocate (grad(3,nat),source=0.0_wp)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lat,npbc=3,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,lattice=lat, &
      &                    iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    !> the whole range the matrix is supposed to take
    call check(error,int(minval(calc%neigh%bpair)),0)
    if (allocated(error)) return
    call check(error,int(maxval(calc%neigh%bpair)),5)
    if (allocated(error)) return

    !> water_box lays out each molecule as O, H, H, so the first molecule's
    !> bond counts are known
    call check(error,int(calc%neigh%bpair(2,1,1)),1)   !> O-H, bonded
    if (allocated(error)) return
    call check(error,int(calc%neigh%bpair(3,2,1)),2)   !> H-H, two bonds apart
    if (allocated(error)) return
    !> and the oxygens of two different molecules are beyond the count
    call check(error,int(calc%neigh%bpair(4,1,1)),5)
    if (allocated(error)) return

    call calc%deallocate()
  end subroutine test_bpair_values

  subroutine test_frag_charge(error)
    !***********************************
    !* Water and NH4+ in a box, the +1 pinned to the second fragment via
    !* reference charges. The Ewald matrix arrives with a total-charge row in
    !* the first constraint slot; left in place it puts -1 on the water and
    !* makes the cell neutral. The energy must also not depend on which
    !* fragment comes first.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc,calc2
    integer,parameter :: nat = 8
    integer,parameter :: swap(nat) = [4,5,6,7,8,1,2,3]
    integer :: at(nat),io,i
    real(wp) :: xyz(3,nat),lat(3,3),e1,e2,grad(3,nat),refq(nat),qf(2)
    real(wp),parameter :: aa = 1.0_wp/0.52917721067_wp

    at = [8,1,1,7,1,1,1,1]
    xyz(:,1) = [0.00_wp,0.00_wp,0.00_wp]
    xyz(:,2) = [0.76_wp,0.59_wp,0.00_wp]
    xyz(:,3) = [-0.76_wp,0.59_wp,0.00_wp]
    xyz(:,4) = [4.00_wp,4.00_wp,4.00_wp]
    xyz(:,5) = [4.59_wp,4.59_wp,4.59_wp]
    xyz(:,6) = [3.41_wp,3.41_wp,4.59_wp]
    xyz(:,7) = [3.41_wp,4.59_wp,3.41_wp]
    xyz(:,8) = [4.59_wp,3.41_wp,3.41_wp]
    xyz = xyz*aa
    refq = [-0.6_wp,0.3_wp,0.3_wp,-0.6_wp,0.4_wp,0.4_wp,0.4_wp,0.4_wp]
    lat = 0.0_wp
    do i = 1,3
      lat(i,i) = 8.0_wp*aa
    end do

    allocate (calc%userinput)
    calc%userinput%refq = refq
    call gfnff_initialize(nat,at,xyz,calc,ichrg=1,lattice=lat,npbc=3,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,e1,grad,lattice=lat,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    qf = 0.0_wp
    do i = 1,nat
      qf(calc%topo%fraglist(i)) = qf(calc%topo%fraglist(i))+calc%nlist%q(i)
    end do
    call check(error,qf(1),0.0_wp,thr=1.0e-8_wp)
    if (allocated(error)) return
    call check(error,qf(2),1.0_wp,thr=1.0e-8_wp)
    if (allocated(error)) return
    call calc%deallocate()

    !>-- NH4+ first: the charged fragment is now the first one
    allocate (calc2%userinput)
    calc2%userinput%refq = refq(swap)
    call gfnff_initialize(nat,at(swap),xyz(:,swap),calc2,ichrg=1,lattice=lat,npbc=3, &
      &                   iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at(swap),xyz(:,swap),calc2,e2,grad,lattice=lat, &
      &                    iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call check(error,e1,e2,thr=1.0e-10_wp)
    if (allocated(error)) return
    call calc2%deallocate()
  end subroutine test_frag_charge

  subroutine test_hb_linear_n(error)
    !***********************************
    !* Water donating to the central N of an exactly linear, symmetric azide.
    !* The two bond vectors at N cancel, so the lone pair of abhgfnff_eg2_rnr
    !* has no direction; the gradient used to come out as NaN. It must be
    !* finite and, for the water atoms whose displacement keeps N3 linear,
    !* match finite differences. Molecular, but kept next to the other
    !* kernel guards.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer,parameter :: nat = 6
    integer :: at(nat),io,i,ic
    real(wp) :: xyz(3,nat),x0(3,nat),grad(3,nat),gdum(3,nat),e,ep,em
    real(wp),allocatable :: hess(:,:)
    real(wp),parameter :: aa = 1.0_wp/0.52917721067_wp,h = 1.0e-4_wp

    at = [7,7,7,8,1,1]
    xyz(:,1) = [-1.18_wp,0.0_wp,0.0_wp]
    xyz(:,2) = [0.0_wp,0.0_wp,0.0_wp]
    xyz(:,3) = [1.18_wp,0.0_wp,0.0_wp]
    xyz(:,4) = [0.10_wp,2.95_wp,0.05_wp]
    xyz(:,5) = [0.05_wp,1.98_wp,0.02_wp]
    xyz(:,6) = [0.95_wp,3.25_wp,0.30_wp]
    xyz = xyz*aa
    x0 = xyz

    call gfnff_initialize(nat,at,xyz,calc,ichrg=-1,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,e,grad,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    !>-- the premise: the central N is an acceptor on the unbound list
    if (.not.any(calc%nlist%hblist2(2,1:calc%nlist%nhb2) == 2)) then
      call test_failed(error,"central azide N is not a hydrogen bond acceptor")
      return
    end if
    if (any(grad /= grad)) then
      call test_failed(error,"NaN in the gradient for a linear two-neighbour acceptor")
      return
    end if

    do i = 4,nat
      do ic = 1,3
        xyz = x0
        xyz(ic,i) = x0(ic,i)+h
        call gfnff_singlepoint(nat,at,xyz,calc,ep,gdum,iostat=io,printlevel=0)
        xyz(ic,i) = x0(ic,i)-h
        call gfnff_singlepoint(nat,at,xyz,calc,em,gdum,iostat=io,printlevel=0)
        call check(error,grad(ic,i),(ep-em)/(2.0_wp*h),thr=1.0e-7_wp)
        if (allocated(error)) return
      end do
    end do

    allocate (hess(3*nat,3*nat))
    call calc%hessian(nat,at,x0,hess,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    if (any(hess /= hess)) then
      call test_failed(error,"NaN in the Hessian for a linear two-neighbour acceptor")
      return
    end if
    call calc%deallocate()
  end subroutine test_hb_linear_n

  subroutine test_hb_list_growth(error)
    !***********************************
    !* A water pushed from 30 bohr outside a 226 atom cluster to contact, on one
    !* reused calculator. With this many atoms the displacement rule of
    !* gfnff_hbset never fires, so only the candidate-count check in
    !* setup_hb_lists can bring the new hydrogen bonds into the lists; without
    !* it the energy is off by more than 1 mEh at contact. The reference has
    !* its lists forced fresh at every step. Molecular, kept with the other
    !* list guards.
    !***********************************
    use supermol,only:snat => testnat,sat => testat,sxyz => testxyz
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: walk,ref
    integer,parameter :: nat = snat+3
    integer :: at(nat),io,istep
    real(wp) :: xyz(3,nat),gw(3,nat),gr(3,nat),ew,er,dmin

    at(1:snat) = sat
    at(snat+1:) = [8,1,1]
    xyz(:,1:snat) = sxyz
    xyz(:,snat+1) = [maxval(sxyz(1,:))+30.0_wp,sum(sxyz(2,:))/snat,sum(sxyz(3,:))/snat]
    xyz(:,snat+2) = xyz(:,snat+1)+[-1.81_wp,0.0_wp,0.0_wp]
    xyz(:,snat+3) = xyz(:,snat+1)+[0.45_wp,1.75_wp,0.0_wp]

    call gfnff_initialize(nat,at,xyz,walk,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_initialize(nat,at,xyz,ref,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return

    do istep = 0,40
      if (istep > 0) xyz(1,snat+1:) = xyz(1,snat+1:)-1.0_wp
      dmin = sqrt(minval(sum((xyz(:,1:snat)-spread(xyz(:,snat+2),2,snat))**2,dim=1)))
      if (dmin < 3.5_wp) exit
      call walk%singlepoint(nat,at,xyz,ew,gw,printlevel=0)
      ref%nlist%hbrefgeo = xyz
      call ref%singlepoint(nat,at,xyz,er,gr,printlevel=0)
      call check(error,ew,er,thr=1.0e-8_wp)
      if (allocated(error)) return
    end do

    !>-- the premise: the walk went far enough for the miss to matter
    if (dmin > 6.0_wp) then
      call test_failed(error,"incoming water never reached the cluster")
      return
    end if
    call walk%deallocate()
    call ref%deallocate()
  end subroutine test_hb_list_growth

end module test_pbc_kernels
