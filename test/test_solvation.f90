module test_solvation
!> Unit tests for the ALPB solvation path: thread invariance, list-builder
!> agreement, and gradient checked against finite differences.
  use testdrive,only:new_unittest,unittest_type,error_type,check
  use iso_fortran_env,only:wp => real64
  use gfnff_interface
  use gfnff_solv_gbsa,only:TBorn,update_nnlist_gbsa
  use gfnff_hess_solv,only:born_weighted_hessian,sasa_weighted_hessian
  use gfnff_hess_driver,only:hess_analytic_terms,hess_term_reference, &
    &                           hess_available_terms
  use gfnff_hess_analysis,only:hess_transl_sumrule,hess_asymmetry, &
    &                          hess_symmetrize
  use gfnff_eg_driver,only:gff_term_es
  use gfnff_math_wrapper,only:gemm
  implicit none
  private

  public :: collect_solvation

contains

  subroutine collect_solvation(testsuite)
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
    new_unittest("Solvation thread invariance    ",test_solv_threads), &
    new_unittest("Solvation neighbour list paths ",test_solv_nnlist), &
    new_unittest("Solvation gradient O(h^2) scan ",test_solv_grad_scan), &
    new_unittest("Solvation energy is bracketed  ",test_solv_energy), &
    new_unittest("Born radius Hessian O(h^2) scan",test_solv_hess_born), &
    new_unittest("SASA Hessian O(h^2) scan       ",test_solv_hess_sasa), &
    new_unittest("Solvated ES Hessian vs FD      ",test_solv_hess_full), &
    new_unittest("Solvated ES Hessian, anion     ",test_solv_hess_ion) &
    ]
!&>
  end subroutine collect_solvation

  subroutine watercluster(nat,at,xyz)
    !***********************************************************************
    !* Four irregular waters; no two atoms are symmetry-equivalent.
    !***********************************************************************
    integer,intent(out) :: nat
    integer,allocatable,intent(out) :: at(:)
    real(wp),allocatable,intent(out) :: xyz(:,:)
    real(wp),parameter :: aatoau = 1.0_wp/0.52917726_wp
    real(wp),parameter :: raw(3,12) = reshape([ &
       &  0.000_wp,0.000_wp,0.000_wp, &
       &  0.958_wp,0.000_wp,0.000_wp, &
       & -0.239_wp,0.927_wp,0.000_wp, &
       &  2.780_wp,0.180_wp,0.310_wp, &
       &  3.120_wp,-0.700_wp,0.060_wp, &
       &  3.410_wp,0.830_wp,-0.040_wp, &
       &  0.410_wp,2.610_wp,1.180_wp, &
       &  0.170_wp,3.400_wp,0.690_wp, &
       &  1.310_wp,2.760_wp,1.470_wp, &
       & -1.120_wp,-1.430_wp,2.240_wp, &
       & -0.640_wp,-0.940_wp,2.920_wp, &
       & -1.980_wp,-1.020_wp,2.360_wp],[3,12])

    nat = 12
    allocate (at(nat),xyz(3,nat))
    at = [8,1,1,8,1,1,8,1,1,8,1,1]
    xyz = raw*aatoau
  end subroutine watercluster

  subroutine test_solv_threads(error)
    !***********************************************************************
    !* Bit-identical at any thread count; parallelised without reordering sums.
    !***********************************************************************
    !$ use omp_lib,only:omp_get_max_threads,omp_set_num_threads
    use coffeine
    type(error_type),allocatable,intent(out) :: error

    integer :: nat,io,k,nthr(3),it,keep
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    real(wp),allocatable :: brad0(:),sasa0(:),bmat0(:,:),brdr0(:,:,:)
    real(wp) :: energy
    type(gfnff_data) :: calc

    nat = testnat
    allocate (at(nat),xyz(3,nat),grad(3,nat))
    at = testat
    xyz = testxyz

    call calc%init(nat,at,xyz,ichrg=0,solvent='water',iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call check(error,allocated(calc%solvation))
    if (allocated(error)) return

    keep = 1
    !$ keep = omp_get_max_threads()
    nthr = [1,2,4]

    do it = 1,size(nthr)
      !$ call omp_set_num_threads(min(nthr(it),keep))
      call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,printlevel=0)
      if (it .eq. 1) then
        brad0 = calc%solvation%brad
        sasa0 = calc%solvation%sasa
        bmat0 = calc%solvation%bornMat
        brdr0 = calc%solvation%brdr
      else
        do k = 1,nat
          call check(error,calc%solvation%brad(k),brad0(k),thr=0.0_wp)
          if (allocated(error)) exit
          call check(error,calc%solvation%sasa(k),sasa0(k),thr=0.0_wp)
          if (allocated(error)) exit
        end do
        if (allocated(error)) exit
        call check(error,maxval(abs(calc%solvation%bornMat-bmat0)),0.0_wp, &
           & thr=0.0_wp)
        if (allocated(error)) exit
        call check(error,maxval(abs(calc%solvation%brdr-brdr0)),0.0_wp, &
           & thr=0.0_wp)
        if (allocated(error)) exit
      end if
    end do

    !$ call omp_set_num_threads(keep)

  end subroutine test_solv_threads

  subroutine test_solv_nnlist(error)
    !***********************************************************************
    !* Sequential and parallel builders must agree entry by entry, in order.
    !***********************************************************************
    use coffeine
    type(error_type),allocatable,intent(out) :: error

    integer :: nat,io,i,j,nrad_s,nrad_p
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: nnsas_s(:),nnlists_s(:,:),nnlistr_s(:,:)
    real(wp),allocatable :: ddpair_s(:,:)
    real(wp) :: energy
    type(gfnff_data) :: calc

    nat = testnat
    allocate (at(nat),xyz(3,nat),grad(3,nat))
    at = testat
    xyz = testxyz

    call calc%init(nat,at,xyz,ichrg=0,solvent='water',iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,printlevel=0)

    associate (g => calc%solvation)
      allocate (nnsas_s(g%nat),nnlists_s(size(g%nnlists,1),g%nat))
      allocate (nnlistr_s(3,g%ntpair),ddpair_s(4,g%ntpair))

      call update_nnlist_gbsa(g%nat,g%ntpair,g%ppind,xyz,g%lrcut,g%srcut, &
         & nnsas_s,nnlists_s,nrad_s,nnlistr_s,ddpair_s,.false.)
      call update_nnlist_gbsa(g%nat,g%ntpair,g%ppind,xyz,g%lrcut,g%srcut, &
         & g%nnsas,g%nnlists,nrad_p,g%nnlistr,g%ddpair,.true.)

      call check(error,nrad_p,nrad_s)
      if (allocated(error)) return
      call check(error,nrad_s > 0)
      if (allocated(error)) return

      do i = 1,g%nat
        call check(error,g%nnsas(i),nnsas_s(i))
        if (allocated(error)) return
      end do

      do i = 1,g%nat
        do j = 1,nnsas_s(i)
          call check(error,g%nnlists(j,i),nnlists_s(j,i))
          if (allocated(error)) return
        end do
      end do

      do i = 1,nrad_s
        do j = 1,3
          call check(error,g%nnlistr(j,i),nnlistr_s(j,i))
          if (allocated(error)) return
        end do
      end do

      call check(error,maxval(abs(g%ddpair-ddpair_s)),0.0_wp,thr=0.0_wp)
    end associate

  end subroutine test_solv_nnlist

  subroutine test_solv_grad_scan(error)
    !***********************************************************************
    !* Halving h should quarter the deviation (O(h^2) truncation of the FD).
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error

    integer,parameter :: nscan = 4
    integer :: nat,io,i,k,s
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),xp(:,:),g0(:,:),gtmp(:,:)
    real(wp) :: energy,ep,em,h,dev(0:nscan),ratio
    type(gfnff_data) :: calc

    call watercluster(nat,at,xyz)
    allocate (xp(3,nat),g0(3,nat),gtmp(3,nat))

    call calc%init(nat,at,xyz,ichrg=0,solvent='water',iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,energy,g0,printlevel=0)

    do s = 0,nscan
      h = 5.0e-3_wp/2.0_wp**s
      dev(s) = 0.0_wp
      do i = 1,nat
        do k = 1,3
          xp = xyz
          xp(k,i) = xp(k,i)+h
          call gfnff_singlepoint(nat,at,xp,calc,ep,gtmp,printlevel=0)
          xp = xyz
          xp(k,i) = xp(k,i)-h
          call gfnff_singlepoint(nat,at,xp,calc,em,gtmp,printlevel=0)
          dev(s) = max(dev(s),abs((ep-em)/(2.0_wp*h)-g0(k,i)))
        end do
      end do
    end do

    do s = 1,nscan
      ratio = dev(s-1)/max(dev(s),1.0e-30_wp)
      call check(error,ratio > 3.5_wp.and.ratio < 4.5_wp)
      if (allocated(error)) return
    end do

  end subroutine test_solv_grad_scan

  subroutine test_solv_energy(error)
    !***********************************************************************
    !* Solvation shifts the energy by more than g_solv alone via EEQ/Born charge repolarisation.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error

    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    real(wp) :: evac,esol,dg
    type(gfnff_data) :: vac,sol
    real(wp),parameter :: autokcal = 627.50947428_wp

    call watercluster(nat,at,xyz)
    allocate (grad(3,nat))

    call vac%init(nat,at,xyz,ichrg=0,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,vac,evac,grad,printlevel=0)

    call sol%init(nat,at,xyz,ichrg=0,solvent='water',iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,sol,esol,grad,printlevel=0)

    dg = (esol-evac)*autokcal
    call check(error,dg < 0.0_wp)
    if (allocated(error)) return
    call check(error,dg > -200.0_wp)  !> four waters in water: tens of kcal/mol, stabilising
    if (allocated(error)) return

    call check(error,sol%res%g_solv, &
       & sol%res%g_born+sol%res%g_sasa+sol%res%g_hb+sol%res%g_shift, &
       & thr=1.0e-10_wp)
    if (allocated(error)) return

    call check(error,maxval(abs(sol%nlist%q-vac%nlist%q)) > 1.0e-3_wp)
    if (allocated(error)) return

    call check(error,esol-evac, &
       & sol%res%g_solv+(sol%res%e_es-vac%res%e_es),thr=1.0e-8_wp)

  end subroutine test_solv_energy

  subroutine test_solv_hess_born(error)
    !***********************************************************************
    !* Born-radii Hessian vs. FD of brdr; weights avoid a dE/db cancellation.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error

    integer,parameter :: nscan = 3
    integer :: nat,io,i,k,s,ndof
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),xp(:,:),grad(:,:),w(:)
    real(wp),allocatable :: ha(:,:),hn(:,:),kmat(:,:),tmp(:,:),dbdr(:,:)
    real(wp),allocatable :: gp(:),gm(:)
    real(wp) :: energy,h,dev(0:nscan),ratio
    type(gfnff_data) :: calc

    call watercluster(nat,at,xyz)
    ndof = 3*nat
    allocate (xp(3,nat),grad(3,nat),w(nat))
    allocate (ha(ndof,ndof),hn(ndof,ndof),kmat(nat,nat),tmp(ndof,nat))
    allocate (dbdr(ndof,nat),gp(ndof),gm(ndof))

    call calc%init(nat,at,xyz,ichrg=0,solvent='water',iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,printlevel=0)

    do i = 1,nat
      w(i) = sin(1.7_wp*i)+0.3_wp*cos(0.9_wp*i)
    end do

    call calc%solvation%update(at,xyz)
    ha = 0.0_wp
    kmat = 0.0_wp
    call born_weighted_hessian(nat,xyz,calc%solvation,w,kmat,ha)
    !> contract the rank-one half the routine leaves on kmat
    dbdr = reshape(calc%solvation%brdr,[ndof,nat])
    call gemm(dbdr,kmat,tmp)
    call gemm(tmp,dbdr,ha,transb='T',alpha=1.0_wp,beta=1.0_wp)

    do s = 0,nscan
      h = 4.0e-3_wp/2.0_wp**s
      do i = 1,nat
        do k = 1,3
          xp = xyz
          xp(k,i) = xp(k,i)+h
          call calc%solvation%update(at,xp)
          call wgrad(nat,calc%solvation%brdr,w,gp)
          xp = xyz
          xp(k,i) = xp(k,i)-h
          call calc%solvation%update(at,xp)
          call wgrad(nat,calc%solvation%brdr,w,gm)
          hn(:,3*(i-1)+k) = (gp-gm)/(2.0_wp*h)
        end do
      end do
      hn = 0.5_wp*(hn+transpose(hn))
      dev(s) = maxval(abs(ha-hn))
    end do

    do s = 1,nscan
      ratio = dev(s-1)/max(dev(s),1.0e-30_wp)
      call check(error,ratio > 3.5_wp.and.ratio < 4.5_wp)
      if (allocated(error)) return
    end do

  end subroutine test_solv_hess_born

  subroutine test_solv_hess_sasa(error)
    !***********************************************************************
    !* Two atoms only: C1-but-not-C2 cubic breaks FD convergence past that.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error

    integer,parameter :: nat = 2,nscan = 3
    integer :: io,i,k,s
    integer :: at(nat)
    real(wp) :: xyz(3,nat),xp(3,nat),grad(3,nat),w(nat)
    real(wp) :: ha(6,6),hn(6,6),gp(6),gm(6)
    real(wp) :: energy,h,dev(0:nscan),ratio
    type(gfnff_data) :: calc

    at = [8,8]
    w = [1.0_wp,0.7_wp]
    xyz = 0.0_wp
    xyz(1,2) = 5.35_wp

    call calc%init(nat,at,xyz,ichrg=0,solvent='water',iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,printlevel=0)

    call calc%solvation%update(at,xyz)
    ha = 0.0_wp
    call sasa_weighted_hessian(nat,xyz,calc%solvation,w,ha)

    !> guard against a zero Hessian passing the O(h^2) scan trivially
    call check(error,maxval(abs(ha)) > 1.0e-3_wp)
    if (allocated(error)) return
    call check(error,maxval(abs(ha-transpose(ha))) < 1.0e-12_wp)
    if (allocated(error)) return

    do s = 0,nscan
      h = 4.0e-3_wp/2.0_wp**s
      do i = 1,nat
        do k = 1,3
          xp = xyz
          xp(k,i) = xp(k,i)+h
          call calc%solvation%update(at,xp)
          call wgrad(nat,calc%solvation%dsdrt,w,gp)
          xp = xyz
          xp(k,i) = xp(k,i)-h
          call calc%solvation%update(at,xp)
          call wgrad(nat,calc%solvation%dsdrt,w,gm)
          hn(:,3*(i-1)+k) = (gp-gm)/(2.0_wp*h)
        end do
      end do
      hn = 0.5_wp*(hn+transpose(hn))
      dev(s) = maxval(abs(ha-hn))
    end do

    do s = 1,nscan
      ratio = dev(s-1)/max(dev(s),1.0e-30_wp)
      call check(error,ratio > 3.5_wp.and.ratio < 4.5_wp)
      if (allocated(error)) return
    end do

  end subroutine test_solv_hess_sasa

  subroutine test_solv_hess_full(error)
    !***********************************************************************
    !* Solvated ES Hessian on the neutral water cluster; see hess_full_check.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error
    call hess_full_check(0,error)
  end subroutine test_solv_hess_full

  subroutine test_solv_hess_ion(error)
    !***********************************************************************
    !* Charged solute exercises the ALPB correction (scales with charge^2).
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error
    call hess_full_check(-1,error)
  end subroutine test_solv_hess_ion

  subroutine hess_full_check(ichrg,error)
    !***********************************************************************
    !* Full solvated ES Hessian vs. FD where step halvings agree to 1e-9.
    !***********************************************************************
    integer,intent(in) :: ichrg
    type(error_type),allocatable,intent(out) :: error

    integer :: nat,io,ndof,a,b,nconv
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    real(wp),allocatable :: ha(:,:),h1(:,:),h2(:,:),h3(:,:)
    real(wp) :: energy,efield(3),d1,d2,scale,worst
    type(gfnff_data) :: calc

    call watercluster(nat,at,xyz)
    ndof = 3*nat
    allocate (grad(3,nat),ha(ndof,ndof),h1(ndof,ndof))
    allocate (h2(ndof,ndof),h3(ndof,ndof))
    efield = 0.0_wp

    call calc%init(nat,at,xyz,ichrg=ichrg,solvent='water',iostat=io, &
       & printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,printlevel=0)

    call check(error,iand(hess_available_terms(calc%cell,calc%neigh, &
       & calc%topo,calc%solvation,efield),gff_term_es) .ne. 0)
    if (allocated(error)) return

    ha = 0.0_wp
    call hess_analytic_terms(nat,at,xyz,calc%accuracy,calc%version, &
       & gff_term_es,calc%param,calc%topo,calc%neigh,calc%nlist,ha, &
       & calc%solvation)
    call hess_symmetrize(ha)

    call fdref(1.0e-2_wp,h1)
    call fdref(5.0e-3_wp,h2)
    call fdref(2.5e-3_wp,h3)

    scale = maxval(abs(ha))
    call check(error,scale > 1.0e-4_wp)
    if (allocated(error)) return

    call check(error,hess_transl_sumrule(ha) < 1.0e-10_wp)
    if (allocated(error)) return
    call check(error,hess_asymmetry(ha) < 1.0e-12_wp)
    if (allocated(error)) return

    nconv = 0
    worst = 0.0_wp
    do b = 1,ndof
      do a = 1,ndof
        d1 = abs(h1(a,b)-h2(a,b))
        d2 = abs(h2(a,b)-h3(a,b))
        if (d1 < 1.0e-6_wp*scale.and.d2 < 1.0e-9_wp*scale) then
          nconv = nconv+1
          worst = max(worst,abs(ha(a,b)-h3(a,b)))
        end if
      end do
    end do

    call check(error,nconv > ndof*ndof/10)  !> guards against a vacuous pass
    if (allocated(error)) return
    call check(error,worst < 1.0e-7_wp*max(scale,1.0_wp))

  contains

    subroutine fdref(h,hh)
      real(wp),intent(in) :: h
      real(wp),intent(out) :: hh(ndof,ndof)
      call hess_term_reference(0,nat,at,xyz,calc%cell,calc%ichrg,calc%param, &
         & calc%topo,calc%neigh,calc%nlist,efield,calc%solvation, &
         & calc%version,calc%accuracy,gff_term_es,hh,step=h)
    end subroutine fdref

  end subroutine hess_full_check

  subroutine wgrad(nat,d,w,gv)
    !***********************************************************************
    !* Contracts a (3,nat,nat) array with weights into a flat 3nat gradient.
    !***********************************************************************
    integer,intent(in) :: nat
    real(wp),intent(in) :: d(3,nat,nat)
    real(wp),intent(in) :: w(nat)
    real(wp),intent(out) :: gv(3*nat)
    integer :: i
    gv = 0.0_wp
    do i = 1,nat
      gv = gv+w(i)*reshape(d(:,:,i),[3*nat])
    end do
  end subroutine wgrad

end module test_solvation
