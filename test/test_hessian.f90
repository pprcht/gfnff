module test_hessian
!> Unit tests for the GFN-FF Cartesian Hessian.
!> Two checks matter: the translational sum rule and raw asymmetry hold for
!> an exact Hessian regardless of the reference, catching a term that forgot
!> a partner atom or scattered inconsistently. Per-term scans difference one
!> closed-form term, isolated via gfnff_eg's gff_term_* mask, against finite
!> differences at shrinking steps; halving the step must quarter the
!> deviation for an exact analytic Hessian, while a wrong one leaves a
!> residual that stops shrinking.
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_interface
  use gfnff_hess_driver
  use gfnff_hess_analysis
  use gfnff_hess_rep
  use gfnff_hess_bonded
  use gfnff_hess_es
  use gfnff_hess_disp
  use gfnff_hess_hb
  use gfnff_eg_driver,only:gff_term_rep,gff_term_angl,gff_term_batm, &
    &                          gff_term_bond,gff_term_tors,gff_term_es, &
    &                          gff_term_disp,gff_term_hb,gff_term_xb, &
    &                          gff_term_all
  use gfnff_param,only:gfnff_thresholds
  use gfnff_cn,only:gfnff_dlogcoord
  use gfnff_eg_hb,only:dncoord_erf
  implicit none
  private

  public :: collect_hessian

contains  !> Unit tests for the analytic Hessian

!> Collect all exported unit tests
  subroutine collect_hessian(testsuite)
    !***********************************
    !* Collection of tests
    !***********************************
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
    new_unittest("Hessian repulsion vs FD        ",test_hess_rep_fd), &
    new_unittest("Hessian repulsion O(h^2) scan  ",test_hess_rep_scan), &
    new_unittest("Hessian angles O(h^2) scan     ",test_hess_angl_scan), &
    new_unittest("Hessian ATM O(h^2) scan        ",test_hess_batm_scan), &
    new_unittest("Hessian bonds O(h^2) scan      ",test_hess_bond_scan), &
    new_unittest("Hessian torsions O(h^2) scan   ",test_hess_tors_scan), &
    new_unittest("Hessian EEQ electrostatics scan",test_hess_es_scan), &
    new_unittest("Hessian D3 dispersion scan     ",test_hess_disp_scan), &
    new_unittest("Hessian hydrogen bonds scan    ",test_hess_hb_scan), &
    new_unittest("Hessian HB acceptor forms      ",test_hess_hb_forms), &
    new_unittest("Hessian halogen bonds          ",test_hess_xb), &
    new_unittest("Hessian EEQ two-fragment border",test_hess_es_frag), &
    new_unittest("Hessian special torsions       ",test_hess_stors), &
    new_unittest("Hessian HB-corrected bonds     ",test_hess_hbbond), &
    new_unittest("Hessian sum rule and symmetry  ",test_hess_sumrule), &
    new_unittest("Hessian full vs pure FD        ",test_hess_full), &
    new_unittest("Hessian translational modes    ",test_hess_freq), &
    new_unittest("Hessian conformer2020 bonds    ",test_hess_bond_conformer) &
    ]
!&>
  end subroutine collect_hessian

  subroutine setup(nat,at,xyz,calc,error)
    !***********************************
    !* Shared caffeine setup: geometry plus an initialised calculator.
    !***********************************
    use coffeine
    integer,intent(out) :: nat
    integer,allocatable,intent(out) :: at(:)
    real(wp),allocatable,intent(out) :: xyz(:,:)
    type(gfnff_data),intent(out) :: calc
    type(error_type),allocatable,intent(out) :: error
    integer :: io

    nat = testnat
    allocate (at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz
    call gfnff_initialize(nat,at,xyz,calc,ichrg=0,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call prime(nat,at,xyz,calc,error)
  end subroutine setup

  subroutine prime(nat,at,xyz,calc,error)
    !***********************************
    !* One full single point, so nlist%q holds the EEQ charges of this
    !* geometry; gfnff_initialize alone does not produce them.
    !***********************************
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(gfnff_data),intent(inout) :: calc
    type(error_type),allocatable,intent(out) :: error
    integer :: io
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)

    allocate (grad(3,nat))
    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,printlevel=0,iostat=io)
    call check(error,io,0)
  end subroutine prime

  subroutine analytic_terms(nat,at,xyz,calc,terms,hess)
    !***********************************
    !* Closed-form Hessian of the selected terms; same kernels the driver dispatches to.
    !***********************************
    integer,intent(in) :: nat,at(nat),terms
    real(wp),intent(in) :: xyz(3,nat)
    type(gfnff_data),intent(inout) :: calc
    real(wp),intent(out) :: hess(3*nat,3*nat)

    integer :: i,j,k,ncov
    real(wp) :: dispthr,cnthr,repthr,hbthr1,hbthr2
    real(wp),allocatable :: sqrab(:),srab(:),cn(:),dcn(:,:,:)
    real(wp),allocatable :: hb_cn(:),hb_dcn(:,:,:),dhbcndL(:,:,:)

    call gfnff_thresholds(calc%accuracy,dispthr,cnthr,repthr,hbthr1,hbthr2)
    allocate (sqrab(nat*(nat+1)/2),srab(nat*(nat+1)/2))
    do i = 1,nat
      k = i*(i-1)/2
      do j = 1,i
        sqrab(k+j) = sum((xyz(:,i)-xyz(:,j))**2)
        srab(k+j) = sqrt(sqrab(k+j))
      end do
    end do

    hess = 0.0_wp
    if (iand(terms,gff_term_rep) .ne. 0) then
      call hess_repulsion_nb(nat,at,xyz,sqrab,repthr,1.0_wp,calc%param,calc%topo, &
         & calc%neigh,hess)
      call hess_repulsion_bonded(nat,at,xyz,calc%param,calc%neigh,hess)
    end if
    if (iand(terms,gff_term_angl) .ne. 0) then
      call hess_angles(nat,at,xyz,calc%param,calc%topo,calc%neigh,hess)
    end if
    if (iand(terms,gff_term_batm) .ne. 0) then
      call hess_batm(nat,at,xyz,calc%param,calc%topo,calc%neigh,hess)
    end if
    if (iand(terms,gff_term_tors) .ne. 0) then
      call hess_torsions(nat,at,xyz,calc%param,calc%topo,calc%neigh,hess)
      call hess_storsions(nat,xyz,calc%topo,hess)
    end if
    if (iand(terms,gff_term_xb) .ne. 0) then
      call hess_xbonds(nat,at,xyz,calc%param,calc%topo,calc%neigh,calc%nlist,hess)
    end if
    if (iand(terms,gff_term_hb) .ne. 0) then
      call hess_hbonds_bound(nat,at,xyz,1.0_wp,calc%param,calc%topo, &
         & calc%neigh,calc%nlist,hess)
      call hess_hbonds_unbound(nat,at,xyz,1.0_wp,calc%param,calc%topo, &
         & calc%neigh,calc%nlist,hess,ncov)
    end if
    if (iand(terms,gff_term_es+gff_term_bond+gff_term_disp) .ne. 0) then
      allocate (cn(nat),dcn(3,nat,nat),hb_cn(nat),hb_dcn(3,nat,nat), &
         & dhbcndL(3,3,nat))
      call gfnff_dlogcoord(nat,at,xyz,srab,cn,dcn,cnthr,calc%param)
      if (iand(terms,gff_term_es) .ne. 0) then
        call hess_electrostatics(nat,at,xyz,srab,cnthr,cn,dcn,calc%nlist%q, &
           & calc%param,calc%topo,hess)
      end if
      if (iand(terms,gff_term_disp) .ne. 0) then
        call hess_dispersion(nat,at,xyz,sqrab,srab,dispthr,cnthr,cn,dcn, &
           & calc%param,calc%topo,hess)
      end if
      if (iand(terms,gff_term_bond) .ne. 0) then
        hb_cn = 0.0_wp
        hb_dcn = 0.0_wp
        if (sum(calc%neigh%nr_hb) > 0) call dncoord_erf(nat,at,xyz,calc%param%rcov, &
           & hb_cn,hb_dcn,900.0_wp,calc%topo,calc%neigh,dhbcndL)
        call hess_bonds(nat,at,xyz,srab,cnthr,cn,dcn,hb_cn,hb_dcn,calc%param, &
           & calc%topo,calc%neigh,calc%version,hess)
      end if
    end if

  end subroutine analytic_terms

  subroutine scan_term(terms,error)
    !***********************************
    !* Halving the finite-difference step must quarter the deviation from the
    !* closed form; a wrong Hessian instead leaves a constant offset.
    !***********************************
    integer,intent(in) :: terms
    type(error_type),allocatable,intent(out) :: error
    integer :: nat,i
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),hana(:,:),hfd(:,:)
    real(wp) :: dev(3),h
    type(gfnff_data) :: calc

    call setup(nat,at,xyz,calc,error)
    if (allocated(error)) return
    allocate (hana(3*nat,3*nat),hfd(3*nat,3*nat))

    call analytic_terms(nat,at,xyz,calc,terms,hana)
    call check(error,hess_transl_sumrule(hana) < 1.0e-10_wp)
    if (allocated(error)) return

    do i = 1,3
      h = 1.0e-2_wp/2.0_wp**(i-1)
      call fd_reference(nat,at,xyz,calc,terms,h,hfd)
      dev(i) = maxval(abs(hana-hfd))
    end do

    do i = 1,2
      call check(error,dev(i)/dev(i+1) > 3.7_wp)
      if (allocated(error)) return
      call check(error,dev(i)/dev(i+1) < 4.3_wp)
      if (allocated(error)) return
    end do

  end subroutine scan_term

  subroutine test_hess_angl_scan(error)
    type(error_type),allocatable,intent(out) :: error
    call scan_term(gff_term_angl,error)
  end subroutine test_hess_angl_scan

  subroutine test_hess_batm_scan(error)
    type(error_type),allocatable,intent(out) :: error
    call scan_term(gff_term_batm,error)
  end subroutine test_hess_batm_scan

  subroutine test_hess_bond_scan(error)
    type(error_type),allocatable,intent(out) :: error
    call scan_term(gff_term_bond,error)
  end subroutine test_hess_bond_scan

  subroutine test_hess_tors_scan(error)
    type(error_type),allocatable,intent(out) :: error
    call scan_term(gff_term_tors,error)
  end subroutine test_hess_tors_scan

  subroutine test_hess_disp_scan(error)
    !***********************************
    !* C6 in D3(BJ) is nonlinear in both coordination numbers of a pair; the
    !* scan covers the reference weights differentiated twice.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    call scan_term(gff_term_disp,error)
  end subroutine test_hess_disp_scan

  subroutine test_hess_hb_scan(error)
    !***********************************
    !* Hydrogen bonds on caffeine, which populates the bound list (abhgfnff_eg1)
    !* heavily and the carbonyl form through its two C=O groups.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    call scan_term(gff_term_hb,error)
  end subroutine test_hess_hb_scan

  subroutine test_hess_es_scan(error)
    !***********************************
    !* Electrostatics is the only term whose energy depends on geometry
    !* implicitly, through the EEQ charges; the scan also tests the response solve.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    call scan_term(gff_term_es,error)
  end subroutine test_hess_es_scan

  subroutine test_hess_hbbond(error)
    !***********************************
    !* Water dimer: every bond takes the hydrogen-bond-corrected egbond_hb
    !* path, untested by caffeine, which has none.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: nat = 6
    integer,parameter :: at(nat) = [8,1,1,8,1,1]
    real(wp),parameter :: xyz(3,nat) = reshape([ &
      & -2.7855_wp,0.0000_wp, 0.1200_wp, -0.9525_wp,0.0000_wp,-0.2400_wp, &
      & -3.4020_wp,0.0000_wp,-1.5120_wp,  2.7855_wp,0.0000_wp,-0.1200_wp, &
      &  3.4020_wp,1.4460_wp, 0.8400_wp,  3.4020_wp,-1.4460_wp,0.8400_wp],[3,nat])
    real(wp),allocatable :: hana(:,:),hfd(:,:)
    real(wp) :: dev(3),h
    integer :: io,i
    type(gfnff_data) :: calc

    call gfnff_initialize(nat,at,xyz,calc,ichrg=0,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call prime(nat,at,xyz,calc,error)
    if (allocated(error)) return

    !>-- every bond must actually be HB-corrected, or the test proves nothing
    call check(error,count(calc%neigh%nr_hb(1:calc%neigh%nbond) >= 1), &
       & calc%neigh%nbond)
    if (allocated(error)) return

    allocate (hana(3*nat,3*nat),hfd(3*nat,3*nat))
    call analytic_terms(nat,at,xyz,calc,gff_term_bond,hana)
    call check(error,hess_transl_sumrule(hana) < 1.0e-10_wp)
    if (allocated(error)) return

    do i = 1,3
      h = 1.0e-2_wp/2.0_wp**(i-1)
      call fd_reference(nat,at,xyz,calc,gff_term_bond,h,hfd)
      dev(i) = maxval(abs(hana-hfd))
    end do
    do i = 1,2
      call check(error,dev(i)/dev(i+1) > 3.7_wp)
      if (allocated(error)) return
      call check(error,dev(i)/dev(i+1) < 4.3_wp)
      if (allocated(error)) return
    end do

  end subroutine test_hess_hbbond

  subroutine test_hess_bond_conformer(error)
    !***********************************
    !* Bond Hessian of conformer2020 with one bond pulled past the inflection
    !* point, exercising the linear continuation's own closed-form partials. The
    !* energy is checked against the published version first: on an undisturbed
    !* geometry the branches coincide and the scan never enters the code it targets.
    !***********************************
    use coffeine
    use gfnff_interface,only:gffVersion
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: iat = 15
    real(wp),parameter :: stretch = 3.0_wp
    integer :: nat,io,i
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),hana(:,:),hfd(:,:),grad(:,:)
    real(wp) :: dev(3),h,dir(3),eref,econ
    type(gfnff_data) :: calc,calc_ref

    nat = testnat
    allocate (at(nat),xyz(3,nat),grad(3,nat))
    at = testat
    xyz = testxyz

    !> topologies are built at the reference geometry and then held, so the
    !> bond survives being pulled apart
    call gfnff_initialize(nat,at,xyz,calc,ichrg=0,iostat=io,printlevel=0, &
       & version=gffVersion%conformer2020)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_initialize(nat,at,xyz,calc_ref,ichrg=0,iostat=io,printlevel=0, &
       & version=gffVersion%angewChem2020_2)
    call check(error,io,0)
    if (allocated(error)) return

    dir = xyz(:,iat)-xyz(:,1)
    dir = dir/sqrt(sum(dir*dir))
    xyz(:,iat) = xyz(:,iat)+dir*stretch

    call gfnff_singlepoint(nat,at,xyz,calc_ref,eref,grad,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,econ,grad,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    if (abs(econ-eref) < 0.1_wp) then
      call test_failed(error,"continuation branch is not active; test is void")
      write (*,'(a,2f14.6)') "  E =",eref,econ
      return
    end if

    allocate (hana(3*nat,3*nat),hfd(3*nat,3*nat))
    call analytic_terms(nat,at,xyz,calc,gff_term_bond,hana)
    call check(error,hess_transl_sumrule(hana) < 1.0e-10_wp)
    if (allocated(error)) return

    do i = 1,3
      h = 1.0e-2_wp/2.0_wp**(i-1)
      call fd_reference(nat,at,xyz,calc,gff_term_bond,h,hfd)
      dev(i) = maxval(abs(hana-hfd))
    end do
    do i = 1,2
      call check(error,dev(i)/dev(i+1) > 3.7_wp)
      if (allocated(error)) return
      call check(error,dev(i)/dev(i+1) < 4.3_wp)
      if (allocated(error)) return
    end do

  end subroutine test_hess_bond_conformer

  subroutine dpa_geometry(alpha,xyz)
    !***********************************
    !* Diphenylacetylene, parametrisation molecule for the special torsion
    !* around a triple-bonded carbon; alpha twists the second ring about the axis.
    !***********************************
    real(wp),intent(in) :: alpha
    real(wp),intent(out) :: xyz(3,24)
    real(wp),parameter :: aatoau = 1.0_wp/0.52917726_wp
    real(wp),parameter :: flat(3,24) = reshape([ &
      &  0.000000_wp, 0.000000_wp,0.000000_wp,  1.200000_wp, 0.000000_wp,0.000000_wp, &
      & -1.430000_wp, 0.000000_wp,0.000000_wp, -2.125000_wp, 1.203775_wp,0.000000_wp, &
      & -3.515000_wp, 1.203775_wp,0.000000_wp, -4.210000_wp, 0.000000_wp,0.000000_wp, &
      & -3.515000_wp,-1.203775_wp,0.000000_wp, -2.125000_wp,-1.203775_wp,0.000000_wp, &
      &  2.630000_wp, 0.000000_wp,0.000000_wp,  3.325000_wp, 1.203775_wp,0.000000_wp, &
      &  4.715000_wp, 1.203775_wp,0.000000_wp,  5.410000_wp, 0.000000_wp,0.000000_wp, &
      &  4.715000_wp,-1.203775_wp,0.000000_wp,  3.325000_wp,-1.203775_wp,0.000000_wp, &
      & -1.585000_wp, 2.139083_wp,0.000000_wp, -4.055000_wp, 2.139083_wp,0.000000_wp, &
      & -5.290000_wp, 0.000000_wp,0.000000_wp, -4.055000_wp,-2.139083_wp,0.000000_wp, &
      & -1.585000_wp,-2.139083_wp,0.000000_wp,  2.785000_wp, 2.139083_wp,0.000000_wp, &
      &  5.255000_wp, 2.139083_wp,0.000000_wp,  6.490000_wp, 0.000000_wp,0.000000_wp, &
      &  5.255000_wp,-2.139083_wp,0.000000_wp,  2.785000_wp,-2.139083_wp,0.000000_wp], &
      & [3,24])
    integer :: i
    real(wp) :: ca,sa,y,z

    xyz = flat
    ca = cos(alpha)
    sa = sin(alpha)
    do i = 1,24
      if ((i >= 9 .and. i <= 14) .or. i >= 20) then
        y = xyz(2,i)
        z = xyz(3,i)
        xyz(2,i) = ca*y-sa*z
        xyz(3,i) = sa*y+ca*z
      end if
    end do
    xyz = xyz*aatoau

  end subroutine dpa_geometry

  subroutine test_hess_stors(error)
    !***********************************
    !* Special torsion around a triple-bonded carbon, E = e0(1 - cos 2phi).
    !* Checks a twisted geometry and the planar reference (phi = 0), where a
    !* formulation differentiating through phi rather than cos(phi) would fail.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: nat = 24
    integer,parameter :: at(nat) = [6,6,6,6,6,6,6,6,6,6,6,6,6,6, &
       &                            1,1,1,1,1,1,1,1,1,1]
    real(wp) :: xyz(3,nat)
    real(wp),allocatable :: hana(:,:),hfd(:,:),hs(:,:)
    real(wp) :: dev(3),h
    integer :: io,i
    type(gfnff_data) :: calc

    call dpa_geometry(0.7_wp,xyz)
    call gfnff_initialize(nat,at,xyz,calc,ichrg=0,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call prime(nat,at,xyz,calc,error)
    if (allocated(error)) return

    !>-- special torsion list must be non-empty, or the test proves nothing
    call check(error,calc%topo%nstors > 0)
    if (allocated(error)) return

    allocate (hana(3*nat,3*nat),hfd(3*nat,3*nat),hs(3*nat,3*nat))

    !>-- its contribution must be large enough for the scan to see it
    hs = 0.0_wp
    call hess_storsions(nat,xyz,calc%topo,hs)
    call check(error,maxval(abs(hs)) > 1.0e-5_wp)
    if (allocated(error)) return
    call check(error,hess_transl_sumrule(hs) < 1.0e-12_wp)
    if (allocated(error)) return

    call analytic_terms(nat,at,xyz,calc,gff_term_tors,hana)
    do i = 1,3
      h = 1.0e-2_wp/2.0_wp**(i-1)
      call fd_reference(nat,at,xyz,calc,gff_term_tors,h,hfd)
      dev(i) = maxval(abs(hana-hfd))
    end do
    do i = 1,2
      call check(error,dev(i)/dev(i+1) > 3.7_wp)
      if (allocated(error)) return
      call check(error,dev(i)/dev(i+1) < 4.3_wp)
      if (allocated(error)) return
    end do

    !>-- planar reference, phi = 0 exactly
    call dpa_geometry(0.0_wp,xyz)
    call prime(nat,at,xyz,calc,error)
    if (allocated(error)) return
    call analytic_terms(nat,at,xyz,calc,gff_term_tors,hana)
    call check(error,hess_transl_sumrule(hana) < 1.0e-10_wp)
    if (allocated(error)) return
    call fd_reference(nat,at,xyz,calc,gff_term_tors,2.5e-3_wp,hfd)
    call check(error,maxval(abs(hana-hfd)) < 1.0e-5_wp)

  end subroutine test_hess_stors

  subroutine test_hess_es_frag(error)
    !***********************************
    !* Water dimer: the topology splits it into two fragments, giving the EEQ
    !* system a constraint border that caffeine's single fragment never exercises.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: nat = 6
    integer,parameter :: at(nat) = [8,1,1,8,1,1]
    real(wp),parameter :: xyz(3,nat) = reshape([ &
      & -2.7855_wp,0.0000_wp, 0.1200_wp, -0.9525_wp,0.0000_wp,-0.2400_wp, &
      & -3.4020_wp,0.0000_wp,-1.5120_wp,  2.7855_wp,0.0000_wp,-0.1200_wp, &
      &  3.4020_wp,1.4460_wp, 0.8400_wp,  3.4020_wp,-1.4460_wp,0.8400_wp],[3,nat])
    real(wp),allocatable :: hana(:,:),hfd(:,:)
    real(wp) :: dev(3),h
    integer :: io,i
    type(gfnff_data) :: calc

    call gfnff_initialize(nat,at,xyz,calc,ichrg=0,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call prime(nat,at,xyz,calc,error)
    if (allocated(error)) return

    call check(error,calc%topo%nfrag,2)
    if (allocated(error)) return

    allocate (hana(3*nat,3*nat),hfd(3*nat,3*nat))
    call analytic_terms(nat,at,xyz,calc,gff_term_es,hana)
    call check(error,hess_transl_sumrule(hana) < 1.0e-10_wp)
    if (allocated(error)) return

    do i = 1,3
      h = 1.0e-2_wp/2.0_wp**(i-1)
      call fd_reference(nat,at,xyz,calc,gff_term_es,h,hfd)
      dev(i) = maxval(abs(hana-hfd))
    end do
    do i = 1,2
      call check(error,dev(i)/dev(i+1) > 3.7_wp)
      if (allocated(error)) return
      call check(error,dev(i)/dev(i+1) < 4.3_wp)
      if (allocated(error)) return
    end do

  end subroutine test_hess_es_frag

  subroutine hb_case(label,nat,at,xyz,terms,error)
    !***********************************
    !* Closed-form Hessian of `terms` vs. a step-scanned finite difference;
    !* `label` only orients a failure message.
    !***********************************
    character(len=*),intent(in) :: label
    integer,intent(in) :: nat,at(nat),terms
    real(wp),intent(in) :: xyz(3,nat)
    type(error_type),allocatable,intent(out) :: error
    real(wp),allocatable :: hana(:,:),hfd(:,:)
    real(wp) :: dev(3),h
    integer :: io,i
    type(gfnff_data) :: calc

    call gfnff_initialize(nat,at,xyz,calc,ichrg=0,iostat=io,printlevel=0)
    call check(error,io,0)
    if (allocated(error)) return
    call prime(nat,at,xyz,calc,error)
    if (allocated(error)) return

    allocate (hana(3*nat,3*nat),hfd(3*nat,3*nat))
    call analytic_terms(nat,at,xyz,calc,terms,hana)
    call check(error,hess_transl_sumrule(hana) < 1.0e-10_wp)
    if (allocated(error)) return
    !>-- term must actually be present, or the scan proves nothing
    call check(error,maxval(abs(hana)) > 1.0e-6_wp)
    if (allocated(error)) return

    do i = 1,3
      h = 1.0e-2_wp/2.0_wp**(i-1)
      call fd_reference(nat,at,xyz,calc,terms,h,hfd)
      dev(i) = maxval(abs(hana-hfd))
    end do
    do i = 1,2
      call check(error,dev(i)/dev(i+1) > 3.7_wp)
      if (allocated(error)) return
      call check(error,dev(i)/dev(i+1) < 4.3_wp)
      if (allocated(error)) return
    end do

  end subroutine hb_case

  subroutine test_hess_hb_forms(error)
    !***********************************
    !* Three unbound hydrogen-bond acceptor forms; caffeine exercises only two,
    !*   water dimer      - DEFAULT form, oriented by the neighbours of B
    !*   acetone + water  - CARBONYL form, with a bend and torsion factor
    !*   pyridine + water - N-heteroaromatic form, lone pair on the bisector
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: nw = 6
    integer,parameter :: atw(nw) = [8,1,1,8,1,1]
    real(wp),parameter :: xyzw(3,nw) = reshape([ &
      & -2.7855_wp,0.0000_wp, 0.1200_wp, -0.9525_wp,0.0000_wp,-0.2400_wp, &
      & -3.4020_wp,0.0000_wp,-1.5120_wp,  2.7855_wp,0.0000_wp,-0.1200_wp, &
      &  3.4020_wp,1.4460_wp, 0.8400_wp,  3.4020_wp,-1.4460_wp,0.8400_wp],[3,nw])
    integer,parameter :: nc = 13
    integer,parameter :: atc(nc) = [6,8,6,6,1,1,1,1,1,1,8,1,1]
    real(wp),parameter :: xyzc(3,nc) = reshape([ &
      &  0.0000_wp, 0.0000_wp, 0.0000_wp,  0.0000_wp, 0.0000_wp, 1.2200_wp, &
      &  1.2900_wp, 0.0000_wp,-0.7700_wp, -1.2900_wp, 0.0000_wp,-0.7700_wp, &
      &  1.3400_wp, 0.8800_wp,-1.4200_wp,  1.3400_wp,-0.8800_wp,-1.4200_wp, &
      &  2.1400_wp, 0.0000_wp,-0.0900_wp, -1.3400_wp, 0.8800_wp,-1.4200_wp, &
      & -1.3400_wp,-0.8800_wp,-1.4200_wp, -2.1400_wp, 0.0000_wp,-0.0900_wp, &
      &  1.6500_wp, 0.9500_wp, 3.4000_wp,  1.1000_wp, 0.6300_wp, 2.6300_wp, &
      &  2.5300_wp, 0.6100_wp, 3.2200_wp],[3,nc])
    integer,parameter :: np = 14
    integer,parameter :: atp(np) = [7,6,6,6,6,6,1,1,1,1,1,8,1,1]
    real(wp),parameter :: xyzp(3,np) = reshape([ &
      &  0.0000_wp, 1.4180_wp, 0.0000_wp,  1.1400_wp, 0.7100_wp, 0.0000_wp, &
      &  1.1900_wp,-0.6800_wp, 0.0000_wp,  0.0000_wp,-1.3900_wp, 0.0000_wp, &
      & -1.1900_wp,-0.6800_wp, 0.0000_wp, -1.1400_wp, 0.7100_wp, 0.0000_wp, &
      &  2.0600_wp, 1.2900_wp, 0.0000_wp,  2.1400_wp,-1.2100_wp, 0.0000_wp, &
      &  0.0000_wp,-2.4800_wp, 0.0000_wp, -2.1400_wp,-1.2100_wp, 0.0000_wp, &
      & -2.0600_wp, 1.2900_wp, 0.0000_wp,  0.0000_wp, 4.2000_wp, 0.0000_wp, &
      &  0.0000_wp, 3.2400_wp, 0.0000_wp,  0.9500_wp, 4.4500_wp, 0.0000_wp],[3,np])
    real(wp),parameter :: aatoau = 1.0_wp/0.52917726_wp

    call hb_case('water dimer',nw,atw,xyzw,gff_term_hb,error)
    if (allocated(error)) return
    call hb_case('acetone-water',nc,atc,xyzc*aatoau,gff_term_hb,error)
    if (allocated(error)) return
    call hb_case('pyridine-water',np,atp,xyzp*aatoau,gff_term_hb,error)

  end subroutine test_hess_hb_forms

  subroutine test_hess_xb(error)
    !***********************************
    !* Halogen bond N...Br-C; the out-of-line factor measures bending via
    !* (r_AX + r_BX)/r_AB rather than an explicit angle.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: nat = 9
    integer,parameter :: at(nat) = [35,6,1,1,1,7,1,1,1]
    real(wp),parameter :: xyz(3,nat) = reshape([ &
      &  0.0000_wp, 0.0000_wp, 0.0000_wp,  0.0000_wp, 0.0000_wp, 1.9700_wp, &
      &  1.0270_wp, 0.0000_wp, 2.3300_wp, -0.5135_wp, 0.8894_wp, 2.3300_wp, &
      & -0.5135_wp,-0.8894_wp, 2.3300_wp,  0.0000_wp, 0.0000_wp,-3.0500_wp, &
      &  0.4776_wp, 0.8272_wp,-3.4300_wp,  0.4776_wp,-0.8272_wp,-3.4300_wp, &
      & -0.9552_wp, 0.0000_wp,-3.4300_wp],[3,nat])
    real(wp),parameter :: aatoau = 1.0_wp/0.52917726_wp

    call hb_case('CH3Br-NH3',nat,at,xyz*aatoau,gff_term_xb,error)

  end subroutine test_hess_xb

  subroutine fd_reference(nat,at,xyz,calc,terms,step,hess)
    !***********************************
    !* Finite-difference Hessian of the selected terms only.
    !***********************************
    integer,intent(in) :: nat,at(nat),terms
    real(wp),intent(in) :: xyz(3,nat),step
    type(gfnff_data),intent(inout) :: calc
    real(wp),intent(out) :: hess(3*nat,3*nat)
    real(wp) :: efield(3)

    efield = 0.0_wp
    call hess_term_reference(0,nat,at,xyz,calc%cell,calc%ichrg,calc%param, &
       & calc%topo,calc%neigh,calc%nlist,efield,calc%solvation,calc%version, &
       & calc%accuracy,terms,hess,step=step)

  end subroutine fd_reference

  subroutine test_hess_rep_fd(error)
    !***********************************
    !* Closed-form repulsion vs. finite difference; tolerance is set by the
    !* O(h^2) truncation error of the reference step, not the analytic side.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),hana(:,:),hfd(:,:)
    type(gfnff_data) :: calc

    call setup(nat,at,xyz,calc,error)
    if (allocated(error)) return
    allocate (hana(3*nat,3*nat),hfd(3*nat,3*nat))

    call analytic_terms(nat,at,xyz,calc,gff_term_rep,hana)
    call fd_reference(nat,at,xyz,calc,gff_term_rep,2.5e-3_wp,hfd)

    call check(error,maxval(abs(hana-hfd)) < 1.0e-4_wp)
    if (allocated(error)) return
    call check(error,sqrt(sum((hana-hfd)**2)/size(hana)) < 1.0e-6_wp)

  end subroutine test_hess_rep_fd

  subroutine test_hess_rep_scan(error)
    !***********************************
    !* Halving the step must quarter the deviation; a wrong Hessian instead
    !* leaves a constant offset no refinement removes.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    integer :: nat,i
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),hana(:,:),hfd(:,:)
    real(wp) :: dev(3),h
    type(gfnff_data) :: calc

    call setup(nat,at,xyz,calc,error)
    if (allocated(error)) return
    allocate (hana(3*nat,3*nat),hfd(3*nat,3*nat))

    call analytic_terms(nat,at,xyz,calc,gff_term_rep,hana)

    do i = 1,3
      h = 1.0e-2_wp/2.0_wp**(i-1)
      call fd_reference(nat,at,xyz,calc,gff_term_rep,h,hfd)
      dev(i) = maxval(abs(hana-hfd))
    end do

    do i = 1,2
      call check(error,dev(i)/dev(i+1) > 3.7_wp)
      if (allocated(error)) return
      call check(error,dev(i)/dev(i+1) < 4.3_wp)
      if (allocated(error)) return
    end do

  end subroutine test_hess_rep_scan

  subroutine test_hess_sumrule(error)
    !***********************************
    !* Repulsion Hessian must satisfy the translational sum rule and come out
    !* symmetric to machine precision, exercising the own-row parallel scatter.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),hana(:,:)
    type(gfnff_data) :: calc

    call setup(nat,at,xyz,calc,error)
    if (allocated(error)) return
    allocate (hana(3*nat,3*nat))

    call analytic_terms(nat,at,xyz,calc,gff_term_rep,hana)

    call check(error,hess_transl_sumrule(hana) < 1.0e-12_wp)
    if (allocated(error)) return
    call check(error,hess_asymmetry(hana) < 1.0e-14_wp)

  end subroutine test_hess_sumrule

  subroutine test_hess_full(error)
    !***********************************
    !* The public driver, which mixes closed-form and finite-difference terms,
    !* must agree with a Hessian obtained by differencing every term.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),hhyb(:,:),hall(:,:)
    real(wp) :: energy
    type(gfnff_data) :: calc

    call setup(nat,at,xyz,calc,error)
    if (allocated(error)) return
    allocate (hhyb(3*nat,3*nat),hall(3*nat,3*nat))

    call gfnff_hessian(nat,at,xyz,calc,hhyb,energy=energy,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    call fd_reference(nat,at,xyz,calc,gff_term_all,5.0e-3_wp,hall)

    call check(error,maxval(abs(hhyb-hall)) < 1.0e-4_wp)

  end subroutine test_hess_full

  subroutine test_hess_freq(error)
    !***********************************
    !* Diagonalisation must return three vanishing translational modes; this
    !* geometry is not a stationary point, so rotational modes need not vanish.
    !***********************************
    type(error_type),allocatable,intent(out) :: error
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),hess(:,:),freq(:)
    type(gfnff_data) :: calc

    call setup(nat,at,xyz,calc,error)
    if (allocated(error)) return
    allocate (hess(3*nat,3*nat))

    call gfnff_hessian(nat,at,xyz,calc,hess)
    call hess_frequencies(hess,nat,at,freq,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,size(freq),3*nat)
    if (allocated(error)) return
    call check(error,maxval(abs(freq(1:3))) < 1.0_wp)

  end subroutine test_hess_freq

end module test_hessian
