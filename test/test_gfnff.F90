module test_gfnff
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_interface
  use gfnff_type_timer
  implicit none
  private

  public :: collect_gfnff

  real(wp),parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
  real(wp),parameter :: thr2 = 1e-6

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for PV calculations
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_gfnff(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
    new_unittest("GFN-FF singlepoint calculation ",test_gfnff_sp), &
    new_unittest("GFN-FF singlepoint with ALPB   ",test_gfnff_alpb), &
!    new_unittest("GFN-FF OpenMP parallel SP      ",test_gfnff_openmp), &
    new_unittest("GFN-FF numerical gradient      ",test_gfnff_numgrad), &
    new_unittest("GFN-FF net force vanishes      ",test_gfnff_netforce), &
    new_unittest("GFN-FF translation invariance  ",test_gfnff_translation), &
    new_unittest("GFN-FF energy decomposition    ",test_gfnff_edecomp), &
    new_unittest("GFN-FF energy components       ",test_gfnff_components), &
    new_unittest("GFN-FF energy components ALPB  ",test_gfnff_components_alpb), &
    new_unittest("GFN-FF supermol singlepoint    ",test_gfnff_supermol) &
    ]
!&>
  end subroutine collect_gfnff

!========================================================================================!

  subroutine test_gfnff_sp(error)
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io,ichrg
    type(gfnff_data) :: calculator

!&<
    real(wp),parameter :: e_ref = -4.672792533926004_wp
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    &   0.005301570264175_wp,   0.000273970046453_wp,   0.000002235966967_wp, &
    &   0.008166037109104_wp,  -0.008220839180901_wp,  -0.000025577434354_wp, &
    &  -0.003078325363552_wp,  -0.009432996299921_wp,   0.000033248959973_wp, &
    &   0.009919832601386_wp,   0.008086633755534_wp,  -0.000022035642360_wp, &
    &  -0.015632596323341_wp,  -0.026672391134961_wp,   0.000004837606473_wp, &
    &   0.014525642097464_wp,  -0.001976846297509_wp,   0.000067168700586_wp, &
    &   0.006146643879669_wp,   0.009561075520487_wp,  -0.000017505347402_wp, &
    &  -0.008820848986042_wp,  -0.001068632415840_wp,  -0.000078000868871_wp, &
    &  -0.000983352664777_wp,   0.014873585269955_wp,   0.000032976017459_wp, &
    &  -0.006683041231125_wp,   0.007422826993429_wp,  -0.000019221295612_wp, &
    &   0.012839290909399_wp,  -0.012743003179261_wp,  -0.000039643202527_wp, &
    &  -0.023422681331404_wp,   0.021005865685095_wp,  -0.000002459581560_wp, &
    &  -0.001884040385407_wp,  -0.003906626891817_wp,  -0.000013746286938_wp, &
    &  -0.003754972778577_wp,   0.003730224340046_wp,  -0.000073759269033_wp, &
    &   0.000742833683906_wp,   0.003621866120529_wp,   0.000003807413478_wp, &
    &   0.001069304816898_wp,  -0.000350576122870_wp,   0.003705269737750_wp, &
    &   0.001070927784981_wp,  -0.000349786548726_wp,  -0.003711459088386_wp, &
    &  -0.002984451042725_wp,   0.000421241696501_wp,   0.000013800909435_wp, &
    &   0.004499278381845_wp,   0.000660466765138_wp,   0.000002343067289_wp, &
    &   0.000371386732209_wp,   0.001498980424092_wp,   0.003776579480469_wp, &
    &   0.000381320666177_wp,   0.001507608943116_wp,  -0.003766961246878_wp, &
    &  -0.001010476530386_wp,  -0.004606702150652_wp,   0.000057155057056_wp, &
    &   0.001617340754014_wp,  -0.001636910544025_wp,   0.003219103002614_wp, &
    &   0.001603376956109_wp,  -0.001699034793892_wp,  -0.003148156655631_wp  &
    & ], shape(g_ref))
!&>

    !> setup
    nat = testnat
    allocate (at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz
    ichrg = 0 !> mol. charge
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)

    !> calculation
    call gfnff_initialize(nat,at,xyz,calculator,ichrg=ichrg,iostat=io)
    call gfnff_singlepoint(nat,at,xyz,calculator,energy,grad,iostat=io)
    !write (*,'(F25.15)') energy
    !write (*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=5e-4_wp)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > thr2)) then
      call test_failed(error,"Gradient does not match reference")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
      print'("---")'
      print'(3es21.14)',grad-g_ref
    end if

    deallocate (grad)
  end subroutine test_gfnff_sp

!========================================================================================!

  subroutine test_gfnff_numgrad(error)
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io,i,j,ichrg
    real(wp) :: step,bw,bw2,fw,fw2
    type(gfnff_data) :: calculator
    real(wp),allocatable :: gradient(:,:),g_ref(:,:),stencil(:,:)
!&<
    real(wp),parameter :: e_ref = -4.672792533926004_wp
!&>

    !> setup
    nat = testnat
    allocate (at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz
    ichrg = 0
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)
    allocate (gradient(3,nat),g_ref(3,nat),stencil(3,nat),source=0.0_wp)

    !> calculation
    call gfnff_initialize(nat,at,xyz,calculator,ichrg=ichrg,iostat=io)
    call gfnff_singlepoint(nat,at,xyz,calculator,energy,grad,iostat=io)
    ! write (*,'(F25.15)') energy
    ! write (*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    stencil = xyz
    step = 0.001_wp
    do i = 1,nat
      do j = 1,3
        !write (*,*) 'Numerical gradient dimension ', (i-1)*3+j
        stencil(j,i) = stencil(j,i)-2.0_wp*step
        call gfnff_singlepoint(nat,at,stencil,calculator,bw2,gradient,iostat=io)
        stencil(j,i) = xyz(j,i)
        stencil(j,i) = stencil(j,i)-1.0_wp*step
        call gfnff_singlepoint(nat,at,stencil,calculator,bw,gradient,iostat=io)
        stencil(j,i) = xyz(j,i)
        stencil(j,i) = stencil(j,i)+1.0_wp*step
        call gfnff_singlepoint(nat,at,stencil,calculator,fw,gradient,iostat=io)
        stencil(j,i) = xyz(j,i)
        stencil(j,i) = stencil(j,i)+2.0_wp*step
        call gfnff_singlepoint(nat,at,stencil,calculator,fw2,gradient,iostat=io)
        stencil(j,i) = xyz(j,i)
        g_ref(j,i) = (bw2/12.0_wp-8.0_wp*bw/12.0_wp+8.0_wp*fw/12.0_wp-fw2/12.0_wp)/step
      end do
    end do

    call check(error,energy,e_ref,thr=5e-4_wp)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > 5e-4_wp)) then
      call test_failed(error,"Gradient does not match reference")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
      print'("---")'
      print'(3es21.14)',grad-g_ref
    end if

    deallocate (grad)
  end subroutine test_gfnff_numgrad

!========================================================================================!

  subroutine test_gfnff_alpb(error)
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io,ichrg
    character(len=:),allocatable :: alpbsolvent
    type(gfnff_data) :: calculator

!&<
    real(wp),parameter :: e_ref = -4.689906356924923_wp
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    &   0.005946311679280_wp,   0.000224845497002_wp,   0.000002359259482_wp, &
    &   0.008276016514049_wp,  -0.008052841401439_wp,  -0.000024497567079_wp, &
    &  -0.002933432054651_wp,  -0.008965786999383_wp,   0.000032565118339_wp, &
    &   0.009327449621835_wp,   0.007163090250006_wp,  -0.000020616237641_wp, &
    &  -0.015612387108414_wp,  -0.026601053781601_wp,   0.000003393879688_wp, &
    &   0.014224598433093_wp,  -0.002328349901387_wp,   0.000066592783236_wp, &
    &   0.006359112133353_wp,   0.009728706323651_wp,  -0.000012692555021_wp, &
    &  -0.008060303655593_wp,  -0.001017324386475_wp,  -0.000067764705364_wp, &
    &  -0.000928875347148_wp,   0.014721272004310_wp,   0.000064282338341_wp, &
    &  -0.007032628120571_wp,   0.007686457118708_wp,  -0.000018779457567_wp, &
    &   0.012172269830222_wp,  -0.012198147898115_wp,  -0.000031535161167_wp, &
    &  -0.023075279911968_wp,   0.020590487741886_wp,   0.000000950277558_wp, &
    &  -0.002527745786103_wp,  -0.004378684545449_wp,  -0.000014408181911_wp, &
    &  -0.003644476831810_wp,   0.004533746951955_wp,  -0.000098967022506_wp, &
    &   0.000763588512030_wp,   0.003493542876475_wp,   0.000003659962064_wp, &
    &   0.001177972011645_wp,  -0.000489794359700_wp,   0.003518470896582_wp, &
    &   0.001179450328452_wp,  -0.000489265524270_wp,  -0.003525395593443_wp, &
    &  -0.002858206649569_wp,  -0.000053699956870_wp,   0.000013389521786_wp, &
    &   0.004179439229098_wp,   0.000474452680083_wp,   0.000002349455888_wp, &
    &   0.000262934291862_wp,   0.001522204493760_wp,   0.003711717192812_wp, &
    &   0.000272333136041_wp,   0.001531298019793_wp,  -0.003702180450665_wp, &
    &  -0.001005504791715_wp,  -0.004218548582919_wp,   0.000052754706566_wp, &
    &   0.001779352856960_wp,  -0.001420259046787_wp,   0.003106520639578_wp, &
    &   0.001758011679622_wp,  -0.001456347573236_wp,  -0.003062169099557_wp  &
    & ], shape(g_ref))
!&>

    !> setup
    nat = testnat
    allocate (at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz
    ichrg = 0 !> mol. charge
    alpbsolvent = 'h2o'
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)

    !> calculation
    call calculator%init(nat,at,xyz,ichrg=ichrg,iostat=io,solvent=alpbsolvent)
    call gfnff_singlepoint(nat,at,xyz,calculator,energy,grad,iostat=io)
    !write (*,'(F25.15)') energy
    !write (*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=5e-4_wp)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > thr2)) then
      call test_failed(error,"Gradient does not match reference")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
      print'("---")'
      print'(3es21.14)',grad-g_ref
    end if

    deallocate (grad)
  end subroutine test_gfnff_alpb

!========================================================================================!

  subroutine test_gfnff_netforce(error)
    !***********************************************
    !* Net force (sum of gradients) must vanish.   *
    !* No reference data needed; purely physical.  *
    !***********************************************
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io,k
    real(wp) :: fnet(3)
    type(gfnff_data) :: calculator

    nat = testnat
    allocate(at(nat),xyz(3,nat),grad(3,nat))
    at = testat
    xyz = testxyz
    energy = 0.0_wp
    grad = 0.0_wp

    call gfnff_initialize(nat,at,xyz,calculator,ichrg=0,iostat=io)
    call check(error,io,0); if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calculator,energy,grad,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    fnet = 0.0_wp
    do k = 1,nat
      fnet(:) = fnet(:) + grad(:,k)
    end do

    if (any(abs(fnet) > 1.0e-10_wp)) then
      call test_failed(error,"Net force is not zero")
      write(*,'(a,3es12.4)') "  F_net =",fnet
    end if
  end subroutine test_gfnff_netforce

!========================================================================================!

  subroutine test_gfnff_translation(error)
    !***********************************************
    !* Translating the molecule must leave energy  *
    !* and gradient magnitudes unchanged.          *
    !***********************************************
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: e0,e1
    real(wp),allocatable :: xyz(:,:),xyz_shifted(:,:),grad0(:,:),grad1(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io,i
    type(gfnff_data) :: calc0,calc1
    real(wp),parameter :: shift(3) = [10.0_wp, -7.3_wp, 4.1_wp]

    nat = testnat
    allocate(at(nat),xyz(3,nat),xyz_shifted(3,nat),grad0(3,nat),grad1(3,nat))
    at = testat
    xyz = testxyz
    xyz_shifted = xyz
    do i = 1,nat
      xyz_shifted(:,i) = xyz_shifted(:,i) + shift(:)
    end do

    e0 = 0.0_wp; grad0 = 0.0_wp
    call gfnff_initialize(nat,at,xyz,calc0,ichrg=0,iostat=io)
    call check(error,io,0); if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc0,e0,grad0,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    e1 = 0.0_wp; grad1 = 0.0_wp
    call gfnff_initialize(nat,at,xyz_shifted,calc1,ichrg=0,iostat=io)
    call check(error,io,0); if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz_shifted,calc1,e1,grad1,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    call check(error,e1,e0,thr=1.0e-10_wp)
    if (allocated(error)) then
      call test_failed(error,"Energy changed after translation")
      return
    end if

    if (any(abs(grad1 - grad0) > 1.0e-10_wp)) then
      call test_failed(error,"Gradient changed after translation")
    end if
  end subroutine test_gfnff_translation

!========================================================================================!

  subroutine test_gfnff_edecomp(error)
    !***********************************************
    !* The sum of all energy components must equal *
    !* the total energy stored in res%e_total.     *
    !***********************************************
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy,esum
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io
    type(gfnff_data) :: calculator

    nat = testnat
    allocate(at(nat),xyz(3,nat),grad(3,nat))
    at = testat
    xyz = testxyz

    call gfnff_initialize(nat,at,xyz,calculator,ichrg=0,iostat=io)
    call check(error,io,0); if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calculator,energy,grad,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    associate(r => calculator%res)
      esum = r%e_bond + r%e_angl + r%e_tors + r%e_batm &
           + r%e_rep  + r%e_es   + r%e_disp + r%e_hb   &
           + r%e_xb   + r%e_ext

      call check(error,esum,r%e_total,thr=1.0e-12_wp)
      if (allocated(error)) then
        call test_failed(error,"Energy component sum does not match e_total")
        write(*,'(a,es20.12)') "  component sum =",esum
        write(*,'(a,es20.12)') "  e_total       =",r%e_total
      end if
    end associate
  end subroutine test_gfnff_edecomp

!========================================================================================!

  subroutine test_gfnff_components(error)
    !*****************************************************
    !* Regression test for all individual energy terms  *
    !* of caffeine (gas phase, no solvation).            *
    !*****************************************************
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io
    type(gfnff_data) :: calculator
!&<
    real(wp),parameter :: e_bond_ref = -4.798457754311154_wp
    real(wp),parameter :: e_angl_ref =  0.018023083376745_wp
    real(wp),parameter :: e_tors_ref =  0.000891368462549_wp
    real(wp),parameter :: e_batm_ref = -0.000967109393591_wp
    real(wp),parameter :: e_rep_ref  =  0.300182688960152_wp
    real(wp),parameter :: e_es_ref   = -0.174351370463650_wp
    real(wp),parameter :: e_disp_ref = -0.018113406607768_wp
    real(wp),parameter :: e_hb_ref   = -0.000000033949287_wp
    real(wp),parameter :: e_xb_ref   =  0.000000000000000_wp
!&>
    real(wp),parameter :: cthr = 1.0e-7_wp

    nat = testnat
    allocate(at(nat),xyz(3,nat),grad(3,nat))
    at = testat; xyz = testxyz
    call gfnff_initialize(nat,at,xyz,calculator,ichrg=0,iostat=io)
    call check(error,io,0); if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calculator,energy,grad,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    associate(r => calculator%res)
      call check(error,r%e_bond,e_bond_ref,thr=cthr); if (allocated(error)) return
      call check(error,r%e_angl,e_angl_ref,thr=cthr); if (allocated(error)) return
      call check(error,r%e_tors,e_tors_ref,thr=cthr); if (allocated(error)) return
      call check(error,r%e_batm,e_batm_ref,thr=cthr); if (allocated(error)) return
      call check(error,r%e_rep, e_rep_ref, thr=cthr); if (allocated(error)) return
      call check(error,r%e_es,  e_es_ref,  thr=cthr); if (allocated(error)) return
      call check(error,r%e_disp,e_disp_ref,thr=cthr); if (allocated(error)) return
      call check(error,r%e_hb,  e_hb_ref,  thr=cthr); if (allocated(error)) return
      call check(error,r%e_xb,  e_xb_ref,  thr=cthr)
    end associate
  end subroutine test_gfnff_components

!========================================================================================!

  subroutine test_gfnff_components_alpb(error)
    !*****************************************************
    !* Regression test for all individual energy terms  *
    !* of caffeine with ALPB implicit solvation (water). *
    !*****************************************************
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io
    type(gfnff_data) :: calculator
!&<
    real(wp),parameter :: e_bond_ref  = -4.798457754311155_wp
    real(wp),parameter :: e_angl_ref  =  0.018023083376745_wp
    real(wp),parameter :: e_tors_ref  =  0.000891368462549_wp
    real(wp),parameter :: e_batm_ref  = -0.000967109393591_wp
    real(wp),parameter :: e_rep_ref   =  0.300182688960152_wp
    real(wp),parameter :: e_es_ref    = -0.171030958335199_wp
    real(wp),parameter :: e_disp_ref  = -0.018113406607768_wp
    real(wp),parameter :: e_hb_ref    = -0.000000033949287_wp
    real(wp),parameter :: e_xb_ref    =  0.000000000000000_wp
    real(wp),parameter :: g_born_ref  = -0.015702779304596_wp
    real(wp),parameter :: g_sasa_ref  =  0.000932775805422_wp
    real(wp),parameter :: g_hb_ref    = -0.006029170883118_wp
    real(wp),parameter :: g_shift_ref =  0.000364939254922_wp
    real(wp),parameter :: g_solv_ref  = -0.020434235127370_wp
!&>
    real(wp),parameter :: cthr = 1.0e-7_wp

    nat = testnat
    allocate(at(nat),xyz(3,nat),grad(3,nat))
    at = testat; xyz = testxyz
    call calculator%init(nat,at,xyz,ichrg=0,iostat=io,solvent='h2o')
    call check(error,io,0); if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calculator,energy,grad,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    associate(r => calculator%res)
      call check(error,r%e_bond, e_bond_ref, thr=cthr); if (allocated(error)) return
      call check(error,r%e_angl, e_angl_ref, thr=cthr); if (allocated(error)) return
      call check(error,r%e_tors, e_tors_ref, thr=cthr); if (allocated(error)) return
      call check(error,r%e_batm, e_batm_ref, thr=cthr); if (allocated(error)) return
      call check(error,r%e_rep,  e_rep_ref,  thr=cthr); if (allocated(error)) return
      call check(error,r%e_es,   e_es_ref,   thr=cthr); if (allocated(error)) return
      call check(error,r%e_disp, e_disp_ref, thr=cthr); if (allocated(error)) return
      call check(error,r%e_hb,   e_hb_ref,   thr=cthr); if (allocated(error)) return
      call check(error,r%e_xb,   e_xb_ref,   thr=cthr); if (allocated(error)) return
      call check(error,r%g_born, g_born_ref, thr=cthr); if (allocated(error)) return
      call check(error,r%g_sasa, g_sasa_ref, thr=cthr); if (allocated(error)) return
      call check(error,r%g_hb,   g_hb_ref,   thr=cthr); if (allocated(error)) return
      call check(error,r%g_shift,g_shift_ref,thr=cthr); if (allocated(error)) return
      call check(error,r%g_solv, g_solv_ref, thr=cthr)
    end associate
  end subroutine test_gfnff_components_alpb

!========================================================================================!

  subroutine test_gfnff_supermol(error)
    !***********************************************
    !* Regression test for the 226-atom supermol  *
    !* geometry. Tests scalability of the neighbor *
    !* list and topology routines.                 *
    !* Input: supermol.f90 (226 atoms, diverse     *
    !*        elements incl. halogens, P, B, Li)   *
    !***********************************************
    use supermol
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io,k
    real(wp) :: fnet(3)
    type(gfnff_data) :: calculator
!&<
    real(wp),parameter :: e_ref = -29.644066967343573_wp
!&>

    nat = testnat
    allocate(at(nat),xyz(3,nat),grad(3,nat))
    at = testat
    xyz = testxyz
    energy = 0.0_wp; grad = 0.0_wp

    call gfnff_initialize(nat,at,xyz,calculator,ichrg=0,iostat=io)
    call check(error,io,0); if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calculator,energy,grad,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    ! net force check
    fnet = 0.0_wp
    do k = 1,nat
      fnet(:) = fnet(:) + grad(:,k)
    end do
    if (any(abs(fnet) > 1.0e-10_wp)) then
      call test_failed(error,"Net force is not zero for supermol")
      write(*,'(a,3es12.4)') "  F_net =",fnet
      return
    end if

    call check(error,energy,e_ref,thr=1.0e-3_wp)
    if (allocated(error)) then
      call test_failed(error,"Supermol energy does not match reference")
      write(*,'(a,f25.15)') "  energy =",energy
    end if
  end subroutine test_gfnff_supermol

!========================================================================================!
!========================================================================================!
end module test_gfnff
