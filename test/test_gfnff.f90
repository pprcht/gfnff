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
    new_unittest("EEQ atom-order invariance",test_eeq_atom_order), &
    new_unittest("GFN-FF singlepoint with ALPB   ",test_gfnff_alpb), &
!    new_unittest("GFN-FF OpenMP parallel SP      ",test_gfnff_openmp), &
    new_unittest("GFN-FF numerical gradient      ",test_gfnff_numgrad), &
    new_unittest("GFN-FF net force vanishes      ",test_gfnff_netforce), &
    new_unittest("GFN-FF translation invariance  ",test_gfnff_translation), &
    new_unittest("GFN-FF energy decomposition    ",test_gfnff_edecomp), &
    new_unittest("GFN-FF energy components       ",test_gfnff_components), &
    new_unittest("GFN-FF energy components ALPB  ",test_gfnff_components_alpb), &
    new_unittest("GFN-FF supermol singlepoint    ",test_gfnff_supermol), &
    new_unittest("conformer2020 matches at min   ",test_conformer_at_minimum), &
    new_unittest("conformer2020 bonds cannot break",test_conformer_no_dissoc), &
    new_unittest("conformer2020 numerical gradient",test_conformer_numgrad), &
    new_unittest("harmonic2020 bond term is live ",test_harmonic_bonds), &
    new_unittest("supplied graph reproduces setup",test_graph_roundtrip), &
    new_unittest("supplied graph rebuilds a soup ",test_graph_soup), &
    new_unittest("malformed graphs are rejected  ",test_graph_validation) &
    ]
!&>
  end subroutine collect_gfnff

!========================================================================================!


!> N-methylacetamide: parameter generation must follow atoms, not list positions.
subroutine test_eeq_atom_order(error)
   use gfnff_interface, only : gfnff_data, gfnff_initialize
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 12
   integer, parameter :: at(nat) = [6,6,8,7,6,1,1,1,1,1,1,1]
   ! Coordinates in bohr; the same molecule in each ordering.
   real(wp), parameter :: xyz(3,nat) = reshape([ &
     & 3.622312140781_wp, 0.584150709785_wp, 0.403823853883_wp, &
     & 0.841690132958_wp, 0.709612631882_wp, 0.875537937257_wp, &
     & 0.001690216911_wp, 2.103743880637_wp, 2.615851423658_wp, &
     & -0.857946784445_wp, -0.712426316667_wp, -0.678850260895_wp, &
     & -3.584295528389_wp, -0.636030709585_wp, -0.295866267463_wp, &
     & 4.066614256909_wp, 1.536668293196_wp, -1.411023959546_wp, &
     & 4.677338535329_wp, 1.531874376674_wp, 1.951448115885_wp, &
     & 4.231725963854_wp, -1.420334469809_wp, 0.299405229455_wp, &
     & -0.129625612682_wp, -1.837950193426_wp, -2.128976201128_wp, &
     & -4.051091096309_wp, -1.315129049540_wp, 1.633832596409_wp, &
     & -4.536218537284_wp, -1.872131986606_wp, -1.699163157357_wp, &
     & -4.282193687633_wp, 1.327952833458_wp, -0.537763562819_wp], [3,nat])
   integer, parameter :: orders(nat,3) = reshape([ &
      & 1,2,3,4,5,6,7,8,9,10,11,12, &
      & 12,11,10,9,8,7,6,5,4,3,2,1, &
      & 4,1,2,3,5,6,7,8,9,10,11,12], [nat,3])
   real(wp) :: reference(nat)
   integer :: pass, i, order(nat)
   integer :: io

   do pass = 1,size(orders,2)
      order = orders(:,pass)
      block
         type(gfnff_data) :: calc
         call gfnff_initialize(nat,at(order),xyz(:,order),calc,iostat=io)
         call check(error,io,0)
         if (allocated(error)) return
         if (pass == 1) then
            reference = calc%topo%gameeq
         else
            do i = 1,nat
               call check(error,calc%topo%gameeq(i),reference(order(i)),thr=1e-8_wp)
               if (allocated(error)) return
            end do
         end if
         call calc%deallocate()
      end block
   end do
end subroutine test_eeq_atom_order


  subroutine test_gfnff_sp(error)
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io,ichrg
    type(gfnff_data) :: calculator

!&<
    real(wp),parameter :: e_ref = -4.672769781691198_wp
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    &  0.005301519993516_wp,  0.000273969363910_wp,  0.000002235959529_wp, &
    &  0.008167320169453_wp, -0.008221105459789_wp, -0.000025576788636_wp, &
    & -0.003077691911525_wp, -0.009433778745828_wp,  0.000033248917316_wp, &
    &  0.009917107785739_wp,  0.008085206601650_wp, -0.000022035760389_wp, &
    & -0.015619889331091_wp, -0.026674405953360_wp,  0.000004837042250_wp, &
    &  0.014524114405198_wp, -0.001977000218303_wp,  0.000067168492967_wp, &
    &  0.006146941767799_wp,  0.009560778231154_wp, -0.000017509764841_wp, &
    & -0.008819906462314_wp, -0.001068258453272_wp, -0.000077995784271_wp, &
    & -0.000981945066646_wp,  0.014871861295735_wp,  0.000032975014178_wp, &
    & -0.006667316982828_wp,  0.007440566142231_wp, -0.000019279207903_wp, &
    &  0.012829096072163_wp, -0.012742592323866_wp, -0.000039626399637_wp, &
    & -0.023432160358615_wp,  0.020998275104046_wp, -0.000002417368967_wp, &
    & -0.001888398357184_wp, -0.003908177584031_wp, -0.000013745132793_wp, &
    & -0.003754858173439_wp,  0.003730236914171_wp, -0.000073759340703_wp, &
    &  0.000742842266713_wp,  0.003621905000371_wp,  0.000003807410713_wp, &
    &  0.001069342816459_wp, -0.000350659204274_wp,  0.003705358718107_wp, &
    &  0.001070966037791_wp, -0.000349869735091_wp, -0.003711548398877_wp, &
    & -0.002984386787343_wp,  0.000421365885265_wp,  0.000013800916223_wp, &
    &  0.004498452939135_wp,  0.000657641500367_wp,  0.000002341165750_wp, &
    &  0.000369369105619_wp,  0.001498995665630_wp,  0.003778129741047_wp, &
    &  0.000379302397946_wp,  0.001507626030698_wp, -0.003768510774416_wp, &
    & -0.001010609628770_wp, -0.004606598765027_wp,  0.000057154528590_wp, &
    &  0.001617376038978_wp, -0.001636926923752_wp,  0.003219282060532_wp, &
    &  0.001603411263247_wp, -0.001699054368637_wp, -0.003148335245769_wp &
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
    real(wp),parameter :: e_ref = -4.672769781691198_wp
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
    real(wp),parameter :: e_ref = -4.689881004521161_wp
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    &  0.005946208652649_wp,  0.000224811865074_wp,  0.000002359240085_wp, &
    &  0.008277407150046_wp, -0.008053219839783_wp, -0.000024497032410_wp, &
    & -0.002932707597764_wp, -0.008966714689328_wp,  0.000032565076289_wp, &
    &  0.009324379007395_wp,  0.007161906063066_wp, -0.000020616694963_wp, &
    & -0.015598834435094_wp, -0.026603386357037_wp,  0.000003395693783_wp, &
    &  0.014222288516666_wp, -0.002328604228673_wp,  0.000066592710130_wp, &
    &  0.006359466557208_wp,  0.009728337263790_wp, -0.000012695623544_wp, &
    & -0.008059111850361_wp, -0.001016829683581_wp, -0.000067760196060_wp, &
    & -0.000927818757908_wp,  0.014719058505107_wp,  0.000064283147700_wp, &
    & -0.007015977194191_wp,  0.007705142260456_wp, -0.000018834371539_wp, &
    &  0.012160923833761_wp, -0.012197711134565_wp, -0.000031516965546_wp, &
    & -0.023085262509599_wp,  0.020582146940539_wp,  0.000000986513042_wp, &
    & -0.002531399274647_wp, -0.004380039449959_wp, -0.000014407099030_wp, &
    & -0.003644260248161_wp,  0.004533871376230_wp, -0.000098969090192_wp, &
    &  0.000763559429034_wp,  0.003493655990505_wp,  0.000003660039402_wp, &
    &  0.001177998829209_wp, -0.000489910550508_wp,  0.003518627656981_wp, &
    &  0.001179477328088_wp, -0.000489381615418_wp, -0.003525552550843_wp, &
    & -0.002858162870343_wp, -0.000053440747795_wp,  0.000013389639796_wp, &
    &  0.004178983498707_wp,  0.000471382477721_wp,  0.000002347767546_wp, &
    &  0.000260877587770_wp,  0.001522493609407_wp,  0.003713766358189_wp, &
    &  0.000270276289200_wp,  0.001531589701024_wp, -0.003704228231022_wp, &
    & -0.001005734198409_wp, -0.004218433170956_wp,  0.000052754010299_wp, &
    &  0.001779381146527_wp, -0.001420317681741_wp,  0.003106702346878_wp, &
    &  0.001758041110219_wp, -0.001456406903575_wp, -0.003062352344973_wp &
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
    real(wp),parameter :: e_es_ref   = -0.174328618228843_wp
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
    real(wp),parameter :: e_es_ref    = -0.171008843183502_wp
    real(wp),parameter :: e_disp_ref  = -0.018113406607768_wp
    real(wp),parameter :: e_hb_ref    = -0.000000033949287_wp
    real(wp),parameter :: e_xb_ref    =  0.000000000000000_wp
    real(wp),parameter :: g_born_ref  = -0.015701428618838_wp
    real(wp),parameter :: g_sasa_ref  =  0.000932775805422_wp
    real(wp),parameter :: g_hb_ref    = -0.006027284316810_wp
    real(wp),parameter :: g_shift_ref =  0.000364939254922_wp
    real(wp),parameter :: g_solv_ref  = -0.020430997875304_wp
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

  subroutine conformer_setup(nat,at,xyz,calc_ref,calc_con,error)
    !> Caffeine with two calculators over the same geometry, one on the
    !> published parametrisation and one on conformer2020. Both topologies
    !> are built here, at the reference geometry, so that later tests can
    !> pull a bond apart without the bond list being redetermined underneath
    !> them -- which is exactly what an optimiser or an MD run sees.
    use coffeine
    integer,intent(out) :: nat
    integer,allocatable,intent(out) :: at(:)
    real(wp),allocatable,intent(out) :: xyz(:,:)
    type(gfnff_data),intent(out) :: calc_ref,calc_con
    type(error_type),allocatable,intent(out) :: error
    integer :: io

    nat = testnat
    allocate (at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz

    call gfnff_initialize(nat,at,xyz,calc_ref,ichrg=0,iostat=io,printlevel=0, &
       & version=gffVersion%angewChem2020_2)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_initialize(nat,at,xyz,calc_con,ichrg=0,iostat=io,printlevel=0, &
       & version=gffVersion%conformer2020)
    call check(error,io,0)
  end subroutine conformer_setup

!========================================================================================!

  subroutine test_conformer_at_minimum(error)
    !> Near equilibrium every bond sits well inside the region where the
    !> Gaussian well is convex, so conformer2020 has to reproduce the
    !> published force field. On one thread it does so bit for bit; with the
    !> OpenMP reductions summing in thread completion order either energy
    !> moves by a couple of ULP between runs, so the assertion is made at
    !> 1e-12 -- twelve orders below what the continuation branch is worth
    !> once it engages, and three above that noise.
    type(error_type),allocatable,intent(out) :: error
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),gref(:,:),gcon(:,:)
    real(wp) :: eref,econ
    type(gfnff_data) :: calc_ref,calc_con

    call conformer_setup(nat,at,xyz,calc_ref,calc_con,error)
    if (allocated(error)) return
    allocate (gref(3,nat),gcon(3,nat),source=0.0_wp)

    call gfnff_singlepoint(nat,at,xyz,calc_ref,eref,gref,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc_con,econ,gcon,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    if (abs(econ-eref) > 1.0e-12_wp) then
      call test_failed(error,"conformer2020 energy differs at the minimum")
      write (*,'(a,2es25.16)') "  E =",eref,econ
      return
    end if
    if (maxval(abs(gcon-gref)) > 1.0e-12_wp) then
      call test_failed(error,"conformer2020 gradient differs at the minimum")
      write (*,'(a,es12.4)') "  max|dg| =",maxval(abs(gcon-gref))
      return
    end if
  end subroutine test_conformer_at_minimum

!========================================================================================!

  subroutine test_conformer_no_dissoc(error)
    !> Pull one C-H well past the inflection point of the well. The Gaussian
    !> flattens out -- successive equal steps cost almost nothing, which is
    !> the dissociation this variant exists to forbid -- while conformer2020
    !> keeps climbing at a constant force, so equal steps cost equal energy.
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: iat = 15    !> a hydrogen bonded to atom 1
    real(wp),parameter :: step = 2.0_wp   !> bohr, per stretch
    integer :: nat,io,k
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),disp(:,:),grad(:,:)
    real(wp) :: eref(0:3),econ(0:3),dir(3),dcon,dref
    type(gfnff_data) :: calc_ref,calc_con

    call conformer_setup(nat,at,xyz,calc_ref,calc_con,error)
    if (allocated(error)) return
    allocate (grad(3,nat),disp(3,nat))

    dir = xyz(:,iat)-xyz(:,1)
    dir = dir/sqrt(sum(dir*dir))

    do k = 0,3
      disp = xyz
      disp(:,iat) = xyz(:,iat)+dir*step*real(k,wp)
      call gfnff_singlepoint(nat,at,disp,calc_ref,eref(k),grad,printlevel=0,iostat=io)
      call check(error,io,0)
      if (allocated(error)) return
      call gfnff_singlepoint(nat,at,disp,calc_con,econ(k),grad,printlevel=0,iostat=io)
      call check(error,io,0)
      if (allocated(error)) return
    end do

    !> the Gaussian has run out of restoring force: the last step is a small
    !> fraction of the one before it
    dref = (eref(3)-eref(2))/(eref(2)-eref(1))
    if (dref > 0.1_wp) then
      call test_failed(error,"reference bond term did not saturate; test is void")
      write (*,'(a,4f14.6)') "  E_ref  =",eref
      return
    end if

    !> the continuation is linear, so equal steps cost equal energy
    dcon = (econ(3)-econ(2))/(econ(2)-econ(1))
    if (dcon < 0.9_wp.or.dcon > 1.1_wp) then
      call test_failed(error,"conformer2020 bond term is not rising linearly")
      write (*,'(a,4f14.6)') "  E_con  =",econ
      write (*,'(a,f14.6)') "  ratio  =",dcon
      return
    end if

    !> and by then it is far above where the Gaussian gave up (Eh)
    if (econ(3)-eref(3) < 0.1_wp) then
      call test_failed(error,"conformer2020 did not outclimb the dissociating well")
      write (*,'(a,2f14.6)') "  E(3) =",eref(3),econ(3)
      return
    end if
  end subroutine test_conformer_no_dissoc

!========================================================================================!

  subroutine test_conformer_numgrad(error)
    !> The continuation branch carries its own derivative expressions, so it
    !> needs its own gradient check. Taken at a geometry where the branch is
    !> demonstrably active, which the energy difference against the published
    !> version establishes before anything else is asserted.
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: iat = 15
    real(wp),parameter :: stretch = 3.0_wp,step = 0.001_wp
    integer :: nat,io,i,j
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),disp(:,:),stencil(:,:),grad(:,:),gnum(:,:)
    real(wp),allocatable :: gtmp(:,:)
    real(wp) :: eref,econ,bw2,bw,fw,fw2,dir(3),dum
    type(gfnff_data) :: calc_ref,calc_con

    call conformer_setup(nat,at,xyz,calc_ref,calc_con,error)
    if (allocated(error)) return
    allocate (grad(3,nat),gnum(3,nat),gtmp(3,nat),disp(3,nat),stencil(3,nat))

    dir = xyz(:,iat)-xyz(:,1)
    dir = dir/sqrt(sum(dir*dir))
    disp = xyz
    disp(:,iat) = xyz(:,iat)+dir*stretch

    call gfnff_singlepoint(nat,at,disp,calc_ref,eref,grad,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,disp,calc_con,econ,grad,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    !> guard: without this the scan below would pass on the Gaussian branch
    if (abs(econ-eref) < 0.1_wp) then
      call test_failed(error,"continuation branch is not active; test is void")
      write (*,'(a,2f14.6)') "  E =",eref,econ
      return
    end if

    !> same five-point stencil as test_gfnff_numgrad; a plain central
    !> difference leaves a truncation error of its own that is larger than
    !> the agreement being asserted
    stencil = disp
    do i = 1,nat
      do j = 1,3
        stencil(j,i) = disp(j,i)-2.0_wp*step
        call gfnff_singlepoint(nat,at,stencil,calc_con,bw2,gtmp,printlevel=0,iostat=io)
        stencil(j,i) = disp(j,i)-step
        call gfnff_singlepoint(nat,at,stencil,calc_con,bw,gtmp,printlevel=0,iostat=io)
        stencil(j,i) = disp(j,i)+step
        call gfnff_singlepoint(nat,at,stencil,calc_con,fw,gtmp,printlevel=0,iostat=io)
        stencil(j,i) = disp(j,i)+2.0_wp*step
        call gfnff_singlepoint(nat,at,stencil,calc_con,fw2,gtmp,printlevel=0,iostat=io)
        stencil(j,i) = disp(j,i)
        gnum(j,i) = (bw2/12.0_wp-8.0_wp*bw/12.0_wp &
           &        +8.0_wp*fw/12.0_wp-fw2/12.0_wp)/step
      end do
    end do

    dum = maxval(abs(grad-gnum))
    if (dum > 1.0e-6_wp) then
      call test_failed(error,"conformer2020 analytic gradient disagrees with FD")
      write (*,'(a,es12.4)') "  max|dg| =",dum
      return
    end if
  end subroutine test_conformer_numgrad

!========================================================================================!

  subroutine perceived_graph(nat,calc,bmat)
    !> The adjacency matrix GFN-FF derived from the geometry, as a graph that
    !> can be handed straight back in.
    integer,intent(in) :: nat
    type(gfnff_data),intent(in) :: calc
    integer,intent(out) :: bmat(nat,nat)
    integer :: i,j,k

    bmat = 0
    do i = 1,nat
      do k = 1,calc%neigh%nb(calc%neigh%numnb,i,1)
        j = calc%neigh%nb(k,i,1)
        bmat(j,i) = 1
        bmat(i,j) = 1
      end do
    end do
  end subroutine perceived_graph

!========================================================================================!

  subroutine test_harmonic_bonds(error)
    !> harmonic2020 must actually evaluate its bond term.
    !>
    !> It read the bond list from topo%blist, which nothing has filled since
    !> the list moved onto the neighbour object. The loop ran zero times and
    !> the version silently degraded to its repulsion alone -- an energy that
    !> still differs from the default parametrisation, so every test that only
    !> asserted "harmonic2020 is different" passed throughout.
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    real(wp) :: energy
    type(gfnff_data) :: calc

    nat = testnat
    allocate (at(nat),xyz(3,nat),grad(3,nat))
    at = testat
    xyz = testxyz

    call gfnff_initialize(nat,at,xyz,calc,ichrg=0,printlevel=0,iostat=io, &
       & version=gffVersion%harmonic2020)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,calc,energy,grad,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    if (calc%neigh%nbond <= 0) then
      call test_failed(error,"harmonic2020 setup produced no bonds at all")
      return
    end if
    if (abs(calc%res%e_bond) < 1.0e-8_wp) then
      call test_failed(error,"harmonic2020 bond term is zero; the bond list is not reaching it")
      write (*,'(a,i0,a,es12.4)') "  nbond = ",calc%neigh%nbond,"  e_bond = ",calc%res%e_bond
      return
    end if
  end subroutine test_harmonic_bonds

!========================================================================================!

  subroutine test_graph_roundtrip(error)
    !> Handing GFN-FF back the graph it perceived must change nothing. This is
    !> what pins the supplied-graph path to the perception path: any drift in
    !> how the neighbour lists are packed shows up here as an energy shift.
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    integer :: nat,io
    integer,allocatable :: at(:),bmat(:,:)
    real(wp),allocatable :: xyz(:,:),g1(:,:),g2(:,:)
    real(wp) :: e1,e2
    type(gfnff_data) :: c1,c2

    nat = testnat
    allocate (at(nat),xyz(3,nat),g1(3,nat),g2(3,nat),bmat(nat,nat))
    at = testat
    xyz = testxyz

    call gfnff_initialize(nat,at,xyz,c1,ichrg=0,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,c1,e1,g1,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    call perceived_graph(nat,c1,bmat)

    allocate (c2%userinput)
    c2%userinput%bondmat = bmat
    call gfnff_initialize(nat,at,xyz,c2,ichrg=0,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return
    call gfnff_singlepoint(nat,at,xyz,c2,e2,g2,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,c2%neigh%nbond,c1%neigh%nbond)
    if (allocated(error)) then
      call test_failed(error,"supplied graph gave a different bond count")
      return
    end if
    !> 1e-12 rather than exact: the OpenMP reductions sum in thread
    !> completion order, which moves either energy by a couple of ULP
    if (abs(e2-e1) > 1.0e-12_wp) then
      call test_failed(error,"supplied graph changed the energy")
      write (*,'(a,2es25.16)') "  E =",e1,e2
      return
    end if
    if (maxval(abs(g2-g1)) > 1.0e-12_wp) then
      call test_failed(error,"supplied graph changed the gradient")
      write (*,'(a,es12.4)') "  max|dg| =",maxval(abs(g2-g1))
      return
    end if
  end subroutine test_graph_roundtrip

!========================================================================================!

  subroutine test_graph_soup(error)
    !> The application this exists for: coordinates that carry no structure at
    !> all, plus a graph, under harmonic2020. Bond perception would find a
    !> different molecule every time here, so the bond count coming out equal
    !> to the graph's edge count is the whole point.
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    integer :: nat,io,i,j,k
    integer,allocatable :: at(:),bmat(:,:)
    real(wp),allocatable :: xyz(:,:),soup(:,:),grad(:,:),gnum(:,:),gscr(:,:)
    real(wp),allocatable :: stencil(:,:)
    real(wp) :: energy,ep,em,dum
    real(wp),parameter :: step = 0.005_wp
    type(gfnff_data) :: cref,calc

    nat = testnat
    allocate (at(nat),xyz(3,nat),soup(3,nat),grad(3,nat),gnum(3,nat),gscr(3,nat))
    allocate (stencil(3,nat),bmat(nat,nat))
    at = testat
    xyz = testxyz

    call gfnff_initialize(nat,at,xyz,cref,ichrg=0,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return
    call perceived_graph(nat,cref,bmat)

    !> A deterministic soup: atoms on a coarse grid with a jitter, so no bond
    !> is anywhere near its length and the connectivity is unrecoverable from
    !> the geometry, but no two atoms sit on top of each other either -- a
    !> near-collision would make the finite differences below meaningless
    !> without saying anything about the graph path. No seed, so no flake.
    k = 0
    do i = 1,nat
      soup(1,i) = 5.0_wp*real(mod(k,3),wp)+0.4_wp*sin(real(i,wp))
      soup(2,i) = 5.0_wp*real(mod(k/3,3),wp)+0.4_wp*sin(real(2*i,wp))
      soup(3,i) = 5.0_wp*real(k/9,wp)+0.4_wp*sin(real(3*i,wp))
      k = k+1
    end do

    allocate (calc%userinput)
    calc%userinput%bondmat = bmat
    call gfnff_initialize(nat,at,soup,calc,ichrg=0,printlevel=0,iostat=io, &
       & version=gffVersion%harmonic2020)
    call check(error,io,0)
    if (allocated(error)) then
      call test_failed(error,"setup from a graph over scrambled coordinates failed")
      return
    end if

    call check(error,calc%neigh%nbond,count(bmat /= 0)/2)
    if (allocated(error)) then
      call test_failed(error,"bond count does not match the supplied graph")
      return
    end if

    call gfnff_singlepoint(nat,at,soup,calc,energy,grad,printlevel=0,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return
    if (abs(calc%res%e_bond) < 1.0e-8_wp) then
      call test_failed(error,"no bond energy on the scrambled geometry")
      return
    end if

    !> the skipped-perception path has its own arrangement of the setup, so
    !> its gradient is checked rather than assumed
    stencil = soup
    do i = 1,nat
      do j = 1,3
        stencil(j,i) = soup(j,i)+step
        call gfnff_singlepoint(nat,at,stencil,calc,ep,gscr,printlevel=0,iostat=io)
        stencil(j,i) = soup(j,i)-step
        call gfnff_singlepoint(nat,at,stencil,calc,em,gscr,printlevel=0,iostat=io)
        stencil(j,i) = soup(j,i)
        gnum(j,i) = 0.5_wp*(ep-em)/step
      end do
    end do
    dum = maxval(abs(grad-gnum))
    if (dum > 1.0e-6_wp) then
      call test_failed(error,"graph/harmonic gradient disagrees with finite differences")
      write (*,'(a,es12.4)') "  max|dg| =",dum
      return
    end if
  end subroutine test_graph_soup

!========================================================================================!

  subroutine test_graph_validation(error)
    !> Every rejected case is a graph the setup would otherwise consume
    !> silently and turn into a strange force field.
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    integer :: nat
    integer,allocatable :: at(:),bmat(:,:)
    real(wp),allocatable :: xyz(:,:)

    nat = testnat
    allocate (at(nat),xyz(3,nat),bmat(nat,nat))
    at = testat
    xyz = testxyz

    !> asymmetric: the two atoms disagree about the bond
    bmat = 0
    bmat(2,1) = 1
    call try(bmat,'asymmetric',error)
    if (allocated(error)) return

    !> an atom bonded to itself
    bmat = 0
    bmat(1,1) = 1
    call try(bmat,'self-bonded',error)
    if (allocated(error)) return

    !> negative bond order
    bmat = 0
    bmat(2,1) = -1
    bmat(1,2) = -1
    call try(bmat,'negative',error)
    if (allocated(error)) return

  contains
    subroutine try(bad,label,err)
      integer,intent(in) :: bad(:,:)
      character(len=*),intent(in) :: label
      type(error_type),allocatable,intent(out) :: err
      type(gfnff_data) :: c
      integer :: ios

      allocate (c%userinput)
      c%userinput%bondmat = bad
      call gfnff_initialize(nat,at,xyz,c,ichrg=0,printlevel=0,iostat=ios)
      if (ios == 0) then
        call test_failed(err,'a '//label//' molecular graph was accepted')
      end if
    end subroutine try
  end subroutine test_graph_validation

!========================================================================================!
!========================================================================================!
end module test_gfnff
