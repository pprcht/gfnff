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
#ifdef WITH_GBSA
    new_unittest("GFN-FF singlepoint with ALPB   ",test_gfnff_alpb), & 
#else
    new_unittest("GFN-FF singlepoint with ALPB   ",test_gfnff_alpb,should_fail=.true.), &   
#endif
!    new_unittest("GFN-FF OpenMP parallel SP      ",test_gfnff_openmp), &
    new_unittest("GFN-FF numerical gradient      ",test_gfnff_numgrad) &
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
    logical :: pr
    character(len=:),allocatable :: alpbsolvent
    type(gfnff_data) :: ffdata

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
    pr = .false.
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)

    !> calculation
    call gfnff_initialize(nat,at,xyz,ffdata,print=pr,ichrg=ichrg,iostat=io)
    call gfnff_singlepoint(nat,at,xyz,ffdata,energy,grad,pr,iostat=io)
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
    logical :: pr
    type(gfnff_data) :: ffdata
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
    pr = .false.
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)
    allocate (gradient(3,nat),g_ref(3,nat),stencil(3,nat),source=0.0_wp)

    !> calculation
    call gfnff_initialize(nat,at,xyz,ffdata,print=pr,ichrg=ichrg,iostat=io)
    call gfnff_singlepoint(nat,at,xyz,ffdata,energy,grad,pr,iostat=io)
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
        call gfnff_singlepoint(nat,at,stencil,ffdata,bw2,gradient,pr,iostat=io)
        stencil(j,i) = xyz(j,i)
        stencil(j,i) = stencil(j,i)-1.0_wp*step
        call gfnff_singlepoint(nat,at,stencil,ffdata,bw,gradient,pr,iostat=io)
        stencil(j,i) = xyz(j,i)
        stencil(j,i) = stencil(j,i)+1.0_wp*step
        call gfnff_singlepoint(nat,at,stencil,ffdata,fw,gradient,pr,iostat=io)
        stencil(j,i) = xyz(j,i)
        stencil(j,i) = stencil(j,i)+2.0_wp*step
        call gfnff_singlepoint(nat,at,stencil,ffdata,fw2,gradient,pr,iostat=io)
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
    logical :: pr
    character(len=:),allocatable :: alpbsolvent
    type(gfnff_data) :: ffdata

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
    alpbsolvent='h2o'
    pr = .false.
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)

    !> calculation
    call ffdata%init(nat,at,xyz,print=pr,ichrg=ichrg,iostat=io,solvent=alpbsolvent)
    call gfnff_singlepoint(nat,at,xyz,ffdata,energy,grad,pr,iostat=io)
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

  subroutine test_gfnff_openmp(error)
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io,i,j,ntimes,ichrg
    logical :: pr
    type(gfnff_data) :: ffdata
    type(gfnff_timer) :: timer
    character(len=40) :: atmp
    logical :: speedup
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
    ichrg=0
    pr=.false.
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)

#ifdef WITH_OpenMP

    ntimes = 2
    call timer%new(ntimes,.true.)
    do i = 1,ntimes
      call OMP_Set_Num_Threads(i)
#ifdef WITH_MKL
      call MKL_Set_Num_Threads(i)
#endif
      call ompprint_intern(atmp)

      call timer%measure(i,atmp)

      !> calculation
      grad = 0.0_wp
      do j = 1,20
        call gfnff_initialize(nat,at,xyz,ffdata,print=pr,ichrg=ichrg,iostat=io)
        call gfnff_singlepoint(nat,at,xyz,ffdata,energy,grad,pr,iostat=io)
        call check(error,io,0)
        if (allocated(error)) return

        call check(error,energy,e_ref,thr=1e-7_wp)
        if (allocated(error)) return

        if (any(abs(grad-g_ref) > thr2)) then
          call test_failed(error,"Gradient does not match reference")
          print'(3es21.14)',grad
          print'("---")'
          print'(3es21.14)',g_ref
          print'("---")'
          print'(3es21.14)',grad-g_ref
        end if
      end do

      call timer%measure(i)
    end do

    call OMP_Set_Num_Threads(1)
#ifdef WITH_MKL
    call MKL_Set_Num_Threads(1)
#endif

    !> check for any noticable speedup
    speedup = 1.1_wp < (timer%get(1)/timer%get(2))
#else /* WITH_OpenMP */
    speedup = .true.
#endif /* WITH_OpenMP */
    call check(error,speedup,.true.)
    if (allocated(error)) return

    deallocate (grad)
  end subroutine test_gfnff_openmp

  subroutine ompprint_intern(str)
    use omp_lib
    implicit none
    integer :: nproc,TID
    character(len=*) :: str
!$OMP PARALLEL PRIVATE(TID)
    TID = OMP_GET_THREAD_NUM()
    IF (TID .EQ. 0) THEN
      nproc = OMP_GET_NUM_THREADS()
#ifdef WITH_MKL
      write (str,'(a,i3)') 'OMP/MKL threads = ',nproc
#else
      write (str,'(a,i3)') 'OMP threads = ',nproc
#endif
    END IF
!$OMP END PARALLEL
  end subroutine ompprint_intern

!========================================================================================!
!========================================================================================!
end module test_gfnff
