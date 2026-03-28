module test_pbc_sp
  !> Singlepoint energy and gradient tests for periodic GFN-FF calculations.
  !> Uses the SiO2 (alpha-quartz-like) unit cell as reference geometry,
  !> caffeine in a periodic box, and SiO2 treated as a molecular system.
  use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
  use iso_fortran_env, only: wp => real64
  use gfnff_interface
  use sio2
  use coffeine, only: testnat, testat, testxyz
  implicit none
  private

  public :: collect_pbc_sp

!========================================================================================!
contains
!========================================================================================!

  subroutine collect_pbc_sp(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("PBC SiO2 singlepoint              ", test_pbc_sio2_sp),         &
      new_unittest("PBC SiO2 net force vanishes       ", test_pbc_sio2_netforce),    &
      new_unittest("PBC SiO2 numerical stress         ", test_pbc_sio2_stress),      &
      new_unittest("PBC caffeine-in-box singlepoint   ", test_pbc_caffeine_sp),      &
      new_unittest("PBC caffeine-in-box net force     ", test_pbc_caffeine_netforce), &
      new_unittest("SiO2 molecular singlepoint        ", test_sio2_molecular_sp)     &
    ]
  end subroutine collect_pbc_sp

!========================================================================================!

  subroutine test_pbc_sio2_sp(error)
    !***********************************************************
    !* GFN-FF singlepoint on the SiO2 unit cell.
    !* Checks that the calculation completes without error and
    !* that energy and gradient match stored reference values.
    !*
    !* INPUT: sio2nat, sio2at, sio2xyz, sio2lattice (from sio2.f90)
    !***********************************************************
    type(error_type), allocatable, intent(out) :: error
    type(gfnff_data) :: calc
    real(wp) :: energy
    real(wp), allocatable :: grad(:,:)
    integer :: io
    real(wp), parameter :: e_ref   = -1.429601370942344_wp
    real(wp), parameter :: thr_e   = 5.0e-4_wp
    real(wp), parameter :: thr_g   = 1.0e-6_wp
    real(wp), parameter :: g_ref(3,sio2nat) = reshape([ &
      & -0.020124023725345_wp,  0.034805181856585_wp,  0.044053404097327_wp, &
      &  0.034842303598706_wp, -0.016298860703557_wp, -0.006525861781396_wp, &
      &  0.005292624932896_wp,  0.023050948861036_wp,  0.002844588646397_wp, &
      & -0.003300742075965_wp, -0.038334744760994_wp,  0.006515122959067_wp, &
      & -0.022620929955122_wp,  0.006952573658589_wp, -0.002822197156296_wp, &
      & -0.020078721125386_wp,  0.034827045716319_wp, -0.044050926778701_wp, &
      &  0.008961279550369_wp, -0.043885191312389_wp,  0.027043358993217_wp, &
      & -0.016503745440319_wp,  0.028573301106312_wp, -0.000001074806692_wp, &
      &  0.033531954240166_wp, -0.029690254421901_wp, -0.027056414172924_wp  &
      & ], shape(g_ref))

    allocate(grad(3, sio2nat), source=0.0_wp)

    call gfnff_initialize(sio2nat, sio2at, sio2xyz, calc, &
      &                   lattice=sio2lattice, npbc=3, iostat=io, printlevel=0)
    call check(error, io, 0)
    if (allocated(error)) return

    call gfnff_singlepoint(sio2nat, sio2at, sio2xyz, calc, energy, grad, &
      &                    lattice=sio2lattice, iostat=io, printlevel=0)
    call check(error, io, 0)
    if (allocated(error)) return

    call check(error, energy, e_ref, thr=thr_e)
    if (allocated(error)) return

    if (any(abs(grad - g_ref) > thr_g)) then
      call test_failed(error, "PBC SiO2 gradient does not match reference")
    end if

  end subroutine test_pbc_sio2_sp

!========================================================================================!

  subroutine test_pbc_sio2_netforce(error)
    !***********************************************************
    !* For a fully periodic system the net force (sum over all
    !* atomic gradients) must vanish: Newton's 3rd law.
    !***********************************************************
    type(error_type), allocatable, intent(out) :: error
    type(gfnff_data) :: calc
    real(wp) :: energy
    real(wp), allocatable :: grad(:,:)
    real(wp) :: netforce(3)
    integer :: io
    real(wp), parameter :: thr = 1.0e-10_wp

    allocate(grad(3, sio2nat), source=0.0_wp)

    call gfnff_initialize(sio2nat, sio2at, sio2xyz, calc, &
      &                   lattice=sio2lattice, npbc=3, iostat=io, printlevel=0)
    call check(error, io, 0)
    if (allocated(error)) return

    call gfnff_singlepoint(sio2nat, sio2at, sio2xyz, calc, energy, grad, &
      &                    lattice=sio2lattice, iostat=io, printlevel=0)
    call check(error, io, 0)
    if (allocated(error)) return

    netforce = sum(grad, dim=2)
    if (any(abs(netforce) > thr)) then
      call test_failed(error, "PBC SiO2 net force does not vanish")
    end if

  end subroutine test_pbc_sio2_netforce

!========================================================================================!

  subroutine test_pbc_sio2_stress(error)
    !***********************************************************
    !* Validate the GFN-FF analytical stress tensor for the
    !* SiO2 unit cell against a central-difference numerical
    !* stress.
    !*
    !* Deformation: F = I + h*e_ij, applied as
    !*   L'   = F · L  (lattice)
    !*   xyz' = F · xyz (positions, keeps fractional coords)
    !* Numerical stress: sigma_num(i,j) = (E(+h) - E(-h)) / (2h)
    !* should equal the analytical sigma(i,j) returned by
    !* gfnff_singlepoint (which equals dE/d_epsilon).
    !***********************************************************
    type(error_type), allocatable, intent(out) :: error
    type(gfnff_data) :: calc, calc_p, calc_m
    real(wp) :: energy, energy_p, energy_m
    real(wp), allocatable :: grad(:,:)
    real(wp) :: sigma(3,3), sigma_num(3,3)
    real(wp) :: lattice_p(3,3), lattice_m(3,3)
    real(wp) :: xyz_p(3,sio2nat), xyz_m(3,sio2nat)
    integer :: io, ii, jj
    real(wp), parameter :: h   = 1.0e-4_wp
    real(wp), parameter :: thr = 5.0e-5_wp

    allocate(grad(3, sio2nat), source=0.0_wp)

    !> Reference init and analytical stress
    call gfnff_initialize(sio2nat, sio2at, sio2xyz, calc, &
      &                   lattice=sio2lattice, npbc=3, iostat=io, printlevel=0)
    call check(error, io, 0)
    if (allocated(error)) return

    call gfnff_singlepoint(sio2nat, sio2at, sio2xyz, calc, energy, grad, &
      &                    lattice=sio2lattice, iostat=io, printlevel=0, sigma=sigma)
    call check(error, io, 0)
    if (allocated(error)) return

    !> Numerical stress: central differences over all 9 strain components
    sigma_num = 0.0_wp
    do ii = 1, 3
      do jj = 1, 3
        !> Forward: L' = (I + h*e_ij)·L, xyz' = (I + h*e_ij)·xyz
        lattice_p = sio2lattice
        lattice_p(ii,:) = sio2lattice(ii,:) + h * sio2lattice(jj,:)
        xyz_p = sio2xyz
        xyz_p(ii,:) = sio2xyz(ii,:) + h * sio2xyz(jj,:)

        call gfnff_initialize(sio2nat, sio2at, xyz_p, calc_p, &
          &                   lattice=lattice_p, npbc=3, iostat=io, printlevel=0)
        call check(error, io, 0)
        if (allocated(error)) return
        call gfnff_singlepoint(sio2nat, sio2at, xyz_p, calc_p, energy_p, grad, &
          &                    lattice=lattice_p, iostat=io, printlevel=0)
        call check(error, io, 0)
        if (allocated(error)) return

        !> Backward: L' = (I - h*e_ij)·L, xyz' = (I - h*e_ij)·xyz
        lattice_m = sio2lattice
        lattice_m(ii,:) = sio2lattice(ii,:) - h * sio2lattice(jj,:)
        xyz_m = sio2xyz
        xyz_m(ii,:) = sio2xyz(ii,:) - h * sio2xyz(jj,:)

        call gfnff_initialize(sio2nat, sio2at, xyz_m, calc_m, &
          &                   lattice=lattice_m, npbc=3, iostat=io, printlevel=0)
        call check(error, io, 0)
        if (allocated(error)) return
        call gfnff_singlepoint(sio2nat, sio2at, xyz_m, calc_m, energy_m, grad, &
          &                    lattice=lattice_m, iostat=io, printlevel=0)
        call check(error, io, 0)
        if (allocated(error)) return

        sigma_num(ii,jj) = (energy_p - energy_m) / (2.0_wp * h)
      end do
    end do

    if (any(abs(sigma - sigma_num) > thr)) then
      call test_failed(error, "PBC SiO2 analytical stress does not match numerical")
    end if

  end subroutine test_pbc_sio2_stress

!========================================================================================!

  subroutine test_pbc_caffeine_sp(error)
    !***********************************************************
    !* GFN-FF singlepoint on caffeine placed in a cubic periodic
    !* box of 30 Bohr side length. Tests the PBC path for a
    !* molecular system with isolated periodic images.
    !*
    !* INPUT: testnat, testat, testxyz (from coffeine.f90)
    !***********************************************************
    type(error_type), allocatable, intent(out) :: error
    type(gfnff_data) :: calc
    real(wp) :: energy
    real(wp), allocatable :: grad(:,:)
    integer :: io
    real(wp), parameter :: e_ref   = -4.672824496568267_wp
    real(wp), parameter :: thr_e   = 5.0e-4_wp
    real(wp), parameter :: thr_g   = 1.0e-6_wp
    real(wp), parameter :: lattice(3,3) = reshape([ &
      & 30.0_wp,  0.0_wp,  0.0_wp, &
      &  0.0_wp, 30.0_wp,  0.0_wp, &
      &  0.0_wp,  0.0_wp, 30.0_wp  &
      & ], shape(lattice))
    real(wp), parameter :: g_ref(3,testnat) = reshape([ &
      &  0.005305429354695_wp,  0.000273396637589_wp,  0.000002229729880_wp, &
      &  0.008155749606608_wp, -0.008218382343138_wp, -0.000025504273643_wp, &
      & -0.003078364060582_wp, -0.009433691578996_wp,  0.000033226239279_wp, &
      &  0.009929561601486_wp,  0.008103015554918_wp, -0.000021941572236_wp, &
      & -0.015628983296790_wp, -0.026669657396164_wp,  0.000004731127436_wp, &
      &  0.014522455162064_wp, -0.001982312220583_wp,  0.000067144858589_wp, &
      &  0.006148113730321_wp,  0.009559864942358_wp, -0.000017842723784_wp, &
      & -0.008827636054556_wp, -0.001056969241862_wp, -0.000077441164087_wp, &
      & -0.000983382358722_wp,  0.014874458482399_wp,  0.000033263280336_wp, &
      & -0.006685192929466_wp,  0.007423923387226_wp, -0.000019596383842_wp, &
      &  0.012840742698299_wp, -0.012721757230174_wp, -0.000039191096717_wp, &
      & -0.023424096559047_wp,  0.021005920324322_wp, -0.000002294976987_wp, &
      & -0.001886034094539_wp, -0.003907249015970_wp, -0.000013743584535_wp, &
      & -0.003755992688866_wp,  0.003724559684310_wp, -0.000073722344180_wp, &
      &  0.000741267612391_wp,  0.003617501633303_wp,  0.000003797334570_wp, &
      &  0.001068019094283_wp, -0.000349341255361_wp,  0.003701790540580_wp, &
      &  0.001069684892227_wp, -0.000348622841591_wp, -0.003708072808293_wp, &
      & -0.002985789172337_wp,  0.000409758991531_wp,  0.000013809377123_wp, &
      &  0.004501205911768_wp,  0.000661479205529_wp,  0.000002293380532_wp, &
      &  0.000375807136583_wp,  0.001497353387914_wp,  0.003775034724670_wp, &
      &  0.000385692252196_wp,  0.001505946007675_wp, -0.003765490708356_wp, &
      & -0.001022525576931_wp, -0.004614285144321_wp,  0.000056935694156_wp, &
      &  0.001624072743271_wp, -0.001646407629381_wp,  0.003223374554169_wp, &
      &  0.001610194995644_wp, -0.001708502341534_wp, -0.003152789204660_wp  &
      & ], shape(g_ref))

    allocate(grad(3, testnat), source=0.0_wp)

    call gfnff_initialize(testnat, testat, testxyz, calc, &
      &                   lattice=lattice, npbc=3, iostat=io, printlevel=0)
    call check(error, io, 0)
    if (allocated(error)) return

    call gfnff_singlepoint(testnat, testat, testxyz, calc, energy, grad, &
      &                    lattice=lattice, iostat=io, printlevel=0)
    call check(error, io, 0)
    if (allocated(error)) return

    call check(error, energy, e_ref, thr=thr_e)
    if (allocated(error)) return

    if (any(abs(grad - g_ref) > thr_g)) then
      call test_failed(error, "PBC caffeine-in-box gradient does not match reference")
    end if

  end subroutine test_pbc_caffeine_sp

!========================================================================================!

  subroutine test_pbc_caffeine_netforce(error)
    !***********************************************************
    !* Net force on caffeine in a periodic box must vanish.
    !***********************************************************
    type(error_type), allocatable, intent(out) :: error
    type(gfnff_data) :: calc
    real(wp) :: energy
    real(wp), allocatable :: grad(:,:)
    real(wp) :: netforce(3)
    integer :: io
    real(wp), parameter :: thr = 1.0e-10_wp
    real(wp), parameter :: lattice(3,3) = reshape([ &
      & 30.0_wp,  0.0_wp,  0.0_wp, &
      &  0.0_wp, 30.0_wp,  0.0_wp, &
      &  0.0_wp,  0.0_wp, 30.0_wp  &
      & ], shape(lattice))

    allocate(grad(3, testnat), source=0.0_wp)

    call gfnff_initialize(testnat, testat, testxyz, calc, &
      &                   lattice=lattice, npbc=3, iostat=io, printlevel=0)
    call check(error, io, 0)
    if (allocated(error)) return

    call gfnff_singlepoint(testnat, testat, testxyz, calc, energy, grad, &
      &                    lattice=lattice, iostat=io, printlevel=0)
    call check(error, io, 0)
    if (allocated(error)) return

    netforce = sum(grad, dim=2)
    if (any(abs(netforce) > thr)) then
      call test_failed(error, "PBC caffeine net force does not vanish")
    end if

  end subroutine test_pbc_caffeine_netforce

!========================================================================================!

  subroutine test_sio2_molecular_sp(error)
    !***********************************************************
    !* GFN-FF singlepoint on the SiO2 unit-cell geometry
    !* treated as a molecular (non-periodic) system.
    !* Same atomic positions as the PBC test; no lattice.
    !* Reference energy ~-0.80221337847 Eh.
    !***********************************************************
    type(error_type), allocatable, intent(out) :: error
    type(gfnff_data) :: calc
    real(wp) :: energy
    real(wp), allocatable :: grad(:,:)
    integer :: io
    real(wp), parameter :: e_ref = -0.80221337847_wp
    real(wp), parameter :: thr_e = 5.0e-4_wp

    allocate(grad(3, sio2nat), source=0.0_wp)

    call gfnff_initialize(sio2nat, sio2at, sio2xyz, calc, &
      &                   iostat=io, printlevel=0)
    call check(error, io, 0)
    if (allocated(error)) return

    call gfnff_singlepoint(sio2nat, sio2at, sio2xyz, calc, energy, grad, &
      &                    iostat=io, printlevel=0)
    call check(error, io, 0)
    if (allocated(error)) return

    call check(error, energy, e_ref, thr=thr_e)

  end subroutine test_sio2_molecular_sp

!========================================================================================!
end module test_pbc_sp
