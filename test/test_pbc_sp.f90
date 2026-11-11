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
    real(wp), parameter :: e_ref   = -4.67280176604433706_wp
    real(wp), parameter :: thr_e   = 5.0e-4_wp
    real(wp), parameter :: thr_g   = 1.0e-6_wp
    real(wp), parameter :: lattice(3,3) = reshape([ &
      & 30.0_wp,  0.0_wp,  0.0_wp, &
      &  0.0_wp, 30.0_wp,  0.0_wp, &
      &  0.0_wp,  0.0_wp, 30.0_wp  &
      & ], shape(lattice))
    real(wp), parameter :: g_ref(3,testnat) = reshape([ &
      &  0.005305377427724_wp,  0.000273395782785_wp,  0.000002229722696_wp, &
      &  0.008157032510193_wp, -0.008218649198900_wp, -0.000025503644124_wp, &
      & -0.003077730337394_wp, -0.009434474669231_wp,  0.000033226205496_wp, &
      &  0.009926846738339_wp,  0.008101598000341_wp, -0.000021941711851_wp, &
      & -0.015616284554338_wp, -0.026671673585643_wp,  0.000004730748402_wp, &
      &  0.014520926307149_wp, -0.001982467157731_wp,  0.000067144649577_wp, &
      &  0.006148410697524_wp,  0.009559568090646_wp, -0.000017847057486_wp, &
      & -0.008826680965292_wp, -0.001056593341812_wp, -0.000077436146745_wp, &
      & -0.000981974778039_wp,  0.014872734384794_wp,  0.000033262218018_wp, &
      & -0.006669473361915_wp,  0.007441652678458_wp, -0.000019653928679_wp, &
      &  0.012830561234821_wp, -0.012721348887888_wp, -0.000039174380689_wp, &
      & -0.023433568654157_wp,  0.020998328564586_wp, -0.000002253415662_wp, &
      & -0.001890394239475_wp, -0.003908800488461_wp, -0.000013742411012_wp, &
      & -0.003755878438731_wp,  0.003724571742672_wp, -0.000073722410850_wp, &
      &  0.000741274365157_wp,  0.003617540667489_wp,  0.000003797336257_wp, &
      &  0.001068054981353_wp, -0.000349428373551_wp,  0.003701882040882_wp, &
      &  0.001069721024351_wp, -0.000348710051849_wp, -0.003708164619712_wp, &
      & -0.002985726178544_wp,  0.000409881155797_wp,  0.000013809384700_wp, &
      &  0.004500375000587_wp,  0.000658663659967_wp,  0.000002291499468_wp, &
      &  0.000373786541839_wp,  0.001497367304434_wp,  0.003776583990462_wp, &
      &  0.000383671083733_wp,  0.001505961774983_wp, -0.003767039159361_wp, &
      & -0.001022661203512_wp, -0.004614177783194_wp,  0.000056935194673_wp, &
      &  0.001624106766305_wp, -0.001646421186259_wp,  0.003223553524210_wp, &
      &  0.001610228032324_wp, -0.001708519082431_wp, -0.003152967628671_wp &
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
