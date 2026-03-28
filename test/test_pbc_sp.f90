module test_pbc_sp
  !> Singlepoint energy and gradient tests for periodic GFN-FF calculations.
  !> Uses the SiO2 (alpha-quartz-like) unit cell as reference geometry.
  use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
  use iso_fortran_env, only: wp => real64
  use gfnff_interface
  use sio2
  implicit none
  private

  public :: collect_pbc_sp

!========================================================================================!
contains
!========================================================================================!

  subroutine collect_pbc_sp(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("PBC SiO2 singlepoint              ", test_pbc_sio2_sp),      &
      new_unittest("PBC SiO2 net force vanishes       ", test_pbc_sio2_netforce)  &
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

    !call calc%resultprint()

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
    !* No reference data needed; purely a physical constraint.
    !***********************************************************
    type(error_type), allocatable, intent(out) :: error
    type(gfnff_data) :: calc
    real(wp) :: energy
    real(wp), allocatable :: grad(:,:)
    real(wp) :: fnet(3)
    integer :: io, k

    allocate(grad(3, sio2nat), source=0.0_wp)

    call gfnff_initialize(sio2nat, sio2at, sio2xyz, calc, &
      &                   lattice=sio2lattice, npbc=3, iostat=io)
    call check(error, io, 0)
    if (allocated(error)) return

    call gfnff_singlepoint(sio2nat, sio2at, sio2xyz, calc, energy, grad, &
      &                    lattice=sio2lattice, iostat=io)
    call check(error, io, 0)
    if (allocated(error)) return

    fnet = 0.0_wp
    do k = 1, sio2nat
      fnet(:) = fnet(:) + grad(:,k)
    end do
    if (any(abs(fnet) > 1.0e-10_wp)) then
      call test_failed(error, "PBC net force is not zero")
      write(*,'(a,3es12.4)') "  F_net =", fnet
    end if

  end subroutine test_pbc_sio2_netforce

end module test_pbc_sp
