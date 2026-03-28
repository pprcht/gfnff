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
    real(wp), parameter :: e_ref   = -1.457828799728746_wp
    real(wp), parameter :: thr_e   = 5.0e-4_wp
    real(wp), parameter :: thr_g   = 1.0e-6_wp
    real(wp), parameter :: g_ref(3,sio2nat) = reshape([ &
      & -0.033206725329150_wp,  0.029375024022352_wp,  0.032656696230473_wp, &
      &  0.041826968098463_wp,  0.001290069729135_wp,  0.001030072820178_wp, &
      &  0.007869272495832_wp,  0.037389211554370_wp, -0.004267315561024_wp, &
      & -0.011535386748665_wp, -0.033419621287930_wp, -0.001037559146122_wp, &
      & -0.037316483923135_wp,  0.009543402978026_wp,  0.005385295118687_wp, &
      & -0.010989440750150_wp,  0.042771819264084_wp, -0.034396977147170_wp, &
      &  0.035301802378919_wp, -0.035648432763542_wp,  0.051608926166113_wp, &
      & -0.013948798232784_wp,  0.001534493760553_wp,  0.001706129125887_wp, &
      &  0.021998792010669_wp, -0.052835967257049_wp, -0.052685267607023_wp  &
      & ], shape(g_ref))

    allocate(grad(3, sio2nat), source=0.0_wp)

    call gfnff_initialize(sio2nat, sio2at, sio2xyz, calc, &
      &                   lattice=sio2lattice, npbc=3, iostat=io, printlevel=3)
    call check(error, io, 0)
    if (allocated(error)) return

    call gfnff_singlepoint(sio2nat, sio2at, sio2xyz, calc, energy, grad, &
      &                    lattice=sio2lattice, iostat=io, printlevel=3)

    call calc%resultprint()

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
