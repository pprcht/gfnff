module test_pbc_ini
  !> Tests for the PBC topology-setup routines wired up in gfnff_ini.f90.
  !> Each test calls gfnff_initialize with npbc=3 and a small periodic
  !> system, then inspects the internal topology / neighbour-list objects
  !> that were populated by the PBC branches.
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use iso_fortran_env,only:wp => real64
  use gfnff_interface
  implicit none
  private

  public :: collect_pbc_ini

!========================================================================================!
contains
!========================================================================================!

  subroutine collect_pbc_ini(testsuite)
    type(unittest_type),allocatable,intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("PBC numctr=27 for 3D box          ",test_pbc_numctr),   &
      new_unittest("PBC allocated                     ",test_pbc_alloc),    &
      new_unittest("PBC EEQ charges charge-conserved  ",test_pbc_charges),  &
      new_unittest("PBC cell volume correct            ",test_pbc_cell_vol), &
      new_unittest("PBC H2 bond detected              ",test_pbc_h2_bond)   &
    ]
  end subroutine collect_pbc_ini

!========================================================================================!
! Shared test geometry helpers
!========================================================================================!

  !> H2 in a 10 Bohr cubic box. Bond along x-axis (r ~ 1.4 Bohr ~ 0.74 A).
  subroutine make_h2_box(nat,at,xyz,lattice)
    integer,intent(out) :: nat
    integer,allocatable,intent(out) :: at(:)
    real(wp),allocatable,intent(out) :: xyz(:,:)
    real(wp),intent(out) :: lattice(3,3)
    nat = 2
    allocate(at(nat),xyz(3,nat))
    at = [1,1]
    xyz(:,1) = [0.0_wp, 0.0_wp, 0.0_wp]
    xyz(:,2) = [1.4_wp, 0.0_wp, 0.0_wp]
    lattice = 0.0_wp
    lattice(1,1) = 10.0_wp
    lattice(2,2) = 10.0_wp
    lattice(3,3) = 10.0_wp
  end subroutine make_h2_box

  !> H2O in a 12 Bohr cubic box.
  !> O-H bond ~ 1.814 Bohr (0.96 A), H-O-H ~ 104.5 deg.
  subroutine make_h2o_box(nat,at,xyz,lattice)
    integer,intent(out) :: nat
    integer,allocatable,intent(out) :: at(:)
    real(wp),allocatable,intent(out) :: xyz(:,:)
    real(wp),intent(out) :: lattice(3,3)
    real(wp),parameter :: rOH = 1.814_wp
    real(wp),parameter :: hoh = 104.5_wp * (3.14159265358979_wp/180.0_wp)
    nat = 3
    allocate(at(nat),xyz(3,nat))
    at = [8,1,1]
    xyz(:,1) = [0.0_wp, 0.0_wp, 0.0_wp]                 ! O
    xyz(:,2) = [rOH, 0.0_wp, 0.0_wp]                    ! H1
    xyz(:,3) = [rOH*cos(hoh), rOH*sin(hoh), 0.0_wp]     ! H2
    lattice = 0.0_wp
    lattice(1,1) = 12.0_wp
    lattice(2,2) = 12.0_wp
    lattice(3,3) = 12.0_wp
  end subroutine make_h2o_box

!========================================================================================!

  subroutine test_pbc_numctr(error)
    !***********************************************************
    !* After gfnff_initialize with npbc=3 and a cubic box the
    !* neighbour type must contain exactly 27 central cells (3^3).
    !***********************************************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: lattice(3,3)

    call make_h2_box(nat,at,xyz,lattice)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lattice,npbc=3,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    !> 3D PBC always includes at least the central 3x3x3 shell of cells
    call check(error,calc%neigh%numctr,27)
  end subroutine test_pbc_numctr

!========================================================================================!

  subroutine test_pbc_alloc(error)
    !***********************************************************
    !* After gfnff_initialize with npbc=3, the key data
    !* structures (topo, neigh, cell) must be allocated and
    !* cell%npbc must equal 3.
    !***********************************************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: lattice(3,3)

    call make_h2_box(nat,at,xyz,lattice)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lattice,npbc=3,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    if (.not.allocated(calc%topo)) then
      call test_failed(error,"topo is not allocated after PBC init")
      return
    end if
    if (.not.allocated(calc%neigh)) then
      call test_failed(error,"neigh is not allocated after PBC init")
      return
    end if
    if (.not.allocated(calc%cell)) then
      call test_failed(error,"cell is not allocated after PBC init")
      return
    end if
    call check(error,calc%cell%npbc,3)
  end subroutine test_pbc_alloc

!========================================================================================!

  subroutine test_pbc_charges(error)
    !***********************************************************
    !* EEQ charges from goedeckera_PBC must be charge-conserved:
    !* sum(topo%qa) must equal the total molecular charge (0).
    !***********************************************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: lattice(3,3),qtot
    real(wp),parameter :: thr = 1.0e-6_wp

    call make_h2o_box(nat,at,xyz,lattice)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lattice,npbc=3,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    if (.not.allocated(calc%topo%qa)) then
      call test_failed(error,"topo%qa is not allocated after PBC init")
      return
    end if
    qtot = sum(calc%topo%qa)
    if (abs(qtot) > thr) then
      call test_failed(error,"EEQ charges not conserved: sum(q) = "//to_str(qtot))
    end if
  end subroutine test_pbc_charges

!========================================================================================!

!========================================================================================!

  subroutine test_pbc_cell_vol(error)
    !***********************************************************
    !* The volume of a cubic 10 Bohr box must be 1000 Bohr^3.
    !* This is a pure math check on cell%init.
    !***********************************************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: lattice(3,3)
    real(wp),parameter :: vol_ref = 1000.0_wp
    real(wp),parameter :: thr = 1.0e-10_wp

    call make_h2_box(nat,at,xyz,lattice)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lattice,npbc=3,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    if (abs(calc%cell%volume - vol_ref) > thr) then
      call test_failed(error,"Cell volume wrong: got "//to_str(calc%cell%volume)// &
        &             ", expected "//to_str(vol_ref))
    end if
  end subroutine test_pbc_cell_vol

!========================================================================================!

  subroutine test_pbc_h2_bond(error)
    !***********************************************************
    !* GFN-FF topology of H2 in a 10 Bohr box must contain
    !* exactly one covalent bond.
    !***********************************************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: lattice(3,3)

    call make_h2_box(nat,at,xyz,lattice)

    call gfnff_initialize(nat,at,xyz,calc,lattice=lattice,npbc=3,iostat=io)
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,calc%neigh%nbond,1)
  end subroutine test_pbc_h2_bond

!========================================================================================!
! Internal helper: real -> character
!========================================================================================!
  pure function to_str(x) result(s)
    real(wp),intent(in) :: x
    character(len=32) :: s
    write(s,'(es13.6)') x
  end function to_str

end module test_pbc_ini
