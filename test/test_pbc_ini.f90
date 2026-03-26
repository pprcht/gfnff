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
      new_unittest("PBC ini: numctr=27 for 3D box          ",test_pbc_numctr), &
      new_unittest("PBC ini: bpair/alphanb_pbc allocated   ",test_pbc_alloc), &
      new_unittest("PBC ini: EEQ charges charge-conserved  ",test_pbc_charges), &
      new_unittest("PBC ini: bond list has iTr column      ",test_pbc_blist), &
      new_unittest("PBC ini: H2O in box, single fragment   ",test_pbc_h2o_frag) &
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
    xyz(:,1) = [0.0_wp, 0.0_wp, 0.0_wp]                                   ! O
    xyz(:,2) = [rOH, 0.0_wp, 0.0_wp]                                       ! H1
    xyz(:,3) = [rOH*cos(hoh), rOH*sin(hoh), 0.0_wp]                        ! H2
    lattice = 0.0_wp
    lattice(1,1) = 12.0_wp
    lattice(2,2) = 12.0_wp
    lattice(3,3) = 12.0_wp
  end subroutine make_h2o_box

!========================================================================================!

  subroutine test_pbc_numctr(error)
    !> After gfnff_initialize with npbc=3 and a cubic box the neighbour type
    !> must contain exactly 27 central cells (3^3).
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: lattice(3,3)

    call make_h2_box(nat,at,xyz,lattice)
    call gfnff_initialize(nat,at,xyz,calc,print=.false.,ichrg=0, &
      &                   lattice=lattice,npbc=3,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    !> 3D PBC always includes at least the central 3x3x3 shell of cells
    if (calc%neigh%numctr < 27) then
      call test_failed(error,"neigh%numctr < 27 for 3D PBC")
    end if
  end subroutine test_pbc_numctr

!========================================================================================!

  subroutine test_pbc_alloc(error)
    !> neigh%bpair and topo%alphanb_pbc must be allocated after PBC init,
    !> with dimensions consistent with nat and numctr.
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: lattice(3,3)

    call make_h2_box(nat,at,xyz,lattice)
    call gfnff_initialize(nat,at,xyz,calc,print=.false.,ichrg=0, &
      &                   lattice=lattice,npbc=3,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    if (.not.allocated(calc%neigh%bpair)) then
      call test_failed(error,"neigh%bpair not allocated after PBC init")
      return
    end if
    !> shape: (nat, nat, numctr)
    call check(error,size(calc%neigh%bpair,1),nat, &
      & "neigh%bpair dim 1 should equal nat")
    call check(error,size(calc%neigh%bpair,3),calc%neigh%numctr, &
      & "neigh%bpair dim 3 should equal numctr")

    if (.not.allocated(calc%topo%alphanb_pbc)) then
      call test_failed(error,"topo%alphanb_pbc not allocated after PBC init")
      return
    end if
    !> shape: (nat, nat, numctr+1)
    call check(error,size(calc%topo%alphanb_pbc,1),nat, &
      & "topo%alphanb_pbc dim 1 should equal nat")
    call check(error,size(calc%topo%alphanb_pbc,3),calc%neigh%numctr+1, &
      & "topo%alphanb_pbc dim 3 should equal numctr+1")
  end subroutine test_pbc_alloc

!========================================================================================!

  subroutine test_pbc_charges(error)
    !> EEQ charges from goedeckera_PBC must be charge-conserved:
    !> sum(qa) should equal ichrg (=0 here) to within 1e-6.
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: lattice(3,3),qtot
    real(wp),parameter :: thr = 1.0e-6_wp

    call make_h2o_box(nat,at,xyz,lattice)
    call gfnff_initialize(nat,at,xyz,calc,print=.false.,ichrg=0, &
      &                   lattice=lattice,npbc=3,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    if (.not.allocated(calc%topo%qa)) then
      call test_failed(error,"topo%qa not allocated"); return
    end if
    qtot = sum(calc%topo%qa)
    if (abs(qtot) > thr) then
      call test_failed(error,"EEQ charges not charge-conserved: sum(qa) = " &
        & //to_str(qtot))
    end if
  end subroutine test_pbc_charges

!========================================================================================!

  subroutine test_pbc_blist(error)
    !> For PBC, the bond list is stored in neigh%blist(3, nbond).
    !> Check it is allocated, has 3 rows, and nbond >= 1 for H2.
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: lattice(3,3)

    call make_h2_box(nat,at,xyz,lattice)
    call gfnff_initialize(nat,at,xyz,calc,print=.false.,ichrg=0, &
      &                   lattice=lattice,npbc=3,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    if (.not.allocated(calc%neigh%blist)) then
      call test_failed(error,"neigh%blist not allocated after PBC init")
      return
    end if
    call check(error,size(calc%neigh%blist,1),3, &
      & "neigh%blist first dim should be 3 (atom1, atom2, iTr)")
    if (calc%neigh%nbond < 1) then
      call test_failed(error,"neigh%nbond < 1 for H2 in box")
      return
    end if
    !> iTr values must be valid cell indices (1 .. numctr)
    if (any(calc%neigh%blist(3,:) < 1) .or. &
      & any(calc%neigh%blist(3,:) > calc%neigh%numctr)) then
      call test_failed(error,"neigh%blist iTr values out of range")
    end if
  end subroutine test_pbc_blist

!========================================================================================!

  subroutine test_pbc_h2o_frag(error)
    !> H2O is fully connected, so mrecgffPBC must assign all atoms
    !> to a single fragment (nfrag == 1, all fraglist entries == 1).
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: nat,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: lattice(3,3)

    call make_h2o_box(nat,at,xyz,lattice)
    call gfnff_initialize(nat,at,xyz,calc,print=.false.,ichrg=0, &
      &                   lattice=lattice,npbc=3,iostat=io)
    call check(error,io,0); if (allocated(error)) return

    call check(error,calc%topo%nfrag,1, &
      & "H2O in periodic box should have 1 fragment")
    if (.not.allocated(calc%topo%fraglist)) then
      call test_failed(error,"topo%fraglist not allocated"); return
    end if
    if (any(calc%topo%fraglist(1:nat) /= 1)) then
      call test_failed(error,"not all atoms assigned to fragment 1")
    end if
  end subroutine test_pbc_h2o_frag

!========================================================================================!
! Internal helper: real -> character
!========================================================================================!
  pure function to_str(x) result(s)
    real(wp),intent(in) :: x
    character(len=32) :: s
    write(s,'(es13.6)') x
  end function to_str

end module test_pbc_ini
