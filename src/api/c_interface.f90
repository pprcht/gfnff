! ──────────────────────────────────────────────────────────────────────────────
! This file is part of gfnff.
!
! Copyright (C) 2023-2026 Philipp Pracht
!
! gfnff is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! gfnff is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with gfnff. If not, see <https://www.gnu.org/licenses/>.
! ──────────────────────────────────────────────────────────────────────────────

!> c bindings for gfnff

module gfnff_interface_c
  use iso_c_binding
  use iso_fortran_env,only:wp => real64,stderr => error_unit
  use gfnff_interface,only:gfnff_data,gfnff_hessian,gfnff_parametrisation_io_available
  use gfnff_param,only:gffVersion
  implicit none
  private

  !> Public C-compatible interface
  public :: c_gfnff_calculator
  public :: c_gfnff_calculator_init
  public :: c_gfnff_calculator_init_pbc
  public :: c_gfnff_calculator_init_ex
  public :: c_gfnff_calculator_deallocate
  public :: c_gfnff_calculator_singlepoint
  public :: c_gfnff_calculator_hessian
  public :: c_gfnff_calculator_charges
  public :: c_gfnff_calculator_results
  public :: c_gfnff_version_from_name
  public :: c_gfnff_toml_available

  !> Returned in iostat when the Hessian is asked for on a periodic system.
  !> Distinct from the library's own codes, which are 0 and +-1.
  integer(c_int),parameter :: GFNFF_C_PBC_UNSUPPORTED = -2_c_int

  !> Returned in iostat when charges are asked for before any singlepoint ran.
  integer(c_int),parameter :: GFNFF_C_NO_CHARGES = -3_c_int

  !> Sentinel meaning "caller did not choose a version"; the library default
  !> then applies. Safe because no gffVersion member is 0.
  integer(c_int),parameter :: GFNFF_C_VERSION_DEFAULT = 0_c_int

  !> C-compatible type containing a pointer to the original Fortran type
  type,bind(C) :: c_gfnff_calculator
    !> C will understand fortran types as pointers
    type(c_ptr) :: ptr
  end type c_gfnff_calculator

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!>--- C-compatible initialization function
  function c_gfnff_calculator_init(c_nat,c_at,c_xyz,c_ichrg,c_printlevel,  &
    &                              c_solvent) &
    &                                  result(calculator) &
    &                                  bind(C,name="c_gfnff_calculator_init")
    implicit none
    type(c_gfnff_calculator) :: calculator
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),target,intent(in) :: c_at(*)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    !> WARNING: row-first vs column-first difference  in Fortran and C!
    real(c_double),target,intent(in) :: c_xyz(3,*)
    !> We assume here that a 3-by-x elements are passed, which in C corresponds
    !> to a vector of length nat for x, y and z coordinates respectively
    !> when xyz[3][nat] was defined.
    !> Hence, it should be defined as xyz[nat][3] in C in order for Fortran
    !> to handle everything correctly in the following!
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    !character(kind=c_char,len=1),intent(in) :: c_model(*)
    integer(c_int),value,intent(in) :: c_ichrg
    !integer(c_int),value,intent(in) :: c_iunit
    integer(c_int),value,intent(in) :: c_printlevel
    character(kind=c_char),intent(in) :: c_solvent(*)
    type(gfnff_data),pointer :: calc

    integer :: nat
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    character(len=:),allocatable :: solvent
    integer :: printlevel,iostatus
    integer :: ichrg

    !> Convert C arguments to Fortran types
    nat = c_nat
    call c_f_pointer(c_loc(c_at),at, [nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,nat]) !> assumes xyz[nat][3] in C
    ichrg = c_ichrg
    printlevel = c_printlevel
    solvent = c_string_to_fortran(c_solvent)
    if (len_trim(solvent) == 0) deallocate (solvent)

    !> Allocate and initialize the Fortran calculator
    allocate (calc)
    call calc%init(nat,at,xyz,ichrg=ichrg, &
    &              printlevel=printlevel,iostat=iostatus,&
    &              solvent=solvent)
    if (iostatus == 0) then
      !> Store the pointer in the C-compatible structure
      calculator%ptr = c_loc(calc)
    else
      write (stderr,'(a,i0)') 'Error initializing GFN-FF calculator. code ',iostatus
      calculator%ptr = c_null_ptr
      deallocate (calc)
    end if
  end function c_gfnff_calculator_init

!========================================================================================!

!>--- C-compatible initialization function with PBC support
  function c_gfnff_calculator_init_pbc(c_nat,c_at,c_xyz,c_ichrg,c_printlevel, &
    &                                  c_lattice,c_npbc) &
    &                                  result(calculator) &
    &                                  bind(C,name="c_gfnff_calculator_init_pbc")
    !***********************************************************
    !* PBC-aware version of c_gfnff_calculator_init.
    !*
    !* INPUT:
    !*   c_nat         - number of atoms
    !*   c_at(c_nat)   - atomic numbers
    !*   c_xyz[nat][3] - Cartesian coordinates (Bohr), C row-major
    !*   c_ichrg       - total molecular charge
    !*   c_printlevel  - verbosity (0=silent)
    !*   c_lattice[3][3] - lattice vectors (Bohr), C row-major
    !*   c_npbc        - number of periodic dimensions (0-3)
    !***********************************************************
    implicit none
    type(c_gfnff_calculator) :: calculator
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),target,intent(in) :: c_at(*)
    real(c_double),target,intent(in) :: c_xyz(3,*)
    integer(c_int),value,intent(in) :: c_ichrg
    integer(c_int),value,intent(in) :: c_printlevel
    !> lattice passed as double[3][3] in C (row-major), maps to Fortran (3,3)
    real(c_double),intent(in) :: c_lattice(3,3)
    integer(c_int),value,intent(in) :: c_npbc
    type(gfnff_data),pointer :: calc

    integer :: nat,npbc
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    integer :: printlevel,iostatus,ichrg

    nat = c_nat
    call c_f_pointer(c_loc(c_at),at, [nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,nat])
    ichrg = c_ichrg
    printlevel = c_printlevel
    npbc = c_npbc

    allocate (calc)
    call calc%init(nat,at,xyz,ichrg=ichrg, &
    &              printlevel=printlevel,iostat=iostatus, &
    &              lattice=c_lattice,npbc=npbc)
    if (iostatus == 0) then
      calculator%ptr = c_loc(calc)
    else
      write (stderr,'(a,i0)') 'Error initializing GFN-FF PBC calculator. code ',iostatus
      calculator%ptr = c_null_ptr
      deallocate (calc)
    end if
  end function c_gfnff_calculator_init_pbc

!========================================================================================!

!>--- Extended C-compatible initializer with optional host-supplied inputs
  function c_gfnff_calculator_init_ex(c_nat,c_at,c_xyz,c_ichrg,c_printlevel, &
    &                                 c_solvent,c_lattice,c_npbc,c_fraglist,c_refq, &
    &                                 c_accuracy,c_version,c_parametrisation, &
    &                                 c_bondmat) &
    &                                 result(calculator) &
    &                                 bind(C,name="c_gfnff_calculator_init_ex")
    !***********************************************************
    !* Superset of c_gfnff_calculator_init / _init_pbc. Every extra input is
    !* optional and skipped when a NULL pointer is passed:
    !*
    !*   c_solvent      - ALPB solvent name (empty string -> none)
    !*   c_lattice[3][3]- lattice vectors (Bohr); NULL -> non-periodic
    !*   c_npbc         - number of periodic dims (used only if c_lattice given)
    !*   c_fraglist[nat]- user-defined fragment index per atom; NULL -> auto.
    !*                    No bonds are formed between atoms of differing fragments.
    !*   c_refq[nat]    - atomic reference charges; NULL -> none. Summed per
    !*                    fragment to define the per-fragment EEQ charge constraint.
    !*   c_accuracy     - cutoff/precision factor; larger is looser and faster,
    !*                    above 1.0 the EEQ system is solved in single
    !*                    precision. Pass <= 0 for the library default
    !*                    (0.1, or 2.0 above 10000 atoms).
    !*   c_version      - gffVersion member selecting the parametrisation
    !*                    version; pass 0 for the library default. Use
    !*                    c_gfnff_version_from_name to resolve a name.
    !*   c_parametrisation - path to a parameter file; NULL or empty -> the
    !*                    internal set for c_version. A ".toml" name is read as
    !*                    an overlay on that set, so a file need only name the
    !*                    keys it changes; any other name is read as the legacy
    !*                    flat format.
    !*   c_bondmat[nat*nat] - molecular graph of integer bond orders; NULL ->
    !*                    GFN-FF perceives the bonds from the geometry. A
    !*                    nonzero element declares a bond. Storage order does
    !*                    not matter because the matrix must be symmetric.
    !*                    Only the zero/nonzero pattern is used today.
    !*                    Molecular systems only.
    !*
    !* The existing _init / _init_pbc entry points are unchanged; this one adds
    !* the host-supplied "bundle" hints for standalone use without an automatic
    !* charge model (the caller simply provides the charge array if desired).
    !***********************************************************
    implicit none
    type(c_gfnff_calculator) :: calculator
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),target,intent(in) :: c_at(*)
    real(c_double),target,intent(in) :: c_xyz(3,*)
    integer(c_int),value,intent(in) :: c_ichrg
    integer(c_int),value,intent(in) :: c_printlevel
    character(kind=c_char),intent(in) :: c_solvent(*)
    type(c_ptr),value,intent(in) :: c_lattice   !> double[3][3] or NULL
    integer(c_int),value,intent(in) :: c_npbc
    type(c_ptr),value,intent(in) :: c_fraglist  !> int[nat]    or NULL
    type(c_ptr),value,intent(in) :: c_refq      !> double[nat] or NULL
    real(c_double),value,intent(in) :: c_accuracy  !> <= 0 -> library default
    integer(c_int),value,intent(in) :: c_version   !> 0 -> library default
    !> Taken by value as a pointer rather than as character(kind=c_char)(*) so
    !> that a NULL is representable; a dummy character array cannot be tested
    !> for NULL and would be dereferenced blindly.
    type(c_ptr),value,intent(in) :: c_parametrisation  !> const char* or NULL
    type(c_ptr),value,intent(in) :: c_bondmat  !> int[nat*nat] or NULL
    type(gfnff_data),pointer :: calc

    integer :: nat,npbc
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    real(wp),pointer :: lattice(:,:)
    integer,pointer :: fraglist(:)
    integer,pointer :: bondmat(:,:)
    real(wp),pointer :: refq(:)
    character(len=:),allocatable :: solvent,pfile
    integer :: printlevel,iostatus,ichrg,version
    real(wp) :: accuracy
    logical :: have_acc,have_ver

    nat = c_nat
    call c_f_pointer(c_loc(c_at),at, [nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,nat])
    ichrg = c_ichrg
    printlevel = c_printlevel
    npbc = c_npbc
    solvent = c_string_to_fortran(c_solvent)
    if (len_trim(solvent) == 0) deallocate (solvent)
    accuracy = c_accuracy
    have_acc = accuracy > 0.0_wp
    version = c_version
    have_ver = c_version /= GFNFF_C_VERSION_DEFAULT

    allocate (calc)

    !> Optional host-supplied bundle hints. Must be set on the data object
    !> BEFORE init runs, since they steer the topology setup. type_init (called
    !> inside init) does not touch the userinput component, so these survive.
    if (c_associated(c_fraglist)) then
      call c_f_pointer(c_fraglist,fraglist, [nat])
      if (.not.allocated(calc%userinput)) allocate (calc%userinput)
      calc%userinput%fraglist = fraglist
    end if
    if (c_associated(c_refq)) then
      call c_f_pointer(c_refq,refq, [nat])
      if (.not.allocated(calc%userinput)) allocate (calc%userinput)
      calc%userinput%refq = refq
    end if
    if (c_associated(c_bondmat)) then
      call c_f_pointer(c_bondmat,bondmat, [nat,nat])
      if (.not.allocated(calc%userinput)) allocate (calc%userinput)
      calc%userinput%bondmat = bondmat
    end if

    !> The parameter file is a data-object field rather than an init argument,
    !> so it has to be set here; init reads it after loading the internal set
    !> for the selected version.
    pfile = c_ptr_to_fortran(c_parametrisation)
    if (len_trim(pfile) > 0) calc%parametrisation = pfile

    if (c_associated(c_lattice)) then
      call c_f_pointer(c_lattice,lattice, [3,3])
      if (have_acc) then
        call calc%init(nat,at,xyz,ichrg=ichrg,printlevel=printlevel, &
        &              iostat=iostatus,solvent=solvent,lattice=lattice,npbc=npbc, &
        &              accuracy=accuracy,version=merge(version,calc%version,have_ver))
      else
        call calc%init(nat,at,xyz,ichrg=ichrg,printlevel=printlevel, &
        &              iostat=iostatus,solvent=solvent,lattice=lattice,npbc=npbc, &
        &              version=merge(version,calc%version,have_ver))
      end if
    else
      if (have_acc) then
        call calc%init(nat,at,xyz,ichrg=ichrg,printlevel=printlevel, &
        &              iostat=iostatus,solvent=solvent,accuracy=accuracy, &
        &              version=merge(version,calc%version,have_ver))
      else
        call calc%init(nat,at,xyz,ichrg=ichrg,printlevel=printlevel, &
        &              iostat=iostatus,solvent=solvent, &
        &              version=merge(version,calc%version,have_ver))
      end if
    end if

    if (iostatus == 0) then
      calculator%ptr = c_loc(calc)
    else
      write (stderr,'(a,i0)') 'Error initializing GFN-FF calculator (ex). code ',iostatus
      calculator%ptr = c_null_ptr
      deallocate (calc)
    end if
  end function c_gfnff_calculator_init_ex

!========================================================================================!

  subroutine c_gfnff_calculator_deallocate(calculator) &
    &     bind(C,name="c_gfnff_calculator_deallocate")
    type(c_gfnff_calculator),intent(inout) :: calculator
    type(gfnff_data),pointer :: calc_ptr

    !> Convert the C pointer to a Fortran pointer
    call c_f_pointer(calculator%ptr,calc_ptr)

    !> Deallocate the Fortran object
    if (associated(calc_ptr)) then
      call calc_ptr%deallocate()
      deallocate (calc_ptr)
    end if

    !> Nullify the C pointer
    calculator%ptr = c_null_ptr
  end subroutine c_gfnff_calculator_deallocate

!========================================================================================!

  subroutine c_gfnff_calculator_singlepoint(c_calculator,c_nat,c_at,c_xyz, &
    &                                       c_energy,c_gradient,c_sigma,c_lattice,c_iostat) &
    &                        bind(C,name="c_gfnff_calculator_singlepoint")
    !***********************************************************
    !* Compute energy, gradient and stress tensor for the
    !* current geometry.
    !*
    !* INPUT:
    !*   c_calculator    - opaque handle (from init)
    !*   c_nat           - number of atoms
    !*   c_at(c_nat)     - atomic numbers
    !*   c_xyz[nat][3]   - Cartesian coordinates (Bohr)
    !*   c_lattice[3][3] - new lattice vectors (Bohr); pass NULL to reuse
    !*                     the stored lattice (non-PBC or unchanged cell)
    !* OUTPUT:
    !*   c_energy        - total energy (Hartree)
    !*   c_gradient[nat][3] - gradient (Eh/Bohr)
    !*   c_sigma[3][3]   - stress tensor (Hartree); zero for non-PBC
    !*   c_iostat        - error status (0 = success)
    !***********************************************************
    implicit none
    !> Input arguments from C
    type(c_gfnff_calculator),intent(inout) :: c_calculator
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),target,intent(in) :: c_at(*)
    real(c_double),target,intent(in) :: c_xyz(3,*) !> NOTE Fortran/C matrix orders

    !> Output arguments to C
    real(c_double),intent(out) :: c_energy
    real(c_double),target,intent(out) :: c_gradient(3,*) !> NOTE Fortran/C matrix orders
    real(c_double),intent(out) :: c_sigma(3,3)            !> stress tensor (Eh); 0 for non-PBC
    type(c_ptr),value,intent(in) :: c_lattice             !> lattice[3][3] or c_null_ptr
    integer(c_int),intent(out) :: c_iostat

    !> Local Fortran variables
    type(gfnff_data),pointer :: calc_ptr
    integer :: nat
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    real(wp),pointer :: grad(:,:)
    real(wp),pointer :: lattice_loc(:,:)
    real(wp) :: energy
    real(wp) :: sigma_loc(3,3)
    integer :: iostat

    !> Convert C pointers to Fortran pointers
    call c_f_pointer(c_calculator%ptr,calc_ptr)
    call c_f_pointer(c_loc(c_at),at, [c_nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,c_nat]) !> Assumes xyz[nat][3] in C
    call c_f_pointer(c_loc(c_gradient),grad, [3,c_nat])  !> Assumes grad[nat][3] in C

    !> Set the integer variable
    nat = c_nat

    !> Use caller-supplied lattice if provided, otherwise reuse stored lattice.
    if (c_associated(c_lattice)) then
      call c_f_pointer(c_lattice,lattice_loc, [3,3])
      call calc_ptr%singlepoint(nat,at,xyz,energy,grad,iostat=iostat, &
      &                         lattice=lattice_loc,sigma=sigma_loc)
    else
      call calc_ptr%singlepoint(nat,at,xyz,energy,grad,iostat=iostat, &
      &                         lattice=calc_ptr%cell%lattice,sigma=sigma_loc)
    end if

    !> Pass back the results to C variables
    c_energy = energy
    c_gradient(1:3,1:nat) = grad(1:3,1:nat)
    !> Zero sigma for non-periodic systems; Fortran accumulates virial-style
    !> contributions regardless of PBC, so we suppress them here.
    if (calc_ptr%cell%npbc > 0) then
      c_sigma(1:3,1:3) = sigma_loc(1:3,1:3)
    else
      c_sigma(1:3,1:3) = 0.0_wp
    end if
    c_iostat = iostat

  end subroutine c_gfnff_calculator_singlepoint

!========================================================================================!

  subroutine c_gfnff_calculator_results(c_calculator,c_iunit) &
    &   bind(C,name="c_gfnff_calculator_results")
    implicit none
    !> Input arguments from C
    type(c_gfnff_calculator),intent(in) :: c_calculator
    integer(c_int),value,intent(in) :: c_iunit
    !> Local Fortran variables
    type(gfnff_data),pointer :: calc_ptr
    integer :: myunit
    !> Convert C pointer to Fortran pointer
    call c_f_pointer(c_calculator%ptr,calc_ptr)
    myunit = c_iunit
    !> Call the Fortran subroutine
    call calc_ptr%resultprint(myunit)
  end subroutine c_gfnff_calculator_results

!========================================================================================!

  function c_string_to_fortran(c_str) result(f_str)
    use iso_c_binding
    implicit none
    character(kind=c_char),intent(in) :: c_str(*)   !> C null-terminated string
    character(len=:),allocatable :: f_str           !> Fortran allocatable string
    integer :: i
    !> Find the null terminator dynamically
    i = 1
    do while (c_str(i) /= c_null_char)
      i = i+1
    end do
    i = i-1  !> Exclude null terminator
    !> Allocate Fortran string with exact length
    allocate (character(len=i) :: f_str)
    !> Copy contents from C string to Fortran string
    f_str = transfer(c_str(1:i),f_str)
  end function c_string_to_fortran

!========================================================================================!

  function c_ptr_to_fortran(c_str) result(f_str)
    !***********************************************************
    !* NULL-tolerant variant of c_string_to_fortran: returns a
    !* zero-length string for a null pointer instead of walking
    !* off the end of it looking for a terminator.
    !***********************************************************
    implicit none
    type(c_ptr),value,intent(in) :: c_str  !> const char* or NULL
    character(len=:),allocatable :: f_str
    character(kind=c_char),pointer :: chars(:)
    integer :: i,n

    if (.not.c_associated(c_str)) then
      f_str = ''
      return
    end if

    !> Length is unknown up front, so map a generous window and scan it. The
    !> bound only limits how long a path may be, it is never dereferenced
    !> beyond the terminator.
    call c_f_pointer(c_str,chars, [4096])
    n = 0
    do i = 1,4096
      if (chars(i) == c_null_char) exit
      n = i
    end do

    allocate (character(len=n) :: f_str)
    do i = 1,n
      f_str(i:i) = chars(i)
    end do
  end function c_ptr_to_fortran

!========================================================================================!

  function c_gfnff_toml_available() result(yes) &
    &   bind(C,name="c_gfnff_toml_available")
    !***********************************************************
    !* Whether this build can read TOML parameter files.
    !*
    !* OUTPUT: 1 if the library was built with toml-f, else 0.
    !*
    !* Lets a host tell "TOML support missing" apart from "file
    !* missing or malformed", which otherwise both surface as a
    !* failed initialisation.
    !***********************************************************
    implicit none
    integer(c_int) :: yes

    if (gfnff_parametrisation_io_available()) then
      yes = 1_c_int
    else
      yes = 0_c_int
    end if
  end function c_gfnff_toml_available

!========================================================================================!

  function c_gfnff_version_from_name(c_name) result(version) &
    &   bind(C,name="c_gfnff_version_from_name")
    !***********************************************************
    !* Resolve a parametrisation version name to its gffVersion
    !* value, so a host never has to hardcode the integers.
    !*
    !* INPUT:
    !*   c_name - version name, e.g. "angewChem2020_2"
    !* OUTPUT:
    !*   the gffVersion member, or 0 if the name is unknown.
    !*   0 is also what _init_ex reads as "use the default", so an
    !*   unrecognised name degrades to the default rather than to
    !*   an invalid version.
    !***********************************************************
    implicit none
    character(kind=c_char),intent(in) :: c_name(*)
    integer(c_int) :: version
    character(len=:),allocatable :: name

    name = c_string_to_fortran(c_name)
    select case (trim(name))
    case ('angewChem2020')
      version = int(gffVersion%angewChem2020,c_int)
    case ('angewChem2020_1')
      version = int(gffVersion%angewChem2020_1,c_int)
    case ('angewChem2020_2')
      version = int(gffVersion%angewChem2020_2,c_int)
    case ('harmonic2020')
      version = int(gffVersion%harmonic2020,c_int)
    case ('mcgfnff2023')
      version = int(gffVersion%mcgfnff2023,c_int)
    case ('conformer2020')
      version = int(gffVersion%conformer2020,c_int)
    case default
      version = GFNFF_C_VERSION_DEFAULT
    end select
  end function c_gfnff_version_from_name

!========================================================================================!

  subroutine c_gfnff_calculator_charges(c_calculator,c_nat,c_charges,c_iostat) &
    &   bind(C,name="c_gfnff_calculator_charges")
    !***********************************************************
    !* Atomic partial charges from the last singlepoint.
    !*
    !* INPUT:
    !*   c_calculator - opaque handle (from init)
    !*   c_nat        - number of atoms
    !* OUTPUT:
    !*   c_charges[nat] - EEQ partial charges (e). Caller-owned: the
    !*                    buffer must already hold nat doubles, and it
    !*                    is written only when c_iostat comes back 0.
    !*   c_iostat     - 0 on success, GFNFF_C_NO_CHARGES (-3) if no
    !*                  singlepoint has run yet, nat disagrees with the
    !*                  stored system, or the version has no charges.
    !*
    !* The charges are a by-product of the energy evaluation, not a
    !* separate model, so they only exist after a singlepoint (or a
    !* Hessian, which runs one).
    !***********************************************************
    implicit none
    !> Input arguments from C
    type(c_gfnff_calculator),intent(in) :: c_calculator
    integer(c_int),value,intent(in) :: c_nat
    !> Output arguments to C
    real(c_double),intent(out) :: c_charges(c_nat)
    integer(c_int),intent(out) :: c_iostat
    !> Local Fortran variables
    type(gfnff_data),pointer :: calc_ptr

    call c_f_pointer(c_calculator%ptr,calc_ptr)

    if (.not.associated(calc_ptr)) then
      c_iostat = GFNFF_C_NO_CHARGES
      return
    end if
    if (.not.allocated(calc_ptr%nlist)) then
      c_iostat = GFNFF_C_NO_CHARGES
      return
    end if
    if (.not.allocated(calc_ptr%nlist%q)) then
      c_iostat = GFNFF_C_NO_CHARGES
      return
    end if
    if (size(calc_ptr%nlist%q) /= c_nat) then
      c_iostat = GFNFF_C_NO_CHARGES
      return
    end if
    !> The harmonic version returns from the energy routine before the EEQ
    !> solve, so q is allocated but never filled. Handing back those zeros
    !> would pass an artefact off as a result.
    if (calc_ptr%version == gffVersion%harmonic2020) then
      c_iostat = GFNFF_C_NO_CHARGES
      return
    end if

    c_charges(1:c_nat) = calc_ptr%nlist%q(1:c_nat)
    c_iostat = 0_c_int
  end subroutine c_gfnff_calculator_charges

!========================================================================================!

  subroutine c_gfnff_calculator_hessian(c_calculator,c_nat,c_at,c_xyz, &
    &                                   c_hessian,c_energy,c_gradient,c_step,c_iostat) &
    &                        bind(C,name="c_gfnff_calculator_hessian")
    !***********************************************************
    !* Cartesian nuclear Hessian for the current geometry.
    !*
    !* INPUT:
    !*   c_calculator    - opaque handle (from init)
    !*   c_nat           - number of atoms
    !*   c_at(c_nat)     - atomic numbers
    !*   c_xyz[nat][3]   - Cartesian coordinates (Bohr)
    !*   c_step          - finite-difference step (Bohr) for the terms not
    !*                     yet available in closed form; pass <= 0 to let
    !*                     the library choose
    !* OUTPUT:
    !*   c_hessian       - (3*nat)x(3*nat) second derivatives (Eh/Bohr^2).
    !*                     Caller-owned: the buffer must already hold
    !*                     9*nat*nat doubles, and it is written only when
    !*                     c_iostat comes back 0.
    !*                     The matrix is symmetrised before it is returned,
    !*                     so it reads the same row-major from C as it does
    !*                     column-major from Fortran and no transpose is
    !*                     needed on either side. Degree of freedom (c,A)
    !*                     sits at index 3*(A-1)+c, one-based.
    !*   c_energy        - total energy (Hartree); pass NULL to skip
    !*   c_gradient[nat][3] - analytic gradient (Eh/Bohr); pass NULL to skip
    !*   c_iostat        - 0 on success, GFNFF_C_PBC_UNSUPPORTED (-2) if the
    !*                     calculator was set up periodic, otherwise the
    !*                     library error code
    !*
    !* Periodic systems are rejected rather than silently finite-differenced:
    !* the analytic terms have no periodic implementation yet, and the
    !* fallback would ignore the images without saying so.
    !***********************************************************
    implicit none
    !> Input arguments from C
    type(c_gfnff_calculator),intent(inout) :: c_calculator
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),target,intent(in) :: c_at(*)
    real(c_double),target,intent(in) :: c_xyz(3,*) !> NOTE Fortran/C matrix orders
    real(c_double),value,intent(in) :: c_step

    !> Output arguments to C
    real(c_double),intent(out) :: c_hessian(3*c_nat,3*c_nat)
    type(c_ptr),value,intent(in) :: c_energy      !> double* or c_null_ptr
    type(c_ptr),value,intent(in) :: c_gradient    !> gradient[nat][3] or c_null_ptr
    integer(c_int),intent(out) :: c_iostat

    !> Local Fortran variables
    type(gfnff_data),pointer :: calc_ptr
    integer :: nat
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    real(wp),pointer :: energy_out
    real(wp),pointer :: grad_out(:,:)
    real(wp),allocatable :: hess(:,:),grad(:,:)
    real(wp) :: energy
    integer :: iostat

    call c_f_pointer(c_calculator%ptr,calc_ptr)

    !> Reject periodic systems before anything is allocated or written, so a
    !> caller that ignores iostat still sees its buffer untouched rather than
    !> filled with a molecular Hessian for a periodic system.
    if (calc_ptr%cell%npbc > 0) then
      c_iostat = GFNFF_C_PBC_UNSUPPORTED
      return
    end if

    call c_f_pointer(c_loc(c_at),at, [c_nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,c_nat]) !> Assumes xyz[nat][3] in C

    nat = c_nat
    allocate (hess(3*nat,3*nat),source=0.0_wp)
    allocate (grad(3,nat),source=0.0_wp)

    if (c_step > 0.0_wp) then
      call gfnff_hessian(nat,at,xyz,calc_ptr,hess,energy=energy,gradient=grad, &
      &                  step=c_step,iostat=iostat)
    else
      call gfnff_hessian(nat,at,xyz,calc_ptr,hess,energy=energy,gradient=grad, &
      &                  iostat=iostat)
    end if

    c_iostat = iostat
    if (iostat /= 0) return

    c_hessian(1:3*nat,1:3*nat) = hess(1:3*nat,1:3*nat)
    if (c_associated(c_energy)) then
      call c_f_pointer(c_energy,energy_out)
      energy_out = energy
    end if
    if (c_associated(c_gradient)) then
      call c_f_pointer(c_gradient,grad_out, [3,nat])
      grad_out(1:3,1:nat) = grad(1:3,1:nat)
    end if

  end subroutine c_gfnff_calculator_hessian

!========================================================================================!
!========================================================================================!
end module gfnff_interface_c
