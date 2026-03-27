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
module gfnff_interface
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_data_types
  use gfnff_neighbor,only:TNeigh
  use gfnff_engrad_module
  use gfnff_gbsa
  use gfnff_param
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: gfnff_data,my_gfnff_data
  public :: gfnff_initialize
  public :: gfnff_singlepoint
  public :: print_gfnff_results
  public :: gfnff_get_fake_wbo

!> this type bundles together most of the
!> data required for a GFN-FF calculation
  type :: gfnff_data
    integer  :: ichrg = 0  !> total molecular charge
    real(wp) :: accuracy = 0.1_wp
    logical  :: make_chrg = .true.
    integer  :: version = gffVersion%angewChem2020_2
    logical  :: update = .true.
    logical  :: write_topo = .true.
    character(len=:),allocatable :: parametrisation
    character(len=:),allocatable :: solvent
    logical :: restart = .false.
    character(len=:),allocatable :: restartfile
    character(len=:),allocatable :: refgeo
    character(len=:),allocatable :: refcharges

    type(TGFFGenerator),allocatable     :: gen
    type(TGFFData),allocatable          :: param
    type(TGFFTopology),allocatable      :: topo
    type(TNeigh),allocatable            :: neigh
    type(TCell),allocatable             :: cell
    type(TGFFNeighbourList),allocatable :: nlist
    type(TBorn),allocatable             :: solvation
    type(gfnff_results),allocatable     :: res
  contains
    procedure :: deallocate => gfnff_data_deallocate
    procedure :: type_reset => gfnff_data_reset_types
    procedure :: type_init => gfnff_data_make_types
    procedure :: singlepoint => gfnff_singlepoint_wrapper
    procedure :: init => gfnff_initialize_wrapper
    procedure :: resultprint => gfnff_print_results_wrapper
  end type gfnff_data

!> This is a semi-global placeholder for a single gfnff_data object.
!> It may be used as storage for running GFN-FF calculations,
!> but ONLY if there is ever a single instance running at the same time!
!> i.e., no parallelization shenanigans !!!
  type(gfnff_data),allocatable :: my_gfnff_data

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine gfnff_singlepoint(nat,at,xyz,dat,energy,gradient,lattice,sigma,printlevel,printunit,iostat)
!**********************************************************************
!* GFN-FF single-point energy and gradient calculation.
!*
!* INPUT:
!*   nat          - number of atoms
!*   at(nat)      - atomic numbers
!*   xyz(3,nat)   - Cartesian coordinates (Bohr)
!*   dat          - bundled GFN-FF data and settings
!*   printlevel   - optional verbosity (0=silent,1=errors,2=info,3=verbose)
!*   printunit    - optional output unit (default: stdout)
!* OUTPUT:
!*   energy       - total energy (Eh)
!*   gradient     - gradient (Eh/Bohr)
!*   iostat       - optional error status
!*   sigma(3,3)   - optional stress tensor (Eh); zero for non-PBC
!*********************************************************************
    implicit none
    !> INPUT
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: at(nat)    !> atom types
    real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    integer,intent(in),optional    :: printlevel  !> verbosity level
    integer,intent(in),optional    :: printunit   !> output unit
    type(gfnff_data),intent(inout) :: dat  !> collection of gfnff datatypes and settings
    real(wp),intent(in),optional :: lattice(3,3)
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    integer,intent(out),optional  :: iostat
    real(wp),intent(out),optional :: sigma(3,3) !> stress tensor (zero for non-PBC)
    !> LOCAL
    integer :: io,mylevel,myunit
    real(wp) :: sigma_loc(3,3),lattice_loc(3,3),efield_loc(3)
    real(wp) :: lthr = sqrt(epsilon(1.0_wp))

    if (present(printlevel)) then
      mylevel = printlevel
    else
      mylevel = 0
    end if
    if (present(printunit)) then
      myunit = printunit
    else
      myunit = stdout
    end if

! ── init datafields ───────────────────────────────────────────────────────────
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    io = 0
    sigma_loc = 0.0_wp
    lattice_loc = 0.0_wp
    efield_loc = 0.0_wp
    if (present(lattice)) lattice_loc(:,:) = lattice(:,:)

! ── update lattice, wsc, ... ──────────────────────────────────────────────────
    if (any(abs(lattice_loc-dat%cell%lattice) .gt. lthr)) then
      call dat%cell%init(lattice_loc)
    end if
    if (dat%cell%npbc > 0) then
      call dat%cell%init_wsc(nat,at,xyz)
    end if

! ── call E+Grd ────────────────────────────────────────────────────────────────
    call gfnff_eg(mylevel,nat,at,xyz,dat%cell,sigma_loc,dat%ichrg,gradient,energy, &
    &            dat%res,dat%param,dat%topo,dat%neigh,dat%nlist,efield_loc,        &
    &            dat%solvation,dat%update,dat%version,dat%accuracy,printunit=myunit)

! ── transfer optional outputs ─────────────────────────────────────────────────
    if (present(sigma)) sigma = sigma_loc

    if (present(iostat)) then
      iostat = io
    end if

  end subroutine gfnff_singlepoint

! ══════════════════════════════════════════════════════════════════════════════
  subroutine gfnff_singlepoint_wrapper(self,nat,at,xyz,energy,gradient, &
  &                                    printlevel,printunit,iostat,lattice,sigma)
!******************************************************************
!* A wrapper to the singlepoint routine, allowing
!* the energy routine to be called with "call dat%singlepoint(...)"
!******************************************************************
    implicit none
    !> INPUT
    class(gfnff_data) :: self
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: at(nat)    !> atom types
    real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    integer,intent(in),optional  :: printlevel  !> verbosity level
    integer,intent(in),optional  :: printunit   !> output unit
    real(wp),intent(in),optional :: lattice(3,3)  !> lattice (optional)
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    real(wp),intent(out),optional :: sigma(3,3) !> stress (optional)
    integer,intent(out),optional  :: iostat
    call gfnff_singlepoint(nat,at,xyz,self,energy,gradient, &
    & printlevel=printlevel,printunit=printunit,iostat=iostat,lattice=lattice,sigma=sigma)
  end subroutine gfnff_singlepoint_wrapper

!========================================================================================!

  subroutine print_gfnff_results(printunit,res_gff,lsolv)
    integer,intent(in) :: printunit ! file handle (usually output_unit=6)
    type(gfnff_results),intent(in) :: res_gff
    logical,intent(in) :: lsolv
    character(len=*),parameter :: outfmt = &
                                  '(2x,a,f23.12,1x,a)'
    write (printunit,outfmt) "total energy      ",res_gff%e_total,"Eh   "
    write (printunit,outfmt) "gradient norm     ",res_gff%gnorm,"Eh/a0"
    write (printunit,'(a)') repeat('-',50)
    write (printunit,outfmt) "bond energy       ",res_gff%e_bond,"Eh   "
    write (printunit,outfmt) "angle energy      ",res_gff%e_angl,"Eh   "
    write (printunit,outfmt) "torsion energy    ",res_gff%e_tors,"Eh   "
    write (printunit,outfmt) "repulsion energy  ",res_gff%e_rep,"Eh   "
    write (printunit,outfmt) "electrostat energy",res_gff%e_es,"Eh   "
    write (printunit,outfmt) "dispersion energy ",res_gff%e_disp,"Eh   "
    write (printunit,outfmt) "HB energy         ",res_gff%e_hb,"Eh   "
    write (printunit,outfmt) "XB energy         ",res_gff%e_xb,"Eh   "
    write (printunit,outfmt) "bonded atm energy ",res_gff%e_batm,"Eh   "
    write (printunit,outfmt) "external energy   ",res_gff%e_ext,"Eh   "
    if (lsolv) then
      write (printunit,'(2x,a)') repeat('-',44)
      write (printunit,outfmt) "-> Gsolv          ",res_gff%g_solv,"Eh   "
      write (printunit,outfmt) "   -> Gborn       ",res_gff%g_born,"Eh   "
      write (printunit,outfmt) "   -> Gsasa       ",res_gff%g_sasa,"Eh   "
      write (printunit,outfmt) "   -> Ghb         ",res_gff%g_hb,"Eh   "
      write (printunit,outfmt) "   -> Gshift      ",res_gff%g_shift,"Eh   "
    end if
    write (printunit,'(a)') repeat('-',50)
  end subroutine print_gfnff_results

  subroutine gfnff_print_results_wrapper(self,printunit)
    implicit none
    class(gfnff_data) :: self
    !> INPUT
    integer,intent(in),optional :: printunit
    !> LOCAL
    integer :: myunit

    if (present(printunit)) then
      myunit = printunit
    else
      myunit = stdout
    end if
    if (allocated(self%res)) then
      call print_gfnff_results(myunit,self%res,allocated(self%solvation))
    end if
  end subroutine gfnff_print_results_wrapper
!========================================================================================!

  subroutine gfnff_initialize(nat,at,xyz,dat, &
  &                 printlevel,printunit,version,iostat,ichrg,lattice,npbc)
    !*************************************************************
    !* Initialize a GFN-FF calculation: load parameters, build
    !* topology, and optionally set up periodic boundary conditions.
    !*
    !* INPUT:
    !*   nat         - number of atoms
    !*   at(nat)     - atomic numbers
    !*   xyz(3,nat)  - Cartesian coordinates (Bohr)
    !*   dat         - bundled GFN-FF data (modified in-place)
    !*   printlevel  - optional verbosity (0=silent,1=errors,2=info,3=verbose)
    !*   printunit   - optional output unit (default: stdout)
    !*   version     - optional: GFN-FF parametrisation version
    !*   ichrg       - optional: total molecular charge (default: 0)
    !*   lattice(3,3)- optional: lattice vectors in Bohr (column-major)
    !*   npbc        - optional: number of periodic dimensions (0-3)
    !* OUTPUT:
    !*   iostat      - optional: error status (0=success)
    !*************************************************************
    use gfnff_param
    use gfnff_setup_mod,only:gfnff_setup
    use gfnff_gdisp0,only:newD3Model
    use gfnff_gbsa
    character(len=*),parameter :: source = 'gfnff_initialize'
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in),optional  :: printlevel  !> verbosity level
    integer,intent(in),optional  :: printunit   !> output unit
    integer,intent(in),optional  :: version
    integer,intent(out),optional :: iostat
    integer,intent(in),optional  :: ichrg
    real(wp),intent(in),optional :: lattice(3,3) !> lattice vectors (Bohr)
    integer,intent(in),optional  :: npbc         !> number of periodic dims (0-3)
    !> OUTPUT
    type(gfnff_data),intent(inout) :: dat
    !> LOCAL
    character(len=:),allocatable :: fname
    integer :: ich,io,mylevel,myunit
    logical :: ex
    logical :: restart

    if (present(printlevel)) then
      mylevel = printlevel
    else
      mylevel = 0
    end if
    if (present(printunit)) then
      myunit = printunit
    else
      myunit = stdout
    end if

!> Reset datatypes
    call dat%type_init()
    if (present(ichrg)) then
      dat%ichrg = ichrg
    end if

!> Periodic boundary conditions setup
    if (present(npbc)) dat%cell%npbc = npbc
    if (present(lattice)) call dat%cell%init(lattice)

!> except restart-related options
    restart = dat%restart
    if (.not.allocated(dat%restartfile)) then
      dat%topo%filename = 'gfnff_topo'
    else
      dat%topo%filename = dat%restartfile
    end if
    if (allocated(dat%refgeo)) restart = .false.
    if (allocated(dat%refcharges)) then
      dat%topo%refcharges = dat%refcharges
    end if

!> Parametrisation version
    if (present(version)) then
      dat%version = version
    else
      dat%version = gffVersion%angewChem2020_2
    end if

    call dat%topo%zero
    dat%update = .true.
!> global accuracy factor similar to acc in xtb used in SCF
    dat%accuracy = 0.1_wp
    if (nat > 10000) then
      dat%accuracy = 2.0_wp
    end if

!> Obtain the parameter file or load internal
    if (allocated(dat%parametrisation)) then
      fname = dat%parametrisation
    else
      fname = 'no file!'
    end if
    inquire (file=fname,exist=ex)
    if (ex) then
      open (newunit=ich,file=fname)
      call gfnff_read_param(ich,dat%param)
      close (ich)
    else !> no parameter file, try to load internal version
      call gfnff_load_param(dat%version,dat%param,ex)
      if (.not.ex.and.mylevel >= 1) then
        write (myunit,'("Parameter file ",a," not found!",a)') fname,source
        return
      end if
    end if

    call newD3Model(dat%topo%dispm,nat,at)

    call gfnff_setup(nat,at,xyz,dat%ichrg,mylevel,restart,dat%write_topo, &
    &        dat%gen,dat%param,dat%topo,dat%neigh,dat%cell,dat%accuracy,dat%version,io, &
    &        printunit=myunit)

    !> Optional, ALPB solvation
    if (allocated(dat%solvent)) then
      if (.not. (allocated(dat%solvation))) allocate (dat%solvation)
      call gfnff_gbsa_init(nat,at,dat%solvent,dat%solvation)
      if (mylevel >= 2) then
        write (myunit,*)
        call gfnff_gbsa_print(dat%solvation,myunit)
      end if
    end if

    if ((io /= 0).and.mylevel >= 1) then
      write (myunit,'("Could not create force field calculator ",a)') source
    end if
    if (present(iostat)) then
      iostat = io
    end if
  end subroutine gfnff_initialize

  subroutine gfnff_initialize_wrapper(self,nat,at,xyz, &
     &                 printlevel,printunit,version,iostat,ichrg,solvent,lattice,npbc)
!******************************************************************
!* A wrapper to the initialize routine, allowing
!* the energy routine to be called with "call dat%init(...)"
!******************************************************************
    implicit none
    class(gfnff_data) :: self
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in),optional  :: printlevel  !> verbosity level
    integer,intent(in),optional  :: printunit   !> output unit
    integer,intent(in),optional  :: version
    integer,intent(out),optional :: iostat
    integer,intent(in),optional  :: ichrg
    character(len=*),intent(in),optional :: solvent
    real(wp),intent(in),optional :: lattice(3,3) !> lattice vectors (Bohr)
    integer,intent(in),optional  :: npbc         !> number of periodic dims (0-3)

    if (present(solvent)) then
      if (solvent .ne. 'none'.and.len_trim(solvent) > 0) self%solvent = solvent
    end if

    call gfnff_initialize(nat,at,xyz,self, &
    &       printlevel=printlevel,printunit=printunit, &
    &       version=version,iostat=iostat,ichrg=ichrg,lattice=lattice,npbc=npbc)
  end subroutine gfnff_initialize_wrapper

!========================================================================================!
  subroutine gfnff_data_deallocate(self)
    implicit none
    class(gfnff_data) :: self
    self%ichrg = 0
    self%accuracy = 0.1_wp
    self%make_chrg = .true.
    self%version = 1
    self%update = .true.
    self%write_topo = .true.
    if (allocated(self%solvent)) deallocate (self%solvent)
    if (allocated(self%gen)) deallocate (self%gen)
    if (allocated(self%param)) deallocate (self%param)
    if (allocated(self%topo)) deallocate (self%topo)
    if (allocated(self%neigh)) deallocate (self%neigh)
    if (allocated(self%cell)) deallocate (self%cell)
    if (allocated(self%nlist)) deallocate (self%nlist)
    if (allocated(self%solvation)) deallocate (self%solvation)
    if (allocated(self%res)) deallocate (self%res)
  end subroutine gfnff_data_deallocate
  subroutine gfnff_data_reset_types(self)
    implicit none
    class(gfnff_data) :: self
    if (allocated(self%gen)) deallocate (self%gen)
    if (allocated(self%param)) deallocate (self%param)
    if (allocated(self%topo)) deallocate (self%topo)
    if (allocated(self%neigh)) deallocate (self%neigh)
    if (allocated(self%cell)) deallocate (self%cell)
    if (allocated(self%nlist)) deallocate (self%nlist)
    if (allocated(self%solvation)) deallocate (self%solvation)
    if (allocated(self%res)) deallocate (self%res)
  end subroutine gfnff_data_reset_types
  subroutine gfnff_data_make_types(self)
    implicit none
    class(gfnff_data) :: self
    call self%type_reset()
    allocate (self%gen)
    allocate (self%param)
    allocate (self%topo)
    allocate (self%neigh)
    allocate (self%cell)
    allocate (self%nlist)
    allocate (self%res)
  end subroutine gfnff_data_make_types

!=========================================================================================!
  subroutine gfnff_get_fake_wbo(ff_dat,nat,wbo)
    implicit none
    type(gfnff_data),intent(in) :: ff_dat
    integer,intent(in) :: nat
    real(wp),intent(out) :: wbo(nat,nat)
    integer :: i,k,l
    wbo = 0.0_wp
    if (allocated(ff_dat%topo)) then
      if (allocated(ff_dat%topo%blist)) then
        do i = 1,ff_dat%topo%nbond
          k = ff_dat%topo%blist(1,i)
          l = ff_dat%topo%blist(2,i)
          wbo(k,l) = 1.0_wp
          wbo(l,k) = wbo(k,l)
        end do
      end if
    end if
  end subroutine gfnff_get_fake_wbo

!========================================================================================!
end module gfnff_interface

