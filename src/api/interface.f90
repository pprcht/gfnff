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
  use gfnff_data_types,only:TCell,TGFFData,TGFFGenerator,TGFFNeighbourList, &
    &                        TGFFTopology,TGFFUserInput,init
  use gfnff_neighbor,only:TNeigh
  use gfnff_eg_driver,only:gfnff_eg,gfnff_results
  use gfnff_alpb,only:TBorn,gfnff_gbsa_init,gfnff_gbsa_print
  use gfnff_param,only:gffVersion,gfnff_load_param,gfnff_read_param, &
    &                  gfnff_set_param
  use gfnff_param_io,only:param_read_toml,param_write_toml,param_io_available
  use gfnff_hess_driver,only:gfnff_hessian_core
  use gfnff_hess_analysis,only:hess_frequencies
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: gfnff_data,my_gfnff_data
  public :: gfnff_initialize
!> the version enumerator, so a host can name what it passes as version=
  public :: gffVersion
  public :: gfnff_singlepoint
  public :: gfnff_hessian
  public :: print_gfnff_results
  public :: gfnff_get_fake_wbo
!> parametrisation I/O: read a TOML parameter set, or write the one in use
  public :: gfnff_write_parametrisation
  public :: gfnff_parametrisation_io_available
!> re-exported so a host can turn a Hessian into frequencies without
!> reaching into the internal modules
  public :: hess_frequencies

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
    !> optional, host-supplied setup hints (e.g. user-defined fragments).
    !> Set this before calling gfnff_initialize to steer the topology setup.
    type(TGFFUserInput),allocatable     :: userinput
  contains
    procedure :: deallocate => gfnff_data_deallocate
    procedure :: type_reset => gfnff_data_reset_types
    procedure :: type_init => gfnff_data_make_types
    procedure :: singlepoint => gfnff_singlepoint_wrapper
    procedure :: hessian => gfnff_hessian_wrapper
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
    if (dat%cell%npbc > 0) then
      if (any(abs(lattice_loc-dat%cell%lattice) .gt. lthr)) then
        call dat%cell%init(lattice_loc)
        call dat%cell%init_wsc(nat,at,xyz)
      end if
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
  subroutine gfnff_hessian(nat,at,xyz,dat,hessian,energy,gradient, &
  &                        step,printlevel,printunit,iostat)
!**********************************************************************
!* GFN-FF Cartesian nuclear Hessian.
!*
!* INPUT:
!*   nat            - number of atoms
!*   at(nat)        - atomic numbers
!*   xyz(3,nat)     - Cartesian coordinates (Bohr)
!*   dat            - bundled GFN-FF data and settings, already initialised
!*   step           - optional finite-difference step (Bohr) for those terms
!*                    that are not yet available in closed form
!*   printlevel     - optional verbosity (0=silent,1=errors,2=info,3=verbose)
!*   printunit      - optional output unit (default: stdout)
!* OUTPUT:
!*   hessian        - (3*nat,3*nat) second derivatives (Eh/Bohr^2), with the
!*                    degree of freedom (c,A) stored at 3*(A-1)+c
!*   energy         - optional total energy at the input geometry (Eh)
!*   gradient       - optional analytic gradient at the input geometry
!*   iostat         - optional error status
!*
!* Terms already derived in closed form are evaluated analytically; the
!* remainder is finite-differenced from the analytic gradient of exactly
!* those remaining terms, so the result is always the complete Hessian.
!*********************************************************************
    implicit none
    !> INPUT
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: at(nat)    !> atom types
    real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    type(gfnff_data),intent(inout) :: dat
    real(wp),intent(in),optional :: step
    integer,intent(in),optional  :: printlevel
    integer,intent(in),optional  :: printunit
    !> OUTPUT
    real(wp),intent(out) :: hessian(3*nat,3*nat)
    real(wp),intent(out),optional :: energy
    real(wp),intent(out),optional :: gradient(3,nat)
    integer,intent(out),optional  :: iostat
    !> LOCAL
    integer :: io,mylevel,myunit
    real(wp) :: e_loc,efield_loc(3)
    real(wp),allocatable :: g_loc(:,:)

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

    hessian(:,:) = 0.0_wp
    e_loc = 0.0_wp
    efield_loc = 0.0_wp
    io = 0
    allocate (g_loc(3,nat),source=0.0_wp)

    call gfnff_hessian_core(mylevel,nat,at,xyz,dat%cell,dat%ichrg,dat%param, &
    &                       dat%topo,dat%neigh,dat%nlist,efield_loc,dat%solvation, &
    &                       dat%version,dat%accuracy,hessian,e_loc,g_loc, &
    &                       step=step,iostat=io,printunit=myunit)

    if (present(energy)) energy = e_loc
    if (present(gradient)) gradient = g_loc
    if (present(iostat)) iostat = io

  end subroutine gfnff_hessian

!========================================================================================!
  subroutine gfnff_hessian_wrapper(self,nat,at,xyz,hessian,energy,gradient, &
  &                                step,printlevel,printunit,iostat)
!******************************************************************
!* A wrapper to the Hessian routine, allowing it to be called
!* with "call dat%hessian(...)"
!******************************************************************
    implicit none
    class(gfnff_data) :: self
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in),optional :: step
    integer,intent(in),optional  :: printlevel
    integer,intent(in),optional  :: printunit
    real(wp),intent(out) :: hessian(3*nat,3*nat)
    real(wp),intent(out),optional :: energy
    real(wp),intent(out),optional :: gradient(3,nat)
    integer,intent(out),optional  :: iostat
    call gfnff_hessian(nat,at,xyz,self,hessian,energy=energy,gradient=gradient, &
    & step=step,printlevel=printlevel,printunit=printunit,iostat=iostat)
  end subroutine gfnff_hessian_wrapper

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
  &                 printlevel,printunit,version,iostat,ichrg,lattice,npbc,accuracy)
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
    !*   accuracy    - optional: cutoff/precision factor. Larger is looser
    !*                 and faster; above 1.0 the EEQ system is solved in
    !*                 single precision. Defaults to 0.1, or 2.0 above
    !*                 10000 atoms.
    !* OUTPUT:
    !*   iostat      - optional: error status (0=success)
    !*************************************************************
    use gfnff_param,only:gfnff_set_param
    use gfnff_api_setup,only:gfnff_setup
    use gfnff_gdisp0,only:newD3Model
    use gfnff_alpb,only:gfnff_gbsa_init
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
    !> loosens the interaction cutoffs and, above 1.0, switches the EEQ
    !> solve to single precision. Passed here rather than set on dat
    !> afterwards so topology and energy use the same thresholds.
    real(wp),intent(in),optional :: accuracy
    !> OUTPUT
    type(gfnff_data),intent(inout) :: dat
    !> LOCAL
    character(len=:),allocatable :: fname
    !> TOML parameter set: label, status, message, and whether one was used
    character(len=:),allocatable :: pname,perr
    integer :: pio
    logical :: toml_param
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

!> Hand any user-defined fragmentation to the neighbor list, so that the
!> topology setup will not create bonds across the supplied fragments.
    if (allocated(dat%userinput)) then
      if (allocated(dat%userinput%fraglist)) then
        if (size(dat%userinput%fraglist) == nat) then
          dat%neigh%user_fraglist = dat%userinput%fraglist
          if (mylevel >= 2) then
            write (myunit,'(10x,a,i0,a)') 'using user-defined fragmentation (', &
            &  maxval(dat%userinput%fraglist)-minval(dat%userinput%fraglist)+1, &
            &  ' groups); no bonds will be formed across them'
          end if
        else if (mylevel >= 1) then
          write (myunit,'("**WARNING** ",a,1x,a)') &
          & 'user fragment list size does not match nat, ignoring it.',source
        end if
      end if
    end if

!> Periodic boundary conditions setup
    if (present(npbc)) dat%cell%npbc = npbc
    if (present(lattice)) then
      call dat%cell%init(lattice)
      call dat%cell%init_wsc(nat,at,xyz)
    end if

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

!> Parametrisation version. Only the explicit argument overrides what the
!> object already holds: this used to reset dat%version to the default
!> whenever the optional argument was absent, so setting the field and then
!> calling init silently got the default instead. The field's own default
!> initialiser supplies angewChem2020_2 for a fresh object.
    if (present(version)) dat%version = version

!> mcGFN-FF is a periodic parametrisation and only makes sense as one. Its
!> four scaling factors were fitted for molecular crystals, and one of them,
!> the dispersion s8, is only ever read on the periodic branch of the energy
!> routine: there the D3 sum is split by fragment, with the intermolecular
!> pairs taking the mc s8 and the intramolecular ones the standard value.
!> A molecular system has no such split, so asking for mcGFN-FF without a
!> cell used to apply three of the four factors and quietly drop the fourth,
!> which is neither GFN-FF nor mcGFN-FF. Refuse the combination instead.
    if (dat%version == gffVersion%mcgfnff2023.and.dat%cell%npbc == 0) then
      if (mylevel >= 1) write (myunit,'("**ERROR** ",a,1x,a)') &
         & 'mcGFN-FF (mcgfnff2023) is a periodic parametrisation and needs a '// &
         & 'lattice; use angewChem2020_2 for molecular systems',source
      if (present(iostat)) iostat = 1
      return
    end if

!> Hand any host-supplied molecular graph to the neighbour list, which then
!> skips the distance criterion entirely. Validated loudly rather than
!> silently repaired: a graph is an assertion about the molecule, and a
!> malformed one would otherwise surface as a strange force field.
    if (allocated(dat%userinput)) then
      if (allocated(dat%userinput%bondmat)) then
        call check_bondmat(nat,dat%userinput%bondmat,dat%cell%npbc, &
           & dat%neigh%numnb,mylevel,myunit,io)
        if (io /= 0) then
          if (present(iostat)) iostat = io
          return
        end if
        dat%neigh%user_bondmat = dat%userinput%bondmat
        if (mylevel >= 2) write (myunit,'(10x,a,i0,a)') &
        & 'using host-supplied molecular graph (', &
        & count(dat%userinput%bondmat /= 0)/2,' bonds); no bond perception'
      end if
    end if

    call dat%topo%zero
    dat%update = .true.

!> Hand any host-supplied atomic reference charges to the topology (must come
!> after topo%zero, which deallocates topo components). These are summed per
!> fragment during setup to define the per-fragment EEQ charge constraint.
    if (allocated(dat%userinput)) then
      if (allocated(dat%userinput%refq)) then
        if (size(dat%userinput%refq) == nat) then
          dat%topo%refq = dat%userinput%refq
          if (mylevel >= 2) write (myunit,'(10x,a)') &
          & 'host-supplied atomic reference charges will set per-fragment charges'
        else if (mylevel >= 1) then
          write (myunit,'("**WARNING** ",a,1x,a)') &
          & 'reference charge array size does not match nat, ignoring it.',source
        end if
      end if
    end if

!> global accuracy factor similar to acc in xtb used in SCF
    if (present(accuracy)) then
      dat%accuracy = accuracy
    else
      dat%accuracy = 0.1_wp
      if (nat > 10000) then
        dat%accuracy = 2.0_wp
      end if
    end if

!> Obtain the parameter file or load internal
    toml_param = .false.
    if (allocated(dat%parametrisation)) then
      fname = dat%parametrisation
    else
      fname = 'no file!'
    end if
    inquire (file=fname,exist=ex)
    if (ex) then
      if (is_toml_name(fname)) then
        !> A TOML set may specify only some keys, so start from the internal
        !> parametrisation and let the file override what it names. That also
        !> means the fixed physical tables are always present.
        call gfnff_load_param(dat%version,dat%param,ex)
        if (.not.allocated(dat%gen)) allocate (dat%gen)
        call gfnff_set_param(nat,dat%gen,dat%param)
        !> gfnff_setup sets this alongside gfnff_set_param; since that call is
        !> skipped for a preset parametrisation, it has to be set here too
        dat%param%dispscale = 1.0_wp
        call param_read_toml(fname,dat%param,dat%gen,dat%version,pname,pio,perr)
        if (pio /= 0) then
          if (mylevel >= 1) write (myunit,'("**ERROR** ",a,1x,a)') perr,source
          if (present(iostat)) iostat = pio
          return
        end if
        if (mylevel >= 2 .and. len_trim(pname) > 0) &
           & write (myunit,'(10x,"parametrisation: ",a)') trim(pname)
        toml_param = .true.
      else
        open (newunit=ich,file=fname)
        call gfnff_read_param(ich,dat%param)
        close (ich)
      end if
    else !> no parameter file, try to load internal version
      call gfnff_load_param(dat%version,dat%param,ex)
      !> the error and the message are separate concerns: this used to be one
      !> condition, so at printlevel 0 an unknown version silently continued
      !> with an unpopulated parameter table instead of failing
      if (.not.ex) then
        if (mylevel >= 1) write (myunit, &
           & '("No internal parametrisation for version ",i0,", and no file ",a,1x,a)') &
           & dat%version,trim(fname),source
        if (present(iostat)) iostat = 1
        return
      end if
    end if

    call newD3Model(dat%topo%dispm,nat,at)

    call gfnff_setup(nat,at,xyz,dat%ichrg,mylevel,restart,dat%write_topo,toml_param, &
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

  subroutine check_bondmat(nat,bondmat,npbc,numnb,printlevel,myunit,io)
    !*************************************************************
    !* Validate a host-supplied molecular graph.
    !*
    !* Every check here is a condition the topology setup would
    !* otherwise violate silently: an asymmetric matrix gives two
    !* atoms different opinions about the same bond, a self-bond
    !* walks into the neighbour packing, and a degree above the
    !* neighbour-list width overruns the count slot.
    !*
    !* INPUT:
    !*   nat       - number of atoms
    !*   bondmat   - the graph, expected (nat,nat)
    !*   npbc      - number of periodic dimensions
    !*   numnb     - neighbour list width; degree must stay below it
    !*   printlevel/myunit - reporting
    !* OUTPUT:
    !*   io        - 0 on success, 1 on any violation
    !*************************************************************
    implicit none
    character(len=*),parameter :: source = 'check_bondmat'
    integer,intent(in) :: nat,npbc,numnb,printlevel,myunit
    integer,intent(in) :: bondmat(:,:)
    integer,intent(out) :: io
    integer :: i,deg

    io = 0
    if (size(bondmat,1) /= nat.or.size(bondmat,2) /= nat) then
      call fail('molecular graph must be (nat,nat)')
      return
    end if
    !> the graph carries no cell index, so an image bond cannot be expressed
    if (npbc /= 0) then
      call fail('a host-supplied molecular graph is molecular-only (npbc = 0)')
      return
    end if
    if (any(bondmat < 0)) then
      call fail('molecular graph has negative bond orders')
      return
    end if
    do i = 1,nat
      if (bondmat(i,i) /= 0) then
        call fail('molecular graph has a nonzero diagonal (atom bonded to itself)')
        return
      end if
    end do
    if (any(bondmat /= transpose(bondmat))) then
      call fail('molecular graph is not symmetric')
      return
    end if
    do i = 1,nat
      deg = count(bondmat(:,i) /= 0)
      if (deg > numnb-1) then
        call fail('molecular graph has an atom with more bonds than the '// &
           & 'neighbour list can hold')
        return
      end if
    end do

  contains
    subroutine fail(msg)
      character(len=*),intent(in) :: msg
      io = 1
      if (printlevel >= 1) write (myunit,'("**ERROR** ",a,1x,a)') msg,source
    end subroutine fail
  end subroutine check_bondmat

  subroutine gfnff_initialize_wrapper(self,nat,at,xyz, &
     &                 printlevel,printunit,version,iostat,ichrg,solvent,lattice,npbc,accuracy)
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
    real(wp),intent(in),optional :: accuracy     !> cutoff/precision factor

    if (present(solvent)) then
      if (solvent .ne. 'none'.and.len_trim(solvent) > 0) self%solvent = solvent
    end if

    call gfnff_initialize(nat,at,xyz,self, &
    &       printlevel=printlevel,printunit=printunit, &
    &       version=version,iostat=iostat,ichrg=ichrg,lattice=lattice,npbc=npbc, &
    &       accuracy=accuracy)
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
    !> neigh%blist, not topo%blist: the latter is never filled, so this used
    !> to hand back an all-zero bond order matrix.
    if (allocated(ff_dat%neigh)) then
      if (allocated(ff_dat%neigh%blist)) then
        do i = 1,ff_dat%neigh%nbond
          k = ff_dat%neigh%blist(1,i)
          l = ff_dat%neigh%blist(2,i)
          wbo(k,l) = 1.0_wp
          wbo(l,k) = wbo(k,l)
        end do
      end if
    end if
  end subroutine gfnff_get_fake_wbo

!========================================================================================!
!========================================================================================!

  !> A parameter file is treated as TOML when its name says so. The legacy
  !> fixed-column format has no marker of its own, so the extension is the
  !> only thing that can distinguish them without reading the file twice.
  pure function is_toml_name(fname) result(yes)
    character(len=*),intent(in) :: fname
    logical :: yes
    integer :: n
    n = len_trim(fname)
    yes = n > 5
    if (yes) yes = fname(n-4:n) == '.toml'
  end function is_toml_name

!========================================================================================!

  subroutine gfnff_write_parametrisation(dat,filename,name,iostat,errmsg)
    !***********************************************************************
    !* Write the parametrisation currently held by a calculator to a TOML
    !* file: the per-element table and the generator constants together.
    !*
    !* The calculator must have been initialised, since that is what fills
    !* param and gen. The file it produces is accepted back by setting
    !* dat%parametrisation to it, which is what makes a fitted or edited
    !* parameter set usable without rebuilding.
    !* Input:
    !*   dat      - an initialised calculator
    !*   filename - destination
    !*   name     - optional label recorded under [meta]
    !* Output:
    !*   iostat   - 0 on success
    !*   errmsg   - set when iostat /= 0
    !***********************************************************************
    type(gfnff_data),intent(in) :: dat
    character(len=*),intent(in) :: filename
    character(len=*),intent(in),optional :: name
    integer,intent(out) :: iostat
    character(len=:),allocatable,intent(out) :: errmsg

    if (.not.allocated(dat%param).or..not.allocated(dat%gen)) then
      iostat = 1
      errmsg = 'calculator holds no parametrisation; initialise it first'
      return
    end if
    call param_write_toml(filename,dat%param,dat%gen,dat%version,name,iostat,errmsg)
  end subroutine gfnff_write_parametrisation

!========================================================================================!

  !> Whether this build can read and write TOML parameter files.
  pure function gfnff_parametrisation_io_available() result(yes)
    logical :: yes
    yes = param_io_available()
  end function gfnff_parametrisation_io_available

end module gfnff_interface

