!================================================================================!
! This file is part of gfnff.
!
! Copyright (C) 2023 Philipp Pracht
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
!================================================================================!
module gfnff_interface
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_data_types
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

!> this type bundles together most of the
!> data required for a GFN-FF calculation
  type :: gfnff_data
    integer  :: ichrg = 0  !> total molecular charge
    real(wp) :: accuracy = 0.1_wp
    logical  :: make_chrg = .true.
    integer  :: version = gffVersion%angewChem2020_2
    logical  :: update = .true.
    logical  :: write_topo = .true.
    character(len=:),allocatable :: solvent

    type(TGFFGenerator),allocatable     :: gen
    type(TGFFData),allocatable          :: param
    type(TGFFTopology),allocatable      :: topo
    type(TGFFNeighbourList),allocatable :: nlist
    type(TBorn),allocatable             :: solvation
    type(gfnff_results),allocatable     :: res
  contains
    procedure deallocate => gfnff_data_deallocate
    procedure type_reset => gfnff_data_reset_types
    procedure type_init => gfnff_data_make_types
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

  subroutine gfnff_singlepoint(nat,at,xyz,dat,energy,gradient,pr,iostat)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: at(nat)    !> atom types
    real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    logical,intent(in)  :: pr         !> printout activation
    type(gfnff_data),intent(inout) :: dat  !> collection of gfnff datatypes and settings
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    integer,intent(out),optional  :: iostat
    !> LOCAL
    integer :: io

    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    io = 0

    call gfnff_eg(pr,nat,dat%ichrg,at,xyz,dat%make_chrg,gradient,energy,   &
    &             dat%res,dat%param,dat%topo,dat%nlist,dat%solvation, &
    &             dat%update,dat%version,dat%accuracy,io)

    if (present(iostat)) then
      iostat = io
    end if

  end subroutine gfnff_singlepoint
!========================================================================================!

  subroutine print_gfnff_results(iunit,res_gff,lsolv)
    integer,intent(in) :: iunit ! file handle (usually output_unit=6)
    type(gfnff_results),intent(in) :: res_gff
    logical,intent(in) :: lsolv
    character(len=*),parameter :: outfmt = &
                                  '(2x,a,f23.12,1x,a)'
    write (iunit,outfmt) "total energy      ",res_gff%e_total,"Eh   "
    write (iunit,outfmt) "gradient norm     ",res_gff%gnorm,"Eh/a0"
    write (iunit,'(a)') repeat('-',50)
    write (iunit,outfmt) "bond energy       ",res_gff%e_bond,"Eh   "
    write (iunit,outfmt) "angle energy      ",res_gff%e_angl,"Eh   "
    write (iunit,outfmt) "torsion energy    ",res_gff%e_tors,"Eh   "
    write (iunit,outfmt) "repulsion energy  ",res_gff%e_rep,"Eh   "
    write (iunit,outfmt) "electrostat energy",res_gff%e_es,"Eh   "
    write (iunit,outfmt) "dispersion energy ",res_gff%e_disp,"Eh   "
    write (iunit,outfmt) "HB energy         ",res_gff%e_hb,"Eh   "
    write (iunit,outfmt) "XB energy         ",res_gff%e_xb,"Eh   "
    write (iunit,outfmt) "bonded atm energy ",res_gff%e_batm,"Eh   "
    write (iunit,outfmt) "external energy   ",res_gff%e_ext,"Eh   "
    if (lsolv) then
      write (iunit,'(2x,a)') repeat('-',44)
      write (iunit,outfmt) "-> Gsolv          ",res_gff%g_solv,"Eh   "
      write (iunit,outfmt) "   -> Gborn       ",res_gff%g_born,"Eh   "
      write (iunit,outfmt) "   -> Gsasa       ",res_gff%g_sasa,"Eh   "
      write (iunit,outfmt) "   -> Ghb         ",res_gff%g_hb,"Eh   "
      write (iunit,outfmt) "   -> Gshift      ",res_gff%g_shift,"Eh   "
    end if
  end subroutine print_gfnff_results
!========================================================================================!

  subroutine gfnff_initialize(nat,at,xyz,dat,fname,restart,pr,version,iostat,ichrg)
    use gfnff_param
    use gfnff_setup_mod,only:gfnff_setup
    use gfnff_gdisp0,only:newD3Model
    use gfnff_gbsa
    character(len=*),parameter :: source = 'gfnff_initialize'
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    character(len=*),intent(in) :: fname
    logical,intent(in) :: restart
    logical,intent(in) :: pr
    integer,intent(in),optional  :: version
    integer,intent(out),optional :: iostat
    integer,intent(in),optional  :: ichrg 
    !> OUTPUT
    type(gfnff_data),intent(inout) :: dat
    !> LOCAL
    integer :: ich,io
    logical :: ex,okbas
    logical :: exitRun

    !> Reset datatypes
    call dat%type_init()
    if(present(ichrg))then
      dat%ichrg = ichrg
    endif

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
    inquire (file=fname,exist=ex)
    if (ex) then
      open (newunit=ich,file=fname)
      call gfnff_read_param(ich,dat%param)
      close (ich)
    else !> no parameter file, try to load internal version
      call gfnff_load_param(dat%version,dat%param,ex)
      if (.not.ex) then
        write (stdout,'("Parameter file ",a," not found!",a)') fname,source
        return
      end if
    end if

    call newD3Model(dat%topo%dispm,nat,at)

    call gfnff_setup(nat,at,xyz,dat%ichrg,pr,restart,dat%write_topo, &
    &         dat%gen,dat%param,dat%topo,dat%accuracy,dat%version,io)

    !> Optional, ALPB solvation
    if (allocated(dat%solvent)) then
      if (.not. (allocated(dat%solvation))) allocate (dat%solvation)
      call gfnff_gbsa_init(nat,at,dat%solvent,dat%solvation)
      if (pr) call gfnff_gbsa_print(dat%solvation,stdout)
    end if

    if ((io /= 0).and.pr) then
      write (stdout,'("Could not create force field calculator ",a)') source
    end if
    if (present(iostat)) then
      iostat = io
    end if
  end subroutine gfnff_initialize

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
    allocate (self%nlist)
    allocate (self%res)
  end subroutine gfnff_data_make_types

!========================================================================================!
end module gfnff_interface

