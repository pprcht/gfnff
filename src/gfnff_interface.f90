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

  subroutine gfnff_singlepoint(nat,at,xyz,dat,energy,gradient,verbose,iostat)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: at(nat)    !> atom types
    real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    logical,intent(in),optional    :: verbose  !> printout activation 
    type(gfnff_data),intent(inout) :: dat  !> collection of gfnff datatypes and settings
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    integer,intent(out),optional  :: iostat
    !> LOCAL
    integer :: io
    logical :: pr

    !> printout activation via verbosity
    if(present(verbose))then
      pr = verbose
    else
      pr =.false. !> (there is close to no printout anyways)
    endif

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

  subroutine gfnff_singlepoint_wrapper(self,nat,at,xyz,energy,gradient,verbose,iostat)
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
    logical,intent(in),optional    :: verbose  !> printout activation 
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    integer,intent(out),optional  :: iostat
    call gfnff_singlepoint(nat,at,xyz,self,energy,gradient,verbose=verbose,iostat=iostat)
  end subroutine gfnff_singlepoint_wrapper

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
    write (iunit,'(a)') repeat('-',50)
  end subroutine print_gfnff_results

  subroutine gfnff_print_results_wrapper(self,iunit)
    implicit none
    class(gfnff_data) :: self
    !> INPUT
    integer,intent(in),optional :: iunit
    !> LOCAL
    integer :: myunit 

    if(present(iunit))then
       myunit = iunit
    else
       myunit = stdout
    endif
    if(allocated(self%res))then
    call print_gfnff_results(myunit,self%res,allocated(self%solvation))
    endif
  end subroutine gfnff_print_results_wrapper
!========================================================================================!

  subroutine gfnff_initialize(nat,at,xyz,dat, &
  &                 print,verbose,iunit,version,iostat,ichrg)
    use gfnff_param
    use gfnff_setup_mod,only:gfnff_setup
    use gfnff_gdisp0,only:newD3Model
    use gfnff_gbsa
    character(len=*),parameter :: source = 'gfnff_initialize'
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    logical,intent(in),optional  :: print
    logical,intent(in),optional  :: verbose
    integer,intent(in),optional  :: iunit
    integer,intent(in),optional  :: version
    integer,intent(out),optional :: iostat
    integer,intent(in),optional  :: ichrg 
    !> OUTPUT
    type(gfnff_data),intent(inout) :: dat
    !> LOCAL
    character(len=:),allocatable :: fname
    integer :: ich,io,myunit
    logical :: ex,okbas,pr,pr2
    logical :: exitRun
    logical :: restart

!> mapping of optional instuctions
    if(present(print))then
      pr = print
    else
      pr = .false.
    endif
    if(present(verbose))then
      pr2 = verbose
    else
      pr2 = .false.
    endif
    if(pr2) pr = pr2    
    if(present(iunit))then
      myunit = iunit
    else
      myunit = stdout
    endif

!> Reset datatypes
    call dat%type_init()
    if(present(ichrg))then
      dat%ichrg = ichrg
    endif
!> except restart-related options
    restart = dat%restart
    if(.not.allocated(dat%restartfile))then
      dat%topo%filename = 'gfnff_topo'
    else
      dat%topo%filename = dat%restartfile
    endif
    if(allocated(dat%refgeo)) restart = .false.
    if(allocated(dat%refcharges))then
      dat%topo%refcharges = dat%refcharges
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
    if(allocated(dat%parametrisation))then
      fname=dat%parametrisation
    else
      fname='no file!'
    endif
    inquire (file=fname,exist=ex)
    if (ex) then
      open (newunit=ich,file=fname)
      call gfnff_read_param(ich,dat%param)
      close (ich)
    else !> no parameter file, try to load internal version
      call gfnff_load_param(dat%version,dat%param,ex)
      if (.not.ex.and.pr) then
        write (myunit,'("Parameter file ",a," not found!",a)') fname,source
        return
      end if
    end if

    call newD3Model(dat%topo%dispm,nat,at)

    call gfnff_setup(nat,at,xyz,dat%ichrg,pr,restart,dat%write_topo, &
    &        dat%gen,dat%param,dat%topo,dat%accuracy,dat%version,io, &
    &        verbose=verbose, iunit=myunit)

    !> Optional, ALPB solvation
    if (allocated(dat%solvent)) then
      if (.not. (allocated(dat%solvation))) allocate (dat%solvation)
      call gfnff_gbsa_init(nat,at,dat%solvent,dat%solvation)
      if (pr)then
        write(myunit,*)
        call gfnff_gbsa_print(dat%solvation,myunit)
      endif
    end if

    if ((io /= 0).and.pr) then
      write (myunit,'("Could not create force field calculator ",a)') source
    end if
    if (present(iostat)) then
      iostat = io
    end if
  end subroutine gfnff_initialize

  subroutine gfnff_initialize_wrapper(self,nat,at,xyz, &
     &                 print,verbose,iunit,version,iostat,ichrg,solvent)
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
    logical,intent(in),optional  :: print
    logical,intent(in),optional  :: verbose
    integer,intent(in),optional  :: iunit
    integer,intent(in),optional  :: version
    integer,intent(out),optional :: iostat
    integer,intent(in),optional  :: ichrg
    character(len=*),intent(in),optional :: solvent

    if(present(solvent))then
      if(solvent.ne.'none'.and.len_trim(solvent)>0) self%solvent=solvent
    endif

    call  gfnff_initialize(nat,at,xyz,self, &
    &       print=print,verbose=verbose,iunit=iunit,&
    &       version=version,iostat=iostat,ichrg=ichrg)
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

!=========================================================================================!
  subroutine gfnff_get_fake_wbo(ff_dat,nat,wbo)
    implicit none
    type(gfnff_data),intent(in) :: ff_dat
    integer,intent(in) :: nat
    real(wp),intent(out) :: wbo(nat,nat)
    integer :: i,k,l
    wbo = 0.0_wp
    if(allocated(ff_dat%topo))then
      if(allocated(ff_dat%topo%blist))then
         do i =1,ff_dat%topo%nbond
            k = ff_dat%topo%blist(1,i)
            l = ff_dat%topo%blist(2,i)
            wbo(k,l) = 1.0_wp
            wbo(l,k) = wbo(k,l)
         enddo
      endif
    endif
  end subroutine gfnff_get_fake_wbo

!========================================================================================!
end module gfnff_interface

