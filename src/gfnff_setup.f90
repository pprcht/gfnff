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
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert, Sebastian Spicher, Stefan Grimme
!> at https://github.com/grimme-lab/xtb
! ──────────────────────────────────────────────────────────────────────────────
module gfnff_setup_mod
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_ini_mod,only:gfnff_ini
  use gfnff_data_types,only:TGFFData,TGFFTopology,TGFFGenerator,TCell
  use gfnff_neighbor,only:TNeigh
  implicit none
  private
  public :: gfnff_setup,gfnff_input

!========================================================================================!
!========================================================================================!
contains   !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine gfnff_setup(nat,at,xyz,ichrg,printlevel,restart,write_topo, &
  &                      gen,param,topo,neigh,cell,accuracy,version,io,printunit)
    !***********************************************************************
    !* Orchestrate GFN-FF topology setup and parameter loading.
    !* Input:
    !*   printlevel  - verbosity (0=silent, 1=errors, 2=info, 3=verbose)
    !*   restart     - attempt to read topology from file
    !*   write_topo  - write topology to file after setup
    !*   printunit   - output unit (optional, default: stdout)
    !* Output:
    !*   io          - error status (0 = success)
    !***********************************************************************
    use gfnff_restart
    use gfnff_param,only:ini,gfnff_set_param
    implicit none
    character(len=*),parameter :: source = 'gfnff_setup'
!> Dummy
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: ichrg

    type(TGFFTopology),intent(inout) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TCell),intent(inout) :: cell
    type(TGFFGenerator),intent(inout) :: gen
    type(TGFFData),intent(inout) :: param
    integer,intent(in) :: version
    logical,intent(in) :: restart
    integer,intent(in) :: printlevel   !< verbosity (0=silent,1=errors,2=info,3=verbose)
    logical,intent(in) :: write_topo
    real(wp) :: efield(3) = 0.0_wp
    real(wp),intent(in) :: accuracy
    integer,intent(out) :: io
    integer,intent(in),optional :: printunit  !< output unit (default: stdout)

!> Stack
    integer :: newichrg,myunit
    logical :: ex,success,exitRun

!> initialize
    io = 0
    if (present(printunit)) then
      myunit = printunit
    else
      myunit = stdout
    end if

    newichrg = ichrg
    call gfnff_input(nat,at,xyz,newichrg,topo)

! ──────────────────────────────────────────────────────────────────────────────
!> Load parameters
! ──────────────────────────────────────────────────────────────────────────────
    call gfnff_set_param(nat,gen,param)
    param%dispscale = 1.0_wp
    if (restart) then
      inquire (file=topo%filename,exist=ex)
      if (ex) then
        call read_restart_gff(topo%filename,nat,version,success,.true.,topo)
        if (success) then
          if (printlevel >= 2) write (myunit,'(/,"> GFN-FF topology read successfully from file ",a," !")') &
          & topo%filename
          return
        else
          if (printlevel >= 1) write (myunit,'("**ERROR** Could not read topology file. ",a)') source
          exitRun = .true.
          if (exitRun) then
            return
          end if
        end if
      end if
    end if

! ──────────────────────────────────────────────────────────────────────────────
!> Run initialization
! ──────────────────────────────────────────────────────────────────────────────
    call gfnff_ini(printlevel,ini,nat,at,xyz,ichrg, &
    &                  gen,param,topo,neigh,cell,efield, &
    &                  accuracy,io,printunit=myunit)

    if (io /= 0) then
      if (printlevel >= 1) write (myunit,'("Failed to generate topology ",a)') source
      return
    end if

    if (write_topo) then
      call write_restart_gff('gfnff_topo',nat,version,topo)
      !call write_gfnff_adjacency('gfnff_adjacency',topo)
    end if

  end subroutine gfnff_setup
!========================================================================================!

  subroutine gfnff_input(nat,at,xyz,ichrg,topo)
    use gfnff_param
    implicit none
    ! Dummy
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(inout)  :: ichrg
    type(TGFFTopology),intent(inout) :: topo
    ! Stack

    if (.false.) write (*,*) at,xyz ! silences -Wunused-dummy-argument

    !if (.not.allocated(topo%nb)) allocate (topo%nb(20,nat),source=0)
    if (.not.allocated(topo%qfrag)) allocate (topo%qfrag(nat),source=0.0d0)
    if (.not.allocated(topo%fraglist)) allocate (topo%fraglist(nat),source=0)

!> We actually do not want to read .CHRG if we have multiple instances of GFN-FF running at once
!> E.g. within an ONIOM setup.
!> Could REALLLY mess things up.....

    ini = .true.
!    inquire (file='.CHRG',exist=ex)
!    if (ex) then
!      open (newunit=ich,file='.CHRG')
!      read (ich,'(a)') atmp
!      close (ich)
!      !call readline(atmp,floats,s,ns,nf)
!      dum1 = -huge(dum1)/2.0_wp
!      floats(:) = dum1
!      read (atmp,*) floats(1:10)
!      nf = 0
!      do i = 1,10
!        if (floats(i) .gt. dum1) nf = nf+1
!      end do
!      topo%qfrag(1:nf) = floats(1:nf)
!      ichrg = int(sum(topo%qfrag(1:nf)))
!      topo%qfrag(nf+1:nat) = 9999
!    else
    topo%qfrag(1) = ichrg
    topo%qfrag(2:nat) = 0
!    end if

  end subroutine gfnff_input
!========================================================================================!

!  subroutine write_gfnff_adjacency(fname,topo)
!    implicit none
!    character(len=*),intent(in) :: fname
!    integer :: ifile ! file handle
!    type(TGFFTopology) :: topo
!    integer :: i,j
!
!    open (newunit=ifile,file=fname)
!    ! looping over topology neighboring list
!    if (ifile .ne. -1) then
!      write (ifile,'(a)') '# indices of neighbouring atoms (max seven)'
!      do i = 1,size(topo%nb,2)
!        write (ifile,'(*(i0:, 1x))') (topo%nb(j,i),j=1,topo%nb(size(topo%nb,1),i))
!      end do
!    end if
!    close (ifile)
!  end subroutine write_gfnff_adjacency

!========================================================================================!
end module gfnff_setup_mod
