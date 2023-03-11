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
!--------------------------------------------------------------------------------!
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert, Sebastian Spicher, Stefan Grimme
!> at https://github.com/grimme-lab/xtb
!================================================================================!
module gfnff_setup_mod
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_ini_mod,only:gfnff_ini
  use gfnff_data_types,only:TGFFData,TGFFTopology,TGFFGenerator
  implicit none
  private
  public :: gfnff_setup,gfnff_input

!========================================================================================!
!========================================================================================!
contains   !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine gfnff_setup(nat,at,xyz,ichrg,verbose,restart,write_topo, &
  &                      gen,param,topo,accuracy,version,io)
    use gfnff_restart
    use gfnff_param,only:ini,gfnff_set_param
    implicit none
    character(len=*),parameter :: source = 'gfnff_setup'
! Dummy
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: ichrg

    !integer,intent(in) :: ich
    type(TGFFTopology),intent(inout) :: topo
    type(TGFFGenerator),intent(inout) :: gen
    type(TGFFData),intent(inout) :: param
    integer,intent(in) :: version
    logical,intent(in) :: restart
    logical,intent(in) :: verbose
    logical,intent(in) :: write_topo
    real(wp),intent(in) :: accuracy
    integer,intent(out) :: io

! Stack
    integer :: newichrg
    logical :: ex
    logical :: success
    logical :: exitRun

    io = 0

    newichrg = ichrg
    call gfnff_input(nat,at,xyz,newichrg,topo)

    call gfnff_set_param(nat,gen,param)
    param%dispscale = 1.0_wp
    if (restart) then
      inquire (file='gfnff_topo',exist=ex)
      if (ex) then
        call read_restart_gff('gfnff_topo',nat,version,success,.true.,topo)
        if (success) then
          write (stdout,'(10x,"GFN-FF topology read from file successfully!")')
          return
        else
          write (stdout,'("Could not read topology file.",a)') source
          exitRun = .true.
          if (exitRun) then
            return
          end if
        end if
      end if
    end if

    call gfnff_ini(verbose,ini,nat,at,xyz,ichrg,gen,param,topo,accuracy,io)
    if (io /= 0) then
      write (stdout,'("Failed to generate topology ",a)') source
      return
    end if

    if (write_topo) then
      call write_restart_gff('gfnff_topo',nat,version,topo)
      call write_gfnff_adjacency('gfnff_adjacency',topo)
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
    integer           :: i,j,k
    integer           :: ni
    integer           :: ns
    integer           :: nf
    integer           :: ich
    integer           :: iatom
    integer           :: iresidue
    integer           :: ifrag
    integer           :: ibond
    integer           :: bond_ij(3)
    real(wp)          :: r
    real(wp)          :: dum1
    real(wp)          :: floats(10)
    logical           :: ex
    character(len=80) :: atmp
    character(len=80) :: s(10)
    integer,allocatable :: rn(:)

    if (.not.allocated(topo%nb)) allocate (topo%nb(20,nat),source=0)
    if (.not.allocated(topo%qfrag)) allocate (topo%qfrag(nat),source=0.0d0)
    if (.not.allocated(topo%fraglist)) allocate (topo%fraglist(nat),source=0)

!  select case(mol%ftype)
!  !--------------------------------------------------------------------
!  ! PDB case
!  case(fileType%pdb)
!    ini = .true.
!    ifrag=0
!    allocate(rn(mol%n))
!    rn(:) = mol%pdb%residue_number
!    do iresidue = minval(rn),maxval(rn)
!      if (any(iresidue .eq. rn)) then
!        ifrag=ifrag+1
!        where(iresidue .eq. rn) topo%fraglist = ifrag
!      end if
!    end do
!    deallocate(rn)
!    topo%nfrag = maxval(topo%fraglist)
!    if (.not.allocated(topo%qpdb)) allocate(topo%qpdb(mol%n))
!    do iatom=1,mol%n
!      topo%qfrag(topo%fraglist(iatom)) = topo%qfrag(topo%fraglist(iatom)) &
!        & + dble(mol%pdb(iatom)%charge)
!      topo%qpdb(iatom) = mol%pdb(iatom)%charge
!    end do
!    ichrg=idint(sum(topo%qfrag(1:topo%nfrag)))
!    write(env%unit,'(10x,"charge from pdb residues: ",i0)') ichrg
!  !--------------------------------------------------------------------
!  ! SDF case
!  case(fileType%sdf,fileType%molfile)
!    ini = .false.
!    topo%nb=0
!    topo%nfrag=0
!    do ibond = 1, len(mol%bonds)
!      call mol%bonds%get_item(ibond,bond_ij)
!      i = bond_ij(1)
!      j = bond_ij(2)
!      ni=topo%nb(20,i)
!      ex=.false.
!      do k=1,ni
!        if(topo%nb(k,i).eq.j) then
!          ex=.true.
!          exit
!        endif
!      enddo
!      if(.not.ex)then
!        topo%nb(20,i)=topo%nb(20,i)+1
!        topo%nb(topo%nb(20,i),i)=j
!        topo%nb(20,j)=topo%nb(20,j)+1
!        topo%nb(topo%nb(20,j),j)=i
!      endif
!    end do
!    do i=1,mol%n
!      if(topo%nb(20,i).eq.0)then
!        dum1=1.d+42
!        do j=1,i
!          r=sqrt(sum((mol%xyz(:,i)-mol%xyz(:,j))**2))
!          if(r.lt.dum1.and.r.gt.0.001)then
!            dum1=r
!            k=j
!          endif
!        enddo
!        topo%nb(20,i)=1
!        topo%nb(1,i)=k
!      endif
!    end do
!    ! initialize qfrag as in the default case
!    topo%qfrag(1)=mol%chrg
!    topo%qfrag(2:mol%n)=0
!  !--------------------------------------------------------------------
!  ! General case: input = xyz or coord
!  case default
!    if (mol%npbc > 0) then
!      call env%error("Input file format not suitable for GFN-FF!")
!      return
!    end if
    ini = .true.
    inquire (file='.CHRG',exist=ex)
    if (ex) then
      open (newunit=ich,file='.CHRG')
      read (ich,'(a)') atmp
      close (ich)
      !call readline(atmp,floats,s,ns,nf)
      dum1 = -huge(dum1)/2.0_wp
      floats(:) = dum1
      read (atmp,*) floats(1:10)
      nf = 0
      do i = 1,10
        if (floats(i) .gt. dum1) nf = nf+1
      end do
      topo%qfrag(1:nf) = floats(1:nf)
      ichrg = int(sum(topo%qfrag(1:nf)))
      topo%qfrag(nf+1:nat) = 9999
    else
      topo%qfrag(1) = ichrg
      topo%qfrag(2:nat) = 0
    end if
!  end select
  end subroutine gfnff_input
!========================================================================================!

  subroutine write_gfnff_adjacency(fname,topo)
    implicit none
    character(len=*),intent(in) :: fname
    integer :: ifile ! file handle
    type(TGFFTopology) :: topo
    integer :: i,j

    open (newunit=ifile,file=fname)
    ! looping over topology neighboring list
    if (ifile .ne. -1) then
      write (ifile,'(a)') '# indices of neighbouring atoms (max seven)'
      do i = 1,size(topo%nb,2)
        write (ifile,'(*(i0:, 1x))') (topo%nb(j,i),j=1,topo%nb(size(topo%nb,1),i))
      end do
    end if
    close (ifile)
  end subroutine write_gfnff_adjacency

!========================================================================================!
end module gfnff_setup_mod
