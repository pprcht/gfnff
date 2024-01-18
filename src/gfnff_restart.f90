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
module gfnff_restart
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_data_types
  use gfnff_param
  implicit none
  private
  public :: write_restart_gff,read_restart_gff

  integer,parameter :: i8 = selected_int_kind(18)

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine write_restart_gff(fname,nat,version,topo)
    implicit none
    type(TGFFTopology),intent(in) :: topo
    character(len=*),intent(in) :: fname
    integer,intent(in)  :: nat
    integer,intent(in)  :: version
    integer :: ich ! file handle
    integer :: err

    open (file=fname,newunit=ich,action='write',form='unformatted',iostat=err)

!>--- Dimensions
    write (ich) int(version,i8),int(nat,i8)
    write (ich) topo%nbond,topo%nangl,topo%ntors,   &
             &  topo%nathbH,topo%nathbAB,topo%natxbAB,topo%nbatm, &
             &  topo%nfrag,topo%nsystem,topo%maxsystem
    write (ich) topo%nbond_blist,topo%nbond_vbond,topo%nangl_alloc, &
             &  topo%ntors_alloc,topo%bond_hb_nr,topo%b_max
!>--- Arrays Integers
    write (ich) topo%nb,topo%bpair,topo%blist,topo%alist, &
             &  topo%tlist,topo%b3list,topo%fraglist,topo%hbatHl,topo%hbatABl, &
             &  topo%xbatABl,topo%ispinsyst,topo%nspinsyst,topo%bond_hb_AH, &
             &  topo%bond_hb_B,topo%bond_hb_Bn,topo%nr_hb
!>--- Arrays Reals
    write (ich) topo%vbond,topo%vangl,topo%vtors,topo%chieeq, &
             &  topo%gameeq,topo%alpeeq,topo%alphanb,topo%qa,  &
             &  topo%xyze0,topo%zetac6, &
             &  topo%qfrag,topo%hbbas,topo%hbaci
    close (ich)
  end subroutine write_restart_gff
!========================================================================================!

  subroutine read_restart_gff(fname,n,version,success,verbose,topo)
    implicit none
    character(len=*),parameter :: source = 'restart_read_restart_gff'
    type(TGFFTopology),intent(inout) :: topo
    character(len=*),intent(in) :: fname
    integer,intent(in)  :: n
    integer,intent(in)  :: version
    logical,intent(out) :: success
    logical,intent(in)  :: verbose

    integer(i8) :: iver8,nat8

    integer :: ich ! file handle
    integer :: err
    logical :: exist

    success = .false.
    open (file=fname,newunit=ich,status='old',action='read', &
         & form='unformatted',iostat=err)
    if (ich .ne. -1) then
!>--- read the first byte, which identify the calculation specs
      read (ich,iostat=err) iver8,nat8
      if (err .eq. 0) then
        if (iver8 .ne. int(version,i8).and.verbose) &
           &  write (stdout,'("Version number missmatch in restart file.",a)') source
        if (nat8 .ne. n.and.verbose) then
          write (stdout,'("Atom number missmatch in restart file.",a)') source
          success = .false.
          close (ich)
          return
        else if (iver8 .eq. int(version)) then
          success = .true.
          read (ich) topo%nbond,topo%nangl,topo%ntors, &
                  &  topo%nathbH,topo%nathbAB,topo%natxbAB,topo%nbatm, &
                  &  topo%nfrag,topo%nsystem,topo%maxsystem
          read (ich) topo%nbond_blist,topo%nbond_vbond,topo%nangl_alloc, &
                  &  topo%ntors_alloc,topo%bond_hb_nr,topo%b_max
!>--- allocate some memory now
          call gfnff_param_alloc(topo,n)
          if (.not.allocated(topo%ispinsyst)) allocate (topo%ispinsyst(n,topo%maxsystem),source=0)
          if (.not.allocated(topo%nspinsyst)) allocate (topo%nspinsyst(topo%maxsystem),source=0)
          read (ich) topo%nb,topo%bpair,topo%blist,topo%alist, &
             & topo%tlist,topo%b3list,topo%fraglist,topo%hbatHl,topo%hbatABl, &
             & topo%xbatABl,topo%ispinsyst,topo%nspinsyst,topo%bond_hb_AH, &
             & topo%bond_hb_B,topo%bond_hb_Bn,topo%nr_hb
          read (ich) topo%vbond,topo%vangl,topo%vtors,topo%chieeq, &
             & topo%gameeq,topo%alpeeq,topo%alphanb,topo%qa, &
             & topo%xyze0,topo%zetac6,&
             & topo%qfrag,topo%hbbas,topo%hbaci
        else
          if (verbose) &
            write (stdout,'("Dimension missmatch in restart file.",a)') source
          success = .false.
        end if
      else
        if (verbose) &
          write (stdout,'("Dimension missmatch in restart file.",a)') source
        success = .false.
      end if
      close (ich)
    end if

  end subroutine read_restart_gff

!========================================================================================!
end module gfnff_restart

