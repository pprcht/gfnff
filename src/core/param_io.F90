! ------------------------------------------------------------------------------
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
! ------------------------------------------------------------------------------
!> Reading and writing a complete GFN-FF parametrisation as TOML: an [elements]
!> table (elements_to_table) and a [generator] table (generator_to_table).
!> param_write_toml/param_read_toml round-trip TGFFData/TGFFGenerator
!> (test/test_param_io.f90). Fixed tables (en, rad, rcov, metal, group, normcn,
!> repz) are physical constants, not fitted, and are not written. Without
!> toml-f both entry points return nonzero.
!>
!> THREAD SAFETY: toml-f's real-to-string conversion is not thread safe; two
!> threads serialising at once can corrupt the file. Every call here runs in one
!> shared critical section, and every write is read back immediately so that a
!> corrupted file fails at write time rather than later. Calling from an OpenMP
!> parallel region remains unsupported; ifx still corrupts despite the lock.
module gfnff_param_io
  use iso_fortran_env,only:wp => real64,sp => real32
  use gfnff_data_types,only:TGFFData,TGFFGenerator
#ifdef WITH_TOMLF
  use tomlf,only:toml_table,toml_array,toml_parse,toml_serialize, &
    &            toml_error,get_value,set_value,add_table,add_array, &
    &            toml_stat, &
    &            len
#endif
  implicit none
  private

  public :: param_write_toml,param_read_toml
  public :: param_io_available

  !> elements covered by the per-element table; the arrays are dimensioned 103
  integer,parameter :: nel = 86

contains  !> MODULE PROCEDURES START HERE

  pure function param_io_available() result(yes)
    !***********************************************************************
    !* Whether the library was built with TOML support.
    !***********************************************************************
    logical :: yes
#ifdef WITH_TOMLF
    yes = .true.
#else
    yes = .false.
#endif
  end function param_io_available

  subroutine param_write_toml(filename,param,gen,version,name,iostat,errmsg)
    !***********************************************************************
    !* Write a complete parametrisation (param, gen) to filename as TOML,
    !* overwriting it if present. name is an optional [meta] label.
    !***********************************************************************
    character(len=*),intent(in) :: filename
    type(TGFFData),intent(in) :: param
    type(TGFFGenerator),intent(in) :: gen
    integer,intent(in) :: version
    character(len=*),intent(in),optional :: name
    integer,intent(out) :: iostat
    character(len=:),allocatable,intent(out) :: errmsg

#ifdef WITH_TOMLF
    type(toml_table) :: table
    type(toml_table),pointer :: child
    character(len=:),allocatable :: serialized
    integer :: unit,err

    iostat = 0
    errmsg = ''
    table = toml_table()

    call add_table(table,'meta',child)
    if (present(name)) then
      call set_value(child,'name',name)
    else
      call set_value(child,'name','gfnff parametrisation')
    end if
    call set_value(child,'version',version)

    call add_table(table,'generator',child)
    call generator_to_table(gen,child)

    call add_table(table,'elements',child)
    call elements_to_table(param,child)

    !>-- serialise, write, and read back to catch corruption (THREAD SAFETY above)
    block
      open (newunit=unit,file=filename,action='write',status='replace',iostat=err)
      if (err /= 0) then
        iostat = err
        errmsg = 'could not open '//filename//' for writing'
        return
      end if
      write (unit,'(a)') '# GFN-FF parametrisation, written by gfnff_param_io.'
      write (unit,'(a)') '# Both halves are here: the per-element table and the'
      write (unit,'(a)') '# global generator constants. Physical constants (en,'
      write (unit,'(a)') '# rad, rcov, metal, group, normcn, repz) are not, they'
      write (unit,'(a)') '# come from gfnff_param_tables.'
      !$omp critical(gfnff_tomlf)
      serialized = toml_serialize(table)
      !$omp end critical(gfnff_tomlf)
      write (unit,'(a)') serialized
      close (unit)

      iostat = 0
      errmsg = ''
      call verify_readable(filename,iostat,errmsg)
    end block
#else
    iostat = 1
    errmsg = 'gfnff was built without toml-f; TOML parameter files unavailable'
    if (.false.) then   !> silences unused-argument warnings
      print *,filename,param%chi(1),gen%linthr,version,present(name)
    end if
#endif
  end subroutine param_write_toml

  subroutine param_read_toml(filename,param,gen,version,name,iostat,errmsg)
    !***********************************************************************
    !* Read a parametrisation written by param_write_toml; omitted keys keep their
    !* current value (overlay). param must be initialised: init(param,103).
    !***********************************************************************
    character(len=*),intent(in) :: filename
    type(TGFFData),intent(inout) :: param
    type(TGFFGenerator),intent(inout) :: gen
    integer,intent(inout) :: version
    character(len=:),allocatable,intent(out) :: name
    integer,intent(out) :: iostat
    character(len=:),allocatable,intent(out) :: errmsg

#ifdef WITH_TOMLF
    type(toml_table),allocatable :: table
    type(toml_table),pointer :: child
    type(toml_error),allocatable :: parse_error
    integer :: unit,err
    logical :: ex

    iostat = 0
    errmsg = ''
    name = ''

    inquire (file=filename,exist=ex)
    if (.not.ex) then
      iostat = 2
      errmsg = 'parameter file not found: '//filename
      return
    end if

    open (newunit=unit,file=filename,action='read',status='old',iostat=err)
    if (err /= 0) then
      iostat = err
      errmsg = 'could not open '//filename//' for reading'
      return
    end if
    !$omp critical(gfnff_tomlf)
    call toml_parse(table,unit,parse_error)
    !$omp end critical(gfnff_tomlf)
    close (unit)
    if (allocated(parse_error)) then
      iostat = 3
      errmsg = 'TOML parse error in '//filename//': '//parse_error%message
      return
    end if
    if (.not.allocated(table)) then
      iostat = 3
      errmsg = 'TOML parse produced no table for '//filename
      return
    end if

    call get_value(table,'meta',child,requested=.false.)
    if (associated(child)) then
      call get_value(child,'name',name,'')
      call get_value(child,'version',version,version)
    end if

    call get_value(table,'generator',child,requested=.false.)
    if (associated(child)) call generator_from_table(gen,child)

    call get_value(table,'elements',child,requested=.false.)
    if (associated(child)) call elements_from_table(param,child)
#else
    iostat = 1
    errmsg = 'gfnff was built without toml-f; TOML parameter files unavailable'
    name = ''
    if (.false.) then   !> silences unused-argument warnings
      print *,filename,param%chi(1),gen%linthr,version
    end if
#endif
  end subroutine param_read_toml

#ifdef WITH_TOMLF

  subroutine verify_readable(filename,iostat,errmsg)
    !***********************************************************************
    !* Re-parse a freshly written file, so a malformed one is caught here.
    !***********************************************************************
    character(len=*),intent(in) :: filename
    integer,intent(inout) :: iostat
    character(len=:),allocatable,intent(inout) :: errmsg
    type(toml_table),allocatable :: check
    type(toml_error),allocatable :: perr
    integer :: unit,err
    open (newunit=unit,file=filename,action='read',status='old',iostat=err)
    if (err /= 0) then
      iostat = err
      errmsg = 'wrote '//filename//' but could not reopen it'
      return
    end if
    !$omp critical(gfnff_tomlf)
    call toml_parse(check,unit,perr)
    !$omp end critical(gfnff_tomlf)
    close (unit)
    if (allocated(perr)) then
      iostat = 4
      errmsg = 'wrote '//filename//' but it does not parse back: '//perr%message
      return
    end if
    if (.not.allocated(check)) then
      iostat = 4
      errmsg = 'wrote '//filename//' but it parsed to nothing'
    end if
  end subroutine verify_readable


  subroutine put_array(table,key,vals)
    !***********************************************************************
    !* Write a real array as a TOML array of the given key.
    !***********************************************************************
    type(toml_table),intent(inout) :: table
    character(len=*),intent(in) :: key
    real(wp),intent(in) :: vals(:)
    type(toml_array),pointer :: arr
    integer :: i
    call add_array(table,key,arr)
    do i = 1,size(vals)
      call set_value(arr,i,vals(i))
    end do
  end subroutine put_array

  subroutine take_array(table,key,vals)
    !***********************************************************************
    !* Read a real array back, left untouched on missing key or length mismatch.
    !***********************************************************************
    type(toml_table),intent(inout) :: table
    character(len=*),intent(in) :: key
    real(wp),intent(inout) :: vals(:)
    type(toml_array),pointer :: arr
    integer :: i,st
    real(wp) :: v
    call get_value(table,key,arr,requested=.false.)
    if (.not.associated(arr)) return
    if (len(arr) /= size(vals)) return
    do i = 1,size(vals)
      call get_value(arr,i,v,stat=st)
      if (st == toml_stat%success) vals(i) = v
    end do
  end subroutine take_array

  subroutine elements_to_table(param,table)
    !***********************************************************************
    !* Write the per-element parameter arrays (elements 1..nel) to table.
    !***********************************************************************
    type(TGFFData),intent(in) :: param
    type(toml_table),intent(inout) :: table
    call put_array(table,'chi',param%chi(1:nel))
    call put_array(table,'gam',param%gam(1:nel))
    call put_array(table,'cnf',param%cnf(1:nel))
    call put_array(table,'alp',param%alp(1:nel))
    call put_array(table,'bond',param%bond(1:nel))
    call put_array(table,'repa',param%repa(1:nel))
    call put_array(table,'repan',param%repan(1:nel))
    call put_array(table,'angl',param%angl(1:nel))
    call put_array(table,'angl2',param%angl2(1:nel))
    call put_array(table,'tors',param%tors(1:nel))
    call put_array(table,'tors2',param%tors2(1:nel))
  end subroutine elements_to_table

  subroutine elements_from_table(param,table)
    !***********************************************************************
    !* Read the per-element parameter arrays (elements 1..nel) from table.
    !***********************************************************************
    type(TGFFData),intent(inout) :: param
    type(toml_table),intent(inout) :: table
    call take_array(table,'chi',param%chi(1:nel))
    call take_array(table,'gam',param%gam(1:nel))
    call take_array(table,'cnf',param%cnf(1:nel))
    call take_array(table,'alp',param%alp(1:nel))
    call take_array(table,'bond',param%bond(1:nel))
    call take_array(table,'repa',param%repa(1:nel))
    call take_array(table,'repan',param%repan(1:nel))
    call take_array(table,'angl',param%angl(1:nel))
    call take_array(table,'angl2',param%angl2(1:nel))
    call take_array(table,'tors',param%tors(1:nel))
    call take_array(table,'tors2',param%tors2(1:nel))
  end subroutine elements_from_table

  subroutine generator_to_table(gen,table)
    !***********************************************************************
    !* Write every TGFFGenerator scalar and array to table, one key each.
    !***********************************************************************
    type(TGFFGenerator),intent(in) :: gen
    type(toml_table),intent(inout) :: table
    real(wp) :: flat(16)
    call set_value(table,'linthr',gen%linthr)
    call set_value(table,'fcthr',gen%fcthr)
    call set_value(table,'tdist_thr',real(gen%tdist_thr,wp))
    call set_value(table,'rthr',gen%rthr)
    call set_value(table,'rthr2',gen%rthr2)
    call set_value(table,'rqshrink',gen%rqshrink)
    call set_value(table,'hqabthr',gen%hqabthr)
    call set_value(table,'qabthr',gen%qabthr)
    call set_value(table,'srb1',gen%srb1)
    call set_value(table,'srb2',gen%srb2)
    call set_value(table,'srb3',gen%srb3)
    call set_value(table,'qrepscal',gen%qrepscal)
    call set_value(table,'nrepscal',gen%nrepscal)
    call set_value(table,'hhfac',gen%hhfac)
    call set_value(table,'hh13rep',gen%hh13rep)
    call set_value(table,'hh14rep',gen%hh14rep)
    call put_array(table,'bstren',gen%bstren)
    call set_value(table,'qfacBEN',gen%qfacBEN)
    call set_value(table,'qfacTOR',gen%qfacTOR)
    call set_value(table,'fr3',gen%fr3)
    call set_value(table,'fr4',gen%fr4)
    call set_value(table,'fr5',gen%fr5)
    call set_value(table,'fr6',gen%fr6)
    call put_array(table,'torsf',gen%torsf)
    call set_value(table,'fbs1',gen%fbs1)
    call set_value(table,'batmscal',gen%batmscal)
    call set_value(table,'mchishift',gen%mchishift)
    call set_value(table,'rabshift',gen%rabshift)
    call set_value(table,'rabshifth',gen%rabshifth)
    call set_value(table,'hyper_shift',gen%hyper_shift)
    call set_value(table,'hshift3',gen%hshift3)
    call set_value(table,'hshift4',gen%hshift4)
    call set_value(table,'hshift5',gen%hshift5)
    call set_value(table,'metal1_shift',gen%metal1_shift)
    call set_value(table,'metal2_shift',gen%metal2_shift)
    call set_value(table,'metal3_shift',gen%metal3_shift)
    call set_value(table,'eta_shift',gen%eta_shift)
    call put_array(table,'qfacbm',gen%qfacbm)
    call set_value(table,'qfacbm0',gen%qfacbm0)
    call set_value(table,'rfgoed1',gen%rfgoed1)
    call set_value(table,'htriple',gen%htriple)
    call set_value(table,'hueckelp2',gen%hueckelp2)
    call set_value(table,'hueckelp3',gen%hueckelp3)
    call put_array(table,'hdiag',gen%hdiag)
    call put_array(table,'hoffdiag',gen%hoffdiag)
    call set_value(table,'hiter',gen%hiter)
    call set_value(table,'hueckelp',gen%hueckelp)
    call set_value(table,'bzref',gen%bzref)
    call set_value(table,'bzref2',gen%bzref2)
    call set_value(table,'pilpf',gen%pilpf)
    call set_value(table,'maxhiter',gen%maxhiter)
    call set_value(table,'d3a1',gen%d3a1)
    call set_value(table,'d3a2',gen%d3a2)
    call set_value(table,'split0',gen%split0)
    call set_value(table,'split1',gen%split1)
    call set_value(table,'fringbo',gen%fringbo)
    call set_value(table,'aheavy3',gen%aheavy3)
    call set_value(table,'aheavy4',gen%aheavy4)
    call set_value(table,'cnmax',gen%cnmax)
    !>-- bsmat is 4x4, flattened in Fortran column order, restored the same way
    flat = reshape(gen%bsmat,[16])
    call put_array(table,'bsmat',flat)
  end subroutine generator_to_table

  subroutine generator_from_table(gen,table)
    !***********************************************************************
    !* Read every TGFFGenerator field from table; missing keys keep gen's value.
    !***********************************************************************
    type(TGFFGenerator),intent(inout) :: gen
    type(toml_table),intent(inout) :: table
    real(wp) :: flat(16),tdum,tmp
    tmp = gen%linthr
    call get_value(table,'linthr',gen%linthr,tmp)
    tmp = gen%fcthr
    call get_value(table,'fcthr',gen%fcthr,tmp)
    tdum = real(gen%tdist_thr,wp)
    call get_value(table,'tdist_thr',tdum,tdum)
    gen%tdist_thr = real(tdum,sp)
    tmp = gen%rthr
    call get_value(table,'rthr',gen%rthr,tmp)
    tmp = gen%rthr2
    call get_value(table,'rthr2',gen%rthr2,tmp)
    tmp = gen%rqshrink
    call get_value(table,'rqshrink',gen%rqshrink,tmp)
    tmp = gen%hqabthr
    call get_value(table,'hqabthr',gen%hqabthr,tmp)
    tmp = gen%qabthr
    call get_value(table,'qabthr',gen%qabthr,tmp)
    tmp = gen%srb1
    call get_value(table,'srb1',gen%srb1,tmp)
    tmp = gen%srb2
    call get_value(table,'srb2',gen%srb2,tmp)
    tmp = gen%srb3
    call get_value(table,'srb3',gen%srb3,tmp)
    tmp = gen%qrepscal
    call get_value(table,'qrepscal',gen%qrepscal,tmp)
    tmp = gen%nrepscal
    call get_value(table,'nrepscal',gen%nrepscal,tmp)
    tmp = gen%hhfac
    call get_value(table,'hhfac',gen%hhfac,tmp)
    tmp = gen%hh13rep
    call get_value(table,'hh13rep',gen%hh13rep,tmp)
    tmp = gen%hh14rep
    call get_value(table,'hh14rep',gen%hh14rep,tmp)
    call take_array(table,'bstren',gen%bstren)
    tmp = gen%qfacBEN
    call get_value(table,'qfacBEN',gen%qfacBEN,tmp)
    tmp = gen%qfacTOR
    call get_value(table,'qfacTOR',gen%qfacTOR,tmp)
    tmp = gen%fr3
    call get_value(table,'fr3',gen%fr3,tmp)
    tmp = gen%fr4
    call get_value(table,'fr4',gen%fr4,tmp)
    tmp = gen%fr5
    call get_value(table,'fr5',gen%fr5,tmp)
    tmp = gen%fr6
    call get_value(table,'fr6',gen%fr6,tmp)
    call take_array(table,'torsf',gen%torsf)
    tmp = gen%fbs1
    call get_value(table,'fbs1',gen%fbs1,tmp)
    tmp = gen%batmscal
    call get_value(table,'batmscal',gen%batmscal,tmp)
    tmp = gen%mchishift
    call get_value(table,'mchishift',gen%mchishift,tmp)
    tmp = gen%rabshift
    call get_value(table,'rabshift',gen%rabshift,tmp)
    tmp = gen%rabshifth
    call get_value(table,'rabshifth',gen%rabshifth,tmp)
    tmp = gen%hyper_shift
    call get_value(table,'hyper_shift',gen%hyper_shift,tmp)
    tmp = gen%hshift3
    call get_value(table,'hshift3',gen%hshift3,tmp)
    tmp = gen%hshift4
    call get_value(table,'hshift4',gen%hshift4,tmp)
    tmp = gen%hshift5
    call get_value(table,'hshift5',gen%hshift5,tmp)
    tmp = gen%metal1_shift
    call get_value(table,'metal1_shift',gen%metal1_shift,tmp)
    tmp = gen%metal2_shift
    call get_value(table,'metal2_shift',gen%metal2_shift,tmp)
    tmp = gen%metal3_shift
    call get_value(table,'metal3_shift',gen%metal3_shift,tmp)
    tmp = gen%eta_shift
    call get_value(table,'eta_shift',gen%eta_shift,tmp)
    call take_array(table,'qfacbm',gen%qfacbm)
    tmp = gen%qfacbm0
    call get_value(table,'qfacbm0',gen%qfacbm0,tmp)
    tmp = gen%rfgoed1
    call get_value(table,'rfgoed1',gen%rfgoed1,tmp)
    tmp = gen%htriple
    call get_value(table,'htriple',gen%htriple,tmp)
    tmp = gen%hueckelp2
    call get_value(table,'hueckelp2',gen%hueckelp2,tmp)
    tmp = gen%hueckelp3
    call get_value(table,'hueckelp3',gen%hueckelp3,tmp)
    call take_array(table,'hdiag',gen%hdiag)
    call take_array(table,'hoffdiag',gen%hoffdiag)
    tmp = gen%hiter
    call get_value(table,'hiter',gen%hiter,tmp)
    tmp = gen%hueckelp
    call get_value(table,'hueckelp',gen%hueckelp,tmp)
    tmp = gen%bzref
    call get_value(table,'bzref',gen%bzref,tmp)
    tmp = gen%bzref2
    call get_value(table,'bzref2',gen%bzref2,tmp)
    tmp = gen%pilpf
    call get_value(table,'pilpf',gen%pilpf,tmp)
    tmp = gen%maxhiter
    call get_value(table,'maxhiter',gen%maxhiter,tmp)
    tmp = gen%d3a1
    call get_value(table,'d3a1',gen%d3a1,tmp)
    tmp = gen%d3a2
    call get_value(table,'d3a2',gen%d3a2,tmp)
    tmp = gen%split0
    call get_value(table,'split0',gen%split0,tmp)
    tmp = gen%split1
    call get_value(table,'split1',gen%split1,tmp)
    tmp = gen%fringbo
    call get_value(table,'fringbo',gen%fringbo,tmp)
    tmp = gen%aheavy3
    call get_value(table,'aheavy3',gen%aheavy3,tmp)
    tmp = gen%aheavy4
    call get_value(table,'aheavy4',gen%aheavy4,tmp)
    tmp = gen%cnmax
    call get_value(table,'cnmax',gen%cnmax,tmp)
    flat = reshape(gen%bsmat,[16])
    call take_array(table,'bsmat',flat)
    gen%bsmat = reshape(flat,[4,4])
  end subroutine generator_from_table
#endif

end module gfnff_param_io
