!================================================================================!
! This file is part of gfnff.
!
! Copyright (C) 2025 Philipp Pracht
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

!> c bindings for gfnff

module gfnff_interface_c
  use iso_c_binding
  use iso_fortran_env,only:wp => real64,stdout => output_unit,stderr => error_unit
  use gfnff_interface
  implicit none
  private

  !> Public C-compatible interface
  public :: c_gfnff_calculator
  public :: c_gfnff_calculator_init
  public :: c_gfnff_calculator_deallocate
  public :: c_gfnff_calculator_singlepoint
  public :: c_gfnff_calculator_results

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
    logical :: verbose,pr
    integer :: iunit,printlevel,iostatus
    integer :: i,ichrg

    !> Convert C arguments to Fortran types
    nat = c_nat
    call c_f_pointer(c_loc(c_at),at, [nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,nat]) !> assumes xyz[nat][3] in C
    ichrg = c_ichrg
    pr = .false.
    verbose = .false.
    iunit = stdout !> use a default here.
    printlevel = c_printlevel
    if (printlevel > 0) pr = .true.
    if (printlevel > 1) verbose = .true.
    solvent = c_string_to_fortran(c_solvent)
    if (len_trim(solvent) == 0) deallocate (solvent)

    !> Allocate and initialize the Fortran calculator
    allocate (calc)
    call calc%init(nat,at,xyz,ichrg=ichrg, &
    &              print=pr,verbose=verbose,iunit=iunit,iostat=iostatus,&
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
    &                                           c_energy,c_gradient,c_iostat) &
    &                        bind(C,name="c_gfnff_calculator_singlepoint")
    implicit none
    !> Input arguments from C
    type(c_gfnff_calculator),intent(inout) :: c_calculator
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),target,intent(in) :: c_at(*)
    real(c_double),target,intent(in) :: c_xyz(3,*) !> NOTE Fortran/C matrix orders

    !> Output arguments to C
    real(c_double),intent(out) :: c_energy
    real(c_double),target,intent(out) :: c_gradient(3,*) !> NOTE Fortran/C matrix orders
    integer(c_int),intent(out) :: c_iostat

    !> Local Fortran variables
    type(gfnff_data),pointer :: calc_ptr
    integer :: nat
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    real(wp),pointer :: grad(:,:)
    real(wp) :: energy
    integer :: iostat,i,j

    !> Convert C pointers to Fortran pointers
    call c_f_pointer(c_calculator%ptr,calc_ptr)
    call c_f_pointer(c_loc(c_at),at, [c_nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,c_nat]) !> Assumes xyz[nat][3] in C
    call c_f_pointer(c_loc(c_gradient),grad, [3,c_nat])  !> Assumes grad[nat][3] in C

    !> Set the integer variable
    nat = c_nat

    !> Call the Fortran subroutine
    call calc_ptr%singlepoint(nat,at,xyz,energy,grad,iostat=iostat)

    !> Pass back the results to C variables
    c_energy = energy
    c_gradient(1:3,1:nat) = grad(1:3,1:nat)
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
!========================================================================================!
end module gfnff_interface_c
