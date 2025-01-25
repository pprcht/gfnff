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
! along with gfnff.  If not, see <https://www.gnu.org/licenses/>.

module xyzreader
  use iso_fortran_env,only:wp => real64
  implicit none
  private

!&<
  !> Element symbols
  character(len=2),parameter,public :: PSE(118) = [ &
   & 'H ',                                                                                'He', &
   & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
   & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
   & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
   & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
   & 'Cs','Ba','La',                                                                            &
   &                'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',      &
   &                'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
   & 'Fr','Ra','Ac',                                                                            &
   &                'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',      &
   &                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' ]
!&>

  public :: readxyz,printxyz
  public :: exists

  real(wp),parameter :: autoaa = 0.529177249_wp

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES BEGIN HERE
!========================================================================================!
!========================================================================================!

  function exists(filename) result(yesno)
!***************************************************
!* one-line wrapper to check for a files' existance
!***************************************************
    character(len=*),intent(in) :: filename
    logical :: yesno
    inquire (file=trim(filename),exist=yesno)
  end function exists

!========================================================================================!

  subroutine readxyz(filename,nat,at,xyz)
!*************************************************
!* A simple reader for xyz files
!* Takes the file name as input and outputs
!* the number of atoms (nat), the atom types (at)
!* and coordinates in BOHR!!!
!*************************************************
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(out) :: nat
    integer,allocatable,intent(out) :: at(:)
    real(wp),allocatable,intent(out) :: xyz(:,:)
    integer :: i,j,k,ich,io
    character(len=80) :: atmp
    character(len=2) :: symb
!>--- reset
    nat = 0
    if (.not.exists(filename)) return
!>--- open and read file
    open (newunit=ich,file=filename)
    read (ich,*,iostat=io) j
    if (io == 0) then
      nat = j
    else
      return
    end if
    allocate (at(nat),source=0)
    allocate (xyz(3,nat),source=0.0_wp)
    read (ich,'(a)') atmp
    do i = 1,nat
      read (ich,*,iostat=io) symb,xyz(1:3,i)
      if (io == 0) then
        at(i) = e2i(symb)
      else
        error stop 'error in readxyz()'
      end if
    end do
    close (ich)
    xyz = xyz/autoaa
  end subroutine readxyz

!========================================================================================!

  subroutine printxyz(chnl,nat,at,xyz)
    implicit none
    integer,intent(in)  :: chnl,nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer :: i
    write (chnl,'(2x,i0)') nat
    write(chnl,*)
    do i = 1,nat
      write (chnl,'(a2,3F20.13)') i2e(at(i),'nc'),xyz(1:3,i)*autoaa
    end do
  end subroutine printxyz

!========================================================================================!

  integer function e2i(cin)
!*********************************************************
!* e2i is used to map the element (as a string) to integer
!*********************************************************
    implicit none
    character(len=*),intent(in) :: cin
    character(len=:),allocatable :: c
    integer :: iout
    integer :: i,j,k,ich,io,Z
    logical :: ex
    c = trim(convertlable(cin))
    read (cin,*,iostat=io) j
    if (io == 0) Z = j
    if (any(PSE(:) .eq. c)) then
      do i = 1,118
        if (trim(PSE(i)) .eq. c) then
          iout = i
          exit
        end if
      end do
    else if (io == 0.and.Z <= 118) then
      iout = Z
    else !> special cases
      select case (trim(c))
      case ('D'); iout = 1
      case ('T'); iout = 1
      case default; iout = 0
      end select
    end if
    e2i = iout
  end function e2i

!========================================================================================!

  character(len=2) function i2e(iin,oformat)
!************************************************************
!* i2e is used to map the element (as a integer) to a string
!************************************************************
    implicit none
    integer,intent(in) :: iin
    character(len=:),allocatable :: c
    character(len=*),optional :: oformat
    if (iin <= 118) then
      c = uppercase(PSE(iin))
    else
      c = 'XX'
    end if
    i2e = trim(c)
    if (present(oformat)) then
      select case (oformat)
      case ('lc','lowercase')
        i2e = lowerCase(trim(c))
      case ('nc','nicecase')
        if (len_trim(c) .gt. 1) then
          c(2:2) = lowerCase(c(2:2))
          i2e = trim(c)
        end if
      case default
        continue
      end select
    end if
  end function i2e

!========================================================================================!

  function upperCase(s)
    implicit none
    character(len=*),intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: upperCase
    integer :: ic,i
    character(26),Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1,LEN_TRIM(s)
      ic = INDEX(low,s(i:i))
      if (ic > 0) sout(i:i) = high(ic:ic)
    end do
    call move_alloc(sout,upperCase)
  end function upperCase

  function lowerCase(s)
    implicit none
    character(len=*),intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: lowerCase
    integer :: ic,i
    character(26),Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1,LEN_TRIM(s)
      ic = INDEX(high,s(i:i))
      if (ic > 0) sout(i:i) = low(ic:ic)
    end do
    call move_alloc(sout,lowerCase)
  end function lowerCase

  function convertlable(s)
    implicit none
    character(len=*),intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: convertlable
    integer :: ic,i
    character(14),parameter :: lab = '0123456789*_+-'
    character(26),parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1,len_trim(s)
      ic = index(lab,s(i:i))
      if (ic > 0) sout(i:i) = ' '
      ic = index(low,s(i:i))
      if (ic > 0) sout(i:i) = high(ic:ic)
    end do
    sout = trim(adjustl(sout))
    if (len_trim(sout) .gt. 1) then
      sout(2:2) = lowerCase(sout(2:2))
    else
      sout = sout//' '
    end if
    call move_alloc(sout,convertlable)
  end function convertlable

!========================================================================================!
!========================================================================================!
end module xyzreader
