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
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert, Sebastian Spicher, Stefan Grimme
!> at https://github.com/grimme-lab/xtb
! ------------------------------------------------------------------------------

!> Free-format number extraction from a text line, used when reading
!> parameter files.
module gfnff_strings
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  public :: readl

contains  !> MODULE PROCEDURES START HERE

  subroutine readl(a1,x,n)
    !***********************************************************************
    !* Reads the numbers on line a1 into x, in order. n counts all numbers
    !* found and can exceed size(x); the surplus is not stored.
    !***********************************************************************
    character(len=*) ::  a1
    real(wp) :: x(:)
    integer,intent(out) :: n
    n = size(x,1) !will be overwritten in getfloats
    call getfloats(a1,x,n)
  end subroutine readl

  subroutine getfloats(line,floats,cf)
    !***********************************************************************
    !* Stores the blank- or tab-separated numeric fields of line in floats.
    !* cf: size of floats on input, number of numeric fields on output.
    !***********************************************************************
    implicit none
    real(wp),intent(inout) :: floats(*)
    character(len=*),intent(in) :: line
    integer,intent(inout) :: cf
    real(wp) :: num
    character(len=128) :: str,stmp
    character(len=80) strings(3)
    character(len=1) digit
    integer :: i,ty,cs,cfmax

    cfmax = cf !on input cf is the dimension of floats()
    stmp = ''
    cs = 0
    cf = 0
    strings(:) = ''
    do i = 1,len(trim(line))
      digit = line(i:i)
      if (digit .ne. ' '.and.digit .ne. char(9)) then  !char(9) is the tab
        stmp = trim(stmp)//trim(digit)
      elseif (stmp .ne. '') then
        call checktype(stmp,num,str,ty)      !ty: 0=number, 1=text
        if (ty .eq. 0) then
          cf = cf+1
          if (cf .le. cfmax) floats(cf) = num
        elseif (ty .eq. 1) then
          cs = cs+1
        else
          write (*,*) 'Problem in checktype, must abort'
          exit
        end if
        stmp = ''
      end if
      if (i .eq. len(trim(line))) then  !flush the last field at end of line
        call checktype(stmp,num,str,ty)
        if (ty .eq. 0) then
          cf = cf+1
          if (cf .le. cfmax) floats(cf) = num
        elseif (ty .eq. 1) then
          cs = cs+1
        else
          write (*,*) 'Problem in checktype, must abort'
          exit
        end if
        stmp = ''
      end if
    end do
  contains
    subroutine checktype(field,num,str,ty)
      !*********************************************************************
      !* Classifies field: ty = 0 with its value in num if it reads as a
      !* number, ty = 1 with the text in str otherwise.
      !*********************************************************************
      implicit none
      character(len=*) :: field,str
      real(wp) :: num
      integer :: e,ty
      logical :: is_num
      ty = 99
      str = ''
      is_num = .false.
      read (field,'(F10.5)',IOSTAT=e) num !e = 0 if field reads as a number
      if (e .eq. 0) is_num = .true.
      if (is_num) then
        if (index(field,'.') .ne. 0) then  !check for integer/real
          read (field,'(F30.16)') num
          ty = 0
        else                       !if integer, add .0 to string; otherwise cast to real does not work
          str = trim(field)//'.0'
          read (str,'(F30.16)') num
          str = ''
          ty = 0
        end if
      else
        str = trim(field)
        ty = 1
      end if
    end subroutine checktype
  end subroutine getfloats
end module gfnff_strings
