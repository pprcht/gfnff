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

!> Bisection search over a sorted array, and an index heap sort that leaves
!> the array untouched and returns the sorting permutation instead.
module gfnff_search
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  public :: indexHeapSort
  public :: bisectSearch
  interface bisectSearch
    module procedure :: bisectSearchReal
    module procedure :: bisectSearchInteger
  end interface bisectSearch

contains  !> MODULE PROCEDURES START HERE

  pure subroutine bisectSearchInteger(j,xx,x)
    !***********************************************************************
    !* Bisection search: on return xx(j) <= x < xx(j+1), xx monotonic
    !* ascending or descending. j = 0 or size(xx) if x is out of range.
    !***********************************************************************
    integer,intent(out) :: j
    integer,intent(in) :: xx(:)
    integer,intent(in) :: x

    integer :: n
    integer :: jlower,jupper,jcurr

    n = size(xx)
    if (n == 0) then
      j = 0
      return
    end if

    if (x < xx(1)) then
      j = 0
    else if (x == xx(1)) then
      j = 1
    else if (x == xx(n)) then
      j = n-1
    else if (x > xx(n)) then
      j = n
    else
      jlower = 0
      jcurr = n+1
      do while ((jcurr-jlower) > 1)
        jupper = (jcurr+jlower)/2
        if ((xx(n) >= xx(1)).eqv.(x >= xx(jupper))) then
          jlower = jupper
        else
          jcurr = jupper
        end if
      end do
      j = jlower
    end if

  end subroutine bisectSearchInteger

  pure subroutine bisectSearchReal(j,xx,x,tol)
    !***********************************************************************
    !* Real-valued variant of bisectSearchInteger; tol is the equality tolerance.
    !***********************************************************************
    integer,intent(out) :: j
    real(wp),intent(in) :: xx(:)
    real(wp),intent(in) :: x
    real(wp),intent(in),optional :: tol

    integer :: n
    integer :: jlower,jupper,jcurr
    real(wp) :: rTol
    logical :: ascending

    n = size(xx)
    if (n == 0) then
      j = 0
      return
    end if

    if (present(tol)) then
      rTol = tol
    else
      rTol = epsilon(0.0_wp)
    end if

    if (x < xx(1)-rTol) then
      j = 0
    else if (abs(x-xx(1)) <= rTol) then
      j = 1
    else if (abs(x-xx(n)) <= rTol) then
      j = n-1
    else if (x > xx(n)+rTol) then
      j = n
    else
      ascending = (xx(n) >= xx(1))
      jlower = 0
      jcurr = n+1
      do while ((jcurr-jlower) > 1)
        jupper = (jcurr+jlower)/2
        if (ascending.eqv.(x >= xx(jupper)+rTol)) then
          jlower = jupper
        else
          jcurr = jupper
        end if
      end do
      j = jlower
    end if

  end subroutine bisectSearchReal

  pure subroutine indexHeapSort(indx,array,tolerance)
    !***********************************************************************
    !* Heap sort of array by index: indx is the permutation such that
    !* array(indx) is ascending; array itself is untouched, and indx must
    !* have the same size. tolerance sets the equality tolerance (default
    !* machine epsilon). Based on Numerical Recipes Software 1986-92.
    !***********************************************************************
    integer,intent(out) :: indx(:)
    real(wp),intent(in) :: array(:)
    real(wp),intent(in),optional :: tolerance

    integer :: n,ir,ij,il,ii
    integer :: indxTmp
    real(wp) :: arrayTmp,tol

    if (present(tolerance)) then
      tol = tolerance
    else
      tol = epsilon(0.0_wp)
    end if

    do ii = 1,size(indx)
      indx(ii) = ii
    end do
    n = size(array)
    if (n <= 1) return
    il = n/2+1
    ir = n
    do
      if (il > 1) then
        il = il-1
        indxTmp = indx(il)
        arrayTmp = array(indxTmp)
      else
        indxTmp = indx(ir)
        arrayTmp = array(indxTmp)
        indx(ir) = indx(1)
        ir = ir-1
        if (ir < 1) then
          indx(1) = indxTmp
          return
        end if
      end if
      ii = il
      ij = 2*il
      do while (ij <= ir)
        if (ij < ir) then
          if (array(indx(ij)) < array(indx(ij+1))-tol) then
            ij = ij+1
          end if
        end if
        if (arrayTmp < array(indx(ij))-tol) then
          indx(ii) = indx(ij)
          ii = ij
          ij = 2*ij
        else
          ij = ir+1
        end if
      end do
      indx(ii) = indxTmp
    end do

  end subroutine indexHeapSort
end module gfnff_search
