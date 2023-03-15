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
! along with gfnff.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!
module solvation_math_wrapper
!> module solvation_math_wrapper
!> contains some interfaces to LAPACK and BLAS routines
!> which must be included via a suitable library.
  use iso_fortran_env,only:wp => real64,sp => real32,stderr => error_unit
  implicit none
  public

  interface dot
    module procedure sdot_wrap
    module procedure ddot_wrap
  end interface dot

  interface symv
    module procedure ssymv_wrap
    module procedure dsymv_wrap
  end interface symv

  interface gemv
    module procedure sgemv_wrap
    module procedure dgemv_wrap
    module procedure dgemv312_wrap
  end interface gemv

contains

!========================================================================================!

!> Determinat of 3×3 matrix
  pure function matDet3x3(a) result(det)

    !> Matrix
    real(wp),intent(in) :: a(3,3)

    !> Determinant
    real(wp) :: det

    det = a(1,1)*a(2,2)*a(3,3)  &
       & -a(1,1)*a(2,3)*a(3,2)  &
       & -a(1,2)*a(2,1)*a(3,3)  &
       & +a(1,2)*a(2,3)*a(3,1)  &
       & +a(1,3)*a(2,1)*a(3,2)  &
       & -a(1,3)*a(2,2)*a(3,1)

  end function matDet3x3

!========================================================================================!
!> The DOT routines perform the dot product of two vectors
!> dot = ∑ᵢ xᵢyᵢ
  function sdot_wrap(xvec,yvec) result(dotp)
    real(sp) :: dotp
    real(sp),intent(in) :: xvec(:)
    real(sp),intent(in) :: yvec(:)
    integer :: incx,incy,n
    !> BLAS
    real(sp),external :: sdot
    incx = 1
    incy = 1
    n = size(xvec)
    dotp = sdot(n,xvec,incx,yvec,incy)
  end function sdot_wrap
  function ddot_wrap(xvec,yvec) result(dotp)
    real(wp) :: dotp
    real(wp),intent(in) :: xvec(:)
    real(wp),intent(in) :: yvec(:)
    integer :: incx,incy,n
    !> BLAS
    real(wp),external :: ddot
    incx = 1
    incy = 1
    n = size(xvec)
    dotp = ddot(n,xvec,incx,yvec,incy)
  end function ddot_wrap
!=======================================================================================!
!> SYMV  performs the matrix-vector  operation
!>
!>   y := α*A*x + β*y,
!>
!> where α and β are scalars, x and y are n element vectors and
!> A is an n by n symmetric matrix.
  subroutine ssymv_wrap(amat,xvec,yvec,uplo,alpha,beta)
    real(sp),intent(in) :: amat(:,:)
    real(sp),intent(in) :: xvec(:)
    real(sp),intent(inout) :: yvec(:)
    character(len=1),intent(in),optional :: uplo
    real(sp),intent(in),optional :: alpha
    real(sp),intent(in),optional :: beta
    character(len=1) :: ula
    real(sp) :: a,b
    integer :: incx,incy,n,lda
    !> BLAS
    external :: ssymv
    if (present(alpha)) then
      a = alpha
    else
      a = 1.0_sp
    end if
    if (present(beta)) then
      b = beta
    else
      b = 0
    end if
    if (present(uplo)) then
      ula = uplo
    else
      ula = 'u'
    end if
    incx = 1
    incy = 1
    lda = max(1,size(amat,1))
    n = size(amat,2)
    call ssymv(ula,n,a,amat,lda,xvec,incx,b,yvec,incy)
  end subroutine ssymv_wrap
  subroutine dsymv_wrap(amat,xvec,yvec,uplo,alpha,beta)
    real(wp),intent(in) :: amat(:,:)
    real(wp),intent(in) :: xvec(:)
    real(wp),intent(inout) :: yvec(:)
    character(len=1),intent(in),optional :: uplo
    real(wp),intent(in),optional :: alpha
    real(wp),intent(in),optional :: beta
    character(len=1) :: ula
    real(wp) :: a,b
    integer :: incx,incy,n,lda
    !> BLAS
    external :: dsymv
    if (present(alpha)) then
      a = alpha
    else
      a = 1.0_wp
    end if
    if (present(beta)) then
      b = beta
    else
      b = 0
    end if
    if (present(uplo)) then
      ula = uplo
    else
      ula = 'u'
    end if
    incx = 1
    incy = 1
    lda = max(1,size(amat,1))
    n = size(amat,2)
    call dsymv(ula,n,a,amat,lda,xvec,incx,b,yvec,incy)
  end subroutine dsymv_wrap

!=======================================================================================!
!> DGEMV  performs one of the matrix-vector operations
!>
!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n matrix.
  subroutine sgemv_wrap(amat,xvec,yvec,alpha,beta,trans)
    real(sp),intent(in) :: amat(:,:)
    real(sp),intent(in) :: xvec(:)
    real(sp),intent(inout) :: yvec(:)
    real(sp),intent(in),optional :: alpha
    real(sp),intent(in),optional :: beta
    character(len=1),intent(in),optional :: trans
    real(sp) :: a,b
    character(len=1) :: tra
    integer :: incx,incy,m,n,lda
    !> BLAS
    external :: sgemv
    if (present(alpha)) then
      a = alpha
    else
      a = 1.0_sp
    end if
    if (present(beta)) then
      b = beta
    else
      b = 0
    end if
    if (present(trans)) then
      tra = trans
    else
      tra = 'n'
    end if
    incx = 1
    incy = 1
    lda = max(1,size(amat,1))
    m = size(amat,1)
    n = size(amat,2)
    call sgemv(tra,m,n,a,amat,lda,xvec,incx,b,yvec,incy)
  end subroutine sgemv_wrap
  subroutine dgemv_wrap(amat,xvec,yvec,alpha,beta,trans)
    real(wp),intent(in) :: amat(:,:)
    real(wp),intent(in) :: xvec(:)
    real(wp),intent(inout) :: yvec(:)
    real(wp),intent(in),optional :: alpha
    real(wp),intent(in),optional :: beta
    character(len=1),intent(in),optional :: trans
    real(wp) :: a,b
    character(len=1) :: tra
    integer :: incx,incy,m,n,lda
    !> BLAS
    external :: dgemv
    if (present(alpha)) then
      a = alpha
    else
      a = 1.0_wp
    end if
    if (present(beta)) then
      b = beta
    else
      b = 0
    end if
    if (present(trans)) then
      tra = trans
    else
      tra = 'n'
    end if
    incx = 1
    incy = 1
    lda = max(1,size(amat,1))
    m = size(amat,1)
    n = size(amat,2)
    call dgemv(tra,m,n,a,amat,lda,xvec,incx,b,yvec,incy)
  end subroutine dgemv_wrap
  subroutine dgemv312_wrap(amat,xvec,yvec,alpha,beta,trans)
    real(wp),intent(in),contiguous,target :: amat(:,:,:)
    real(wp),intent(in) :: xvec(:)
    real(wp),intent(inout),contiguous,target :: yvec(:,:)
    real(wp),intent(in),optional :: alpha
    real(wp),intent(in),optional :: beta
    character(len=1),intent(in),optional :: trans
    real(wp),pointer :: aptr(:,:),yptr(:)
    character(len=1) :: tra
    if (present(trans)) then
      tra = trans
    else
      tra = 'n'
    end if
    if (any(tra == ['n','N'])) then
      aptr(1:size(amat,1)*size(amat,2),1:size(amat,3)) => amat
      yptr(1:size(yvec,1)*size(yvec,2)) => yvec
    else
      aptr(1:size(amat,1),1:size(amat,2)*size(amat,3)) => amat
      yptr(1:size(yvec,1)*size(yvec,2)) => yvec
    end if
    call dgemv_wrap(aptr,xvec,yptr,alpha,beta,tra)
  end subroutine dgemv312_wrap
!=======================================================================================!
end module solvation_math_wrapper
