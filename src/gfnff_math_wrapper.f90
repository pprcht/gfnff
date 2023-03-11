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
module gfnff_math_wrapper
!> module gfnff_math_wrapper
!> contains some interfaces to LAPACK and BLAS routines
!> which must be included via a suitable library.
  use iso_fortran_env,only:wp => real64,sp => real32,stderr => error_unit
  implicit none
  public

  interface contract
    module procedure contract312
    module procedure contract323
  end interface contract

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

  interface gemm
    module procedure sgemm_wrap
    module procedure dgemm_wrap
  end interface gemm

  interface sytrf_wrap
    module procedure ssytrf_wrap
    module procedure dsytrf_wrap
  end interface sytrf_wrap
  interface lapack_sytrf
    pure subroutine ssytrf(uplo,n,a,lda,ipiv,work,lwork,info)
      import :: sp
      real(sp),intent(inout) :: a(lda,*)
      character(len=1),intent(in) :: uplo
      integer,intent(out) :: ipiv(*)
      integer,intent(out) :: info
      integer,intent(in) :: n
      integer,intent(in) :: lda
      real(sp),intent(inout) :: work(*)
      integer,intent(in) :: lwork
    end subroutine ssytrf
    pure subroutine dsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
      import :: wp
      real(wp),intent(inout) :: a(lda,*)
      character(len=1),intent(in) :: uplo
      integer,intent(out) :: ipiv(*)
      integer,intent(out) :: info
      integer,intent(in) :: n
      integer,intent(in) :: lda
      real(wp),intent(inout) :: work(*)
      integer,intent(in) :: lwork
    end subroutine dsytrf
  end interface lapack_sytrf

  interface lapack_sytri
    module procedure ssytri_wrap
    module procedure dsytri_wrap
  end interface lapack_sytri

  interface sytrs_wrap
    module procedure ssytrs_wrap
    module procedure ssytrs1_wrap
    module procedure ssytrs3_wrap
    module procedure dsytrs_wrap
    module procedure dsytrs1_wrap
    module procedure dsytrs3_wrap
  end interface sytrs_wrap
  interface lapack_sytrs
    pure subroutine ssytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
      import :: sp
      real(sp),intent(in) :: a(lda,*)
      real(sp),intent(inout) :: b(ldb,*)
      integer,intent(in) :: ipiv(*)
      character(len=1),intent(in) :: uplo
      integer,intent(out) :: info
      integer,intent(in) :: n
      integer,intent(in) :: nrhs
      integer,intent(in) :: lda
      integer,intent(in) :: ldb
    end subroutine ssytrs
    pure subroutine dsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
      import :: wp
      real(wp),intent(in) :: a(lda,*)
      real(wp),intent(inout) :: b(ldb,*)
      integer,intent(in) :: ipiv(*)
      character(len=1),intent(in) :: uplo
      integer,intent(out) :: info
      integer,intent(in) :: n
      integer,intent(in) :: nrhs
      integer,intent(in) :: lda
      integer,intent(in) :: ldb
    end subroutine dsytrs
  end interface lapack_sytrs

  interface lapack_sygvd
    module procedure ssygvd_wrap
    module procedure dsygvd_wrap
  end interface lapack_sygvd

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
  subroutine contract312(amat,bvec,cvec,alpha,beta)
    !> note: this routine does the same as the gemv312 routines
    real(wp),intent(in),contiguous,target :: amat(:,:,:)
    real(wp),intent(in),contiguous :: bvec(:)
    real(wp),intent(inout),contiguous,target :: cvec(:,:)
    real(wp),intent(in),optional :: alpha
    real(wp),intent(in),optional :: beta
    real(wp),pointer :: aptr(:,:)
    real(wp),pointer :: cptr(:)

    integer :: nn,mm
    real(wp) :: aa,bb
    !> BLAS
    external :: dgemv

    aa = 1.0_wp
    if (present(alpha)) aa = alpha
    bb = 0.0_wp
    if (present(beta)) bb = beta

    nn = size(amat,dim=1)*size(amat,dim=2)
    mm = size(amat,dim=3)
    aptr(1:size(amat,dim=1)*size(amat,dim=2),1:size(amat,dim=3)) => amat
    cptr(1:size(cvec,dim=1)*size(cvec,dim=2)) => cvec

    call dgemv('n',nn,mm,aa,aptr,nn,bvec,1,bb,cptr,1)

  end subroutine contract312
!========================================================================================!
  subroutine contract323(amat,bmat,cmat,alpha,beta)

    real(wp),intent(in),contiguous,target :: amat(:,:,:)
    real(wp),intent(in),contiguous :: bmat(:,:)
    real(wp),intent(inout),contiguous,target :: cmat(:,:,:)
    real(wp),intent(in),optional :: alpha
    real(wp),intent(in),optional :: beta
    real(wp),pointer :: aptr(:,:)
    real(wp),pointer :: cptr(:,:)

    integer :: nn,mm,kk
    real(wp) :: aa,bb
    !> BLAS
    external :: dgemm

    aa = 1.0_wp
    if (present(alpha)) aa = alpha
    bb = 0.0_wp
    if (present(beta)) bb = beta

    nn = size(amat,dim=1)*size(amat,dim=2)
    mm = size(amat,dim=3)
    kk = size(cmat,dim=3)
    aptr(1:size(amat,dim=1)*size(amat,dim=2),1:size(amat,dim=3)) => amat
    cptr(1:size(cmat,dim=1)*size(cmat,dim=2),1:size(cmat,dim=3)) => cmat

    call dgemm('n','n',nn,kk,mm,aa,aptr,nn,bmat,mm,bb,cptr,nn)

  end subroutine contract323
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
!>  SYTRF computes the factorization of a real symmetric matrix A using
!>  the Bunch-Kaufman diagonal pivoting method.  The form of the
!>  factorization is
!>
!>     A = U**T*D*U  or  A = L*D*L**T
!>
!>  where U (or L) is a product of permutation and unit upper (lower)
!>  triangular matrices, and D is symmetric and block diagonal with
!>  1-by-1 and 2-by-2 diagonal blocks.

!  subroutine ssytrf_wrap(uplo,n,a,lda,ipiv,work,lwork,info)
!    real(sp),intent(inout) :: a(lda,*)
!    character(len=1),intent(in) :: uplo
!    integer,intent(out) :: ipiv(*)
!    integer,intent(out) :: info
!    integer,intent(in) :: n
!    integer,intent(in) :: lda
!    real(sp),intent(inout) :: work(*)
!    integer,intent(in) :: lwork
!    !> LAPACK
!    external :: ssytrf
!    call ssytrf(uplo,n,a,lda,ipiv,work,lwork,info)
!  end subroutine ssytrf_wrap
  subroutine ssytrf_wrap(amat,ipiv,info,uplo)
    character(len=*),parameter :: source = 'lapack_sytrf'
    real(sp),intent(inout) :: amat(:,:)
    integer,intent(out) :: ipiv(:)
    integer,intent(out) :: info
    character(len=1),intent(in),optional :: uplo
    character(len=1) :: ula
    integer :: n,lda,lwork,stat_alloc,stat_dealloc
    real(sp),allocatable :: work(:)
    real(sp) :: test(1)
    if (present(uplo)) then
      ula = uplo
    else
      ula = 'u'
    end if
    lda = max(1,size(amat,1))
    n = size(amat,2)
    stat_alloc = 0
    lwork = -1
    call lapack_sytrf(ula,n,amat,lda,ipiv,test,lwork,info)
    if (info == 0) then
      lwork = nint(test(1))
      if (stat_alloc == 0) then
        allocate (work(lwork),stat=stat_alloc)
      end if
      if (stat_alloc == 0) then
        call lapack_sytrf(ula,n,amat,lda,ipiv,work,lwork,info)
      else
        info = -1000
      end if
      deallocate (work,stat=stat_dealloc)
    end if
    if (info /= 0) then
      write (stderr,'("Factorisation of matrix failed ",a)') source
    end if
  end subroutine ssytrf_wrap

!  subroutine dsytrf_wrap(uplo,n,a,lda,ipiv,work,lwork,info)
!    real(wp),intent(inout) :: a(lda,*)
!    character(len=1),intent(in) :: uplo
!    integer,intent(out) :: ipiv(*)
!    integer,intent(out) :: info
!    integer,intent(in) :: n
!    integer,intent(in) :: lda
!    real(wp),intent(inout) :: work(*)
!    integer,intent(in) :: lwork
!    !> LAPACK
!    external :: dsytrf
!    call dsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
!  end subroutine dsytrf_wrap
  subroutine dsytrf_wrap(amat,ipiv,info,uplo)
    character(len=*),parameter :: source = 'lapack_sytrf'
    real(wp),intent(inout) :: amat(:,:)
    integer,intent(out) :: ipiv(:)
    integer,intent(out) :: info
    character(len=1),intent(in),optional :: uplo
    character(len=1) :: ula
    integer :: n,lda,lwork,stat_alloc,stat_dealloc
    real(wp),allocatable :: work(:)
    real(wp) :: test(1)
    if (present(uplo)) then
      ula = uplo
    else
      ula = 'u'
    end if
    lda = max(1,size(amat,1))
    n = size(amat,2)
    stat_alloc = 0
    lwork = -1
    call lapack_sytrf(ula,n,amat,lda,ipiv,test,lwork,info)
    if (info == 0) then
      lwork = nint(test(1))
      if (stat_alloc == 0) then
        allocate (work(lwork),stat=stat_alloc)
      end if
      if (stat_alloc == 0) then
        call lapack_sytrf(ula,n,amat,lda,ipiv,work,lwork,info)
      else
        info = -1000
      end if
      deallocate (work,stat=stat_dealloc)
    end if
    if (info /= 0) then
      write (stderr,'("Factorisation of matrix failed ",a)') source
    end if
  end subroutine dsytrf_wrap

!=======================================================================================!
!> SYTRI computes the inverse of a real symmetric indefinite matrix
!> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
!> SYTRF.
  subroutine ssytri_wrap(uplo,n,a,lda,ipiv,work,info)
    real(sp),intent(inout) :: a(lda,*)
    integer,intent(in) :: ipiv(*)
    character(len=1),intent(in) :: uplo
    integer,intent(out) :: info
    integer,intent(in) :: n
    integer,intent(in) :: lda
    real(sp),intent(in) :: work(*)
    !> LAPACK
    external :: ssytri
    call ssytri(uplo,n,a,lda,ipiv,work,info)
  end subroutine ssytri_wrap
  subroutine dsytri_wrap(uplo,n,a,lda,ipiv,work,info)
    real(wp),intent(inout) :: a(lda,*)
    integer,intent(in) :: ipiv(*)
    character(len=1),intent(in) :: uplo
    integer,intent(out) :: info
    integer,intent(in) :: n
    integer,intent(in) :: lda
    real(wp),intent(in) :: work(*)
    !> LAPACK
    external :: dsytri
    call dsytri(uplo,n,a,lda,ipiv,work,info)
  end subroutine dsytri_wrap

!=======================================================================================!
!> DSYTRS solves a system of linear equations A*X = B with a real
!> symmetric matrix A using the factorization A = U*D*U**T or
!>  A = L*D*L**T computed by DSYTRF.
  subroutine ssytrs_wrap(amat,bmat,ipiv,info,uplo)
    character(len=*),parameter :: source = 'lapack_sytrs'
    real(sp),intent(in) :: amat(:,:)
    real(sp),intent(inout) :: bmat(:,:)
    integer,intent(in) :: ipiv(:)
    integer,intent(out) :: info
    character(len=1),intent(in),optional :: uplo
    character(len=1) :: ula
    integer :: n,nrhs,lda,ldb
    if (present(uplo)) then
      ula = uplo
    else
      ula = 'u'
    end if
    lda = max(1,size(amat,1))
    ldb = max(1,size(bmat,1))
    n = size(amat,2)
    nrhs = size(bmat,2)
    call lapack_sytrs(ula,n,nrhs,amat,lda,ipiv,bmat,ldb,info)
    if (info /= 0) then
      write (stderr,'("Solving linear system failed ",a)') source
    end if
  end subroutine ssytrs_wrap
  subroutine ssytrs1_wrap(amat,bvec,ipiv,info,uplo)
    real(sp),intent(in) :: amat(:,:)
    real(sp),intent(inout),target :: bvec(:)
    integer,intent(in) :: ipiv(:)
    integer,intent(out) :: info
    character(len=1),intent(in),optional :: uplo
    real(sp),pointer :: bptr(:,:)
    bptr(1:size(bvec),1:1) => bvec
    call ssytrs_wrap(amat,bptr,ipiv,info,uplo)
  end subroutine ssytrs1_wrap
  subroutine ssytrs3_wrap(amat,bmat,ipiv,info,uplo)
    real(sp),intent(in) :: amat(:,:)
    real(sp),intent(inout),contiguous,target :: bmat(:,:,:)
    integer,intent(in) :: ipiv(:)
    integer,intent(out) :: info
    character(len=1),intent(in),optional :: uplo
    real(sp),pointer :: bptr(:,:)
    bptr(1:size(bmat,1),1:size(bmat,2)*size(bmat,3)) => bmat
    call ssytrs_wrap(amat,bptr,ipiv,info,uplo)
  end subroutine ssytrs3_wrap

  subroutine dsytrs_wrap(amat,bmat,ipiv,info,uplo)
    character(len=*),parameter :: source = 'lapack_sytrs'
    real(wp),intent(in) :: amat(:,:)
    real(wp),intent(inout) :: bmat(:,:)
    integer,intent(in) :: ipiv(:)
    integer,intent(out) :: info
    character(len=1),intent(in),optional :: uplo
    character(len=1) :: ula
    integer :: n,nrhs,lda,ldb
    if (present(uplo)) then
      ula = uplo
    else
      ula = 'u'
    end if
    lda = max(1,size(amat,1))
    ldb = max(1,size(bmat,1))
    n = size(amat,2)
    nrhs = size(bmat,2)
    call lapack_sytrs(ula,n,nrhs,amat,lda,ipiv,bmat,ldb,info)
    if (info /= 0) then
      write (stderr,'("Solving linear system failed ",a)') source
    end if
  end subroutine dsytrs_wrap
  subroutine dsytrs1_wrap(amat,bvec,ipiv,info,uplo)
    real(wp),intent(in) :: amat(:,:)
    real(wp),intent(inout),target :: bvec(:)
    integer,intent(in) :: ipiv(:)
    integer,intent(out) :: info
    character(len=1),intent(in),optional :: uplo
    real(wp),pointer :: bptr(:,:)
    bptr(1:size(bvec),1:1) => bvec
    call dsytrs_wrap(amat,bptr,ipiv,info,uplo)
  end subroutine dsytrs1_wrap
  subroutine dsytrs3_wrap(amat,bmat,ipiv,info,uplo)
    real(wp),intent(in) :: amat(:,:)
    real(wp),intent(inout),contiguous,target :: bmat(:,:,:)
    integer,intent(in) :: ipiv(:)
    integer,intent(out) :: info
    character(len=1),intent(in),optional :: uplo
    real(wp),pointer :: bptr(:,:)
    bptr(1:size(bmat,1),1:size(bmat,2)*size(bmat,3)) => bmat
    call dsytrs_wrap(amat,bptr,ipiv,info,uplo)
  end subroutine dsytrs3_wrap

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
!     DGEMM  performs one of the matrix-matrix operations
!
!        C := alpha*op( A )*op( B ) + beta*C,
!
!     where  op( X ) is one of
!
!        op( X ) = X   or   op( X ) = X**T,
!
!     alpha and beta are scalars, and A, B and C are matrices, with op( A )
!     an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  pure subroutine dgemm_wrap(amat,bmat,cmat,transa,transb,alpha,beta)
    real(wp),intent(in) :: amat(:,:)
    real(wp),intent(in) :: bmat(:,:)
    real(wp),intent(inout) :: cmat(:,:)
    character(len=1),intent(in),optional :: transa
    character(len=1),intent(in),optional :: transb
    real(wp),intent(in),optional :: alpha
    real(wp),intent(in),optional :: beta
    character(len=1) :: tra,trb
    real(wp) :: a,b
    integer :: m,n,k,lda,ldb,ldc
    !> BLAS
    interface
      pure subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb, &
            & beta,c,ldc)
        import :: wp
        real(wp),intent(in) :: a(lda,*)
        real(wp),intent(in) :: b(ldb,*)
        real(wp),intent(inout) :: c(ldc,*)
        character(len=1),intent(in) :: transa
        character(len=1),intent(in) :: transb
        real(wp),intent(in) :: alpha
        real(wp),intent(in) :: beta
        integer,intent(in) :: m
        integer,intent(in) :: n
        integer,intent(in) :: k
        integer,intent(in) :: lda
        integer,intent(in) :: ldb
        integer,intent(in) :: ldc
      end subroutine dgemm
    end interface
    if (present(alpha)) then
      a = alpha
    else
      a = 1.0_wp
    end if
    if (present(beta)) then
      b = beta
    else
      b = 0.0_wp
    end if
    if (present(transa)) then
      tra = transa
    else
      tra = 'n'
    end if
    if (present(transb)) then
      trb = transb
    else
      trb = 'n'
    end if
    if ((tra .eq. 'n'.or.tra .eq. 'N')) then
      k = size(amat,2)
    else
      k = size(amat,1)
    end if
    lda = max(1,size(amat,1))
    ldb = max(1,size(bmat,1))
    ldc = max(1,size(cmat,1))
    m = size(cmat,1)
    n = size(cmat,2)

    call dgemm(tra,trb,m,n,k,a,amat,lda,bmat,ldb,b,cmat,ldc)
  end subroutine dgemm_wrap
  pure subroutine sgemm_wrap(amat,bmat,cmat,transa,transb,alpha,beta)
    real(sp),intent(in) :: amat(:,:)
    real(sp),intent(in) :: bmat(:,:)
    real(sp),intent(inout) :: cmat(:,:)
    character(len=1),intent(in),optional :: transa
    character(len=1),intent(in),optional :: transb
    real(sp),intent(in),optional :: alpha
    real(sp),intent(in),optional :: beta
    character(len=1) :: tra,trb
    real(sp) :: a,b
    integer :: m,n,k,lda,ldb,ldc
    !> BLAS
    interface
      pure subroutine sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb, &
            & beta,c,ldc)
        import :: sp
        real(sp),intent(in) :: a(lda,*)
        real(sp),intent(in) :: b(ldb,*)
        real(sp),intent(inout) :: c(ldc,*)
        character(len=1),intent(in) :: transa
        character(len=1),intent(in) :: transb
        real(sp),intent(in) :: alpha
        real(sp),intent(in) :: beta
        integer,intent(in) :: m
        integer,intent(in) :: n
        integer,intent(in) :: k
        integer,intent(in) :: lda
        integer,intent(in) :: ldb
        integer,intent(in) :: ldc
      end subroutine sgemm
    end interface
    if (present(alpha)) then
      a = alpha
    else
      a = 1.0_sp
    end if
    if (present(beta)) then
      b = beta
    else
      b = 0.0_sp
    end if
    if (present(transa)) then
      tra = transa
    else
      tra = 'n'
    end if
    if (present(transb)) then
      trb = transb
    else
      trb = 'n'
    end if
    if ((tra .eq. 'n'.or.tra .eq. 'N')) then
      k = size(amat,2)
    else
      k = size(amat,1)
    end if
    lda = max(1,size(amat,1))
    ldb = max(1,size(bmat,1))
    ldc = max(1,size(cmat,1))
    m = size(cmat,1)
    n = size(cmat,2)
    call sgemm(tra,trb,m,n,k,a,amat,lda,bmat,ldb,b,cmat,ldc)
  end subroutine sgemm_wrap
!=======================================================================================!
! SYGVD computes all the eigenvalues, and optionally, the eigenvectors
! of a real generalized symmetric-definite eigenproblem, of the form
! A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
! B are assumed to be symmetric and B is also positive definite.
! If eigenvectors are desired, it uses a divide and conquer algorithm.
  subroutine ssygvd_wrap(itype,jobz,uplo,n,a,lda,b,ldb,w,info)
    real(sp),intent(inout) :: a(lda,*)
    real(sp),intent(inout) :: b(ldb,*)
    real(sp),intent(out) :: w(*)
    integer,intent(in) :: itype
    character(len=1),intent(in) :: jobz
    character(len=1),intent(in) :: uplo
    integer,intent(out) :: info
    integer,intent(in) :: n
    integer,intent(in) :: lda
    integer,intent(in) :: ldb
    real(sp),allocatable :: work(:)
    integer              :: lwork
    integer,allocatable  :: iwork(:)
    integer              :: liwork
    !> LAPACK
    interface
      pure subroutine ssygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork, &
            & iwork,liwork,info)
        import :: sp
        real(sp),intent(inout) :: a(lda,*)
        real(sp),intent(inout) :: b(ldb,*)
        real(sp),intent(out) :: w(*)
        integer,intent(in) :: itype
        character(len=1),intent(in) :: jobz
        character(len=1),intent(in) :: uplo
        integer,intent(out) :: info
        integer,intent(in) :: n
        integer,intent(in) :: lda
        integer,intent(in) :: ldb
        real(sp),intent(inout) :: work(*)
        integer,intent(in) :: lwork
        integer,intent(inout) :: iwork(*)
        integer,intent(in) :: liwork
      end subroutine ssygvd
    end interface
    allocate (work(1),iwork(1))
    !>--- workspace query
    call ssygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work, &
   &           -1,iwork,liwork,info)
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate (work,iwork)
    allocate (work(lwork),iwork(liwork))
    call ssygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work, &
   &           lwork,iwork,liwork,info)
    deallocate (work,iwork)
  end subroutine ssygvd_wrap

  subroutine dsygvd_wrap(itype,jobz,uplo,n,a,lda,b,ldb,w,info)
    real(wp),intent(inout) :: a(lda,*)
    real(wp),intent(inout) :: b(ldb,*)
    real(wp),intent(out) :: w(*)
    integer,intent(in) :: itype
    character(len=1),intent(in) :: jobz
    character(len=1),intent(in) :: uplo
    integer,intent(out) :: info
    integer,intent(in) :: n
    integer,intent(in)   :: lda
    integer,intent(in)   :: ldb
    real(wp),allocatable :: work(:)
    integer              :: lwork
    integer,allocatable  :: iwork(:)
    integer              :: liwork
    !> LAPACK
    interface
      pure subroutine dsygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork, &
            & iwork,liwork,info)
        import :: wp
        real(wp),intent(inout) :: a(lda,*)
        real(wp),intent(inout) :: b(ldb,*)
        real(wp),intent(out) :: w(*)
        integer,intent(in) :: itype
        character(len=1),intent(in) :: jobz
        character(len=1),intent(in) :: uplo
        integer,intent(out) :: info
        integer,intent(in) :: n
        integer,intent(in) :: lda
        integer,intent(in) :: ldb
        real(wp),intent(inout) :: work(*)
        integer,intent(in) :: lwork
        integer,intent(inout) :: iwork(*)
        integer,intent(in) :: liwork
      end subroutine dsygvd
    end interface
    allocate (work(1),iwork(1))
    !>--- workspace query
    call dsygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work, &
   &           -1,iwork,liwork,info)
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate (work,iwork)
    allocate (work(lwork),iwork(liwork))
    call dsygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work, &
   &           lwork,iwork,liwork,info)
    deallocate (work,iwork)
  end subroutine dsygvd_wrap

!=======================================================================================!
end module gfnff_math_wrapper
