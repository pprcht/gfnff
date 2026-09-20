! ------------------------------------------------------------------------------
! This file is part of gfnff.
!
! Copyright (C) 2026 Philipp Pracht
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
! ------------------------------------------------------------------------------
!> Quality metrics and harmonic frequencies for a Cartesian nuclear Hessian.
module gfnff_hess_analysis

  use iso_fortran_env,only:wp => real64
  implicit none
  private

  public :: hess_asymmetry,hess_transl_sumrule
  public :: hess_symmetrize,hess_frequencies
  public :: atomic_mass

  !> amu -> atomic units of mass (electron masses)
  real(wp),parameter :: amu2au = 1822.888486209_wp
  !> sqrt(Eh / (bohr^2 m_e)) -> cm^-1
  real(wp),parameter :: au2rcm = 219474.6313705_wp

contains  !> MODULE PROCEDURES START HERE

  pure function hess_asymmetry(hess) result(amax)
    !***********************************************************************
    !* max|H - H^T|. Identically zero after symmetrisation, so call it before.
    !***********************************************************************
    real(wp),intent(in) :: hess(:,:)
    real(wp) :: amax
    amax = maxval(abs(hess-transpose(hess)))
  end function hess_asymmetry

  pure function hess_transl_sumrule(hess) result(smax)
    !***********************************************************************
    !* Largest violation of the translational sum rule
    !*   sum_B H(a,(c,B)) = 0   for every row a and direction c.
    !***********************************************************************
    real(wp),intent(in) :: hess(:,:)
    real(wp) :: smax

    integer :: ndof,nat,ia,c,iat,ib
    real(wp) :: s

    ndof = size(hess,1)
    nat = ndof/3
    smax = 0.0_wp
    do ia = 1,ndof
      do c = 1,3
        s = 0.0_wp
        do iat = 1,nat
          ib = 3*(iat-1)+c
          s = s+hess(ia,ib)
        end do
        smax = max(smax,abs(s))
      end do
    end do

  end function hess_transl_sumrule

  subroutine hess_symmetrize(hess)
    !***********************************************************************
    !* Replace H by (H + H^T)/2 in place. The array form with transpose() would
    !* materialise a second full matrix and not thread. Iteration j owns column
    !* j above the diagonal and row j to its left, so the loop is race-free.
    !***********************************************************************
    real(wp),intent(inout) :: hess(:,:)
    integer :: i,j,ndof
    real(wp) :: s

    ndof = size(hess,1)
    !>-- triangular work per iteration, hence guided rather than static
    !$omp parallel do default(none) shared(hess,ndof) private(i,j,s) &
    !$omp schedule(guided)
    do j = 1,ndof
      do i = 1,j-1
        s = 0.5_wp*(hess(i,j)+hess(j,i))
        hess(i,j) = s
        hess(j,i) = s
      end do
    end do
    !$omp end parallel do
  end subroutine hess_symmetrize

  subroutine hess_frequencies(hess,nat,at,freq,modes,iostat)
    !***********************************************************************
    !* Harmonic frequencies from the mass-weighted Cartesian Hessian.
    !*   hess   - (3*nat,3*nat) Cartesian Hessian in Eh/bohr^2
    !*   freq   - cm^-1 by ascending |freq| (trans/rot first), imaginary negative
    !*   modes  - optional mass-weighted eigenvectors, same column order
    !*   iostat - optional, non-zero if the diagonalisation failed
    !***********************************************************************
    real(wp),intent(in) :: hess(:,:)
    integer,intent(in) :: nat,at(nat)
    real(wp),allocatable,intent(out) :: freq(:)
    real(wp),allocatable,intent(out),optional :: modes(:,:)
    integer,intent(out),optional :: iostat

    real(wp),allocatable :: fmw(:,:),eval(:),mem(:)
    integer,allocatable :: idx(:)
    integer :: ndof,iat,ic,ib,jb,i,info

    ndof = 3*nat
    allocate (fmw(ndof,ndof),eval(ndof),mem(ndof),freq(ndof))

    do iat = 1,nat
      do ic = 1,3
        ib = 3*(iat-1)+ic
        mem(ib) = atomic_mass(at(iat))*amu2au
      end do
    end do

    do jb = 1,ndof
      do ib = 1,ndof
        fmw(ib,jb) = hess(ib,jb)/sqrt(mem(ib)*mem(jb))
      end do
    end do

    call symeig(fmw,eval,info)
    if (present(iostat)) iostat = info
    if (info /= 0) then
      freq = 0.0_wp
      if (present(modes)) then
        allocate (modes(ndof,ndof),source=0.0_wp)
      end if
      return
    end if

    do i = 1,ndof
      freq(i) = sign(sqrt(abs(eval(i))),eval(i))*au2rcm
    end do

    idx = argsort_abs(freq)
    freq = freq(idx)
    if (present(modes)) then
      allocate (modes(ndof,ndof))
      do i = 1,ndof
        modes(:,i) = fmw(:,idx(i))
      end do
    end if

  end subroutine hess_frequencies

  subroutine symeig(a,w,info)
    !***********************************************************************
    !* Eigenvalues w and in-place eigenvectors of symmetric a via LAPACK dsyev.
    !***********************************************************************
    real(wp),intent(inout) :: a(:,:)
    real(wp),intent(out) :: w(:)
    integer,intent(out) :: info

    integer :: n,lwork
    real(wp) :: wq(1)
    real(wp),allocatable :: work(:)
    external :: dsyev

    n = size(a,1)
    call dsyev('V','U',n,a,n,w,wq,-1,info)
    if (info /= 0) return
    lwork = max(1,nint(wq(1)))
    allocate (work(lwork))
    call dsyev('V','U',n,a,n,w,work,lwork,info)

  end subroutine symeig

  pure function argsort_abs(x) result(idx)
    !***********************************************************************
    !* Indices sorting x by ascending |x|; insertion sort, arrays are 3N.
    !***********************************************************************
    real(wp),intent(in) :: x(:)
    integer,allocatable :: idx(:)
    integer :: n,i,j,k

    n = size(x)
    allocate (idx(n))
    idx = [(i,i=1,n)]
    do i = 2,n
      k = idx(i)
      j = i-1
      do while (j >= 1)
        if (abs(x(idx(j))) <= abs(x(k))) exit
        idx(j+1) = idx(j)
        j = j-1
      end do
      idx(j+1) = k
    end do

  end function argsort_abs

  pure function atomic_mass(z) result(m)
    !***********************************************************************
    !* Standard atomic weight in amu for Z = 1..118 (IUPAC 2021; most stable
    !* isotope where no standard weight exists). Returns 1 outside that range.
    !***********************************************************************
    integer,intent(in) :: z
    real(wp) :: m
!&<
    real(wp),parameter :: w(118) = [ &
       &   1.008_wp,   4.0026_wp,  6.94_wp,    9.0122_wp, 10.81_wp,   12.011_wp, &
       &  14.007_wp,  15.999_wp,  18.998_wp,  20.180_wp, 22.990_wp,  24.305_wp, &
       &  26.982_wp,  28.085_wp,  30.974_wp,  32.06_wp,  35.45_wp,   39.948_wp, &
       &  39.098_wp,  40.078_wp,  44.956_wp,  47.867_wp, 50.942_wp,  51.996_wp, &
       &  54.938_wp,  55.845_wp,  58.933_wp,  58.693_wp, 63.546_wp,  65.38_wp,  &
       &  69.723_wp,  72.630_wp,  74.922_wp,  78.971_wp, 79.904_wp,  83.798_wp, &
       &  85.468_wp,  87.62_wp,   88.906_wp,  91.224_wp, 92.906_wp,  95.95_wp,  &
       &  97.0_wp,   101.07_wp,  102.91_wp,  106.42_wp, 107.87_wp,  112.41_wp,  &
       & 114.82_wp,  118.71_wp,  121.76_wp,  127.60_wp, 126.90_wp,  131.29_wp,  &
       & 132.91_wp,  137.33_wp,  138.91_wp,  140.12_wp, 140.91_wp,  144.24_wp,  &
       & 145.0_wp,   150.36_wp,  151.96_wp,  157.25_wp, 158.93_wp,  162.50_wp,  &
       & 164.93_wp,  167.26_wp,  168.93_wp,  173.05_wp, 174.97_wp,  178.49_wp,  &
       & 180.95_wp,  183.84_wp,  186.21_wp,  190.23_wp, 192.22_wp,  195.08_wp,  &
       & 196.97_wp,  200.59_wp,  204.38_wp,  207.2_wp,  208.98_wp,  209.0_wp,   &
       & 210.0_wp,   222.0_wp,   223.0_wp,   226.0_wp,  227.0_wp,   232.04_wp,  &
       & 231.04_wp,  238.03_wp,  237.0_wp,   244.0_wp,  243.0_wp,   247.0_wp,   &
       & 247.0_wp,   251.0_wp,   252.0_wp,   257.0_wp,  258.0_wp,   259.0_wp,   &
       & 262.0_wp,   267.0_wp,   270.0_wp,   269.0_wp,  270.0_wp,   270.0_wp,   &
       & 278.0_wp,   281.0_wp,   281.0_wp,   285.0_wp,  286.0_wp,   289.0_wp,   &
       & 289.0_wp,   293.0_wp,   293.0_wp,   294.0_wp]
!&>
    if (z >= 1.and.z <= 118) then
      m = w(z)
    else
      m = 1.0_wp
    end if

  end function atomic_mass

end module gfnff_hess_analysis
