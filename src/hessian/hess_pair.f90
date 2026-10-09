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
!> Building blocks for nuclear Hessians of two-body radial energies f(r):
!> the 3x3 pair block and its scatter into the 3N x 3N matrix.
module gfnff_hess_pair

  use iso_fortran_env,only:wp => real64
  implicit none
  private

  public :: pair_hess_block,scatter_pair_hessian,scatter_pair_row

contains  !> MODULE PROCEDURES START HERE

  pure subroutine pair_hess_block(fp,fpp,r1,ehat,blk)
    !***********************************************************************
    !* 3x3 block d2f/dR_i dR_i of a radial pair energy f(r), fp = f', fpp = f'':
    !*   blk = fpp e e^T + (fp/r1)(I - e e^T),  e = ehat = (R_i - R_j)/r1
    !* The cross blocks d2f/dR_i dR_j are -blk.
    !***********************************************************************
    real(wp),intent(in) :: fp,fpp,r1
    real(wp),intent(in) :: ehat(3)
    real(wp),intent(out) :: blk(3,3)

    integer :: a,b
    real(wp) :: fpr,proj

    fpr = fp/r1
    do b = 1,3
      do a = 1,3
        proj = merge(1.0_wp,0.0_wp,a == b)-ehat(a)*ehat(b)
        blk(a,b) = fpp*ehat(a)*ehat(b)+fpr*proj
      end do
    end do

  end subroutine pair_hess_block

  pure subroutine scatter_pair_hessian(hess,iat,jat,blk)
    !***********************************************************************
    !* Add pair block B to hess for a loop that visits each unordered pair
    !* once: H_ii += B, H_jj += B, H_ij -= B, H_ji -= B.
    !***********************************************************************
    real(wp),intent(inout) :: hess(:,:)
    integer,intent(in) :: iat,jat
    real(wp),intent(in) :: blk(3,3)

    integer :: a,b,ia,ja

    ia = 3*(iat-1)
    ja = 3*(jat-1)
    do b = 1,3
      do a = 1,3
        hess(ia+a,ia+b) = hess(ia+a,ia+b)+blk(a,b)
        hess(ja+a,ja+b) = hess(ja+a,ja+b)+blk(a,b)
        hess(ia+a,ja+b) = hess(ia+a,ja+b)-blk(a,b)
        hess(ja+a,ia+b) = hess(ja+a,ia+b)-blk(a,b)
      end do
    end do

  end subroutine scatter_pair_hessian

  pure subroutine scatter_pair_row(hess,iat,jat,blk)
    !***********************************************************************
    !* Own-row half of scatter_pair_hessian: H_ii += B, H_ij -= B only. A loop
    !* over all ordered pairs gives the same matrix, race-free per owned iat.
    !***********************************************************************
    real(wp),intent(inout) :: hess(:,:)
    integer,intent(in) :: iat,jat
    real(wp),intent(in) :: blk(3,3)

    integer :: a,b,ia,ja

    ia = 3*(iat-1)
    ja = 3*(jat-1)
    do b = 1,3
      do a = 1,3
        hess(ia+a,ia+b) = hess(ia+a,ia+b)+blk(a,b)
        hess(ia+a,ja+b) = hess(ia+a,ja+b)-blk(a,b)
      end do
    end do

  end subroutine scatter_pair_row

end module gfnff_hess_pair
