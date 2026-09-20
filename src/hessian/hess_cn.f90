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
!> Weighted second derivative of the GFN-FF coordination number, shared by the
!> Hessians of all CN-dependent terms. logCN_i = L(cn_i) is a cut logarithm of
!> the erf pair sum cn_i (gfnff_dlogcoord), hence
!>   d2 logCN_i = L' d2 cn_i + L'' (d cn_i)(x)(d cn_i),  L'' = -L'(1 - L').
!> The translational sum rule cannot detect a missing rank-one term.
module gfnff_hess_cn

  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData,TGFFTopology
  use gfnff_neighbor,only:TNeigh
  use gfnff_hess_pair,only:pair_hess_block,scatter_pair_row, &
    &                      scatter_pair_hessian
  use gfnff_math_wrapper,only:gemm
  implicit none
  private

  public :: logcn_weighted_hessian,logcn_raw,hbcn_weighted_hessian

  real(wp),parameter :: sqrtpi = 1.77245385091_wp
  !> erf steepness of gfnff_dlogcoord
  real(wp),parameter :: kn = -7.5_wp
  !> hydrogen-bond CN, must match dncoord_erf: steepness, radius scaling, squared cutoff
  real(wp),parameter :: kn_hb = 27.5_wp
  real(wp),parameter :: rcov_scal_hb = 1.78_wp
  real(wp),parameter :: thr_hb = 900.0_wp

contains  !> MODULE PROCEDURES START HERE

  subroutine logcn_raw(n,at,xyz,srab,cnthr,param,cnr)
    !***********************************************************************
    !* Raw pair-sum coordination number cnr(n), before the cut logarithm.
    !* gfnff_dlogcoord returns only logCN, but L' and L'' need the raw value.
    !* srab: packed distances; cnthr: squared CN cutoff as in gfnff_dlogcoord.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: srab(n*(n+1)/2)
    real(wp),intent(in) :: cnthr
    type(TGFFData),intent(in) :: param
    real(wp),intent(out) :: cnr(n)

    integer :: i,j,ij,ii
    real(wp) :: r,r0,dr,thr,fc

    thr = sqrt(cnthr)
    cnr = 0.0_wp
    do i = 2,n
      ii = i*(i-1)/2
      do j = 1,i-1
        ij = ii+j
        r = srab(ij)
        if (r .gt. thr) cycle
        r0 = param%rcov(at(i))+param%rcov(at(j))
        dr = (r-r0)/r0
        fc = 0.5_wp*(1.0_wp+erf(kn*dr))
        cnr(i) = cnr(i)+fc
        cnr(j) = cnr(j)+fc
      end do
    end do

  end subroutine logcn_raw

  subroutine logcn_weighted_hessian(n,at,xyz,srab,cnthr,param,dcn,w,hess)
    !***********************************************************************
    !* Add sum_a w_a d2 logCN_a / dR dR to hess (3n,3n). Molecular only.
    !*   srab, cnthr - as in logcn_raw
    !*   dcn - dlogCN of gfnff_dlogcoord, dcn(:,m,a) = d logCN_a / d R_m,
    !*         (3,n,n) taken as the contiguous (3n,n) matrix
    !*   w   - (n) weights, typically dE/dlogCN_a of the consuming term
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: srab(n*(n+1)/2)
    real(wp),intent(in) :: cnthr
    type(TGFFData),intent(in) :: param
    real(wp),intent(in) :: dcn(3*n,n)
    real(wp),intent(in) :: w(n)
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: iat,jat,ij,ihi,ilo,ndof,a
    real(wp) :: r,r0,dr,thr,fp,fpp,ehat(3),vec(3),wgt
    real(wp) :: blk(3,3)
    real(wp),allocatable :: cnr(:),lp(:),wlp(:),ro(:),bmat(:,:)

    ndof = 3*n
    thr = sqrt(cnthr)

    allocate (cnr(n),lp(n),wlp(n),ro(n))
    call logcn_raw(n,at,xyz,srab,cnthr,param,cnr)

    do a = 1,n
      !>-- L'(cnr), written so that no exponential overflows
      lp(a) = 1.0_wp/(1.0_wp+exp(cnr(a)-param%cnmax))
      wlp(a) = w(a)*lp(a)
      !>-- w L'' (d cn)(x)(d cn) = ro (d logCN)(x)(d logCN), since dcn = L' d cn
      ro(a) = -w(a)*(1.0_wp-lp(a))/lp(a)
    end do

    !>-- pair part sum_a w_a L'_a d2 cn_a; a pair feeds both cn_iat and cn_jat
    !$omp parallel do default(none) shared(n, at, xyz, srab, thr, param, wlp, hess) &
    !$omp private(iat, jat, ij, ihi, ilo, r, r0, dr, fp, fpp, ehat, vec, wgt, blk) &
    !$omp schedule(dynamic)
    do iat = 1,n
      do jat = 1,n
        if (jat .eq. iat) cycle
        ihi = max(iat,jat)
        ilo = min(iat,jat)
        ij = ihi*(ihi-1)/2+ilo
        r = srab(ij)
        if (r .gt. thr.or.r .lt. 1.0e-6_wp) cycle
        r0 = param%rcov(at(iat))+param%rcov(at(jat))
        dr = (r-r0)/r0
        fp = kn/(sqrtpi*r0)*exp(-kn*kn*dr*dr)
        fpp = fp*(-2.0_wp*kn*kn*dr/r0)
        wgt = wlp(iat)+wlp(jat)
        vec = xyz(:,iat)-xyz(:,jat)
        ehat = vec/r
        call pair_hess_block(fp,fpp,r,ehat,blk)
        blk = wgt*blk
        call scatter_pair_row(hess,iat,jat,blk)
      end do
    end do
    !$omp end parallel do

    !>-- rank-one part sum_a ro_a (d logCN_a)(x)(d logCN_a) as one GEMM
    allocate (bmat(ndof,n))
    do a = 1,n
      bmat(:,a) = ro(a)*dcn(:,a)
    end do
    call gemm(bmat,dcn,hess,transb='T',alpha=1.0_wp,beta=1.0_wp)

  end subroutine logcn_weighted_hessian

  subroutine hbcn_weighted_hessian(n,at,xyz,param,topo,neigh,w,hess)
    !***********************************************************************
    !* Add sum_a w_a d2 hbcn_a / dR dR to hess (3n,3n) for the hydrogen-bond
    !* CN of dncoord_erf: a plain pair sum (no cut logarithm, so no rank-one
    !* term) over the H...B pairs in topo%bond_hb_AH / topo%bond_hb_B. Each
    !* pair counts toward both partners. Loop and cutoff mirror dncoord_erf.
    !*   w - (n) weights, dE/dhbcn_a of the consuming term
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(in) :: w(n)
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: i,j,iat,jat,iTrH,iTrB
    real(wp) :: rij(3),r2,r,rc,dr,fp,fpp,wgt,ehat(3)
    real(wp) :: blk(3,3)

    do i = 1,topo%bond_hb_nr
      iat = topo%bond_hb_AH(2,i)          !> the H atom
      iTrH = topo%bond_hb_AH(4,i)
      do j = 1,topo%bond_hb_Bn(i)
        jat = topo%bond_hb_B(1,j,i)       !> the B atom
        iTrB = topo%bond_hb_B(2,j,i)
        if (iTrB .gt. neigh%nTrans.or.iTrH .gt. neigh%nTrans) cycle
        wgt = w(iat)+w(jat)
        if (wgt .eq. 0.0_wp) cycle
        rij = (xyz(:,jat)+neigh%transVec(:,iTrB)) &
           & -(xyz(:,iat)+neigh%transVec(:,iTrH))
        r2 = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
        if (r2 .gt. thr_hb.or.r2 .lt. 1.0e-12_wp) cycle
        r = sqrt(r2)
        rc = rcov_scal_hb*(param%rcov(at(iat))+param%rcov(at(jat)))
        dr = r-rc
        !>-- fp is dncoord_erf's dtmp, fpp its derivative
        fp = -kn_hb/(sqrtpi*rc)*exp(-kn_hb*kn_hb*dr*dr/(rc*rc))
        fpp = fp*(-2.0_wp*kn_hb*kn_hb*dr/(rc*rc))
        ehat = rij/r
        call pair_hess_block(fp,fpp,r,ehat,blk)
        blk = wgt*blk
        call scatter_pair_hessian(hess,iat,jat,blk)
      end do
    end do

  end subroutine hbcn_weighted_hessian

end module gfnff_hess_cn
