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
!> Closed-form Hessian of the D3(BJ) dispersion energy of GFN-FF.
!> C6 is nonlinear in both CNs, so the CN-carrying parts are accumulated as
!>   K(a,b) = sum over pairs of d2E/dCN_a dCN_b        (n x n)
!>   V(:,a) = sum over pairs of (d2E/dx dCN_a) grad x  (3n x n)
!> and added as H += (dcn K + V) dcn^T + dcn V^T, not as O(N^4) rank-one updates.
module gfnff_hess_disp

  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData,TGFFTopology,TDispersionData
  use gfnff_hess_pair,only:pair_hess_block,scatter_pair_row
  use gfnff_hess_cn,only:logcn_weighted_hessian
  use gfnff_geometry,only:lin
  use gfnff_param,only:sqrtZr4r2
  use gfnff_math_wrapper,only:gemm
  implicit none
  private

  public :: hess_dispersion

  !> the Gaussian weighting exponent, hard-coded at the d3_gradient call site
  real(wp),parameter :: wf_d3 = 4.0_wp

contains  !> MODULE PROCEDURES START HERE

  subroutine hess_dispersion(n,at,xyz,sqrab,srab,dispthr,cnthr,cn,dcn, &
     &                       param,topo,hess)
    !***********************************************************************
    !* Add the closed-form D3(BJ) dispersion Hessian to hess, the
    !* counterpart of d3_gradient. Molecular systems only.
    !*   sqrab/srab    - packed squared and plain interatomic distances
    !*   dispthr/cnthr - squared dispersion (same as d3list) and CN cutoffs
    !*   cn            - (n) logarithmic coordination numbers
    !*   dcn           - dcn(:,m,a) = d logCN_a / d R_m as a (3n,n) matrix
    !* Each unordered pair is visited twice (own-row convention), so a thread
    !* writes only the rows and the K/V column of the atom it owns.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: sqrab(n*(n+1)/2),srab(n*(n+1)/2)
    real(wp),intent(in) :: dispthr,cnthr
    real(wp),intent(in) :: cn(n)
    real(wp),intent(in) :: dcn(3*n,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: iat,jat,ati,atj,ihi,ilo,ij,iref,jref,ndof,maxref,ia,ja
    real(wp) :: x,r1,r0,t6,t8,pp,pp1,pp2,sc,r4r2ij
    real(wp) :: c6,dc6i,dc6j,d2c6ii,d2c6ij,d2c6jj,refc6
    real(wp) :: ex,exx,eni,exni,ennii,ennij
    real(wp) :: fp,fpp,vec(3),ehat(3),blk(3,3),gx(3)
    real(wp),allocatable :: gw(:,:),dgw(:,:),d2gw(:,:)
    real(wp),allocatable :: w(:),kmat(:,:),vmat(:,:),wmat(:,:)

    ndof = 3*n
    maxref = maxval(topo%dispm%nref(at))
    allocate (gw(maxref,n),dgw(maxref,n),d2gw(maxref,n))
    call weight_references_d2(topo%dispm,n,at,wf_d3,cn,gw,dgw,d2gw)

    allocate (w(n),source=0.0_wp)
    allocate (kmat(n,n),source=0.0_wp)
    allocate (vmat(ndof,n),source=0.0_wp)

    !$omp parallel do default(none) schedule(dynamic) &
    !$omp shared(n, at, xyz, sqrab, srab, dispthr, param, topo, gw, dgw, d2gw, &
    !$omp&       hess, w, kmat, vmat) &
    !$omp private(iat, jat, ati, atj, ihi, ilo, ij, iref, jref, ia, ja, &
    !$omp&        x, r1, r0, t6, t8, pp, pp1, pp2, sc, r4r2ij, &
    !$omp&        c6, dc6i, dc6j, d2c6ii, d2c6ij, d2c6jj, refc6, &
    !$omp&        ex, exx, eni, exni, ennii, ennij, &
    !$omp&        fp, fpp, vec, ehat, blk, gx)
    do iat = 1,n
      ia = 3*(iat-1)
      ati = at(iat)
      do jat = 1,n
        if (jat .eq. iat) cycle
        ihi = max(iat,jat)
        ilo = min(iat,jat)
        ij = ihi*(ihi-1)/2+ilo
        x = sqrab(ij)
        if (x .ge. dispthr) cycle
        ja = 3*(jat-1)
        atj = at(jat)

        c6 = 0.0_wp
        dc6i = 0.0_wp
        dc6j = 0.0_wp
        d2c6ii = 0.0_wp
        d2c6ij = 0.0_wp
        d2c6jj = 0.0_wp
        do iref = 1,topo%dispm%nref(ati)
          do jref = 1,topo%dispm%nref(atj)
            refc6 = topo%dispm%c6(iref,jref,ati,atj)
            c6 = c6+gw(iref,iat)*gw(jref,jat)*refc6
            dc6i = dc6i+dgw(iref,iat)*gw(jref,jat)*refc6
            dc6j = dc6j+gw(iref,iat)*dgw(jref,jat)*refc6
            d2c6ii = d2c6ii+d2gw(iref,iat)*gw(jref,jat)*refc6
            d2c6ij = d2c6ij+dgw(iref,iat)*dgw(jref,jat)*refc6
            d2c6jj = d2c6jj+gw(iref,iat)*d2gw(jref,jat)*refc6
          end do
        end do

        !>-- radial factor P(x) and its two x derivatives
        r4r2ij = 3.0_wp*sqrtZr4r2(ati)*sqrtZr4r2(atj)
        r0 = param%d3r0(lin(ati,atj))
        t6 = 1.0_wp/(x**3+r0**3)
        t8 = 1.0_wp/(x**4+r0**4)
        pp = t6+2.0_wp*r4r2ij*t8
        pp1 = -3.0_wp*x*x*t6*t6-8.0_wp*r4r2ij*x**3*t8*t8
        pp2 = -6.0_wp*x*t6*t6+18.0_wp*x**4*t6**3 &
           & +2.0_wp*r4r2ij*(-12.0_wp*x*x*t8*t8+32.0_wp*x**6*t8**3)
        sc = topo%zetac6(ij)*param%dispscale

        !>-- partials of E = -C6(n_i,n_j) * sc * P(x)
        ex = -c6*sc*pp1
        exx = -c6*sc*pp2
        eni = -dc6i*sc*pp
        exni = -dc6i*sc*pp1
        ennii = -d2c6ii*sc*pp
        ennij = -d2c6ij*sc*pp

        !>-- radial pair block, converted from x = r^2 back to r
        r1 = srab(ij)
        fp = 2.0_wp*r1*ex
        fpp = 4.0_wp*x*exx+2.0_wp*ex
        vec = xyz(:,iat)-xyz(:,jat)
        ehat = vec/r1
        call pair_hess_block(fp,fpp,r1,ehat,blk)
        call scatter_pair_row(hess,iat,jat,blk)

        !>-- CN weight, identical to the gradient's dEdcn
        w(iat) = w(iat)+eni

        kmat(iat,iat) = kmat(iat,iat)+ennii
        kmat(iat,jat) = kmat(iat,jat)+ennij

        !>-- distance-CN cross term; grad x = +2 vec at iat, -2 vec at jat
        gx = 2.0_wp*exni*vec
        vmat(ia+1:ia+3,iat) = vmat(ia+1:ia+3,iat)+gx
        vmat(ja+1:ja+3,iat) = vmat(ja+1:ja+3,iat)-gx
      end do
    end do
    !$omp end parallel do

    call logcn_weighted_hessian(n,at,xyz,srab,cnthr,param,dcn,w,hess)

    !>-- H += (dcn K + V) dcn^T + dcn V^T
    allocate (wmat(ndof,n))
    wmat = vmat
    call gemm(dcn,kmat,wmat,alpha=1.0_wp,beta=1.0_wp)
    call gemm(wmat,dcn,hess,transb='T',alpha=1.0_wp,beta=1.0_wp)
    call gemm(dcn,vmat,hess,transb='T',alpha=1.0_wp,beta=1.0_wp)

  end subroutine hess_dispersion

  subroutine weight_references_d2(dispm,n,at,wf,cn,gw,dgw,d2gw)
    !***********************************************************************
    !* Gaussian reference weights w_k = e_k/Z, e_k = exp(-wf (CN - CN_k)^2),
    !* and their first two CN derivatives, each (maxref,n). Value and first
    !* derivative match weight_references_d4 exactly, NaN fallbacks included.
    !***********************************************************************
    implicit none
    type(TDispersionData),intent(in) :: dispm
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: wf
    real(wp),intent(in) :: cn(n)
    real(wp),intent(out) :: gw(:,:),dgw(:,:),d2gw(:,:)

    integer :: iat,ati,iref
    real(wp) :: zz,z1,z2,zinv,aa,ee,e1,e2,gwk,dgwk,d2gwk

    gw = 0.0_wp
    dgw = 0.0_wp
    d2gw = 0.0_wp

    do iat = 1,n
      ati = at(iat)
      zz = 0.0_wp
      z1 = 0.0_wp
      z2 = 0.0_wp
      do iref = 1,dispm%nref(ati)
        aa = 2.0_wp*wf*(dispm%cn(iref,ati)-cn(iat))
        ee = exp(-wf*(cn(iat)-dispm%cn(iref,ati))**2)
        zz = zz+ee
        z1 = z1+aa*ee
        z2 = z2+(aa*aa-2.0_wp*wf)*ee
      end do
      zinv = 1.0_wp/zz

      do iref = 1,dispm%nref(ati)
        aa = 2.0_wp*wf*(dispm%cn(iref,ati)-cn(iat))
        ee = exp(-wf*(cn(iat)-dispm%cn(iref,ati))**2)
        e1 = aa*ee
        e2 = (aa*aa-2.0_wp*wf)*ee

        gwk = ee*zinv
        if (gwk /= gwk) then
          !>-- weight_references_d4 fallback for an under- or overflowing Z
          if (maxval(dispm%cn(:dispm%nref(ati),ati)) == dispm%cn(iref,ati)) then
            gwk = 1.0_wp
          else
            gwk = 0.0_wp
          end if
        end if
        gw(iref,iat) = gwk

        dgwk = e1*zinv-ee*z1*zinv*zinv
        if (dgwk /= dgwk) dgwk = 0.0_wp
        dgw(iref,iat) = dgwk

        d2gwk = e2*zinv-2.0_wp*e1*z1*zinv*zinv &
           & -ee*z2*zinv*zinv+2.0_wp*ee*z1*z1*zinv*zinv*zinv
        !>-- NaN also marks a weight replaced by a constant, whose derivatives vanish
        if (d2gwk /= d2gwk) d2gwk = 0.0_wp
        d2gw(iref,iat) = d2gwk
      end do
    end do

  end subroutine weight_references_d2

end module gfnff_hess_disp
