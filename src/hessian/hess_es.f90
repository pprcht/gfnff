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
!> Closed-form Hessian of the GFN-FF EEQ electrostatics, the only term whose
!> energy depends on the geometry implicitly, via the charges q that solve the
!> bordered system M [q; lam] = [x; Q]:
!>   E = 1/2 q^T A q - q^T x,  A_ij = erf(gam_ij r_ij)/r_ij,  x_i = chi_i + cnf_i sqrt(CN_i)
!>   H_ab = 1/2 q^T (d_a d_b A) q - q^T (d_a d_b x) - rho_a^T M^-1 rho_b,
!>   rho_b = (d_b A) q - d_b x
module gfnff_hess_es

  use iso_fortran_env,only:wp => real64,stderr => error_unit
  use gfnff_data_types,only:TGFFData,TGFFTopology
  use gfnff_hess_pair,only:pair_hess_block,scatter_pair_row
  use gfnff_hess_cn,only:logcn_weighted_hessian
  use gfnff_math_wrapper,only:gemm,sytrf_wrap,sytrs_wrap
  use gfnff_solv_gbsa,only:TBorn
  use gfnff_hess_solv,only:hess_solvation
  implicit none
  private

  public :: hess_electrostatics

  real(wp),parameter :: sqrtpi = 1.77245385091_wp
  real(wp),parameter :: tsqrt2pi = 0.797884560802866_wp
  !> same regulariser as the energy code's 1/(2 sqrt(CN)), so the Hessian matches that gradient
  real(wp),parameter :: cnreg = 1.0e-16_wp

contains  !> MODULE PROCEDURES START HERE

  subroutine hess_electrostatics(n,at,xyz,srab,cnthr,cn,dcn,q,param,topo,hess, &
     &                           iostat,gbsa)
    !***********************************************************************
    !* Add the EEQ Hessian incl. charge response to hess (3n,3n); molecular only.
    !*   srab   - packed interatomic distances; cnthr - squared CN cutoff
    !*   cn,dcn - logCN and dcn(:,m,a) = d logCN_a / d R_m as the (3n,n) matrix
    !*   q      - EEQ charges of goed_gfnff, assumed exact; with accuracy > 1
    !*            they come from a single-precision solve, which limits the
    !*            agreement with finite differences to about 1e-7
    !*   gbsa   - optional implicit solvation, adds to hess, rho and M
    !*   iostat - optional, non-zero if the EEQ solve failed
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: srab(n*(n+1)/2)
    real(wp),intent(in) :: cnthr
    real(wp),intent(in) :: cn(n)
    real(wp),intent(in) :: dcn(3*n,n)
    real(wp),intent(in) :: q(n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    real(wp),intent(inout) :: hess(3*n,3*n)
    integer,intent(out),optional :: iostat
    type(TBorn),intent(in),optional :: gbsa

    integer :: iat,jat,ihi,ilo,ij,ndof,m,ia,ja,io1,io2
    real(wp) :: r1,gam,ee,erfv,fp,fpp,qq,vec(3),ehat(3),tmp(3)
    real(wp) :: blk(3,3),dxdcn,ai
    real(wp),allocatable :: w(:),rone(:),amat(:,:),rmat(:,:),ymat(:,:),cmat(:,:)
    integer,allocatable :: ipiv(:)

    ndof = 3*n
    m = n+topo%nfrag
    if (present(iostat)) iostat = 0

    allocate (rmat(m,ndof),source=0.0_wp)

    !>-- 1. pair part 1/2 q^T (d2 A) q and pair half of rho; iteration iat owns row iat of both
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp shared(n, xyz, srab, q, topo, hess, rmat) &
    !$omp private(iat, jat, ihi, ilo, ij, r1, gam, ee, erfv, fp, fpp, qq, &
    !$omp&        vec, ehat, tmp, blk, ia, ja)
    do iat = 1,n
      ia = 3*(iat-1)
      do jat = 1,n
        if (jat .eq. iat) cycle
        ihi = max(iat,jat)
        ilo = min(iat,jat)
        ij = ihi*(ihi-1)/2+ilo
        r1 = srab(ij)
        if (r1 .lt. 1.0e-6_wp) cycle
        ja = 3*(jat-1)
        gam = 1.0_wp/sqrt(topo%alpeeq(iat)+topo%alpeeq(jat))
        ee = exp(-gam*gam*r1*r1)
        erfv = erf(gam*r1)
        fp = 2.0_wp*gam*ee/(sqrtpi*r1)-erfv/(r1*r1)
        fpp = -4.0_wp*gam*gam*gam*ee/sqrtpi &
           & -4.0_wp*gam*ee/(sqrtpi*r1*r1)+2.0_wp*erfv/(r1*r1*r1)
        vec = xyz(:,iat)-xyz(:,jat)
        ehat = vec/r1
        qq = q(iat)*q(jat)
        call pair_hess_block(qq*fp,qq*fpp,r1,ehat,blk)
        call scatter_pair_row(hess,iat,jat,blk)
        !>-- pair term of rho_i = sum_j q_j dA_ij/dR, opposite signs on iat and jat
        tmp = q(jat)*(fp/r1)*vec
        rmat(iat,ia+1:ia+3) = rmat(iat,ia+1:ia+3)+tmp
        rmat(iat,ja+1:ja+3) = rmat(iat,ja+1:ja+3)-tmp
      end do
    end do
    !$omp end parallel do

    !>-- 2. CN part -q^T d2 x and CN half of rho. x_i is nonlinear in CN:
    !>   d2 x_i = cnf_i/(2 sqrt CN_i) d2 CN_i - cnf_i/(4 CN_i^3/2) dCN (x) dCN
    allocate (w(n),rone(n))
    do iat = 1,n
      dxdcn = param%cnf(at(iat))/(2.0_wp*sqrt(cn(iat))+cnreg)
      ai = q(iat)*dxdcn
      w(iat) = -ai
      rone(iat) = ai/(2.0_wp*cn(iat)+cnreg)
      rmat(iat,:) = rmat(iat,:)-dxdcn*dcn(:,iat)
    end do

    call logcn_weighted_hessian(n,at,xyz,srab,cnthr,param,dcn,w,hess)

    !>-- rank-one piece sum_i rone_i (d logCN_i)(x)(d logCN_i) as one GEMM
    allocate (cmat(ndof,n))
    do iat = 1,n
      cmat(:,iat) = rone(iat)*dcn(:,iat)
    end do
    call gemm(cmat,dcn,hess,transb='T',alpha=1.0_wp,beta=1.0_wp)
    deallocate (cmat)

    if (present(gbsa)) then
      call hess_solvation(n,xyz,q,gbsa,hess,rmat)
    end if

    !>-- 3. response -rho^T M^-1 rho, all 3N right-hand sides in one solve. M is rebuilt
    !>   here to leave the energy path untouched; the border rows of rho are zero.
    allocate (ymat(m,ndof))
    ymat = rmat
    call build_eeq_matrix(n,srab,param,topo,amat)
    !>-- the Born matrix screens the same charges, so it enters M
    if (present(gbsa)) then
      amat(1:n,1:n) = amat(1:n,1:n)+gbsa%bornMat(1:n,1:n)
    end if
    allocate (ipiv(m))
    io2 = 0
    call sytrf_wrap(amat,ipiv,io1)
    if (io1 .eq. 0) call sytrs_wrap(amat,ymat,ipiv,io2)
    if (io1 .ne. 0.or.io2 .ne. 0) then
      write (stderr,'("**ERROR** EEQ response solve failed in hess_electrostatics")')
      if (present(iostat)) iostat = 1
      return
    end if

    call gemm(rmat,ymat,hess,transa='T',alpha=-1.0_wp,beta=1.0_wp)

  end subroutine hess_electrostatics

  subroutine build_eeq_matrix(n,srab,param,topo,amat)
    !***********************************************************************
    !* Bordered EEQ matrix M, identical to the one goed_gfnff factorises for
    !* the charges. amat is allocated here as (n+nfrag, n+nfrag).
    !***********************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: srab(n*(n+1)/2)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    real(wp),allocatable,intent(out) :: amat(:,:)

    integer :: i,j,k,ij,m
    real(wp) :: gam

    m = n+topo%nfrag
    allocate (amat(m,m),source=0.0_wp)

    !$omp parallel do default(none) schedule(dynamic) &
    !$omp shared(n, srab, topo, amat) private(i, j, k, ij, gam)
    do i = 1,n
      amat(i,i) = tsqrt2pi/sqrt(topo%alpeeq(i))+topo%gameeq(i)
      k = i*(i-1)/2
      do j = 1,i-1
        ij = k+j
        gam = 1.0_wp/sqrt(topo%alpeeq(i)+topo%alpeeq(j))
        amat(j,i) = erf(gam*srab(ij))/srab(ij)
        amat(i,j) = amat(j,i)
      end do
    end do
    !$omp end parallel do

    do i = 1,topo%nfrag
      do j = 1,n
        if (topo%fraglist(j) .eq. i) then
          amat(n+i,j) = 1.0_wp
          amat(j,n+i) = 1.0_wp
        end if
      end do
    end do

  end subroutine build_eeq_matrix

end module gfnff_hess_es
