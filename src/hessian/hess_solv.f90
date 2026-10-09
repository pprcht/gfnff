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
!> Closed-form second derivatives of the ALPB implicit solvation energy.
!> The model is C1 but not C2: the SASA switching cubic and the descreening
!> overlap branch jump in curvature. The result is the exact one-sided second
!> derivative; a finite difference across such a surface does not converge.
module gfnff_hess_solv

  use iso_fortran_env,only:wp => real64
  use gfnff_solv_gbsa,only:TBorn
  use gfnff_solv_sasa,only:ah0,ah1,ah3,tolsesp
  use gfnff_hess_pair,only:pair_hess_block,scatter_pair_row
  use gfnff_math_wrapper,only:gemm
  implicit none
  private

  public :: hess_solvation
  public :: born_weighted_hessian
  public :: born_radius_derivs
  public :: p16_kernel_derivs
  public :: sasa_weighted_hessian

  !> GBOBCII rescaling constants, as in the Born radii integrator
  real(wp),parameter :: alp = 1.0_wp
  real(wp),parameter :: bet = 0.8_wp
  real(wp),parameter :: gam = 4.85_wp

  !> P16 zeta parameter, as in the interaction kernel
  real(wp),parameter :: zetaP16 = 1.028_wp
  real(wp),parameter :: zetaP16o16 = zetaP16/16.0_wp

contains  !> MODULE PROCEDURES START HERE

  subroutine hess_solvation(n,xyz,q,gbsa,hess,rmat)
    !***********************************************************************
    !* Add the solvent part of the second derivative: explicit terms into
    !* hess, the extra response source into rows 1..n of rmat. Radius-coupled
    !* sums go into dense kmat, vmat, cmat and are contracted with dbdr by GEMM.
    !* q: EEQ charges of the energy code. gbsa: updated at this geometry.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: q(n)
    type(TBorn),intent(in) :: gbsa
    real(wp),intent(inout) :: hess(3*n,3*n)
    real(wp),intent(inout) :: rmat(:,:)

    integer :: iat,jat,ia,ja,ndof
    real(wp) :: r1,dr2,vec(3),ehat(3),blk(3,3),qq,keps
    real(wp) :: v,d1(3),d2(3,3),bi,bj,tmp(3)
    real(wp),allocatable :: kmat(:,:),cmat(:,:),vmat(:,:),wmat(:,:)
    real(wp),allocatable :: dEdb(:),wsasa(:),dbdr(:,:)

    ndof = 3*n
    keps = gbsa%keps

    allocate (kmat(n,n),source=0.0_wp)
    allocate (cmat(n,n),source=0.0_wp)
    allocate (vmat(ndof,n),source=0.0_wp)
    allocate (dEdb(n),source=0.0_wp)

    allocate (dbdr(ndof,n))
    dbdr = reshape(gbsa%brdr, [ndof,n])

    !>-- 1. P16 pair kernel; thread iat owns row iat (column iat of vmat)
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp shared(n,xyz,q,gbsa,hess,rmat,kmat,cmat,vmat,dEdb,keps,ndof) &
    !$omp private(iat,jat,ia,ja,r1,dr2,vec,ehat,blk,qq,v,d1,d2,bi,bj,tmp)
    do iat = 1,n
      ia = 3*(iat-1)
      bi = gbsa%brad(iat)
      do jat = 1,n
        if (jat .eq. iat) cycle
        ja = 3*(jat-1)
        vec = xyz(:,iat)-xyz(:,jat)
        dr2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
        if (dr2 .lt. 1.0e-12_wp) cycle
        r1 = sqrt(dr2)
        ehat = vec/r1
        bj = gbsa%brad(jat)
        qq = q(iat)*q(jat)

        call p16_kernel_derivs(r1,bi,bj,keps,v,d1,d2)

        call pair_hess_block(qq*d1(1),qq*d2(1,1),r1,ehat,blk)
        call scatter_pair_row(hess,iat,jat,blk)

        tmp = qq*d2(1,2)*ehat
        vmat(ia+1:ia+3,iat) = vmat(ia+1:ia+3,iat)+tmp
        vmat(ja+1:ja+3,iat) = vmat(ja+1:ja+3,iat)-tmp

        kmat(iat,iat) = kmat(iat,iat)+qq*d2(2,2)
        kmat(iat,jat) = kmat(iat,jat)+qq*d2(2,3)

        dEdb(iat) = dEdb(iat)+qq*d1(2)

        !>-- rho: distance half directly, radius half through cmat
        tmp = q(jat)*d1(1)*ehat
        rmat(iat,ia+1:ia+3) = rmat(iat,ia+1:ia+3)+tmp
        rmat(iat,ja+1:ja+3) = rmat(iat,ja+1:ja+3)-tmp
        cmat(iat,iat) = cmat(iat,iat)+q(jat)*d1(2)
        cmat(iat,jat) = cmat(iat,jat)+q(jat)*d1(3)
      end do
    end do
    !$omp end parallel do

    !>-- 2. Born self energy, 1/2 q^2 keps / b
    do iat = 1,n
      bi = gbsa%brad(iat)
      dEdb(iat) = dEdb(iat)-0.5_wp*q(iat)*q(iat)*keps/(bi*bi)
      kmat(iat,iat) = kmat(iat,iat)+q(iat)*q(iat)*keps/(bi*bi*bi)
      cmat(iat,iat) = cmat(iat,iat)-q(iat)*keps/(bi*bi)
    end do

    !>-- 3. surface term and hydrogen bond correction, both linear in sasa
    allocate (wsasa(n))
    do iat = 1,n
      wsasa(iat) = gbsa%gamsasa(iat)
      if (gbsa%lhb) wsasa(iat) = wsasa(iat)+q(iat)*q(iat)*gbsa%dhbdw(iat)
    end do
    call sasa_weighted_hessian(n,xyz,gbsa,wsasa,hess)

    if (gbsa%lhb) then
      !>-- B_ii carries 2 hbw_i, so (Bq)_i picks up 2 q_i dhbdw_i sasa_i
      do iat = 1,n
        rmat(iat,1:ndof) = rmat(iat,1:ndof) &
           & +2.0_wp*q(iat)*gbsa%dhbdw(iat)*reshape(gbsa%dsdrt(:,:,iat), [ndof])
      end do
    end if

    call adet_shape_hessian(n,xyz,sum(q),gbsa,hess,rmat)

    !>-- 4. adds the rank-one part to diag(kmat), so it must precede step 5
    call born_weighted_hessian(n,xyz,gbsa,dEdb,kmat,hess)

    !>-- 5. contract everything coupled through the Born radii
    allocate (wmat(ndof,n))
    wmat = vmat
    call gemm(dbdr,kmat,wmat,alpha=1.0_wp,beta=1.0_wp)
    call gemm(wmat,dbdr,hess,transb='T',alpha=1.0_wp,beta=1.0_wp)
    call gemm(dbdr,vmat,hess,transb='T',alpha=1.0_wp,beta=1.0_wp)

    call gemm(cmat,dbdr,rmat(1:n,1:ndof),transb='T',alpha=1.0_wp,beta=1.0_wp)

  end subroutine hess_solvation

  subroutine adet_shape_hessian(n,xyz,qtot,gbsa,hess,rmat)
    !***********************************************************************
    !* ALPB shape correction E = 1/2 keps beta Q^2 / A_det, added to hess and
    !* to rows 1..n of rmat; returns early for a neutral solute. E = K S^(-1/6)
    !* with S = det(I) of the radius-weighted inertia tensor, and
    !*   d2S = C:d2I + dS dS / S - S tr(I^-1 dI I^-1 dI),  C the cofactors.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: qtot
    type(TBorn),intent(in) :: gbsa
    real(wp),intent(inout) :: hess(3*n,3*n)
    real(wp),intent(inout) :: rmat(:,:)

    real(wp),parameter :: tof = 2.0_wp/5.0_wp
    real(wp),parameter :: unity(3,3) = reshape( &
       & [1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp], [3,3])

    integer :: iat,jat,ndof,p,q,ia,ja,m,k
    real(wp) :: rad2,rad3,totRad3,center(3),vec(3),r2
    real(wp) :: inert(3,3),cof(3,3),tmat(3,3),iinv(3,3),mm(3,3),amat3(3,3)
    real(wp) :: sdet,kk,c1,c2,pref,trc,fac
    real(wp),allocatable :: rad3v(:),gvec(:),pm(:,:),qm(:,:)

    ndof = 3*n
    if (abs(qtot) .lt. 1.0e-12_wp) return
    if (gbsa%alpbet .le. 0.0_wp) return

    allocate (rad3v(n),gvec(ndof))

    totRad3 = 0.0_wp
    center = 0.0_wp
    do iat = 1,n
      rad3 = gbsa%vdwr(iat)**3
      rad3v(iat) = rad3
      totRad3 = totRad3+rad3
      center = center+xyz(:,iat)*rad3
    end do
    center = center/totRad3

    inert = 0.0_wp
    do iat = 1,n
      rad2 = gbsa%vdwr(iat)*gbsa%vdwr(iat)
      vec = xyz(:,iat)-center
      r2 = sum(vec**2)
      inert = inert+rad3v(iat)*((r2+tof*rad2)*unity &
         & -spread(vec,1,3)*spread(vec,2,3))
    end do

    cof(1,1) = inert(2,2)*inert(3,3)-inert(2,3)*inert(3,2)
    cof(2,2) = inert(1,1)*inert(3,3)-inert(1,3)*inert(3,1)
    cof(3,3) = inert(1,1)*inert(2,2)-inert(1,2)*inert(2,1)
    cof(1,2) = inert(1,3)*inert(3,2)-inert(1,2)*inert(3,3)
    cof(1,3) = inert(1,2)*inert(2,3)-inert(1,3)*inert(2,2)
    cof(2,3) = inert(1,3)*inert(2,1)-inert(1,1)*inert(2,3)
    cof(2,1) = cof(1,2)
    cof(3,1) = cof(1,3)
    cof(3,2) = cof(2,3)
    sdet = inert(1,1)*cof(1,1)+inert(1,2)*cof(2,1)+inert(1,3)*cof(3,1)
    if (abs(sdet) .lt. 1.0e-30_wp) return
    iinv = cof/sdet

    trc = cof(1,1)+cof(2,2)+cof(3,3)
    tmat = trc*unity-cof

    !>-- E = kk S^(-1/6), c1 = dE/dS, c2 = d2E/dS2
    kk = 0.5_wp*gbsa%keps*gbsa%alpbet*qtot*qtot*sqrt(tof*totRad3)
    c1 = -kk/6.0_wp*sdet**(-7.0_wp/6.0_wp)
    c2 = 7.0_wp*kk/36.0_wp*sdet**(-13.0_wp/6.0_wp)

    !>-- gvec = dS/dx; barycentre motion cancels since sum_a rad3_a v_a = 0
    do iat = 1,n
      vec = xyz(:,iat)-center
      ia = 3*(iat-1)
      gvec(ia+1:ia+3) = 2.0_wp*rad3v(iat)*matmul(tmat,vec)
    end do

    !>-- both rank-one dS (x) dS pieces, from the chain rule and from d2S
    fac = c2+c1/sdet
    do q = 1,ndof
      do p = 1,ndof
        hess(p,q) = hess(p,q)+fac*gvec(p)*gvec(q)
      end do
    end do

    !>-- C : d2I, block (b,c) = 2 rad3_b (delta_bc - rad3_c/totRad3) T
    do jat = 1,n
      ja = 3*(jat-1)
      do iat = 1,n
        ia = 3*(iat-1)
        pref = -2.0_wp*rad3v(iat)*rad3v(jat)/totRad3
        if (iat .eq. jat) pref = pref+2.0_wp*rad3v(iat)
        hess(ia+1:ia+3,ja+1:ja+3) = hess(ia+1:ia+3,ja+1:ja+3) &
           & +c1*pref*tmat
      end do
    end do

    !>-- -S tr(I^-1 dI I^-1 dI) as a nine-component contraction, one GEMM
    allocate (pm(ndof,9),qm(ndof,9))
    do iat = 1,n
      vec = xyz(:,iat)-center
      ia = 3*(iat-1)
      do p = 1,3
        mm = 2.0_wp*vec(p)*unity
        do k = 1,3
          mm(p,k) = mm(p,k)-vec(k)
          mm(k,p) = mm(k,p)-vec(k)
        end do
        amat3 = matmul(iinv,mm)*rad3v(iat)
        m = 0
        do q = 1,3
          do k = 1,3
            m = m+1
            qm(ia+p,m) = amat3(k,q)
            pm(ia+p,m) = amat3(q,k)
          end do
        end do
      end do
    end do
    call gemm(pm,qm,hess,transb='T',alpha=-c1*sdet,beta=1.0_wp)

    !>-- rho: (Bq)_i picks up keps*beta*Q/A_det, the same for every atom
    pref = -gbsa%keps*gbsa%alpbet*qtot*sqrt(tof*totRad3)/6.0_wp &
          & *sdet**(-7.0_wp/6.0_wp)
    do iat = 1,n
      rmat(iat,1:ndof) = rmat(iat,1:ndof)+pref*gvec(1:ndof)
    end do

  end subroutine adet_shape_hessian

  pure subroutine born_radius_derivs(psi,svdwi,vdwri,c1,fp,fpp)
    !***********************************************************************
    !* db/dpsi (fp) and d2b/dpsi2 (fpp) of the GBOBCII-rescaled Born radius
    !*   b = c1/(1/svdw - tanh(A)/vdw), A = alp u - bet u^2 + gam u^3,
    !* u = psi*svdw/2, psi the descreening integral, svdw the offset vdW radius.
    !***********************************************************************
    real(wp),intent(in) :: psi,svdwi,vdwri,c1
    real(wp),intent(out) :: fp,fpp

    real(wp) :: s1,v1,s2,u,a1,a2,argv,th,ch,den,ch2den2

    s1 = 1.0_wp/svdwi
    v1 = 1.0_wp/vdwri
    s2 = 0.5_wp*svdwi

    u = psi*s2
    argv = u*(alp+u*(gam*u-bet))
    a1 = alp-2.0_wp*bet*u+3.0_wp*gam*u*u
    a2 = -2.0_wp*bet+6.0_wp*gam*u

    th = tanh(argv)
    ch = cosh(argv)
    den = s1-v1*th
    ch2den2 = ch*ch*den*den

    fp = c1*s2*v1*a1/ch2den2
    fpp = c1*s2*s2*v1/ch2den2 &
         & *(a2-2.0_wp*a1*a1*th+2.0_wp*a1*a1*v1/(ch*ch*den))

  end subroutine born_radius_derivs

  pure subroutine psi_side_derivs(r,rhoo,rvdws,g1,g2,have)
    !***********************************************************************
    !* dg/dr (g1) and d2g/dr2 (g2) of the descreening an atom of vdW radius
    !* rvdws receives from a neighbour of descreening radius rhoo, branches as
    !* in compute_psi. have is false if that sphere lies entirely inside rvdws.
    !***********************************************************************
    real(wp),intent(in) :: r,rhoo,rvdws
    real(wp),intent(out) :: g1,g2
    logical,intent(out) :: have

    real(wp) :: ap,am,ab,rhab,lnab,pp,p1,p2,r2i

    g1 = 0.0_wp
    g2 = 0.0_wp
    have = .false.
    ap = r+rhoo
    am = r-rhoo

    if (r .lt. rvdws+rhoo) then
      if (ap .le. rvdws) return
      pp = am/(2.0_wp*ap)-(r*r-rhoo*rhoo)/(2.0_wp*rvdws*rvdws)-log(ap/rvdws)
      p1 = rhoo/(ap*ap)-r/(rvdws*rvdws)-1.0_wp/ap
      p2 = -2.0_wp*rhoo/(ap**3)-1.0_wp/(rvdws*rvdws)+1.0_wp/(ap*ap)
      r2i = 1.0_wp/(r*r)
      g1 = 1.0_wp/(ap*ap)+p1/(2.0_wp*r)-0.5_wp*pp*r2i
      g2 = -2.0_wp/(ap**3)+p2/(2.0_wp*r)-p1*r2i+pp/(r**3)
    else
      ab = ap*am
      rhab = rhoo/ab
      lnab = 0.5_wp*log(am/ap)/r
      g1 = -2.0_wp*r*rhab/ab+(rhab-lnab)/r
      g2 = -4.0_wp*rhab/ab+8.0_wp*r*r*rhab/(ab*ab)-2.0_wp*(rhab-lnab)/(r*r)
    end if
    have = .true.

  end subroutine psi_side_derivs

  pure subroutine p16_kernel_derivs(r,bi,bj,keps,v,d1,d2)
    !***********************************************************************
    !* P16 Born interaction kernel V = keps/f, f = r + t (t/(t + c r))^16,
    !* t = sqrt(b_i b_j), c = zetaP16/16. d1 and d2 (symmetric) hold the first
    !* and second derivatives in the order (r, b_i, b_j).
    !***********************************************************************
    real(wp),intent(in) :: r,bi,bj,keps
    real(wp),intent(out) :: v,d1(3),d2(3,3)

    real(wp) :: t,dd,ww,w15,w16,w17,f,f1,fr,ft,frr,frt,ftt
    real(wp) :: fi2,fi3,vr,vt,vrr,vrt,vtt
    real(wp) :: ti,tj,tii,tjj,tij

    t = sqrt(bi*bj)
    dd = t+zetaP16o16*r
    ww = t/dd
    w15 = ww**15
    w16 = w15*ww
    w17 = w16*ww

    f = r+t*w16
    f1 = 1.0_wp/f
    fi2 = f1*f1
    fi3 = fi2*f1

    fr = 1.0_wp-zetaP16*w17
    ft = w16*(17.0_wp-16.0_wp*ww)
    frr = 272.0_wp*zetaP16o16*zetaP16o16*w17/dd
    frt = -272.0_wp*zetaP16o16*w16*(1.0_wp-ww)/dd
    ftt = 272.0_wp*w15*(1.0_wp-ww)*(1.0_wp-ww)/dd

    vr = -keps*fr*fi2
    vt = -keps*ft*fi2
    vrr = keps*(-frr*fi2+2.0_wp*fr*fr*fi3)
    vrt = keps*(-frt*fi2+2.0_wp*fr*ft*fi3)
    vtt = keps*(-ftt*fi2+2.0_wp*ft*ft*fi3)

    ti = 0.5_wp*t/bi
    tj = 0.5_wp*t/bj
    tii = -0.25_wp*t/(bi*bi)
    tjj = -0.25_wp*t/(bj*bj)
    tij = 0.25_wp*t/(bi*bj)

    v = keps*f1
    d1(1) = vr
    d1(2) = vt*ti
    d1(3) = vt*tj

    d2(1,1) = vrr
    d2(1,2) = vrt*ti
    d2(1,3) = vrt*tj
    d2(2,1) = d2(1,2)
    d2(3,1) = d2(1,3)
    d2(2,2) = vtt*ti*ti+vt*tii
    d2(3,3) = vtt*tj*tj+vt*tjj
    d2(2,3) = vtt*ti*tj+vt*tij
    d2(3,2) = d2(2,3)

  end subroutine p16_kernel_derivs

  subroutine born_weighted_hessian(n,xyz,gbsa,w,kmat,hess)
    !***********************************************************************
    !* Add sum_i w_i d2 b_i / dR dR to hess, w_i = dE/db_i. With b_i = F(psi_i),
    !*   d2 b_i = F'_i d2 psi_i + F''_i/F'_i^2 (d b_i)(x)(d b_i).
    !* Only the pair part enters hess here. The rank-one coefficient is added to
    !* diag(kmat), which the caller contracts with the Born radius gradients.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(TBorn),intent(in) :: gbsa
    real(wp),intent(in) :: w(n)
    real(wp),intent(inout) :: kmat(n,n)
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: iat,jat
    real(wp) :: r,lrcut2,dr2,vec(3),ehat(3),blk(3,3)
    real(wp) :: gi1,gi2,gj1,gj2,wgt1,wgt2
    logical :: havei,havej
    real(wp),allocatable :: fp(:),fpp(:),wf(:)

    allocate (fp(n),fpp(n),wf(n))

    do iat = 1,n
      call born_radius_derivs(gbsa%psi(iat),gbsa%svdw(iat),gbsa%vdwr(iat), &
         & gbsa%bornScale,fp(iat),fpp(iat))
      wf(iat) = w(iat)*fp(iat)
      if (abs(fp(iat)) > 1.0e-300_wp) then
        kmat(iat,iat) = kmat(iat,iat)+w(iat)*fpp(iat)/(fp(iat)*fp(iat))
      end if
    end do

    lrcut2 = gbsa%lrcut*gbsa%lrcut

    !>-- pair part, sum_i (w_i F'_i) d2 psi_i; thread iat writes only its rows
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp shared(n,xyz,gbsa,wf,hess,lrcut2) &
    !$omp private(iat,jat,r,dr2,vec,ehat,blk,gi1,gi2,gj1,gj2,havei,havej, &
    !$omp&        wgt1,wgt2)
    do iat = 1,n
      do jat = 1,n
        if (jat .eq. iat) cycle
        vec = xyz(:,iat)-xyz(:,jat)
        dr2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
        if (dr2 .ge. lrcut2.or.dr2 .lt. 1.0e-12_wp) cycle
        r = sqrt(dr2)
        ehat = vec/r

        !>-- what iat receives from jat, and what jat receives from iat
        call psi_side_derivs(r,gbsa%rho(jat),gbsa%vdwr(iat),gi1,gi2,havei)
        call psi_side_derivs(r,gbsa%rho(iat),gbsa%vdwr(jat),gj1,gj2,havej)

        wgt1 = 0.0_wp
        wgt2 = 0.0_wp
        if (havei) then
          wgt1 = wgt1+wf(iat)*gi1
          wgt2 = wgt2+wf(iat)*gi2
        end if
        if (havej) then
          wgt1 = wgt1+wf(jat)*gj1
          wgt2 = wgt2+wf(jat)*gj2
        end if
        if (wgt1 .eq. 0.0_wp.and.wgt2 .eq. 0.0_wp) cycle

        call pair_hess_block(wgt1,wgt2,r,ehat,blk)
        call scatter_pair_row(hess,iat,jat,blk)
      end do
    end do
    !$omp end parallel do

  end subroutine born_weighted_hessian

  pure subroutine sasa_point_derivs(nno,nnlist,trj2,vdwsa,xyz,xyzp, &
        & sasap,nni,idx,ehat,dist,lp,lpp)
    !***********************************************************************
    !* One Lebedev point of an atom's surface integral, mirroring compute_w_sp.
    !* Returns the switch product sasap (0 if buried) and, for the nni neighbours
    !* inside the smoothing window: idx, unit vector ehat and distance dist from
    !* neighbour to grid point, lp = d(log s)/d(dist), lpp = d2(log s)/d(dist)2.
    !***********************************************************************
    integer,intent(in) :: nno
    integer,intent(in) :: nnlist(nno)
    real(wp),intent(in) :: trj2(:,:),vdwsa(:),xyz(:,:),xyzp(3)
    real(wp),intent(out) :: sasap
    integer,intent(out) :: nni
    integer,intent(out) :: idx(nno)
    real(wp),intent(out) :: ehat(3,nno),dist(nno),lp(nno),lpp(nno)

    integer :: i,ia
    real(wp) :: tj(3),tj2,sqtj,uj,ah3uj2,sasaij,ds,dss,li

    nni = 0
    sasap = 1.0_wp
    do i = 1,nno
      ia = nnlist(i)
      tj(:) = xyzp(:)-xyz(:,ia)
      tj2 = tj(1)*tj(1)+tj(2)*tj(2)+tj(3)*tj(3)
      if (tj2 .lt. trj2(2,ia)) then
        if (tj2 .le. trj2(1,ia)) then
          sasap = 0.0_wp
          nni = 0
          return
        else
          sqtj = sqrt(tj2)
          uj = sqtj-vdwsa(ia)
          ah3uj2 = ah3*uj*uj
          sasaij = ah0+(ah1+ah3uj2)*uj
          ds = ah1+3.0_wp*ah3uj2
          dss = 6.0_wp*ah3*uj

          sasap = sasap*sasaij
          li = ds/sasaij

          nni = nni+1
          idx(nni) = ia
          dist(nni) = sqtj
          ehat(:,nni) = tj(:)/sqtj
          lp(nni) = li
          lpp(nni) = dss/sasaij-li*li
        end if
      end if
    end do

  end subroutine sasa_point_derivs

  subroutine sasa_weighted_hessian(n,xyz,gbsa,w,hess)
    !***********************************************************************
    !* Add sum_i w_i d2 s_i / dR dR to hess, w_i = dE/ds_i. Per Lebedev point
    !* P = prod_j s_j and, with L_j = log s_j,
    !*   d2 P = P [ (sum_j dL_j) (x) (sum_j dL_j) + sum_j d2 L_j ].
    !* Only the weighted sum is accumulated; no per-atom 3N x 3N object is formed.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(TBorn),intent(in) :: gbsa
    real(wp),intent(in) :: w(n)
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: iat,ip,jj,kk,nno,nni,nloc,a,b,ja,ka,p,q
    integer :: maxnn
    real(wp) :: rsas,wr,sasap,wsa,xyzp(3),blkm(3,3)
    real(wp) :: gsum(3)
    integer,allocatable :: idx(:),slotof(:),loclist(:),slot(:)
    real(wp),allocatable :: ehat(:,:),dist(:),lp(:),lpp(:)
    real(wp),allocatable :: blk(:,:),gam(:,:)

    maxnn = maxval(gbsa%nnsas)
    if (maxnn .le. 0) return

    !$omp parallel default(none) &
    !$omp shared(n,xyz,gbsa,w,hess,maxnn) &
    !$omp private(iat,ip,jj,kk,nno,nni,nloc,a,b,ja,ka,p,q,rsas,wr,sasap,wsa, &
    !$omp&        xyzp,blkm,gsum,idx,slotof,loclist,slot,ehat,dist,lp,lpp, &
    !$omp&        blk,gam)
    allocate (idx(maxnn),ehat(3,maxnn),dist(maxnn),lp(maxnn),lpp(maxnn))
    allocate (slot(maxnn))
    allocate (slotof(n),source=0)
    allocate (loclist(maxnn+1))
    allocate (blk(3*(maxnn+1),3*(maxnn+1)))
    allocate (gam(3,maxnn+1))

    !$omp do schedule(dynamic)
    do iat = 1,n

      if (w(iat) .eq. 0.0_wp) cycle
      nno = gbsa%nnsas(iat)
      rsas = gbsa%vdwsa(iat)
      wr = gbsa%wrp(iat)*w(iat)

      !>-- slot 1 is the atom itself; neighbours are added as they appear
      nloc = 1
      loclist(1) = iat
      slotof(iat) = 1
      blk(1:3,1:3) = 0.0_wp

      do ip = 1,size(gbsa%angGrid,2)
        xyzp(:) = xyz(:,iat)+rsas*gbsa%angGrid(1:3,ip)
        call sasa_point_derivs(nno,gbsa%nnlists(:nno,iat),gbsa%trj2, &
           & gbsa%vdwsa,xyz,xyzp,sasap,nni,idx,ehat,dist,lp,lpp)
        if (sasap .le. tolsesp) cycle
        wsa = gbsa%angWeight(ip)*wr*sasap

        do jj = 1,nni
          if (slotof(idx(jj)) .eq. 0) then
            nloc = nloc+1
            loclist(nloc) = idx(jj)
            slotof(idx(jj)) = nloc
            a = 3*(nloc-1)
            blk(1:3*nloc,a+1:a+3) = 0.0_wp
            blk(a+1:a+3,1:3*nloc) = 0.0_wp
          end if
          slot(jj) = slotof(idx(jj))
        end do

        !>-- gradient of log P over the block, and its outer product
        gsum = 0.0_wp
        gam(:,1:nloc) = 0.0_wp
        do jj = 1,nni
          gsum = gsum+lp(jj)*ehat(:,jj)
          gam(:,slot(jj)) = gam(:,slot(jj))-lp(jj)*ehat(:,jj)
        end do
        gam(:,1) = gsum

        do q = 1,nloc
          ka = 3*(q-1)
          do p = 1,nloc
            ja = 3*(p-1)
            do b = 1,3
              do a = 1,3
                blk(ja+a,ka+b) = blk(ja+a,ka+b)+wsa*gam(a,p)*gam(b,q)
              end do
            end do
          end do
        end do

        !>-- second derivatives of log s, one pair block per neighbour
        do jj = 1,nni
          call pair_hess_block(lp(jj),lpp(jj),dist(jj),ehat(:,jj),blkm)
          blkm = wsa*blkm
          ka = 3*(slot(jj)-1)
          blk(1:3,1:3) = blk(1:3,1:3)+blkm
          blk(1:3,ka+1:ka+3) = blk(1:3,ka+1:ka+3)-blkm
          blk(ka+1:ka+3,1:3) = blk(ka+1:ka+3,1:3)-blkm
          blk(ka+1:ka+3,ka+1:ka+3) = blk(ka+1:ka+3,ka+1:ka+3)+blkm
        end do

      end do

      !>-- one scatter per atom, after the whole quadrature
      do q = 1,nloc
        ka = 3*(loclist(q)-1)
        do p = 1,nloc
          ja = 3*(loclist(p)-1)
          do b = 1,3
            do a = 1,3
              !$omp atomic update
              hess(ja+a,ka+b) = hess(ja+a,ka+b)+blk(3*(p-1)+a,3*(q-1)+b)
            end do
          end do
        end do
      end do

      do q = 1,nloc
        slotof(loclist(q)) = 0
      end do

    end do
    !$omp end do
    deallocate (idx,ehat,dist,lp,lpp,slot,slotof,loclist,blk,gam)
    !$omp end parallel

  end subroutine sasa_weighted_hessian

end module gfnff_hess_solv
