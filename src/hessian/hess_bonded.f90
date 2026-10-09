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
!> Closed-form Cartesian Hessians of the bonded terms in gfnff_eg_bonded.
!> Each term is E = f(s_1,...,s_k) in the gfnff_hess_prim primitives, so that
!>   H = sum_ab f_ab grad(s_a) (x) grad(s_b) + sum_a f_a hess(s_a)
module gfnff_hess_bonded

  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData,TGFFTopology
  use gfnff_neighbor,only:TNeigh
  use gfnff_hess_prim,only:prim_cos_angle,rsq_in_block,psi_from_cos, &
    &                      scatter_block,prim_cos_dihedral,prim_sin_oop, &
    &                      cheby_t
  use gfnff_hess_pair,only:pair_hess_block,scatter_pair_hessian
  use gfnff_hess_cn,only:logcn_weighted_hessian,hbcn_weighted_hessian
  use gfnff_rab,only:gfnffdrab
  use gfnff_bond_potential,only:bond_potential,bond_potential_hb
  use gfnff_math_wrapper,only:gemm
  implicit none
  private

  public :: hess_angles,hess_batm,hess_bonds
  public :: hess_torsions,hess_storsions,hess_torsions_available

  real(wp),parameter :: pi = 3.1415926535897932385_wp
  !> half the barrier of the special torsion around a triple-bonded carbon,
  !> the value sTors_eg uses (DLPNO-CCSD(T)/CBS on diphenylacetylene)
  real(wp),parameter :: sTors_erefhalf = 3.75e-4_wp

contains  !> MODULE PROCEDURES START HERE

  pure subroutine damp_derivs(rcut,x,d0,d1,d2)
    !***********************************************************************
    !* Bending and torsion damping D = 1/(1+(x/rcut)^2) and its true first two
    !* derivatives in x = r^2. gfnffdampa/gfnffdampt return 2*D'(x) instead.
    !***********************************************************************
    real(wp),intent(in) :: rcut,x
    real(wp),intent(out) :: d0,d1,d2

    real(wp) :: u,up,upp,den

    u = (x/rcut)**2
    den = 1.0_wp+u
    up = 2.0_wp*u/x
    upp = 2.0_wp*u/(x*x)
    d0 = 1.0_wp/den
    d1 = -up/(den*den)
    d2 = -upp/(den*den)+2.0_wp*up*up/(den*den*den)

  end subroutine damp_derivs

  subroutine hess_angles(n,at,xyz,param,topo,neigh,hess)
    !***********************************************************************
    !* Angle bending Hessian (cf. eg_angles/egbend), added to hess(3n,3n):
    !*   E = ea(c) D(x1) D(x2), c = cos(theta), x = squared vertex distances
    !* Vertex is alist(1,m), cf. eg_angles; the sum rule cannot detect a swap.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: m,iv,ja,ka,iTrj,iTrk,idx(3),a,b
    real(wp) :: va(3),vb(3),c0,kijk,cosv,c0cos
    real(wp) :: x1,x2,rcut1,rcut2,d10,d11,d12,d20,d21,d22
    real(wp) :: ea,eac,eacc,eps,gg,gg1,gg2,psi,ps1,ps2
    real(wp) :: fc,f1,f2,fcc,fc1,fc2,f11,f12,f22
    real(wp) :: gc(9),g1(9),g2(9)
    real(wp) :: hc(9,9),h1(9,9),h2(9,9),blk(9,9)

    if (topo%nangl .le. 0) return

    do m = 1,topo%nangl
      iv = topo%alist(1,m)
      ja = topo%alist(2,m)
      ka = topo%alist(3,m)
      iTrj = topo%alist(4,m)
      iTrk = topo%alist(5,m)
      c0 = topo%vangl(1,m)
      kijk = topo%vangl(2,m)

      va = xyz(:,ja)+neigh%transVec(:,iTrj)-xyz(:,iv)
      vb = xyz(:,ka)+neigh%transVec(:,iTrk)-xyz(:,iv)

      call prim_cos_angle(va,vb,cosv,gc,hc)
      call rsq_in_block(3,2,1,va,x1,g1,h1)
      call rsq_in_block(3,3,1,vb,x2,g2,h2)

      cosv = min(1.0_wp,max(-1.0_wp,cosv))

      rcut1 = param%atcuta*(param%rcov(at(ja))+param%rcov(at(iv)))**2
      rcut2 = param%atcuta*(param%rcov(at(ka))+param%rcov(at(iv)))**2
      call damp_derivs(rcut1,x1,d10,d11,d12)
      call damp_derivs(rcut2,x2,d20,d21,d22)

      if (pi-c0 .lt. 1.0e-6_wp) then
        !>-- near-linear: ea = k (psi - eps)^2, psi = pi - theta, eps = pi - c0.
        !>   The eps terms are singular at psi = 0; eps = 0 in all known cases.
        eps = pi-c0
        call psi_from_cos(cosv,gg,gg1,gg2,psi,ps1,ps2)
        ea = kijk*(gg-2.0_wp*eps*psi+eps*eps)
        eac = kijk*(gg1-2.0_wp*eps*ps1)
        eacc = kijk*(gg2-2.0_wp*eps*ps2)
      else
        c0cos = cos(c0)
        ea = kijk*(cosv-c0cos)**2
        eac = 2.0_wp*kijk*(cosv-c0cos)
        eacc = 2.0_wp*kijk
      end if

      fc = eac*d10*d20
      f1 = ea*d11*d20
      f2 = ea*d10*d21
      fcc = eacc*d10*d20
      fc1 = eac*d11*d20
      fc2 = eac*d10*d21
      f11 = ea*d12*d20
      f12 = ea*d11*d21
      f22 = ea*d10*d22

      do b = 1,9
        do a = 1,9
          blk(a,b) = fc*hc(a,b)+f1*h1(a,b)+f2*h2(a,b) &
             & +fcc*gc(a)*gc(b)+f11*g1(a)*g1(b)+f22*g2(a)*g2(b) &
             & +fc1*(gc(a)*g1(b)+g1(a)*gc(b)) &
             & +fc2*(gc(a)*g2(b)+g2(a)*gc(b)) &
             & +f12*(g1(a)*g2(b)+g2(a)*g1(b))
        end do
      end do

      idx = [iv,ja,ka]
      call scatter_block(hess,3,idx,blk)
    end do

  end subroutine hess_angles

  pure function hess_torsions_available(topo) result(ok)
    !***********************************************************************
    !* True if every proper torsion has sin(rn*phi0) = 0 (Chebyshev form exact).
    !***********************************************************************
    type(TGFFTopology),intent(in) :: topo
    logical :: ok
    integer :: m,rn

    ok = .true.
    do m = 1,topo%ntors
      rn = topo%tlist(5,m)
      if (rn .le. 0) cycle
      if (abs(sin(real(rn,wp)*topo%vtors(1,m))) .gt. 1.0e-8_wp) then
        ok = .false.
        return
      end if
    end do

  end function hess_torsions_available

  subroutine hess_torsions(n,at,xyz,param,topo,neigh,hess)
    !***********************************************************************
    !* Torsion Hessian (cf. egtors), added to hess(3n,3n). E = ea D1 D2 D3 with
    !*   rn > 0:  ea = V(1 - sigma T_rn(c)), c = cos(phi), sigma = cos(rn phi0),
    !*            polynomial in c, as phi has a branch point at planarity
    !*   rn == 0: ea = V(1 - cos(w - phi0)), out-of-plane angle w, s = sin(w)
    !*   rn < 0:  ea = V(cos(w) - cos(phi0))^2, double minimum at +-phi0
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: m,i,j,k,l,iTrl,iTrj,iTrk,rn,idx(4),a,b,p,q
    real(wp) :: vTrl(3),vTrj(3),vTrk(3),ra(3),rb(3),rc(3),re(3),rd(3),rv(3)
    real(wp) :: cv,phi0,vv,sigma,t0,t1,t2,ea,eac,eacc
    real(wp) :: co,so,rt,c0cos,s0sin
    real(wp) :: x(3),rcut,dv(3,3)
    real(wp) :: gc(12),gx(3,12),hc(12,12),hx(3,12,12),blk(12,12)
    real(wp) :: fd(4),fdd(4,4),dprod

    if (topo%ntors .le. 0) return

    do m = 1,topo%ntors
      i = topo%tlist(1,m)
      j = topo%tlist(2,m)
      k = topo%tlist(3,m)
      l = topo%tlist(4,m)
      rn = topo%tlist(5,m)
      iTrl = topo%tlist(6,m)
      iTrj = topo%tlist(7,m)
      iTrk = topo%tlist(8,m)
      if (iTrj .gt. neigh%nTrans.or.iTrk .gt. neigh%nTrans &
         & .or.iTrl .gt. neigh%nTrans) cycle
      vTrl = neigh%transVec(:,iTrl)
      vTrj = neigh%transVec(:,iTrj)
      vTrk = neigh%transVec(:,iTrk)
      phi0 = topo%vtors(1,m)
      vv = topo%vtors(2,m)

      if (rn .gt. 0) then
        !>-- proper torsion, bond vectors as in torsPBC(mo=1)
        ra = xyz(:,j)-(xyz(:,i)+vTrl)
        rb = (xyz(:,k)+vTrj)-xyz(:,j)
        rc = (xyz(:,l)+vTrk)-(xyz(:,k)+vTrj)
        call prim_cos_dihedral(ra,rb,rc,cv,gc,hc)
        cv = min(1.0_wp,max(-1.0_wp,cv))
        sigma = cos(real(rn,wp)*phi0)
        call cheby_t(rn,cv,t0,t1,t2)
        ea = vv*(1.0_wp-sigma*t0)
        eac = -vv*sigma*t1
        eacc = -vv*sigma*t2
        call rsq_in_block(4,1,2,-ra,x(1),gx(1,:),hx(1,:,:))
        call rsq_in_block(4,2,3,-rb,x(2),gx(2,:),hx(2,:,:))
        call rsq_in_block(4,3,4,-rc,x(3),gx(3,:),hx(3,:,:))
        rcut = param%atcutt*(param%rcov(at(i))+param%rcov(at(j)))**2
        call damp_derivs(rcut,x(1),dv(1,1),dv(1,2),dv(1,3))
        rcut = param%atcutt*(param%rcov(at(k))+param%rcov(at(j)))**2
        call damp_derivs(rcut,x(2),dv(2,1),dv(2,2),dv(2,3))
        rcut = param%atcutt*(param%rcov(at(k))+param%rcov(at(l)))**2
        call damp_derivs(rcut,x(3),dv(3,1),dv(3,2),dv(3,3))
      else
        !>-- out-of-plane, vectors as in omegaPBC
        re = xyz(:,i)-(xyz(:,j)+vTrj)
        rd = (xyz(:,k)+vTrk)-(xyz(:,j)+vTrj)
        rv = (xyz(:,l)+vTrl)-xyz(:,i)
        call prim_sin_oop(re,rd,rv,cv,gc,hc)
        cv = min(1.0_wp,max(-1.0_wp,cv))
        rt = sqrt(max(1.0_wp-cv*cv,1.0e-14_wp))
        co = rt
        so = -cv/rt                      !> d cos(omega) / ds
        rt = -1.0_wp/(rt*rt*rt)          !> d2 cos(omega) / ds2
        if (rn .eq. 0) then
          c0cos = cos(phi0)
          s0sin = sin(phi0)
          ea = vv*(1.0_wp-c0cos*co-s0sin*cv)
          eac = vv*(-c0cos*so-s0sin)
          eacc = -vv*c0cos*rt
        else
          c0cos = cos(phi0)
          ea = vv*(co-c0cos)**2
          eac = 2.0_wp*vv*(co-c0cos)*so
          eacc = 2.0_wp*vv*(so*so+(co-c0cos)*rt)
        end if
        call rsq_in_block(4,1,2,re,x(1),gx(1,:),hx(1,:,:))
        call rsq_in_block(4,3,2,rd,x(2),gx(2,:),hx(2,:,:))
        !>-- third damped pair is j-l: (R_j+vTrj) - (R_l+vTrl) = -(re+rv)
        call rsq_in_block(4,2,4,-(re+rv),x(3),gx(3,:),hx(3,:,:))
        rcut = param%atcutt*(param%rcov(at(i))+param%rcov(at(j)))**2
        call damp_derivs(rcut,x(1),dv(1,1),dv(1,2),dv(1,3))
        rcut = param%atcutt*(param%rcov(at(k))+param%rcov(at(j)))**2
        call damp_derivs(rcut,x(2),dv(2,1),dv(2,2),dv(2,3))
        rcut = param%atcutt*(param%rcov(at(j))+param%rcov(at(l)))**2
        call damp_derivs(rcut,x(3),dv(3,1),dv(3,2),dv(3,3))
      end if

      !>-- fd/fdd: partials of E in the primitives (c, x1, x2, x3)
      dprod = dv(1,1)*dv(2,1)*dv(3,1)
      fd(1) = eac*dprod
      fd(2) = ea*dv(1,2)*dv(2,1)*dv(3,1)
      fd(3) = ea*dv(1,1)*dv(2,2)*dv(3,1)
      fd(4) = ea*dv(1,1)*dv(2,1)*dv(3,2)
      fdd(1,1) = eacc*dprod
      fdd(1,2) = eac*dv(1,2)*dv(2,1)*dv(3,1)
      fdd(1,3) = eac*dv(1,1)*dv(2,2)*dv(3,1)
      fdd(1,4) = eac*dv(1,1)*dv(2,1)*dv(3,2)
      fdd(2,2) = ea*dv(1,3)*dv(2,1)*dv(3,1)
      fdd(3,3) = ea*dv(1,1)*dv(2,3)*dv(3,1)
      fdd(4,4) = ea*dv(1,1)*dv(2,1)*dv(3,3)
      fdd(2,3) = ea*dv(1,2)*dv(2,2)*dv(3,1)
      fdd(2,4) = ea*dv(1,2)*dv(2,1)*dv(3,2)
      fdd(3,4) = ea*dv(1,1)*dv(2,2)*dv(3,2)
      do q = 1,4
        do p = q+1,4
          fdd(p,q) = fdd(q,p)
        end do
      end do

      do b = 1,12
        do a = 1,12
          blk(a,b) = fd(1)*hc(a,b)+fdd(1,1)*gc(a)*gc(b)
        end do
      end do
      do p = 1,3
        do b = 1,12
          do a = 1,12
            blk(a,b) = blk(a,b)+fd(1+p)*hx(p,a,b) &
               & +fdd(1,1+p)*(gc(a)*gx(p,b)+gx(p,a)*gc(b))
          end do
        end do
      end do
      do q = 1,3
        do p = 1,3
          do b = 1,12
            do a = 1,12
              blk(a,b) = blk(a,b)+fdd(1+p,1+q)*gx(p,a)*gx(q,b)
            end do
          end do
        end do
      end do

      idx = [i,j,k,l]
      call scatter_block(hess,4,idx,blk)
    end do

  end subroutine hess_torsions

  subroutine hess_storsions(n,xyz,topo,hess)
    !***********************************************************************
    !* Hessian of the special torsion around a triple-bonded carbon (sTors_eg),
    !* added to hess(3n,3n); part of the gff_term_tors mask. Undamped, no phase:
    !*   E = e0 (1 - cos(2 phi)) = 2 e0 (1 - c^2),  c = cos(phi)
    !***********************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFTopology),intent(in) :: topo
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: m,a,b,c1,c2,c3,c4,idx(4)
    real(wp) :: ra(3),rb(3),rc(3),cv,ec,ecc
    real(wp) :: gc(12),hc(12,12),blk(12,12)

    if (.not.allocated(topo%sTorsl)) return
    if (topo%nstors .le. 0) return

    ecc = -4.0_wp*sTors_erefhalf
    do m = 1,topo%nstors
      !>-- same completeness guard as sTors_eg
      if (any(topo%sTorsl(:,m) .eq. 0)) cycle
      c1 = topo%sTorsl(1,m)
      c2 = topo%sTorsl(2,m)
      c3 = topo%sTorsl(5,m)
      c4 = topo%sTorsl(6,m)

      ra = xyz(:,c2)-xyz(:,c1)
      rb = xyz(:,c3)-xyz(:,c2)
      rc = xyz(:,c4)-xyz(:,c3)
      call prim_cos_dihedral(ra,rb,rc,cv,gc,hc)
      cv = min(1.0_wp,max(-1.0_wp,cv))
      ec = -4.0_wp*sTors_erefhalf*cv

      do b = 1,12
        do a = 1,12
          blk(a,b) = ec*hc(a,b)+ecc*gc(a)*gc(b)
        end do
      end do

      idx = [c1,c2,c3,c4]
      call scatter_block(hess,4,idx,blk)
    end do

  end subroutine hess_storsions

  subroutine hess_bonds(n,at,xyz,srab,cnthr,cn,dcn,hb_cn,hb_dcn, &
        & param,topo,neigh,version,hess)
    !***********************************************************************
    !* Bond stretch Hessian (cf. egbond), added to hess(3n,3n), for molecular
    !* systems only. With d = r - r0 and r0 linear in the logCN (gfnffdrab),
    !*   H = sum_b E_d hess(r_b) + sum_a w_a hess(logCN_a)
    !*     + sum_b E_dd grad(d_b) (x) grad(d_b), done as one GEMM,
    !* where w_a = -sum_b E_d dr0_b/dlogCN_a is the gradient code's dEdcn.
    !* Bonds with nr_hb >= 1 scale alpha by (1 - t1*hbcn_H), a second primitive.
    !*   srab   - packed distances;  cnthr - squared CN cutoff
    !*   dcn    - dcn(:,m,a) = d logCN_a / d R_m, cn being logCN
    !*   hb_dcn - opposite index order, hb_dcn(:,a,m) = d hbcn_a / d R_m
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: srab(n*(n+1)/2)
    real(wp),intent(in) :: cnthr
    real(wp),intent(in) :: cn(n)
    real(wp),intent(in) :: dcn(3*n,n)
    real(wp),intent(in) :: hb_cn(n)
    real(wp),intent(in) :: hb_dcn(3,n,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    integer,intent(in) :: version
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: b,iat,jat,iTr,nb,ndof,a,m,hbH,ncol,cd,cn2
    real(wp) :: vec(3),ehat(3),r,d,amp,ee,ed,edd,k1,k2
    real(wp) :: t1,al0,en,enn,edn
    real(wp) :: blk(3,3)
    real(wp),allocatable :: rab0(:),gfac(:,:),rabdcn(:,:)
    real(wp),allocatable :: w(:),wh(:),gmat(:,:),bmat(:,:)
    real(wp),allocatable :: eddv(:),ennv(:),ednv(:)
    logical :: ishb

    nb = neigh%nbond
    if (nb .le. 0) return
    ndof = 3*n

    allocate (rab0(nb),gfac(3,nb),rabdcn(2,nb))
    rab0(:) = neigh%vbond(1,:)
    call gfnffdrab(n,at,cn,nb,neigh%blist,rab0,gfac,rabdcn)

    t1 = 1.0_wp-param%vbond_scale

    !>-- gmat: two columns per bond, grad d and grad n (zero for a plain bond)
    ncol = 2*nb
    allocate (w(n),source=0.0_wp)
    allocate (wh(n),source=0.0_wp)
    allocate (gmat(ndof,ncol),source=0.0_wp)
    allocate (eddv(nb),ennv(nb),ednv(nb),source=0.0_wp)

    do b = 1,nb
      jat = neigh%blist(1,b)
      iat = neigh%blist(2,b)
      iTr = neigh%blist(3,b)
      if (iTr .gt. neigh%nTrans) cycle

      vec = xyz(:,iat)-xyz(:,jat)-neigh%transVec(:,iTr)
      r = sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))
      ehat = vec/r
      al0 = neigh%vbond(2,b)
      amp = neigh%vbond(3,b)
      d = r-rab0(b)

      ishb = neigh%nr_hb(b) .ge. 1
      hbH = 0
      if (ishb) then
        if (at(iat) .eq. 1) then
          hbH = iat
        else if (at(jat) .eq. 1) then
          hbH = jat
        else
          ishb = .false.
        end if
      end if

      if (ishb) then
        call bond_potential_hb(version,al0,t1,hb_cn(hbH),amp,d,ee,ed,en, &
           & edd=edd,enn=enn,edn=edn)
      else
        call bond_potential(version,al0,amp,d,ee,ed,edd=edd)
      end if
      k1 = rabdcn(1,b)          !> dr0 / dlogCN_iat
      k2 = rabdcn(2,b)          !> dr0 / dlogCN_jat

      call pair_hess_block(ed,0.0_wp,r,ehat,blk)
      call scatter_pair_hessian(hess,iat,jat,blk)

      w(iat) = w(iat)-ed*k1
      w(jat) = w(jat)-ed*k2

      cd = 2*b-1
      cn2 = 2*b
      gmat(:,cd) = -k1*dcn(:,iat)-k2*dcn(:,jat)
      do a = 1,3
        gmat(3*(iat-1)+a,cd) = gmat(3*(iat-1)+a,cd)+ehat(a)
        gmat(3*(jat-1)+a,cd) = gmat(3*(jat-1)+a,cd)-ehat(a)
      end do
      eddv(b) = edd

      if (ishb) then
        wh(hbH) = wh(hbH)+en
        do m = 1,n
          gmat(3*(m-1)+1:3*(m-1)+3,cn2) = hb_dcn(:,hbH,m)
        end do
        ennv(b) = enn
        ednv(b) = edn
      end if
    end do

    call logcn_weighted_hessian(n,at,xyz,srab,cnthr,param,dcn,w,hess)
    if (any(wh .ne. 0.0_wp)) then
      call hbcn_weighted_hessian(n,at,xyz,param,topo,neigh,wh,hess)
    end if

    !>-- bmat = gmat times the per-bond 2x2 curvature; one GEMM gives the sum
    allocate (bmat(ndof,ncol))
    do b = 1,nb
      cd = 2*b-1
      cn2 = 2*b
      bmat(:,cd) = eddv(b)*gmat(:,cd)+ednv(b)*gmat(:,cn2)
      bmat(:,cn2) = ednv(b)*gmat(:,cd)+ennv(b)*gmat(:,cn2)
    end do
    call gemm(bmat,gmat,hess,transb='T',alpha=1.0_wp,beta=1.0_wp)

  end subroutine hess_bonds

  subroutine hess_batm(n,at,xyz,param,topo,neigh,hess)
    !***********************************************************************
    !* Bonded Axilrod-Teller-Muto Hessian (cf. eg_batm/batmgfnff_eg), added to
    !* hess(3n,3n). In the squared sides x = r_ij^2, y = r_jk^2, z = r_ik^2,
    !*   E = c9 [ 3/8 P S^(-5/2) + S^(-3/2) ],  S = xyz,
    !*   P = (x+y-z)(x-y+z)(-x+y+z). c9 is constant: topo%qa is frozen at setup.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(inout) :: hess(3*n,3*n)

    real(wp),parameter :: fqq = 3.0_wp
    integer :: i,iat,jat,kat,iTrj,iTrk,idx(3),a,b,p,q
    real(wp) :: fi,fj,fk,c9,x,y,z,s,pp
    real(wp) :: pd(3),pdd(3,3),sd(3),sdd(3,3)
    real(wp) :: um(3),umm(3,3),vm(3),vmm(3,3),fd(3),fdd(3,3)
    real(wp) :: sm25,sm35,sm45
    real(wp) :: gp(3,9),hp(3,9,9),blk(9,9)

    if (topo%nbatm .le. 0) return

    do i = 1,topo%nbatm
      iat = topo%b3list(1,i)
      jat = topo%b3list(2,i)
      kat = topo%b3list(3,i)
      iTrj = topo%b3list(4,i)
      iTrk = topo%b3list(5,i)

      fi = min(max(1.0_wp-fqq*topo%qa(iat),-4.0_wp),4.0_wp)
      fj = min(max(1.0_wp-fqq*topo%qa(jat),-4.0_wp),4.0_wp)
      fk = min(max(1.0_wp-fqq*topo%qa(kat),-4.0_wp),4.0_wp)
      c9 = fi*fj*fk*param%zb3atm(at(iat))*param%zb3atm(at(jat)) &
          & *param%zb3atm(at(kat))

      call rsq_in_block(3,1,2,xyz(:,iat)-xyz(:,jat)-neigh%transVec(:,iTrj), &
         & x,gp(1,:),hp(1,:,:))
      call rsq_in_block(3,2,3,xyz(:,jat)+neigh%transVec(:,iTrj) &
         & -xyz(:,kat)-neigh%transVec(:,iTrk),y,gp(2,:),hp(2,:,:))
      call rsq_in_block(3,1,3,xyz(:,iat)-xyz(:,kat)-neigh%transVec(:,iTrk), &
         & z,gp(3,:),hp(3,:,:))

      pp = -x**3-y**3-z**3+x*x*(y+z)+y*y*(x+z)+z*z*(x+y)-2.0_wp*x*y*z
      pd(1) = -3.0_wp*x*x+2.0_wp*x*(y+z)+y*y+z*z-2.0_wp*y*z
      pd(2) = -3.0_wp*y*y+2.0_wp*y*(x+z)+x*x+z*z-2.0_wp*x*z
      pd(3) = -3.0_wp*z*z+2.0_wp*z*(x+y)+x*x+y*y-2.0_wp*x*y
      pdd(1,1) = -6.0_wp*x+2.0_wp*(y+z)
      pdd(2,2) = -6.0_wp*y+2.0_wp*(x+z)
      pdd(3,3) = -6.0_wp*z+2.0_wp*(x+y)
      pdd(1,2) = 2.0_wp*(x+y-z)
      pdd(1,3) = 2.0_wp*(x+z-y)
      pdd(2,3) = 2.0_wp*(y+z-x)
      pdd(2,1) = pdd(1,2)
      pdd(3,1) = pdd(1,3)
      pdd(3,2) = pdd(2,3)

      s = x*y*z
      sd(1) = y*z
      sd(2) = x*z
      sd(3) = x*y
      sdd = 0.0_wp
      sdd(1,2) = z
      sdd(2,1) = z
      sdd(1,3) = y
      sdd(3,1) = y
      sdd(2,3) = x
      sdd(3,2) = x

      sm25 = s**(-2.5_wp)          !> S^-5/2, the factor multiplying P
      sm35 = sm25/s
      sm45 = sm35/s
      do p = 1,3
        um(p) = -2.5_wp*sm35*sd(p)          !> d(S^-5/2)
        vm(p) = -1.5_wp*sm25*sd(p)          !> d(S^-3/2)
      end do
      do q = 1,3
        do p = 1,3
          umm(p,q) = 8.75_wp*sm45*sd(p)*sd(q)-2.5_wp*sm35*sdd(p,q)
          vmm(p,q) = 3.75_wp*sm35*sd(p)*sd(q)-1.5_wp*sm25*sdd(p,q)
        end do
      end do

      do p = 1,3
        fd(p) = c9*(0.375_wp*(pd(p)*sm25+pp*um(p))+vm(p))
      end do
      do q = 1,3
        do p = 1,3
          fdd(p,q) = c9*(0.375_wp*(pdd(p,q)*sm25+pd(p)*um(q)+pd(q)*um(p) &
             & +pp*umm(p,q))+vmm(p,q))
        end do
      end do

      blk = 0.0_wp
      do q = 1,3
        do p = 1,3
          do b = 1,9
            do a = 1,9
              blk(a,b) = blk(a,b)+fdd(p,q)*gp(p,a)*gp(q,b)
            end do
          end do
        end do
      end do
      do p = 1,3
        do b = 1,9
          do a = 1,9
            blk(a,b) = blk(a,b)+fd(p)*hp(p,a,b)
          end do
        end do
      end do

      idx = [iat,jat,kat]
      call scatter_block(hess,3,idx,blk)
    end do

  end subroutine hess_batm

end module gfnff_hess_bonded
