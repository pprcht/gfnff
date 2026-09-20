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
!> Closed-form Hessian of the GFN-FF halogen and hydrogen bond terms. Apart from
!> the carbonyl bend and torsion every primitive is a distance. The charges are
!> topo%qa, frozen at setup, so there is no charge response term.
module gfnff_hess_hb

  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData,TGFFTopology,TGFFNeighbourList
  use gfnff_neighbor,only:TNeigh
  use gfnff_hess_prim,only:dist_in_block,combine_prims,scatter_block, &
    &                      prim_cos_angle,prim_cos_dihedral,embed_block
  use gfnff_topo_hbset,only:hbonds
  implicit none
  private

  public :: hess_xbonds,hess_hbonds_bound,hess_hbonds_unbound

  real(wp),parameter :: eps12 = 1.0e-12_wp
  real(wp),parameter :: pi = 3.1415926535897932385_wp

contains  !> MODULE PROCEDURES START HERE

  subroutine hess_xbonds(n,at,xyz,param,topo,neigh,nlist,hess)
    !***********************************************************************
    !* Halogen bond Hessian over hblist3, counterpart of eg_xbonds and
    !* rbxgfnff_eg. hess is the (3n,3n) Hessian and is incremented.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    type(TGFFNeighbourList),intent(in) :: nlist
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: i,ia,ib,ix,iTrB,iTrX,idx(3)
    logical :: ok
    real(wp) :: drax(3),drbx(3),drab(3)
    real(wp) :: s1,s2,s3,cnst,ex1,ex2,qb,qx
    real(wp) :: gg(3,9),hh(3,9,9),fd(3),fdd(3,3),blk(9,9)

    if (nlist%nxb .le. 0) return

    do i = 1,nlist%nxb
      ia = nlist%hblist3(1,i)     !> A, in the central cell
      ib = nlist%hblist3(2,i)     !> B, the donor
      ix = nlist%hblist3(3,i)     !> X, the halogen
      iTrB = nlist%hblist3(4,i)
      iTrX = nlist%hblist3(5,i)
      if (iTrB .gt. neigh%nTrans.or.iTrX .gt. neigh%nTrans) cycle
      if (ia .eq. 0.or.ib .eq. 0) cycle

      drax = xyz(:,ia)-(xyz(:,ix)+neigh%transVec(:,iTrX))
      drbx = (xyz(:,ib)+neigh%transVec(:,iTrB))-(xyz(:,ix)+neigh%transVec(:,iTrX))
      drab = xyz(:,ia)-(xyz(:,ib)+neigh%transVec(:,iTrB))

      ex1 = exp(param%xbst*topo%qa(ix))
      ex2 = ex1+param%xbsf
      qx = ex1/ex2
      ex1 = exp(-param%xbst*topo%qa(ib))
      ex2 = ex1+param%xbsf
      qb = ex1/ex2
      cnst = qb*param%xbaci(at(ix))*qx

      !>-- 1e-12 offsets as in rbxgfnff_eg, so the Hessian matches its gradient
      call dist_in_block(3,1,3,drax,s1,gg(1,:),hh(1,:,:))
      call dist_in_block(3,2,3,drbx,s2,gg(2,:),hh(2,:,:))
      call dist_in_block(3,1,2,drab,s3,gg(3,:),hh(3,:,:))
      s1 = s1+1.0e-12_wp
      s2 = s2+1.0e-12_wp
      if (s3 .lt. 1.0e-12_wp) cycle

      call xb_partials(param,at(ia),at(ib),s1,s2,s3,cnst,fd,fdd,ok)
      if (.not.ok) cycle

      call combine_prims(3,3,fd,fdd,gg,hh,blk)
      idx = [ia,ib,ix]
      call scatter_block(hess,3,idx,blk)
    end do

  end subroutine hess_xbonds

  pure subroutine xb_partials(param,ata,atb,s1,s2,s3,cnst,fd,fdd,ok)
    !***********************************************************************
    !* First and second partials of E = -cnst * R(s2) * O(s1,s2,s3) in s1 = r_AX,
    !* s2 = r_BX, s3 = r_AB. ok is false when the rbxgfnff_eg overflow guard trips.
    !***********************************************************************
    type(TGFFData),intent(in) :: param
    integer,intent(in) :: ata,atb
    real(wp),intent(in) :: s1,s2,s3,cnst
    real(wp),intent(out) :: fd(3),fdd(3,3)
    logical,intent(out) :: ok

    real(wp) :: aa,pp,shortcut,u1,u3,k1,k3,u1p,u1pp,u3p,u3pp
    real(wp) :: dl0,dl1,dl2,ds0,ds1,ds2,c0,c1,c2
    real(wp) :: r0,r1,r2,zz,ee,den,oo,oz,ozz
    real(wp) :: z1,z2,z3,z13,z23,z33
    real(wp) :: o1,o2,o3,o11,o12,o13,o22,o23,o33

    fd = 0.0_wp
    fdd = 0.0_wp
    ok = .false.

    aa = param%xbacut
    pp = param%hbalp

    zz = aa*((s1+s2)/s3-1.0_wp)
    if (zz .gt. 15.0_wp) return
    ok = .true.

    !>-- R = damp_s damp_l / s2^3; ratios are powers: v' = k v/s, v'' = k(k-1)v/s^2
    k1 = 2.0_wp*pp
    k3 = -2.0_wp*pp
    u1 = (s2*s2/param%hblongcut_xb)**pp
    shortcut = param%xbscut*(param%rad(ata)+param%rad(atb))
    u3 = (shortcut/(s2*s2))**pp
    u1p = k1*u1/s2
    u1pp = k1*(k1-1.0_wp)*u1/(s2*s2)
    u3p = k3*u3/s2
    u3pp = k3*(k3-1.0_wp)*u3/(s2*s2)

    den = 1.0_wp+u1
    dl0 = 1.0_wp/den
    dl1 = -u1p/(den*den)
    dl2 = -u1pp/(den*den)+2.0_wp*u1p*u1p/(den*den*den)

    den = 1.0_wp+u3
    ds0 = 1.0_wp/den
    ds1 = -u3p/(den*den)
    ds2 = -u3pp/(den*den)+2.0_wp*u3p*u3p/(den*den*den)

    c0 = 1.0_wp/(s2*s2*s2)
    c1 = -3.0_wp*c0/s2
    c2 = 12.0_wp*c0/(s2*s2)

    r0 = dl0*ds0*c0
    r1 = dl1*ds0*c0+dl0*ds1*c0+dl0*ds0*c1
    r2 = dl2*ds0*c0+dl0*ds2*c0+dl0*ds0*c2 &
       & +2.0_wp*(dl1*ds1*c0+dl1*ds0*c1+dl0*ds1*c1)

    ee = exp(zz)
    den = 1.0_wp+ee
    oo = 2.0_wp/den
    oz = -2.0_wp*ee/(den*den)
    ozz = 2.0_wp*ee*(ee-1.0_wp)/(den*den*den)

    z1 = aa/s3
    z2 = z1
    z3 = -aa*(s1+s2)/(s3*s3)
    z13 = -aa/(s3*s3)
    z23 = z13
    z33 = 2.0_wp*aa*(s1+s2)/(s3*s3*s3)

    o1 = oz*z1
    o2 = oz*z2
    o3 = oz*z3
    o11 = ozz*z1*z1
    o12 = ozz*z1*z2
    o13 = ozz*z1*z3+oz*z13
    o22 = ozz*z2*z2
    o23 = ozz*z2*z3+oz*z23
    o33 = ozz*z3*z3+oz*z33

    fd(1) = -cnst*r0*o1
    fd(2) = -cnst*(r1*oo+r0*o2)
    fd(3) = -cnst*r0*o3
    fdd(1,1) = -cnst*r0*o11
    fdd(1,2) = -cnst*(r1*o1+r0*o12)
    fdd(1,3) = -cnst*r0*o13
    fdd(2,2) = -cnst*(r2*oo+2.0_wp*r1*o2+r0*o22)
    fdd(2,3) = -cnst*(r1*o3+r0*o23)
    fdd(3,3) = -cnst*r0*o33
    fdd(2,1) = fdd(1,2)
    fdd(3,1) = fdd(1,3)
    fdd(3,2) = fdd(2,3)

  end subroutine xb_partials

  pure subroutine sfun_mul(ns,a0,a1,a2,b0,b1,b2)
    !***********************************************************************
    !* Product rule: multiply b into a (value, gradient, Hessian); a overwritten.
    !***********************************************************************
    integer,intent(in) :: ns
    real(wp),intent(inout) :: a0,a1(ns),a2(ns,ns)
    real(wp),intent(in) :: b0,b1(ns),b2(ns,ns)

    integer :: p,q

    do q = 1,ns
      do p = 1,ns
        a2(p,q) = a2(p,q)*b0+a1(p)*b1(q)+b1(p)*a1(q)+a0*b2(p,q)
      end do
    end do
    do p = 1,ns
      a1(p) = a1(p)*b0+a0*b1(p)
    end do
    a0 = a0*b0

  end subroutine sfun_mul

  pure subroutine outl_factor(ns,ip,iq,ir,aa,sp,sq,sr,flip,f0,f1,f2,ok)
    !***********************************************************************
    !* Out-of-line (non-collinearity) factor in z = aa ((s_p+s_q+1e-12)/s_r - 1):
    !*   flip = .false. :  f = 2/(1 + e^z)        (the A-H...B factor)
    !*   flip = .true.  :  f = 2/(1 + e^-z) - 1   (the A...nb(B)-B factor)
    !* ok is false if the overflow guard of the energy code trips.
    !***********************************************************************
    integer,intent(in) :: ns,ip,iq,ir
    real(wp),intent(in) :: aa,sp,sq,sr
    logical,intent(in) :: flip
    real(wp),intent(out) :: f0,f1(ns),f2(ns,ns)
    logical,intent(out) :: ok

    real(wp) :: sm,zz,ee,den,fz,fzz,z1,z3,z13,z33

    f0 = 0.0_wp
    f1 = 0.0_wp
    f2 = 0.0_wp
    ok = .false.

    sm = sp+sq+eps12
    zz = aa*(sm/sr-1.0_wp)
    !>-- the energy code guards only the unflipped variant; the other saturates
    if (.not.flip.and.zz .gt. 15.0_wp) return
    ok = .true.

    if (flip) then
      !>-- same shape in y = -z, so f_z = -f_y and f_zz = f_yy
      ee = exp(-zz)
      den = 1.0_wp+ee
      f0 = 2.0_wp/den-1.0_wp
      fz = 2.0_wp*ee/(den*den)
      fzz = 2.0_wp*ee*(ee-1.0_wp)/(den*den*den)
    else
      ee = exp(zz)
      den = 1.0_wp+ee
      f0 = 2.0_wp/den
      fz = -2.0_wp*ee/(den*den)
      fzz = 2.0_wp*ee*(ee-1.0_wp)/(den*den*den)
    end if

    z1 = aa/sr                       !> = z_p = z_q
    z3 = -aa*sm/(sr*sr)              !> = z_r
    z13 = -aa/(sr*sr)                !> = z_pr = z_qr
    z33 = 2.0_wp*aa*sm/(sr*sr*sr)    !> = z_rr

    f1(ip) = fz*z1
    f1(iq) = fz*z1
    f1(ir) = fz*z3

    f2(ip,ip) = fzz*z1*z1
    f2(iq,iq) = fzz*z1*z1
    f2(ip,iq) = fzz*z1*z1
    f2(iq,ip) = f2(ip,iq)
    f2(ip,ir) = fzz*z1*z3+fz*z13
    f2(ir,ip) = f2(ip,ir)
    f2(iq,ir) = fzz*z1*z3+fz*z13
    f2(ir,iq) = f2(iq,ir)
    f2(ir,ir) = fzz*z3*z3+fz*z33

  end subroutine outl_factor

  pure subroutine damp_rab(param,radab,s,d0,d1,d2)
    !***********************************************************************
    !* Long- times short-range hydrogen bond damping in the A-B distance s, with
    !* derivatives.
    !***********************************************************************
    type(TGFFData),intent(in) :: param
    real(wp),intent(in) :: radab,s
    real(wp),intent(out) :: d0,d1,d2

    real(wp) :: u1,u3,k1,k3,u1p,u1pp,u3p,u3pp,den
    real(wp) :: l0,l1,l2,s0,s1,s2

    k1 = 2.0_wp*param%hbalp
    k3 = -2.0_wp*param%hbalp
    u1 = (s*s/param%hblongcut)**param%hbalp
    u3 = (param%hbscut*radab/(s*s))**param%hbalp
    u1p = k1*u1/s
    u1pp = k1*(k1-1.0_wp)*u1/(s*s)
    u3p = k3*u3/s
    u3pp = k3*(k3-1.0_wp)*u3/(s*s)

    den = 1.0_wp+u1
    l0 = 1.0_wp/den
    l1 = -u1p/(den*den)
    l2 = -u1pp/(den*den)+2.0_wp*u1p*u1p/(den*den*den)

    den = 1.0_wp+u3
    s0 = 1.0_wp/den
    s1 = -u3p/(den*den)
    s2 = -u3pp/(den*den)+2.0_wp*u3p*u3p/(den*den*den)

    d0 = l0*s0
    d1 = l1*s0+l0*s1
    d2 = l2*s0+2.0_wp*l1*s1+l0*s2

  end subroutine damp_rab

  pure subroutine quartic_mix(ns,ip,iq,sp,sq,cp,cq,f0,f1,f2)
    !***********************************************************************
    !* Bound-HB donor/acceptor mixing f = (cp s_p^4 + cq s_q^4)/(s_p^4 + s_q^4)
    !*   = cq + (cp - cq) G,  G = u/(u+v),  u = s_p^4,  v = s_q^4.
    !***********************************************************************
    integer,intent(in) :: ns,ip,iq
    real(wp),intent(in) :: sp,sq,cp,cq
    real(wp),intent(out) :: f0,f1(ns),f2(ns,ns)

    real(wp) :: uu,vv,dd,gu,gv,guu,guv,gvv,up,upp,vp,vpp,dc,g0

    f1 = 0.0_wp
    f2 = 0.0_wp

    uu = sp**4
    vv = sq**4
    dd = uu+vv
    g0 = uu/dd
    gu = vv/(dd*dd)
    gv = -uu/(dd*dd)
    guu = -2.0_wp*vv/(dd*dd*dd)
    guv = (uu-vv)/(dd*dd*dd)
    gvv = 2.0_wp*uu/(dd*dd*dd)

    up = 4.0_wp*sp**3
    upp = 12.0_wp*sp*sp
    vp = 4.0_wp*sq**3
    vpp = 12.0_wp*sq*sq

    dc = cp-cq
    f0 = cq+dc*g0
    f1(ip) = dc*gu*up
    f1(iq) = dc*gv*vp
    f2(ip,ip) = dc*(guu*up*up+gu*upp)
    f2(iq,iq) = dc*(gvv*vp*vp+gv*vpp)
    f2(ip,iq) = dc*guv*up*vp
    f2(iq,ip) = f2(ip,iq)

  end subroutine quartic_mix

  pure function qscale(st,sf,qq) result(f)
    !***********************************************************************
    !* Charge-scaled prefactor from a frozen topo%qa charge; no derivative.
    !***********************************************************************
    real(wp),intent(in) :: st,sf,qq
    real(wp) :: f,e1
    e1 = exp(st*qq)
    f = e1/(e1+sf)
  end function qscale

  subroutine hess_hbonds_bound(n,at,xyz,mcf_ehb,param,topo,neigh,nlist,hess)
    !***********************************************************************
    !* Hessian of the bound-hydrogen hydrogen bonds (hblist1), counterpart of
    !* eg_hbonds_bound/abhgfnff_eg1: E = -q_H bas aci R(r_AB) O(r_AH,r_BH,r_AB).
    !* mcf_ehb is the mcGFN-FF HB scaling (1.0 for GFN-FF); hess is incremented.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n),mcf_ehb
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    type(TGFFNeighbourList),intent(in) :: nlist
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer,parameter :: ns = 3
    integer :: i,ia,ib,ih,iTrA,iTrB,idx(3),p,q
    logical :: ok
    real(wp) :: drah(3),drbh(3),drab(3),ca(2),cb(2)
    real(wp) :: s1,s2,s3,radab,qh,qa,qb,pre
    real(wp) :: gg(ns,9),hh(ns,9,9),blk(9,9)
    real(wp) :: e0,e1(ns),e2(ns,ns),f0,f1(ns),f2(ns,ns)
    real(wp) :: d0,d1,d2,r0,r1,r2

    if (nlist%nhb1 .le. 0) return

    do i = 1,nlist%nhb1
      ia = nlist%hblist1(1,i)
      ib = nlist%hblist1(2,i)
      ih = nlist%hblist1(3,i)
      iTrA = nlist%hblist1(4,i)
      iTrB = nlist%hblist1(5,i)
      if (iTrA .gt. neigh%nTrans.or.iTrB .gt. neigh%nTrans) cycle

      drah = (xyz(:,ia)+neigh%transVec(:,iTrA))-xyz(:,ih)
      drbh = (xyz(:,ib)+neigh%transVec(:,iTrB))-xyz(:,ih)
      drab = (xyz(:,ia)+neigh%transVec(:,iTrA))-(xyz(:,ib)+neigh%transVec(:,iTrB))

      !>-- slots 1 = A, 2 = B, 3 = H
      call dist_in_block(3,1,3,drah,s1,gg(1,:),hh(1,:,:))
      call dist_in_block(3,2,3,drbh,s2,gg(2,:),hh(2,:,:))
      call dist_in_block(3,1,2,drab,s3,gg(3,:),hh(3,:,:))
      if (s1 .lt. eps12.or.s2 .lt. eps12.or.s3 .lt. eps12) cycle

      call hbonds(ia,ib,ca,cb,param,topo)
      radab = param%rad(at(ia))+param%rad(at(ib))
      qh = qscale(param%hbst,param%hbsf,topo%qa(ih))
      qa = qscale(-param%hbst,param%hbsf,topo%qa(ia))
      qb = qscale(-param%hbst,param%hbsf,topo%qa(ib))

      call outl_factor(ns,1,2,3,param%hbacut/radab,s1,s2,s3,.false.,e0,e1,e2,ok)
      if (.not.ok) cycle

      !>-- bas, then aci: same mixing shape with the acidities swapped
      call quartic_mix(ns,1,2,s1,s2,qa*ca(1),qb*cb(1),f0,f1,f2)
      call sfun_mul(ns,e0,e1,e2,f0,f1,f2)

      call quartic_mix(ns,1,2,s1,s2,cb(2),ca(2),f0,f1,f2)
      call sfun_mul(ns,e0,e1,e2,f0,f1,f2)

      call damp_rab(param,radab,s3,d0,d1,d2)
      r0 = d0/(s3*s3*s3)
      r1 = d1/(s3*s3*s3)-3.0_wp*d0/(s3**4)
      r2 = d2/(s3*s3*s3)-6.0_wp*d1/(s3**4)+12.0_wp*d0/(s3**5)
      f0 = r0
      f1 = 0.0_wp
      f2 = 0.0_wp
      f1(3) = r1
      f2(3,3) = r2
      call sfun_mul(ns,e0,e1,e2,f0,f1,f2)

      pre = -qh*mcf_ehb
      do q = 1,ns
        do p = 1,ns
          e2(p,q) = pre*e2(p,q)
        end do
      end do
      do p = 1,ns
        e1(p) = pre*e1(p)
      end do

      call combine_prims(3,ns,e1,e2,gg,hh,blk)
      idx = [ia,ib,ih]
      call scatter_block(hess,3,idx,blk)
    end do

  end subroutine hess_hbonds_bound

  pure subroutine rdamp_mix(ns,ibh,iab,dd0,dd1,dd2,sbh,sab,pbh,pab,f0,f1,f2)
    !***********************************************************************
    !* Radial factor Rd = damp(r_AB) (p_bh/r_BH^3 + p_ab/r_AB^3) of the unbound
    !* forms; dd0/dd1/dd2 are damp(r_AB) and its r_AB derivatives.
    !***********************************************************************
    integer,intent(in) :: ns,ibh,iab
    real(wp),intent(in) :: dd0,dd1,dd2,sbh,sab,pbh,pab
    real(wp),intent(out) :: f0,f1(ns),f2(ns,ns)

    real(wp) :: pv,p2,p3,p22,p33

    f1 = 0.0_wp
    f2 = 0.0_wp

    pv = pbh/sbh**3+pab/sab**3
    p2 = -3.0_wp*pbh/sbh**4
    p3 = -3.0_wp*pab/sab**4
    p22 = 12.0_wp*pbh/sbh**5
    p33 = 12.0_wp*pab/sab**5

    f0 = dd0*pv
    f1(ibh) = dd0*p2
    f1(iab) = dd1*pv+dd0*p3
    f2(ibh,ibh) = dd0*p22
    f2(ibh,iab) = dd1*p2
    f2(iab,ibh) = f2(ibh,iab)
    f2(iab,iab) = dd2*pv+2.0_wp*dd1*p3+dd0*p33

  end subroutine rdamp_mix

  subroutine hess_hbonds_unbound(n,at,xyz,mcf_ehb,param,topo,neigh,nlist,hess, &
     &                           ncovered)
    !***********************************************************************
    !* Hessian of the unbound-hydrogen hydrogen bonds (hblist2), counterpart of
    !* eg_hbonds_unbound with the same branch selection. Carbonyl and N lone-pair
    !* acceptors are delegated; the default form (abhgfnff_eg2new) is done here.
    !* ncovered returns the number of list entries treated.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n),mcf_ehb
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TGFFNeighbourList),intent(in) :: nlist
    real(wp),intent(inout) :: hess(3*n,3*n)
    integer,intent(out) :: ncovered

    integer :: i,ia,ib,ih,iTrA,iTrB,iTr,inb,nbb,nbk,nbnbk,atnb,iTrDum
    integer :: ns,nc,ii,ip,iq,p,q
    logical :: ok
    integer,allocatable :: idx(:),nbl(:)
    real(wp),allocatable :: gg(:,:),hh(:,:,:),blk(:,:)
    real(wp),allocatable :: e1(:),e2(:,:),f1(:),f2(:,:)
    real(wp),allocatable :: sanb(:),sbnb(:)
    real(wp) :: drah(3),drbh(3),drab(3),dranb(3),drbnb(3),vTrN(3)
    real(wp) :: s1,s2,s3,sa,sb,radab,qh,qa,qb,cnst,pre,e0,f0
    real(wp) :: d0,d1,d2,p_bh,p_ab,hbnbc,ca(2),cb(2)

    ncovered = 0
    if (nlist%nhb2 .le. 0) return

    p_bh = 1.0_wp+param%hbabmix
    p_ab = -param%hbabmix

    do i = 1,nlist%nhb2
      ia = nlist%hblist2(1,i)     !> A, the donor
      ib = nlist%hblist2(2,i)     !> B, the acceptor
      ih = nlist%hblist2(3,i)     !> H, always in the central cell
      iTrA = nlist%hblist2(4,i)
      iTrB = nlist%hblist2(5,i)
      if (iTrA .gt. neigh%nTrans.or.iTrB .gt. neigh%nTrans) cycle

      nbnbk = 0
      atnb = 0
      if (at(ib) .eq. 8.and.sum(neigh%nb(neigh%numnb,ib,:)) .eq. 1) then
        nbk = 0
        iTr = 0
        call neigh%jth_nb(n,xyz,nbk,1,ib,iTr)
        iTrDum = neigh%fTrSum(iTr,iTrB)
        if (iTrDum .eq. -1.or.iTrDum .gt. neigh%nTrans) cycle
        if (nbk .ne. 0) then
          nbnbk = sum(neigh%nb(neigh%numnb,nbk,:))
          atnb = at(nbk)
        end if
      end if
      if (at(ib) .eq. 8.and.sum(neigh%nb(neigh%numnb,ib,:)) .eq. 1 &
         & .and.(atnb .eq. 6.or.atnb .eq. 7).and.nbnbk .gt. 1) then
        call hess_hb_carbonyl(n,at,xyz,mcf_ehb,ia,ib,ih,nbk,iTrA,iTrB,iTrDum, &
           & param,topo,neigh,hess)
        ncovered = ncovered+1
        cycle
      end if
      if (at(ib) .eq. 7.and.sum(neigh%nb(neigh%numnb,ib,:)) .eq. 2) then
        call hess_hb_lonepair(n,at,xyz,mcf_ehb,ia,ib,ih,iTrA,iTrB,param,topo, &
           & neigh,hess,ok)
        if (ok) ncovered = ncovered+1
        cycle
      end if

      nbb = sum(neigh%nb(neigh%numnb,ib,:))
      nc = 3+nbb
      ns = 3+2*nbb
      allocate (gg(ns,3*nc),hh(ns,3*nc,3*nc),blk(3*nc,3*nc))
      allocate (e1(ns),e2(ns,ns),f1(ns),f2(ns,ns),idx(nc),nbl(max(nbb,1)))
      allocate (sanb(max(nbb,1)),sbnb(max(nbb,1)))

      drah = (xyz(:,ia)+neigh%transVec(:,iTrA))-xyz(:,ih)
      drbh = (xyz(:,ib)+neigh%transVec(:,iTrB))-xyz(:,ih)
      drab = (xyz(:,ia)+neigh%transVec(:,iTrA))-(xyz(:,ib)+neigh%transVec(:,iTrB))

      !>-- slots 1 = A, 2 = B, 3 = H, 3+k = k-th neighbour of B
      call dist_in_block(nc,1,3,drah,s1,gg(1,:),hh(1,:,:))
      call dist_in_block(nc,2,3,drbh,s2,gg(2,:),hh(2,:,:))
      call dist_in_block(nc,1,2,drab,s3,gg(3,:),hh(3,:,:))
      idx(1) = ia
      idx(2) = ib
      idx(3) = ih

      ok = s1 .ge. eps12.and.s2 .ge. eps12.and.s3 .ge. eps12
      do ii = 1,nbb
        inb = 0
        iTr = 0
        call neigh%jth_nb(n,xyz,inb,ii,ib,iTr)
        nbl(ii) = inb
        idx(3+ii) = inb
        vTrN = neigh%transVec(:,iTr)+neigh%transVec(:,iTrB)
        dranb = (xyz(:,ia)+neigh%transVec(:,iTrA))-(xyz(:,inb)+vTrN)
        drbnb = (xyz(:,ib)+neigh%transVec(:,iTrB))-(xyz(:,inb)+vTrN)
        ip = 2+2*ii
        iq = 3+2*ii
        call dist_in_block(nc,1,3+ii,dranb,sa,gg(ip,:),hh(ip,:,:))
        call dist_in_block(nc,2,3+ii,drbnb,sb,gg(iq,:),hh(iq,:,:))
        sanb(ii) = sa
        sbnb(ii) = sb
        if (sa .lt. eps12.or.sb .lt. eps12) ok = .false.
      end do
      if (.not.ok) then
        deallocate (gg,hh,blk,e1,e2,f1,f2,idx,nbl,sanb,sbnb)
        cycle
      end if

      call hbonds(ia,ib,ca,cb,param,topo)
      radab = param%rad(at(ia))+param%rad(at(ib))
      qh = qscale(param%hbst,param%hbsf,topo%qa(ih))
      qa = qscale(-param%hbst,param%hbsf,topo%qa(ia))
      qb = qscale(-param%hbst,param%hbsf,topo%qa(ib))
      cnst = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

      call outl_factor(ns,1,2,3,param%hbacut/radab,s1,s2,s3,.false.,e0,e1,e2,ok)
      if (.not.ok) then
        deallocate (gg,hh,blk,e1,e2,f1,f2,idx,nbl,sanb,sbnb)
        cycle
      end if

      !>-- one A...nb(B)-B factor per neighbour of the acceptor
      if (at(ib) .eq. 7.and.nbb .eq. 1) then
        hbnbc = 2.0_wp
      else
        hbnbc = param%hbnbcut
      end if
      do ii = 1,nbb
        ip = 2+2*ii
        iq = 3+2*ii
        call outl_factor(ns,ip,iq,3,hbnbc/radab,sanb(ii),sbnb(ii),s3,.true., &
           & f0,f1,f2,ok)
        call sfun_mul(ns,e0,e1,e2,f0,f1,f2)
      end do

      call damp_rab(param,radab,s3,d0,d1,d2)
      call rdamp_mix(ns,2,3,d0,d1,d2,s2,s3,p_bh,p_ab,f0,f1,f2)
      call sfun_mul(ns,e0,e1,e2,f0,f1,f2)

      pre = -cnst*qh*mcf_ehb
      do q = 1,ns
        do p = 1,ns
          e2(p,q) = pre*e2(p,q)
        end do
      end do
      do p = 1,ns
        e1(p) = pre*e1(p)
      end do

      call combine_prims(nc,ns,e1,e2,gg,hh,blk)
      call scatter_block(hess,nc,idx,blk)
      ncovered = ncovered+1
      deallocate (gg,hh,blk,e1,e2,f1,f2,idx,nbl,sanb,sbnb)
    end do
  end subroutine hess_hbonds_unbound

  subroutine hess_hb_carbonyl(n,at,xyz,mcf_ehb,ia,ib,ih,ic,iTrA,iTrB,iTrC, &
     &                        param,topo,neigh,hess)
    !***********************************************************************
    !* Carbonyl/nitro form R-C=O...H-A of the unbound hydrogen bond (abhgfnff_eg3):
    !* default form times a bend and one torsion factor per other neighbour R_i of C,
    !*   eangl = 1 - k (cos t - cos 120 deg)^2,  etors_i = 2 fc cos^2 phi_i + tshift
    !* (rn = 2, phi0 = pi/2, so cos(rn (phi - phi0) + pi) = 2 cos^2 phi - 1).
    !* Block slots: 1 = A, 2 = B (carbonyl O), 3 = H, 4 = C, 4+i = R_i.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n),ia,ib,ih,ic,iTrA,iTrB,iTrC
    real(wp),intent(in) :: xyz(3,n),mcf_ehb
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: nc,ns,ntors,i,j,iTr,p,q,slots(4)
    integer,allocatable :: idx(:),rlist(:),rtr(:)
    logical :: ok
    real(wp),allocatable :: gg(:,:),hh(:,:,:),blk(:,:)
    real(wp),allocatable :: e1(:),e2(:,:),f1(:),f2(:,:)
    real(wp) :: drah(3),drbh(3),drab(3),dranb(3),drbnb(3)
    real(wp) :: vTrA(3),vTrB(3),vTrC(3),vTrR(3),ra(3),rb(3),rcv(3)
    real(wp) :: s1,s2,s3,s4,s5,cv,radab,qh,qa,qb,cnst,pre,e0,f0
    real(wp) :: d0,d1,d2,p_bh,p_ab,ca(2),cb(2)
    real(wp) :: kijk,c0cos,fcb,fct,tshift
    real(wp) :: d1s9(9),d2s9(9,9),d1s12(12),d2s12(12,12)

    vTrA = neigh%transVec(:,iTrA)
    vTrB = neigh%transVec(:,iTrB)
    vTrC = neigh%transVec(:,iTrC)

    ntors = 0
    allocate (rlist(max(1,sum(neigh%nb(neigh%numnb,ic,:)))))
    allocate (rtr(size(rlist)))
    do iTr = 1,neigh%numctr
      do i = 1,neigh%nb(neigh%numnb,ic,iTr)
        if (neigh%nb(i,ic,iTr) .eq. ib) cycle
        j = neigh%fTrSum(iTr,iTrC)
        if (j .le. 0.or.j .gt. neigh%nTrans) cycle
        ntors = ntors+1
        rlist(ntors) = neigh%nb(i,ic,iTr)
        rtr(ntors) = j
      end do
    end do

    nc = 4+ntors
    ns = 6+ntors
    allocate (gg(ns,3*nc),hh(ns,3*nc,3*nc),blk(3*nc,3*nc))
    allocate (e1(ns),e2(ns,ns),f1(ns),f2(ns,ns),idx(nc))

    idx(1) = ia
    idx(2) = ib
    idx(3) = ih
    idx(4) = ic
    do i = 1,ntors
      idx(4+i) = rlist(i)
    end do

    drah = (xyz(:,ia)+vTrA)-xyz(:,ih)
    drbh = (xyz(:,ib)+vTrB)-xyz(:,ih)
    drab = (xyz(:,ia)+vTrA)-(xyz(:,ib)+vTrB)
    dranb = (xyz(:,ia)+vTrA)-(xyz(:,ic)+vTrC)
    drbnb = (xyz(:,ib)+vTrB)-(xyz(:,ic)+vTrC)

    call dist_in_block(nc,1,3,drah,s1,gg(1,:),hh(1,:,:))
    call dist_in_block(nc,2,3,drbh,s2,gg(2,:),hh(2,:,:))
    call dist_in_block(nc,1,2,drab,s3,gg(3,:),hh(3,:,:))
    call dist_in_block(nc,1,4,dranb,s4,gg(4,:),hh(4,:,:))
    call dist_in_block(nc,2,4,drbnb,s5,gg(5,:),hh(5,:,:))

    !>-- bend about B between (C - B) and (H - B)
    call prim_cos_angle((xyz(:,ic)+vTrC)-(xyz(:,ib)+vTrB), &
       & xyz(:,ih)-(xyz(:,ib)+vTrB),cv,d1s9,d2s9)
    slots(1:3) = [2,4,3]
    call embed_block(nc,3,slots(1:3),d1s9,d2s9,gg(6,:),hh(6,:,:))

    c0cos = cos(2.0_wp*pi/3.0_wp)
    fcb = 1.0_wp-param%bend_hb
    kijk = fcb/(1.0_wp-c0cos)**2
    tshift = param%tors_hb
    fct = (1.0_wp-tshift)/2.0_wp

    radab = param%rad(at(ia))+param%rad(at(ib))
    p_bh = 1.0_wp+param%hbabmix
    p_ab = -param%hbabmix
    call hbonds(ia,ib,ca,cb,param,topo)
    qh = qscale(param%hbst,param%hbsf,topo%qa(ih))
    qa = qscale(-param%hbst,param%hbsf,topo%qa(ia))
    qb = qscale(-param%hbst,param%hbsf,topo%qa(ib))
    cnst = ca(2)*qa*cb(1)*qb*param%xhaci_coh

    call outl_factor(ns,1,2,3,param%hbacut/radab,s1,s2,s3,.false.,e0,e1,e2,ok)
    if (.not.ok) then
      deallocate (gg,hh,blk,e1,e2,f1,f2,idx,rlist,rtr)
      return
    end if

    call outl_factor(ns,4,5,3,param%hbnbcut/radab,s4,s5,s3,.true.,f0,f1,f2,ok)
    call sfun_mul(ns,e0,e1,e2,f0,f1,f2)

    call damp_rab(param,radab,s3,d0,d1,d2)
    call rdamp_mix(ns,2,3,d0,d1,d2,s2,s3,p_bh,p_ab,f0,f1,f2)
    call sfun_mul(ns,e0,e1,e2,f0,f1,f2)

    f1 = 0.0_wp
    f2 = 0.0_wp
    f0 = 1.0_wp-kijk*(cv-c0cos)**2
    f1(6) = -2.0_wp*kijk*(cv-c0cos)
    f2(6,6) = -2.0_wp*kijk
    call sfun_mul(ns,e0,e1,e2,f0,f1,f2)

    do i = 1,ntors
      vTrR = neigh%transVec(:,rtr(i))
      ra = (xyz(:,ib)+vTrB)-(xyz(:,rlist(i))+vTrR)
      rb = (xyz(:,ic)+vTrC)-(xyz(:,ib)+vTrB)
      rcv = xyz(:,ih)-(xyz(:,ic)+vTrC)
      call prim_cos_dihedral(ra,rb,rcv,cv,d1s12,d2s12)
      slots = [4+i,2,4,3]
      call embed_block(nc,4,slots,d1s12,d2s12,gg(6+i,:),hh(6+i,:,:))
      f1 = 0.0_wp
      f2 = 0.0_wp
      f0 = 2.0_wp*fct*cv*cv+tshift
      f1(6+i) = 4.0_wp*fct*cv
      f2(6+i,6+i) = 4.0_wp*fct
      call sfun_mul(ns,e0,e1,e2,f0,f1,f2)
    end do

    pre = -cnst*qh*mcf_ehb
    do q = 1,ns
      do p = 1,ns
        e2(p,q) = pre*e2(p,q)
      end do
    end do
    do p = 1,ns
      e1(p) = pre*e1(p)
    end do

    call combine_prims(nc,ns,e1,e2,gg,hh,blk)
    call scatter_block(hess,nc,idx,blk)
    deallocate (gg,hh,blk,e1,e2,f1,f2,idx,rlist,rtr)

  end subroutine hess_hb_carbonyl

  pure subroutine outl_factor_c(ns,ip,ir,aa,sp,coff,sr,f0,f1,f2)
    !***********************************************************************
    !* Out-of-line factor with one summed distance held constant (the fixed
    !* B-to-lone-pair offset): z = aa ((s_p + coff + 1e-12)/s_r - 1).
    !***********************************************************************
    integer,intent(in) :: ns,ip,ir
    real(wp),intent(in) :: aa,sp,coff,sr
    real(wp),intent(out) :: f0,f1(ns),f2(ns,ns)

    real(wp) :: sm,zz,ee,den,fz,fzz,z1,z3,z13,z33

    f1 = 0.0_wp
    f2 = 0.0_wp

    sm = sp+coff+eps12
    zz = aa*(sm/sr-1.0_wp)
    ee = exp(zz)
    den = 1.0_wp+ee
    f0 = 2.0_wp/den
    fz = -2.0_wp*ee/(den*den)
    fzz = 2.0_wp*ee*(ee-1.0_wp)/(den*den*den)

    z1 = aa/sr
    z3 = -aa*sm/(sr*sr)
    z13 = -aa/(sr*sr)
    z33 = 2.0_wp*aa*sm/(sr*sr*sr)

    !>-- accumulate: the degenerate lone-pair branch calls this with ip == ir
    f1(ip) = f1(ip)+fz*z1
    f1(ir) = f1(ir)+fz*z3
    f2(ip,ip) = f2(ip,ip)+fzz*z1*z1
    f2(ip,ir) = f2(ip,ir)+fzz*z1*z3+fz*z13
    f2(ir,ip) = f2(ir,ip)+fzz*z1*z3+fz*z13
    f2(ir,ir) = f2(ir,ir)+fzz*z3*z3+fz*z33

  end subroutine outl_factor_c

  pure subroutine ralp_prim(nc,sA,sB,snb,nbb,drab,vv,ldist,val,d1,d2)
    !***********************************************************************
    !* A-to-lone-pair distance of the N-heteroaromatic hydrogen bond with
    !* gradient d1 and Hessian d2 over the nc-atom block:
    !*   val = |drab + ldist u|,  u = v/|v|,  v = sum_i nb_i - nbb B.
    !* sA/sB, snb - block slots of A, B and the neighbours of B
    !* drab - (A + shift) - (B + shift);  vv - bisector sum v;  ldist - lp offset
    !***********************************************************************
    integer,intent(in) :: nc,sA,sB,nbb,snb(nbb)
    real(wp),intent(in) :: drab(3),vv(3),ldist
    real(wp),intent(out) :: val
    real(wp),intent(out) :: d1(3*nc),d2(3*nc,3*nc)

    integer :: a,b,i,j,x,y,ox,oy,ns
    real(wp) :: vn,uu(3),jm(3,3),ww(3),what(3),qq(3,3),wk(3,3)
    real(wp) :: pmat(3,3,nc+1),cx(nc+1),tmp
    integer :: slot(nc+1)

    d1 = 0.0_wp
    d2 = 0.0_wp

    vn = sqrt(vv(1)*vv(1)+vv(2)*vv(2)+vv(3)*vv(3))
    uu = vv/vn

    !>-- normalisation Jacobian J = du/dv and the w vector
    do b = 1,3
      do a = 1,3
        jm(a,b) = (merge(1.0_wp,0.0_wp,a == b)-uu(a)*uu(b))/vn
      end do
    end do
    ww = drab+ldist*uu
    val = sqrt(ww(1)*ww(1)+ww(2)*ww(2)+ww(3)*ww(3))
    if (val < eps12) return
    what = ww/val
    do b = 1,3
      do a = 1,3
        qq(a,b) = (merge(1.0_wp,0.0_wp,a == b)-what(a)*what(b))/val
      end do
    end do

    !>-- W(i,j) = ldist sum_a what_a d2u_a/dv_i dv_j
    do j = 1,3
      do i = 1,3
        tmp = -what(i)*uu(j)-what(j)*uu(i) &
           & -merge(1.0_wp,0.0_wp,i == j)*dot_product(what,uu) &
           & +3.0_wp*dot_product(what,uu)*uu(i)*uu(j)
        wk(i,j) = ldist*tmp/(vn*vn)
      end do
    end do

    !>-- dw/dx = L_x + ldist c_x J, with c_x the scalar in dv/dx = c_x I
    ns = 2+nbb
    slot(1) = sA
    cx(1) = 0.0_wp
    slot(2) = sB
    cx(2) = -real(nbb,wp)
    do i = 1,nbb
      slot(2+i) = snb(i)
      cx(2+i) = 1.0_wp
    end do
    do x = 1,ns
      do b = 1,3
        do a = 1,3
          pmat(a,b,x) = ldist*cx(x)*jm(a,b)
        end do
      end do
    end do
    do a = 1,3
      pmat(a,a,1) = pmat(a,a,1)+1.0_wp      !> A enters w with +I
      pmat(a,a,2) = pmat(a,a,2)-1.0_wp      !> B with -I
    end do

    do x = 1,ns
      ox = 3*(slot(x)-1)
      do i = 1,3
        d1(ox+i) = d1(ox+i)+dot_product(what,pmat(:,i,x))
      end do
    end do
    do y = 1,ns
      oy = 3*(slot(y)-1)
      do x = 1,ns
        ox = 3*(slot(x)-1)
        do j = 1,3
          do i = 1,3
            tmp = dot_product(pmat(:,i,x),matmul(qq,pmat(:,j,y))) &
               & +cx(x)*cx(y)*wk(i,j)
            d2(ox+i,oy+j) = d2(ox+i,oy+j)+tmp
          end do
        end do
      end do
    end do

  end subroutine ralp_prim

  subroutine hess_hb_lonepair(n,at,xyz,mcf_ehb,ia,ib,ih,iTrA,iTrB, &
     &                        param,topo,neigh,hess,ok)
    !***********************************************************************
    !* N-heteroaromatic form of the unbound hydrogen bond (abhgfnff_eg2_rnr):
    !* default form times an out-of-line factor Olp(r_A-lp, r_AB) on a lone-pair
    !* site on the bisector of the acceptor's two bonds.
    !* Block slots: 1 = A, 2 = B, 3 = H, 4/5 = the two neighbours of B.
    !* ok is false when the overflow guard trips or a distance vanishes.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n),ia,ib,ih,iTrA,iTrB
    real(wp),intent(in) :: xyz(3,n),mcf_ehb
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    real(wp),intent(inout) :: hess(3*n,3*n)
    logical,intent(out) :: ok

    integer,parameter :: nbb = 2,nc = 3+nbb,ns = 4+2*nbb
    integer :: i,inb,iTr,ip,iq,p,q,idx(nc),snb(nbb)
    real(wp) :: gg(ns,3*nc),hh(ns,3*nc,3*nc),blk(3*nc,3*nc)
    real(wp) :: e1(ns),e2(ns,ns),f1(ns),f2(ns,ns)
    real(wp) :: drah(3),drbh(3),drab(3),dranb(3),drbnb(3),vTrN(3),vv(3)
    real(wp) :: sanb(nbb),sbnb(nbb)
    real(wp) :: s1,s2,s3,slp,sa,sb,radab,qh,qa,qb,cnst,pre,e0,f0
    real(wp) :: d0,d1,d2,p_bh,p_ab,ldist,vn,ca(2),cb(2)
    real(wp),parameter :: hblpcut = 56.0_wp
    logical :: degen

    ok = .false.
    p_bh = 1.0_wp+param%hbabmix
    p_ab = -param%hbabmix
    ldist = 0.50_wp-0.018_wp*param%repz(at(ib))

    drah = (xyz(:,ia)+neigh%transVec(:,iTrA))-xyz(:,ih)
    drbh = (xyz(:,ib)+neigh%transVec(:,iTrB))-xyz(:,ih)
    drab = (xyz(:,ia)+neigh%transVec(:,iTrA))-(xyz(:,ib)+neigh%transVec(:,iTrB))

    call dist_in_block(nc,1,3,drah,s1,gg(1,:),hh(1,:,:))
    call dist_in_block(nc,2,3,drbh,s2,gg(2,:),hh(2,:,:))
    call dist_in_block(nc,1,2,drab,s3,gg(3,:),hh(3,:,:))
    idx(1) = ia
    idx(2) = ib
    idx(3) = ih
    if (s1 .lt. eps12.or.s2 .lt. eps12.or.s3 .lt. eps12) return

    vv = 0.0_wp
    do i = 1,nbb
      inb = 0
      iTr = 0
      call neigh%jth_nb(n,xyz,inb,i,ib,iTr)
      idx(3+i) = inb
      snb(i) = 3+i
      vTrN = neigh%transVec(:,iTr)+neigh%transVec(:,iTrB)
      dranb = (xyz(:,ia)+neigh%transVec(:,iTrA))-(xyz(:,inb)+vTrN)
      drbnb = (xyz(:,ib)+neigh%transVec(:,iTrB))-(xyz(:,inb)+vTrN)
      vv = vv+((xyz(:,inb)+vTrN)-(xyz(:,ib)+neigh%transVec(:,iTrB)))
      ip = 2+2*i
      iq = 3+2*i
      call dist_in_block(nc,1,3+i,dranb,sa,gg(ip,:),hh(ip,:,:))
      call dist_in_block(nc,2,3+i,drbnb,sb,gg(iq,:),hh(iq,:,:))
      sanb(i) = sa
      sbnb(i) = sb
      if (sa .lt. eps12.or.sb .lt. eps12) return
    end do

    vn = sqrt(vv(1)*vv(1)+vv(2)*vv(2)+vv(3)*vv(3))
    degen = vn .le. 1.0e-10_wp

    !>-- lone-pair distance, primitive 8. For collinear bonds at B the energy code
    !>   puts the lone pair on B; Olp then uses r_AB.
    if (.not.degen) then
      call ralp_prim(nc,1,2,snb,nbb,drab,vv,ldist,slp,gg(8,:),hh(8,:,:))
      if (slp .lt. eps12) return
    else
      slp = s3
      gg(8,:) = 0.0_wp
      hh(8,:,:) = 0.0_wp
    end if

    call hbonds(ia,ib,ca,cb,param,topo)
    radab = param%rad(at(ia))+param%rad(at(ib))
    qh = qscale(param%hbst,param%hbsf,topo%qa(ih))
    qa = qscale(-param%hbst,param%hbsf,topo%qa(ia))
    qb = qscale(-param%hbst,param%hbsf,topo%qa(ib))
    cnst = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

    call outl_factor(ns,1,2,3,param%hbacut/radab,s1,s2,s3,.false.,e0,e1,e2,ok)
    if (.not.ok) return

    do i = 1,nbb
      call outl_factor(ns,2+2*i,3+2*i,3,param%hbnbcut/radab,sanb(i),sbnb(i), &
         & s3,.true.,f0,f1,f2,ok)
      call sfun_mul(ns,e0,e1,e2,f0,f1,f2)
    end do
    if (.not.degen) then
      call outl_factor_c(ns,8,3,hblpcut/radab,slp,ldist,s3,f0,f1,f2)
    else
      call outl_factor_c(ns,3,3,hblpcut/radab,s3,0.0_wp,s3,f0,f1,f2)
    end if
    call sfun_mul(ns,e0,e1,e2,f0,f1,f2)

    call damp_rab(param,radab,s3,d0,d1,d2)
    call rdamp_mix(ns,2,3,d0,d1,d2,s2,s3,p_bh,p_ab,f0,f1,f2)
    call sfun_mul(ns,e0,e1,e2,f0,f1,f2)

    pre = -cnst*qh*mcf_ehb
    do q = 1,ns
      do p = 1,ns
        e2(p,q) = pre*e2(p,q)
      end do
    end do
    do p = 1,ns
      e1(p) = pre*e1(p)
    end do

    call combine_prims(nc,ns,e1,e2,gg,hh,blk)
    call scatter_block(hess,nc,idx,blk)
    ok = .true.

  end subroutine hess_hb_lonepair

end module gfnff_hess_hb
