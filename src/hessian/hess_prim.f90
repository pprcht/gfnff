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
!> Internal-coordinate primitives for the closed-form GFN-FF Hessian: values,
!> Cartesian gradients and Cartesian Hessians. hess_bonded combines them.
!> r^2 is used rather than r because gfnffdampa/gfnffdampt take r^2, and
!> cos(theta) rather than theta because egbend's non-linear branch is in cosa.
!> Component (c,A) sits at 3*(A-1)+c, with A counting the atoms of the
!> primitive in the order documented at each routine.
module gfnff_hess_prim

  use iso_fortran_env,only:wp => real64
  implicit none
  private

  public :: prim_rsq,prim_cos_angle,cos_kernel,rsq_in_block,dist_in_block
  public :: psi_from_cos,chain_linear,scatter_block,combine_prims
  public :: embed_block
  public :: prim_cos_dihedral,prim_sin_oop,cheby_t

contains  !> MODULE PROCEDURES START HERE

  pure subroutine prim_rsq(vec,val,d1,d2)
    !***********************************************************************
    !* r^2 = |R_i - R_j|^2 from vec = R_i - R_j, with (6) gradient and (6,6)
    !* Hessian over the atoms ordered (i,j). The Hessian is constant.
    !***********************************************************************
    real(wp),intent(in) :: vec(3)
    real(wp),intent(out) :: val
    real(wp),intent(out) :: d1(6)
    real(wp),intent(out) :: d2(6,6)

    integer :: a

    val = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
    d1(1:3) = 2.0_wp*vec
    d1(4:6) = -2.0_wp*vec

    d2 = 0.0_wp
    do a = 1,3
      d2(a,a) = 2.0_wp
      d2(3+a,3+a) = 2.0_wp
      d2(a,3+a) = -2.0_wp
      d2(3+a,a) = -2.0_wp
    end do

  end subroutine prim_rsq

  pure subroutine rsq_in_block(nc,si,sj,vec,val,d1,d2)
    !***********************************************************************
    !* r^2 = |R_si - R_sj|^2 over an nc-atom block, for combination with the
    !* angle and torsion primitives, whose damping is a function of r^2.
    !*   si,sj - slots of the two atoms within the block, 1..nc
    !*   vec   - R_si - R_sj, including any periodic translation
    !***********************************************************************
    integer,intent(in) :: nc,si,sj
    real(wp),intent(in) :: vec(3)
    real(wp),intent(out) :: val
    real(wp),intent(out) :: d1(3*nc),d2(3*nc,3*nc)

    integer :: a,oi,oj

    val = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
    oi = 3*(si-1)
    oj = 3*(sj-1)

    d1 = 0.0_wp
    d2 = 0.0_wp
    do a = 1,3
      d1(oi+a) = 2.0_wp*vec(a)
      d1(oj+a) = -2.0_wp*vec(a)
      d2(oi+a,oi+a) = 2.0_wp
      d2(oj+a,oj+a) = 2.0_wp
      d2(oi+a,oj+a) = -2.0_wp
      d2(oj+a,oi+a) = -2.0_wp
    end do

  end subroutine rsq_in_block

  pure subroutine dist_in_block(nc,si,sj,vec,val,d1,d2)
    !***********************************************************************
    !* r = |R_si - R_sj| over an nc-atom block, arguments as in rsq_in_block.
    !* The hydrogen and halogen bond potentials are written in r, not r^2.
    !***********************************************************************
    integer,intent(in) :: nc,si,sj
    real(wp),intent(in) :: vec(3)
    real(wp),intent(out) :: val
    real(wp),intent(out) :: d1(3*nc),d2(3*nc,3*nc)

    integer :: a,b,oi,oj
    real(wp) :: e(3),rinv,blk

    val = sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))
    d1 = 0.0_wp
    d2 = 0.0_wp
    if (val < 1.0e-12_wp) return

    rinv = 1.0_wp/val
    e = vec*rinv
    oi = 3*(si-1)
    oj = 3*(sj-1)

    do a = 1,3
      d1(oi+a) = e(a)
      d1(oj+a) = -e(a)
    end do
    do b = 1,3
      do a = 1,3
        blk = (merge(1.0_wp,0.0_wp,a == b)-e(a)*e(b))*rinv
        d2(oi+a,oi+b) = blk
        d2(oj+a,oj+b) = blk
        d2(oi+a,oj+b) = -blk
        d2(oj+a,oi+b) = -blk
      end do
    end do

  end subroutine dist_in_block

  pure subroutine cos_kernel(a,b,val,d1,d2)
    !***********************************************************************
    !* c = (a.b)/(|a||b|) of two free vectors, with first and second
    !* derivatives over the stacked six-vector (a,b); zeros for a null vector.
    !* Used on bond vectors for angles and on plane normals for dihedrals.
    !***********************************************************************
    real(wp),intent(in) :: a(3),b(3)
    real(wp),intent(out) :: val
    real(wp),intent(out) :: d1(6)
    real(wp),intent(out) :: d2(6,6)

    integer :: i,j
    real(wp) :: na,nb,p,q,ab,c,p2,q2,p3q,pq3,pq,dij

    val = 0.0_wp
    d1 = 0.0_wp
    d2 = 0.0_wp

    na = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
    nb = sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))
    if (na < 1.0e-12_wp.or.nb < 1.0e-12_wp) return

    p = 1.0_wp/na
    q = 1.0_wp/nb
    ab = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
    c = ab*p*q
    val = c

    p2 = p*p
    q2 = q*q
    pq = p*q
    p3q = p2*p*q
    pq3 = p*q2*q

    d1(1:3) = pq*b-c*p2*a
    d1(4:6) = pq*a-c*q2*b

    do j = 1,3
      do i = 1,3
        dij = merge(1.0_wp,0.0_wp,i == j)
        d2(i,j) = -p3q*(a(i)*b(j)+b(i)*a(j))+3.0_wp*c*p2*p2*a(i)*a(j)-c*p2*dij
        d2(3+i,3+j) = -pq3*(b(i)*a(j)+a(i)*b(j))+3.0_wp*c*q2*q2*b(i)*b(j)-c*q2*dij
        d2(i,3+j) = -pq3*b(i)*b(j)+pq*dij-p3q*a(i)*a(j)+c*p2*q2*a(i)*b(j)
      end do
    end do
    do j = 1,3
      do i = 1,3
        d2(3+j,i) = d2(i,3+j)
      end do
    end do

  end subroutine cos_kernel

  pure subroutine prim_cos_angle(va,vb,val,d1,d2)
    !***********************************************************************
    !* cos(theta) of the bond angle j-i-k about vertex i from the bond vectors
    !* va = R_j - R_i and vb = R_k - R_i, with (9) gradient and (9,9) Hessian
    !* over the atoms ordered (i,j,k).
    !***********************************************************************
    real(wp),intent(in) :: va(3),vb(3)
    real(wp),intent(out) :: val
    real(wp),intent(out) :: d1(9)
    real(wp),intent(out) :: d2(9,9)

    integer :: a
    real(wp) :: k1(6),k2(6,6),tmat(6,9)

    call cos_kernel(va,vb,val,k1,k2)

    tmat = 0.0_wp
    do a = 1,3
      tmat(a,a) = -1.0_wp        !> d va / d R_i
      tmat(a,3+a) = 1.0_wp       !> d va / d R_j
      tmat(3+a,a) = -1.0_wp      !> d vb / d R_i
      tmat(3+a,6+a) = 1.0_wp     !> d vb / d R_k
    end do

    call chain_linear(6,9,tmat,k1,k2,d1,d2)

  end subroutine prim_cos_angle

  pure function skew(a) result(m)
    !***********************************************************************
    !* Cross-product matrix [a]x, so that [a]x v = a x v.
    !***********************************************************************
    real(wp),intent(in) :: a(3)
    real(wp) :: m(3,3)
    m(1,1) = 0.0_wp; m(1,2) = -a(3); m(1,3) = a(2)
    m(2,1) = a(3); m(2,2) = 0.0_wp; m(2,3) = -a(1)
    m(3,1) = -a(2); m(3,2) = a(1); m(3,3) = 0.0_wp
  end function skew

  pure subroutine prim_cos_dihedral(ra,rb,rc,val,d1,d2)
    !***********************************************************************
    !* cos(phi) of the proper dihedral i-j-k-l, consistent with torsPBC's
    !* snanb, with (12) gradient and (12,12) Hessian over the atoms (i,j,k,l).
    !*   ra,rb,rc - R_j - R_i, R_k - R_j, R_l - R_k; a constant periodic
    !*              translation may be folded in
    !* Evaluated as cos_kernel on the plane normals ra x rb and rb x rc.
    !***********************************************************************
    real(wp),intent(in) :: ra(3),rb(3),rc(3)
    real(wp),intent(out) :: val
    real(wp),intent(out) :: d1(12)
    real(wp),intent(out) :: d2(12,12)

    integer :: a
    real(wp) :: pv(3),qv(3),k1(6),k2(6,6)
    real(wp) :: jac(6,9),g9(9),h9(9,9),tmat(9,12)
    real(wp) :: mp(3,3),mq(3,3)

    pv = cross(ra,rb)
    qv = cross(rb,rc)
    call cos_kernel(pv,qv,val,k1,k2)

    !>-- Jacobian of (p,q) with respect to (ra,rb,rc)
    jac = 0.0_wp
    jac(1:3,1:3) = -skew(rb)
    jac(1:3,4:6) = skew(ra)
    jac(4:6,4:6) = -skew(rc)
    jac(4:6,7:9) = skew(rb)

    call chain_linear(6,9,jac,k1,k2,g9,h9)

    !>-- the cross products are bilinear and add a constant second derivative
    !>   weighted by dcos/dp and dcos/dq
    mp = -skew(k1(1:3))
    mq = -skew(k1(4:6))
    h9(1:3,4:6) = h9(1:3,4:6)+mp
    h9(4:6,1:3) = h9(4:6,1:3)+transpose(mp)
    h9(4:6,7:9) = h9(4:6,7:9)+mq
    h9(7:9,4:6) = h9(7:9,4:6)+transpose(mq)

    tmat = 0.0_wp
    do a = 1,3
      tmat(a,a) = -1.0_wp
      tmat(a,3+a) = 1.0_wp
      tmat(3+a,3+a) = -1.0_wp
      tmat(3+a,6+a) = 1.0_wp
      tmat(6+a,6+a) = -1.0_wp
      tmat(6+a,9+a) = 1.0_wp
    end do

    call chain_linear(9,12,tmat,g9,h9,d1,d2)

  end subroutine prim_cos_dihedral

  pure subroutine prim_sin_oop(re,rd,rv,val,d1,d2)
    !***********************************************************************
    !* sin(omega) of the out-of-plane angle of GFN-FF's inversion term, the
    !* quantity omegaPBC feeds to asin, with (12) gradient and (12,12) Hessian
    !* over the atoms ordered (i,j,k,l), i central.
    !*   re,rd,rv - R_i - R_j, R_k - R_j, R_l - R_i
    !* Evaluated as cos_kernel on the plane normal re x rd and rv. The energy
    !* needs omega only through cos(omega), see hess_torsions.
    !***********************************************************************
    real(wp),intent(in) :: re(3),rd(3),rv(3)
    real(wp),intent(out) :: val
    real(wp),intent(out) :: d1(12)
    real(wp),intent(out) :: d2(12,12)

    integer :: a
    real(wp) :: nv(3),k1(6),k2(6,6)
    real(wp) :: jac(6,9),g9(9),h9(9,9),tmat(9,12),mp(3,3)

    nv = cross(re,rd)
    call cos_kernel(nv,rv,val,k1,k2)

    !>-- Jacobian of (n,rv) with respect to (re,rd,rv)
    jac = 0.0_wp
    jac(1:3,1:3) = -skew(rd)
    jac(1:3,4:6) = skew(re)
    do a = 1,3
      jac(3+a,6+a) = 1.0_wp
    end do

    call chain_linear(6,9,jac,k1,k2,g9,h9)

    mp = -skew(k1(1:3))
    h9(1:3,4:6) = h9(1:3,4:6)+mp
    h9(4:6,1:3) = h9(4:6,1:3)+transpose(mp)

    tmat = 0.0_wp
    do a = 1,3
      tmat(a,a) = 1.0_wp
      tmat(a,3+a) = -1.0_wp
      tmat(3+a,3+a) = -1.0_wp
      tmat(3+a,6+a) = 1.0_wp
      tmat(6+a,a) = -1.0_wp
      tmat(6+a,9+a) = 1.0_wp
    end do

    call chain_linear(9,12,tmat,g9,h9,d1,d2)

  end subroutine prim_sin_oop

  pure function cross(a,b) result(c)
    !***********************************************************************
    !* Vector cross product a x b.
    !***********************************************************************
    real(wp),intent(in) :: a(3),b(3)
    real(wp) :: c(3)
    c(1) = a(2)*b(3)-a(3)*b(2)
    c(2) = a(3)*b(1)-a(1)*b(3)
    c(3) = a(1)*b(2)-a(2)*b(1)
  end function cross

  pure subroutine cheby_t(nn,c,t0,t1,t2)
    !***********************************************************************
    !* Chebyshev polynomial of the first kind T_nn(c) with its first two
    !* derivatives, by the recurrence T_{k+1} = 2 c T_k - T_{k-1}.
    !* cos(n phi) = T_n(cos phi) makes the proper torsion a polynomial in
    !* cos(phi), free of the sqrt(1-c^2) branch point at planar geometries.
    !***********************************************************************
    integer,intent(in) :: nn
    real(wp),intent(in) :: c
    real(wp),intent(out) :: t0,t1,t2

    integer :: k
    real(wp) :: tm0,tm1,tm2,tp0,tp1,tp2

    if (nn .le. 0) then
      t0 = 1.0_wp
      t1 = 0.0_wp
      t2 = 0.0_wp
      return
    end if

    tm0 = 1.0_wp; tm1 = 0.0_wp; tm2 = 0.0_wp       !> T_0
    t0 = c; t1 = 1.0_wp; t2 = 0.0_wp        !> T_1
    do k = 1,nn-1
      tp0 = 2.0_wp*c*t0-tm0
      tp1 = 2.0_wp*t0+2.0_wp*c*t1-tm1
      tp2 = 4.0_wp*t1+2.0_wp*c*t2-tm2
      tm0 = t0; tm1 = t1; tm2 = t2
      t0 = tp0; t1 = tp1; t2 = tp2
    end do

  end subroutine cheby_t

  pure subroutine chain_linear(nv,nc,tmat,d1in,d2in,d1out,d2out)
    !***********************************************************************
    !* Chain rule through a linear map v = T x, exact since T is constant.
    !*   tmat - (nv,nc) Jacobian dv/dx
    !***********************************************************************
    integer,intent(in) :: nv,nc
    real(wp),intent(in) :: tmat(nv,nc)
    real(wp),intent(in) :: d1in(nv),d2in(nv,nv)
    real(wp),intent(out) :: d1out(nc),d2out(nc,nc)

    d1out = matmul(transpose(tmat),d1in)
    d2out = matmul(transpose(tmat),matmul(d2in,tmat))

  end subroutine chain_linear

  pure subroutine psi_from_cos(cosv,gg,gg1,gg2,psi,ps1,ps2)
    !***********************************************************************
    !* Deviation from linearity psi = pi - theta = acos(-c), c = cos(theta).
    !*   gg,gg1,gg2  - G = psi^2, dG/dc, d2G/dc2; smooth at linearity
    !*   psi,ps1,ps2 - psi, dpsi/dc, d2psi/dc2; singular at linearity
    !* With u = 1 + c, psi^2 = 2u + u^2/3 + 4u^3/45 + u^4/35 + ..., so k*psi^2
    !* has a finite Hessian at an exactly linear geometry although psi does
    !* not. Below u = 1e-2 the series replaces the closed form, whose G'' is
    !* a difference of two terms diverging as 1/u. psi is needed when the
    !* equilibrium angle is close to but not exactly pi, see hess_angles.
    !***********************************************************************
    real(wp),intent(in) :: cosv
    real(wp),intent(out) :: gg,gg1,gg2
    real(wp),intent(out) :: psi,ps1,ps2

    real(wp) :: u,w,rw

    u = 1.0_wp+cosv
    u = max(u,0.0_wp)

    !>-- half-angle form, stable for small u
    psi = 2.0_wp*asin(sqrt(0.5_wp*u))

    w = u*(2.0_wp-u)                      !> = 1 - c^2 = sin(theta)^2
    if (w > 1.0e-14_wp) then
      rw = sqrt(w)
      ps1 = 1.0_wp/rw
      ps2 = -(1.0_wp-u)/(w*rw)
    else
      ps1 = 0.0_wp
      ps2 = 0.0_wp
    end if

    if (u > 1.0e-2_wp) then
      gg = psi*psi
      gg1 = 2.0_wp*psi*ps1
      gg2 = 2.0_wp*ps1*ps1+2.0_wp*psi*ps2
    else
      !>-- series for psi^2, relative error below 1e-9 at the switch point
      gg = u*(2.0_wp+u*(1.0_wp/3.0_wp+u*(4.0_wp/45.0_wp &
         & +u*(1.0_wp/35.0_wp+u*1.015873015873016e-2_wp))))
      gg1 = 2.0_wp+u*(2.0_wp/3.0_wp+u*(4.0_wp/15.0_wp &
         & +u*(4.0_wp/35.0_wp+u*5.079365079365079e-2_wp)))
      gg2 = 2.0_wp/3.0_wp+u*(8.0_wp/15.0_wp &
         & +u*(12.0_wp/35.0_wp+u*0.2031746031746032_wp))
    end if

  end subroutine psi_from_cos

  pure subroutine embed_block(nc,nsub,slots,d1s,d2s,d1,d2)
    !***********************************************************************
    !* Place a primitive built on nsub atoms into the nc-atom block of a
    !* larger term; slots(a) is the position of its a-th atom in that block.
    !*   d1,d2 - (3nc) and (3nc,3nc), overwritten, zero outside the slots
    !***********************************************************************
    integer,intent(in) :: nc,nsub,slots(nsub)
    real(wp),intent(in) :: d1s(3*nsub),d2s(3*nsub,3*nsub)
    real(wp),intent(out) :: d1(3*nc),d2(3*nc,3*nc)

    integer :: a,b,ia,ja,p,q

    d1 = 0.0_wp
    d2 = 0.0_wp
    do a = 1,nsub
      ia = 3*(slots(a)-1)
      do p = 1,3
        d1(ia+p) = d1(ia+p)+d1s(3*(a-1)+p)
      end do
      do b = 1,nsub
        ja = 3*(slots(b)-1)
        do q = 1,3
          do p = 1,3
            d2(ia+p,ja+q) = d2(ia+p,ja+q)+d2s(3*(a-1)+p,3*(b-1)+q)
          end do
        end do
      end do
    end do

  end subroutine embed_block

  pure subroutine combine_prims(nc,ns,fd,fdd,g,h,blk)
    !***********************************************************************
    !* Cartesian block of a term E = f(s_1..s_ns) in ns primitive scalars:
    !*   blk = sum_p f_p hess(s_p) + sum_pq f_pq grad(s_p) (x) grad(s_q)
    !*   fd,fdd - (ns) first and (ns,ns) second partials of f
    !*   g,h    - (ns,3nc) primitive gradients and (ns,3nc,3nc) Hessians
    !*   blk    - (3nc,3nc), overwritten
    !***********************************************************************
    integer,intent(in) :: nc,ns
    real(wp),intent(in) :: fd(ns),fdd(ns,ns)
    real(wp),intent(in) :: g(ns,3*nc),h(ns,3*nc,3*nc)
    real(wp),intent(out) :: blk(3*nc,3*nc)

    integer :: a,b,p,q

    blk = 0.0_wp
    do p = 1,ns
      do b = 1,3*nc
        do a = 1,3*nc
          blk(a,b) = blk(a,b)+fd(p)*h(p,a,b)
        end do
      end do
    end do
    do q = 1,ns
      do p = 1,ns
        do b = 1,3*nc
          do a = 1,3*nc
            blk(a,b) = blk(a,b)+fdd(p,q)*g(p,a)*g(q,b)
          end do
        end do
      end do
    end do

  end subroutine combine_prims

  pure subroutine scatter_block(hess,nc,idx,blk)
    !***********************************************************************
    !* Add the (3nc,3nc) block of an nc-atom term into the full Hessian.
    !*   idx - (nc) atom indices, in the order the block was built
    !* Not thread safe: two terms sharing an atom write the same elements, so
    !* a parallel caller must guard the call.
    !***********************************************************************
    real(wp),intent(inout) :: hess(:,:)
    integer,intent(in) :: nc,idx(nc)
    real(wp),intent(in) :: blk(3*nc,3*nc)

    integer :: ia,ja,a,b,ig,jg

    do ja = 1,nc
      jg = 3*(idx(ja)-1)
      do ia = 1,nc
        ig = 3*(idx(ia)-1)
        do b = 1,3
          do a = 1,3
            hess(ig+a,jg+b) = hess(ig+a,jg+b)+blk(3*(ia-1)+a,3*(ja-1)+b)
          end do
        end do
      end do
    end do

  end subroutine scatter_block

end module gfnff_hess_prim
