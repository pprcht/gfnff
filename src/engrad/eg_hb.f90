! ------------------------------------------------------------------------------
! This file is part of gfnff.
!
! Copyright (C) 2019-2020 Stefan Grimme
! Copyright (C) 2023-2026 Philipp Pracht
!
! The energy and gradient routines in this file originate from the GFN-FF
! implementation in the xtb code by S. Spicher and S. Grimme
! (https://github.com/grimme-lab/xtb). They were reorganised into per-term
! modules for this library and are expected to diverge further from the
! upstream implementation over time.
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
!> Hydrogen- and halogen-bond energy terms of GFN-FF.
!> Holds the four A...H...B cases (generic, neighbour-resolved, N-hetero
!> aromatic and carbonyl/nitro), the halogen-bond term, and the ERF
!> coordination number the HB terms are scaled with.
module gfnff_eg_hb

  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData,TGFFNeighbourList,TGFFTopology
  use gfnff_neighbor,only:TNeigh
  use gfnff_topo_hbset,only:hbonds
  use gfnff_geometry,only:crprod,dphidrPBC,impsc,valijklffPBC,vlen,vsub
  implicit none
  private

  !> one driver per HB/XB term; the single-bond kernels (abhgfnff_eg1, ...) stay private
  public :: eg_hbonds_bound,eg_hbonds_unbound,eg_xbonds
  !> the ERF coordination number is also needed by the bond term
  public :: dncoord_erf

  real(wp),private,parameter :: pi = 3.1415926535897932385_wp
  real(wp),private,parameter :: sqrtpi = 1.77245385091_wp

contains  !> MODULE PROCEDURES START HERE

  subroutine eg_hbonds_bound(n,at,xyz,mcf_ehb,param,topo,neigh,nlist,ehb,g,sigma)
    !***********************************************************************
    !* Hydrogen bonds with a bound hydrogen, case A...H...B (hblist1).
    !* Input:  n/at/xyz, param/topo/neigh, nlist (provides hblist1)
    !*         mcf_ehb - mcGFN-FF HB scaling factor (1.0 for standard GFN-FF)
    !* In/out: ehb, g, sigma - HB energy, gradient and stress, incremented
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: mcf_ehb
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TGFFNeighbourList),intent(in) :: nlist
    real(wp),intent(inout) :: ehb
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)

    integer :: i,j,k,l,iTri,iTrj
    real(wp) :: etmp,g3tmp(3,3),sig(3,3)

    if (nlist%nhb1 .le. 0) return

    sig = 0.0_wp

    !$omp parallel do default(none) reduction(+:ehb, g, sig) &
    !$omp shared(topo,nlist, neigh, param, n, at, xyz, mcf_ehb) &
    !$omp private(i, j, k, l, iTri, iTrj, etmp, g3tmp)
    !> Atom indices come out of a list, so two iterations can hit the same
    !> entry of the gradient. Vectorising that scatter needs conflict
    !> detection; ifx omits it above -O1 and silently drops contributions.
    !> A comment to compilers that do not know the directive.
!DIR$ NOVECTOR
    do i = 1,nlist%nhb1
      j = nlist%hblist1(1,i)
      k = nlist%hblist1(2,i)
      l = nlist%hblist1(3,i)
      iTri = nlist%hblist1(4,i)
      iTrj = nlist%hblist1(5,i)
      if (iTri .gt. neigh%nTrans.or.iTrj .gt. neigh%nTrans) cycle
      call abhgfnff_eg1(n,j,k,l,iTri,iTrj,at,xyz,topo%qa,etmp,&
              & g3tmp,param,topo,neigh,sig,mcf_ehb)
      g(1:3,j) = g(1:3,j)+g3tmp(1:3,1)*mcf_ehb
      g(1:3,k) = g(1:3,k)+g3tmp(1:3,2)*mcf_ehb
      g(1:3,l) = g(1:3,l)+g3tmp(1:3,3)*mcf_ehb
      ehb = ehb+etmp*mcf_ehb
    end do
    !$omp end parallel do

    if (neigh%nTrans .ne. 1) sigma = sigma+sig   !> stress only wanted for PBC

  end subroutine eg_hbonds_bound

  subroutine eg_hbonds_unbound(n,at,xyz,sqrab,srab,mcf_ehb, &
        & param,topo,neigh,nlist,ehb,g,sigma)
    !***********************************************************************
    !* Hydrogen bonds with an unbound hydrogen (hblist2). The acceptor
    !* environment decides which of the four A...H...B forms is used:
    !* carbonyl and nitro go through abhgfnff_eg3, N hetero aromats through
    !* abhgfnff_eg2_rnr, everything else through abhgfnff_eg2new.
    !* Input:  n/at/xyz, param/topo/neigh, mcf_ehb (see eg_hbonds_bound)
    !*         sqrab/srab - packed squared and plain interatomic distances
    !* In/out: nlist - hbe2 receives the per-bond energies
    !*         ehb, g, sigma - HB energy, gradient and stress, incremented
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: sqrab(n*(n+1)/2),srab(n*(n+1)/2)
    real(wp),intent(in) :: mcf_ehb
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TGFFNeighbourList),intent(inout) :: nlist
    real(wp),intent(inout) :: ehb
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)

    integer :: i,j,k,l,iTr,iTrj,iTrk,iTrDum,nbb,nbk,nbnbk,atnb
    real(wp) :: etmp
    real(wp) :: g5tmp(3,n)  !> automatic, must not be allocatable: it is OMP private
    real(wp) :: sig(3,3)

    if (nlist%nhb2 .le. 0) return

    sig = 0.0_wp

    !$omp parallel do default(none) reduction(+:ehb, g, sig) &
    !$omp shared(topo,nlist, neigh, param, n, at, xyz, sqrab, srab, mcf_ehb) &
    !$omp private(i, j, k, l,iTrj,iTrk,iTrDum,iTr, nbb, nbk, nbnbk, atnb, etmp, g5tmp)
    !> Atom indices come out of a list, so two iterations can hit the same
    !> entry of the gradient. Vectorising that scatter needs conflict
    !> detection; ifx omits it above -O1 and silently drops contributions.
    !> A comment to compilers that do not know the directive.
!DIR$ NOVECTOR
    do i = 1,nlist%nhb2
      j = nlist%hblist2(1,i)  !   A
      k = nlist%hblist2(2,i)  !   B
      l = nlist%hblist2(3,i)  !   H -> always in central cell
      iTrj = nlist%hblist2(4,i) ! iTrA
      iTrk = nlist%hblist2(5,i) ! iTrB

      if (iTrj .gt. neigh%nTrans.or.iTrk .gt. neigh%nTrans) cycle
      !>-- carbonyl/nitro test: terminal O whose neighbour nbk has more neighbours
      nbnbk = 0
      if (at(k) .eq. 8.and.sum(neigh%nb(neigh%numnb,k,:)) .eq. 1) then
        nbk = 0
        iTr = 0 ! nbk is the first neighbor of k !
        atnb = 0
        call neigh%jth_nb(n,xyz,nbk,1,k,iTr)
        !>-- iTrC; cycle if beyond the cutoff of the last getTransVec call
        iTrDum = neigh%fTrSum(iTr,iTrk)
        if (iTrDum .eq. -1.or.iTrDum .gt. neigh%nTrans) cycle
        if (nbk .ne. 0) then
          nbnbk = sum(neigh%nb(neigh%numnb,nbk,:))
          atnb = at(nbk)
        end if
      end if

      !>-- carbonyl case R-C=O...H_A
      if (at(k) .eq. 8.and.sum(neigh%nb(neigh%numnb,k,:)) .eq. 1.and.atnb .eq. 6 &
            & .and.nbnbk .gt. 1) then
        call abhgfnff_eg3(n,j,k,l,iTrj,iTrk,nbk,iTrDum,at,xyz,topo%qa,sqrab,&
                & srab,etmp,g5tmp,param,topo,neigh,sig,mcf_ehb)

        !>-- nitro case R-N=O...H_A
      else if (at(k) .eq. 8.and.sum(neigh%nb(neigh%numnb,k,:)) .eq. 1.and.atnb .eq. 7 &
            &  .and.nbnbk .gt. 1) then
        call abhgfnff_eg3(n,j,k,l,iTrj,iTrk,nbk,iTrDum,at,xyz,topo%qa,sqrab,&
                & srab,etmp,g5tmp,param,topo,neigh,sig,mcf_ehb)

        !>-- N hetero aromat
      else if (at(k) .eq. 7.and.sum(neigh%nb(neigh%numnb,k,:)) .eq. 2) then
        call abhgfnff_eg2_rnr(n,j,k,l,iTrj,iTrk,at,xyz,topo%qa,sqrab,&
                 & srab,etmp,g5tmp,param,topo,neigh,sig,mcf_ehb)

      else
        nbb = sum(neigh%nb(neigh%numnb,k,:))
        call abhgfnff_eg2new(n,j,k,l,iTrj,iTrk,nbb,at,xyz,topo%qa,sqrab,srab, &
           & etmp,g5tmp,param,topo,neigh,sig,mcf_ehb)
      end if
      g = g+g5tmp*mcf_ehb
      ehb = ehb+etmp*mcf_ehb
      nlist%hbe2(i) = etmp

    end do
    !$omp end parallel do

    if (neigh%nTrans .ne. 1) sigma = sigma+sig   !> stress only wanted for PBC

  end subroutine eg_hbonds_unbound

  subroutine eg_xbonds(n,at,xyz,param,topo,neigh,nlist,exb,g,sigma)
    !***********************************************************************
    !* Halogen bonds A...X-B over hblist3.
    !* Input:  n/at/xyz, param/topo/neigh
    !* In/out: nlist - hbe3 receives the per-bond energies
    !*         exb, g, sigma - XB energy, gradient and stress, incremented
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TGFFNeighbourList),intent(inout) :: nlist
    real(wp),intent(inout) :: exb
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)

    integer :: i,j,k,l,iTrk,iTrl
    real(wp) :: etmp,g3tmp(3,3),sig(3,3)

    if (nlist%nxb .le. 0) return

    sig = 0.0_wp

    !$omp parallel do default(none) reduction(+:exb, g, sig) &
    !$omp shared(topo, neigh, nlist, param, n, at, xyz) &
    !$omp private(i, j, k, l, iTrk, iTrl, etmp, g3tmp)
    !> Atom indices come out of a list, so two iterations can hit the same
    !> entry of the gradient. Vectorising that scatter needs conflict
    !> detection; ifx omits it above -O1 and silently drops contributions.
    !> A comment to compilers that do not know the directive.
!DIR$ NOVECTOR
    do i = 1,nlist%nxb
      j = nlist%hblist3(1,i)   ! A in central cell
      k = nlist%hblist3(2,i)   ! B
      l = nlist%hblist3(3,i)   ! X
      iTrk = nlist%hblist3(4,i) !iTrB
      iTrl = nlist%hblist3(5,i) !iTrX
      if (iTrk .gt. neigh%nTrans.or.iTrl .gt. neigh%nTrans) cycle
      if (j .ne. 0.and.k .ne. 0) then
        call rbxgfnff_eg(n,j,k,l,iTrk,iTrl,at,xyz,topo%qa,etmp,g3tmp,param,neigh,sig)
        g(1:3,j) = g(1:3,j)+g3tmp(1:3,1)
        g(1:3,k) = g(1:3,k)+g3tmp(1:3,2)
        g(1:3,l) = g(1:3,l)+g3tmp(1:3,3)
        exb = exb+etmp
        nlist%hbe3(i) = etmp
      end if
    end do
    !$omp end parallel do

    if (neigh%nTrans .ne. 1) sigma = sigma+sig   !> stress only wanted for PBC

  end subroutine eg_xbonds

  subroutine dncoord_erf(nat,at,xyz,rcov,cn,dcn,thr,topo,neigh,dcndL)
    !***********************************************************************
    !* ERF coordination number between HB hydrogens and their acceptors B
    !* (topo%bond_hb_* pairs only), with Cartesian and strain derivatives.
    !* Input:  rcov - covalent radii by atomic number
    !*         thr  - squared distance cutoff, optional (default 1600)
    !* Output: cn, dcn(3,nat,nat), dcndL (derivative w.r.t. strain)
    !***********************************************************************

    implicit none

    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    integer,intent(in)   :: nat
    integer,intent(in)   :: at(nat)
    real(wp),intent(in)  :: xyz(3,nat)
    real(wp),intent(in)  :: rcov(:)
    real(wp),intent(out) :: cn(nat)
    real(wp),intent(out) :: dcn(3,nat,nat)
    real(wp),intent(out) :: dcndL(:,:,:)
    real(wp),intent(in),optional :: thr
    real(wp) :: cn_thr

    integer  :: i,j,iTrB,iTrH
    integer  :: iat,jat
    integer  :: ati,atj
    real(wp) :: r,r2,rij(3),stress(3,3)
    real(wp) :: rcovij
    real(wp) :: dtmp,tmp
    real(wp),parameter :: hlfosqrtpi = 1.0_wp/1.77245385091_wp
    real(wp),parameter :: kn = 27.5_wp
    real(wp),parameter :: rcov_scal = 1.78

    cn = 0._wp
    dcn = 0._wp
    dcndL = 0.0_wp
    cn_thr = 1600.0_wp
    if (present(thr)) cn_thr = thr

    do i = 1,topo%bond_hb_nr
      iat = topo%bond_hb_AH(2,i) ! H atom
      iTrH = topo%bond_hb_AH(4,i)
      ati = at(iat)
      do j = 1,topo%bond_hb_Bn(i)
        jat = topo%bond_hb_B(1,j,i) ! B atom
        iTrB = topo%bond_hb_B(2,j,i)
        atj = at(jat)
        if (iTrB .gt. neigh%nTrans.or.iTrH .gt. neigh%nTrans) cycle
        rij = (xyz(:,jat)+neigh%transVec(:,iTrB))-(xyz(:,iat)+neigh%transVec(:,iTrH))
        r2 = sum(rij**2)
        if (r2 .gt. cn_thr) cycle
        r = sqrt(r2)
        rcovij = rcov_scal*(rcov(ati)+rcov(atj))
        tmp = 0.5_wp*(1.0_wp+erf(-kn*(r-rcovij)/rcovij))
        dtmp = -hlfosqrtpi*kn*exp(-kn**2*(r-rcovij)**2/rcovij**2)/rcovij
        cn(iat) = cn(iat)+tmp
        cn(jat) = cn(jat)+tmp
        dcn(:,jat,jat) = dtmp*rij/r+dcn(:,jat,jat)
        dcn(:,iat,jat) = dtmp*rij/r+dcn(:,iat,jat)
        dcn(:,jat,iat) = -dtmp*rij/r+dcn(:,jat,iat)
        dcn(:,iat,iat) = -dtmp*rij/r+dcn(:,iat,iat)

        stress(:,1) = rij(1)*dtmp*rij/r
        stress(:,2) = rij(2)*dtmp*rij/r
        stress(:,3) = rij(3)*dtmp*rij/r
        dcndL(:,:,iat) = dcndL(:,:,iat)+stress
        if (iat .ne. jat.or.iTrH .ne. iTrB) then
          dcndL(:,:,jat) = dcndL(:,:,jat)+stress
        end if
      end do
    end do

  end subroutine dncoord_erf

  subroutine abhgfnff_eg1(n,A,B,H,iTrA,iTrB,at,xyz,q,energy,gdr,param,topo,neigh,sigma,mcf_ehb)
    !***********************************************************************
    !* HB case 1, A...H...B: energy and gradient of a single bond, with the
    !* basicity/acidity of A and B mixed by rah^4 and rbh^4 weights.
    !* A and B sit in cells iTrA and iTrB, H in the central cell.
    !* gdr(:,1:3) returns the A, B, H gradients unscaled; sigma is
    !* incremented with the mcf_ehb-scaled stress.
    !***********************************************************************
    implicit none
    type(TGFFData),intent(in)     :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout)    :: neigh
    real(wp),intent(inout)        :: sigma(3,3)
    real(wp),intent(in)           :: mcf_ehb
    integer A,B,H,n,at(n),iTrA,iTrB
    real(wp) xyz(3,n),energy,gdr(3,3)
    real(wp) q(n)

    real(wp) outl,dampl,damps,rdamp,damp,dd24a,dd24b
    real(wp) ratio1,ratio2,ratio3
    real(wp) rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
    real(wp) drah(3),drbh(3),drab(3)
    real(wp) dg(3),dga(3),dgb(3),dgh(3)
    real(wp) ga(3),gb(3),gh(3)
    real(wp) gi,denom,tmp,qhoutl,radab,rahprbh
    real(wp) ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo
    real(wp) bas,aci
    real(wp) aterm,rterm,dterm,sterm
    real(wp) qa,qb,qh
    real(wp) ca(2),cb(2)
    real(wp) caa,cbb
    real(wp) shortcut


    gdr = 0
    energy = 0
    call hbonds(A,B,ca,cb,param,topo) ! get HB strength
    rab = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,B)+neigh%transVec(:,iTrB)))
    rab2 = rab**2
    rah = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,H)))
    rah2 = rah**2
    rbh = NORM2((xyz(:,B)+neigh%transVec(:,iTrB))-xyz(:,H))
    rbh2 = rbh**2

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))

    !>-- out-of-line damping
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    rdamp = damp/rab2/rab

    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

    !>-- donor-acceptor term
    rah4 = rah2*rah2
    rbh4 = rbh2*rbh2
    denom = 1.d0/(rah4+rbh4)

    caa = qa*ca(1)
    cbb = qb*cb(1)
    qhoutl = qh*outl

    bas = (caa*rah4+cbb*rbh4)*denom
    aci = (cb(2)*rah4+ca(2)*rbh4)*denom

    rterm = -aci*rdamp*qhoutl
    energy = bas*rterm

    drah(1:3) = xyz(1:3,A)-xyz(1:3,H)+neigh%transVec(1:3,iTrA)
    drbh(1:3) = xyz(1:3,B)-xyz(1:3,H)+neigh%transVec(1:3,iTrB)
    drab(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))

    aterm = -aci*bas*rdamp*qh
    sterm = -rdamp*bas*qhoutl
    dterm = -aci*bas*qhoutl

    tmp = denom*denom*4.0d0
    dd24a = rah2*rbh4*tmp
    dd24b = rbh2*rah4*tmp

    !>-- donor-acceptor part: bas
    gi = (caa-cbb)*dd24a*rterm
    ga(1:3) = gi*drah(1:3)
    gi = (cbb-caa)*dd24b*rterm
    gb(1:3) = gi*drbh(1:3)
    gh(1:3) = -ga(1:3)-gb(1:3)

    !>-- donor-acceptor part: aci
    gi = (cb(2)-ca(2))*dd24a
    dga(1:3) = gi*drah(1:3)*sterm
    ga(1:3) = ga(1:3)+dga(1:3)

    gi = (ca(2)-cb(2))*dd24b
    dgb(1:3) = gi*drbh(1:3)*sterm
    gb(1:3) = gb(1:3)+dgb(1:3)

    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

    !>-- damping part: rab
   gi = rdamp*(-(2.d0*param%hbalp*ratio1/(1+ratio1))+(2.d0*param%hbalp*ratio3/(1+ratio3))-3.d0)/rab2
    dg(1:3) = gi*drab(1:3)*dterm
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

    !>-- out-of-line term: rab, then rah and rbh
    gi = aterm*2.d0*ratio2*expo*rahprbh/(1+ratio2)**2/(rahprbh-rab)/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

    tmp = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    dga(1:3) = drah(1:3)*tmp/rah
    ga(1:3) = ga(1:3)+dga(1:3)
    dgb(1:3) = drbh(1:3)*tmp/rbh
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)
    sigma(:,1) = sigma(:,1)+mcf_ehb*ga(1)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,2) = sigma(:,2)+mcf_ehb*ga(2)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,3) = sigma(:,3)+mcf_ehb*ga(3)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,1) = sigma(:,1)+mcf_ehb*gb(1)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,2) = sigma(:,2)+mcf_ehb*gb(2)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,3) = sigma(:,3)+mcf_ehb*gb(3)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,1) = sigma(:,1)+mcf_ehb*gh(1)*xyz(:,H)
    sigma(:,2) = sigma(:,2)+mcf_ehb*gh(2)*xyz(:,H)
    sigma(:,3) = sigma(:,3)+mcf_ehb*gh(3)*xyz(:,H)
    gdr(1:3,1) = ga(1:3)
    gdr(1:3,2) = gb(1:3)
    gdr(1:3,3) = gh(1:3)

  end subroutine abhgfnff_eg1

  subroutine abhgfnff_eg2new(n,A,B,H,iTrA,iTrB,nbb,at,xyz,q,sqrab, &
                  & srab,energy,gdr,param,topo,neigh,sigma,mcf_ehb)
    !***********************************************************************
    !* HB case 2, A-H...B, including the orientation of the nbb neighbours
    !* of B via an A...nb(B)-B out-of-line damping per neighbour.
    !* gdr(3,n) returns the unscaled gradient on A, B, H and the neighbours
    !* of B; sigma is incremented with the mcf_ehb-scaled stress.
    !***********************************************************************
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    real(wp),intent(inout)        :: sigma(3,3)
    real(wp),intent(in)           :: mcf_ehb
    integer A,B,H,n,at(n),iTrA,iTrB,nbb
    real(wp) xyz(3,n),energy,gdr(3,n)
    real(wp) q(n)
    real(wp) sqrab(n*(n+1)/2)   ! squared dist
    real(wp) srab(n*(n+1)/2)    ! dist

    real(wp) outl,dampl,damps,rdamp,damp
    real(wp) ddamp,rabdamp,rbhdamp
    real(wp) ratio1,ratio2,ratio2_nb(nbb),ratio3
    real(wp) rab,rah,rbh,rab2,rah2,rbh2
    real(wp) ranb(nbb),ranb2(nbb),rbnb(nbb),rbnb2(nbb)
    real(wp) drah(3),drbh(3),drab(3)
    real(wp) dranb(3,nbb),drbnb(3,nbb)
    real(wp) dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
    real(wp) ga(3),gb(3),gh(3),gnb(3,nbb)
    real(wp) qhoutl,radab
    real(wp) gi,gi_nb(nbb)
    real(wp) tmp1,tmp2(nbb)
    real(wp) rahprbh,ranbprbnb(nbb)
    real(wp) ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb(nbb)
    real(wp) aterm,dterm,nbterm
    real(wp) qa,qb,qh
    real(wp) ca(2),cb(2)
    real(wp) shortcut
    real(wp) const
    real(wp) outl_nb(nbb),outl_nb_tot
    real(wp) hbnbcut_save
    real(wp) vecDum(3)
    logical mask_nb(nbb)

    !>-- proportion between the rbh and rab distance dependencies
    real(wp) :: p_bh
    real(wp) :: p_ab

    integer i,inb,iTr
    p_bh = 1.d0+param%hbabmix
    p_ab = -param%hbabmix

    gdr = 0
    energy = 0

    call hbonds(A,B,ca,cb,param,topo)
    do i = 1,nbb
      inb = 0; iTr = 0 ! jth_nb output
      call neigh%jth_nb(n,xyz,inb,i,B,iTr) ! inb is the i-th nb of B when inb is shifted to iTr
      vecDum = neigh%transVec(:,iTr)+neigh%transVec(:,iTrB)
      dranb(1:3,i) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,inb)+vecDum)
      drbnb(1:3,i) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-(xyz(1:3,inb)+vecDum)
      ranb2(i) = sum(dranb(1:3,i)**2)
      ranb(i) = sqrt(ranb2(i))
      rbnb2(i) = sum(drbnb(1:3,i)**2)
      rbnb(i) = sqrt(rbnb2(i))
    end do

    rab = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,B)+neigh%transVec(:,iTrB)))
    rab2 = rab**2
    rah = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,H)))
    rah2 = rah**2
    rbh = NORM2((xyz(:,B)+neigh%transVec(:,iTrB))-xyz(:,H))
    rbh2 = rbh**2

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))

    !>-- out-of-line damping: A-H...B
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

    !>-- out-of-line damping: A...nb(B)-B
    if (at(B) .eq. 7.and.nbb .eq. 1) then
      hbnbcut_save = 2.0
    else
      hbnbcut_save = param%hbnbcut
    end if
    do i = 1,nbb
      ranbprbnb(i) = ranb(i)+rbnb(i)+1.d-12
      expo_nb(i) = (hbnbcut_save/radab)*(ranbprbnb(i)/rab-1.d0)
      ratio2_nb(i) = exp(-expo_nb(i))**(1.0)
      outl_nb(i) = (2.d0/(1.d0+ratio2_nb(i)))-1.0d0
    end do
    outl_nb_tot = product(outl_nb)

    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
    rbhdamp = damp*((p_bh/rbh2/rbh))
    rabdamp = damp*((p_ab/rab2/rab))
    rdamp = rbhdamp+rabdamp

    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

    qhoutl = qh*outl*outl_nb_tot

    const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh
    energy = -rdamp*qhoutl*const
    drah(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-xyz(1:3,H)
    drbh(1:3) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-xyz(1:3,H)
    drab(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))

    aterm = -rdamp*qh*outl_nb_tot*const
    nbterm = -rdamp*qh*outl*const
    dterm = -qhoutl*const

    !>-- damping part: rab
    gi = ((rabdamp+rbhdamp)*ddamp-3.d0*rabdamp)/rab2
    gi = gi*dterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = dg(1:3)
    gb(1:3) = -dg(1:3)

    !>-- damping part: rbh
    gi = -3.d0*rbhdamp/rbh2
    gi = gi*dterm
    dg(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dg(1:3)
    gh(1:3) = -dg(1:3)

    !>-- angular A-H...B term: rab, then rah and rbh
    tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    gi = -tmp1*rahprbh/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

    gi = tmp1/rah
    dga(1:3) = gi*drah(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp1/rbh
    dgb(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

    !>-- angular A...nb(B)-B term: rab, then ranb and rbnb
    mask_nb = .true.
    do i = 1,nbb
      mask_nb(i) = .false.
      tmp2(i) = 2.d0*nbterm*product(outl_nb,mask_nb)*ratio2_nb(i)*expo_nb(i)/&
               & (1+ratio2_nb(i))**2/(ranbprbnb(i)-rab)
      gi_nb(i) = -tmp2(i)*ranbprbnb(i)/rab2
      dg(1:3) = gi_nb(i)*drab(1:3)
      ga(1:3) = ga(1:3)+dg(1:3)
      gb(1:3) = gb(1:3)-dg(1:3)
      mask_nb = .true.
    end do

    do i = 1,nbb
      gi_nb(i) = tmp2(i)/ranb(i)
      dga(1:3) = gi_nb(i)*dranb(1:3,i)
      ga(1:3) = ga(1:3)+dga(1:3)
      gi_nb(i) = tmp2(i)/rbnb(i)
      dgb(1:3) = gi_nb(i)*drbnb(1:3,i)
      gb(1:3) = gb(1:3)+dgb(1:3)
      dgnb(1:3) = -dga(1:3)-dgb(1:3)
      gnb(1:3,i) = dgnb(1:3)
    end do

    if (nbb .lt. 1) then
      gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
      gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
      gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
      sigma(:,1) = sigma(:,1)+mcf_ehb*ga(1)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
      sigma(:,2) = sigma(:,2)+mcf_ehb*ga(2)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
      sigma(:,3) = sigma(:,3)+mcf_ehb*ga(3)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
      sigma(:,1) = sigma(:,1)+mcf_ehb*gb(1)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
      sigma(:,2) = sigma(:,2)+mcf_ehb*gb(2)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
      sigma(:,3) = sigma(:,3)+mcf_ehb*gb(3)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
      sigma(:,1) = sigma(:,1)+mcf_ehb*gh(1)*xyz(:,H)
      sigma(:,2) = sigma(:,2)+mcf_ehb*gh(2)*xyz(:,H)
      sigma(:,3) = sigma(:,3)+mcf_ehb*gh(3)*xyz(:,H)
      return
    end if

    gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
    gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
    gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
    do i = 1,nbb
      inb = 0; iTr = 0 ! jth_nb output
      call neigh%jth_nb(n,xyz,inb,i,B,iTr) ! inb is the i-th nb of B when shifted to iTr
      gdr(1:3,inb) = gdr(1:3,inb)+gnb(1:3,i)
    end do

    sigma(:,1) = sigma(:,1)+mcf_ehb*ga(1)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,2) = sigma(:,2)+mcf_ehb*ga(2)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,3) = sigma(:,3)+mcf_ehb*ga(3)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,1) = sigma(:,1)+mcf_ehb*gb(1)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,2) = sigma(:,2)+mcf_ehb*gb(2)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,3) = sigma(:,3)+mcf_ehb*gb(3)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,1) = sigma(:,1)+mcf_ehb*gh(1)*xyz(:,H)
    sigma(:,2) = sigma(:,2)+mcf_ehb*gh(2)*xyz(:,H)
    sigma(:,3) = sigma(:,3)+mcf_ehb*gh(3)*xyz(:,H)
    do i = 1,nbb
      inb = 0; iTr = 0 ! jth_nb output
      call neigh%jth_nb(n,xyz,inb,i,B,iTr) ! inb is the i-th nb of B when shifted to iTr
      vecDum = neigh%transVec(:,iTr)+neigh%transVec(:,iTrB)
      sigma(:,1) = sigma(:,1)+mcf_ehb*gnb(1,i)*(xyz(:,inb)+vecDum)
      sigma(:,2) = sigma(:,2)+mcf_ehb*gnb(2,i)*(xyz(:,inb)+vecDum)
      sigma(:,3) = sigma(:,3)+mcf_ehb*gnb(3,i)*(xyz(:,inb)+vecDum)
    end do

  end subroutine abhgfnff_eg2new

  subroutine abhgfnff_eg2_rnr(n,A,B,H,iTrA,iTrB,at,xyz,q,sqrab,srab,energy,gdr,param,topo,neigh,sigma,mcf_ehb)
    !***********************************************************************
    !* HB case 2 for N hetero aromats, A-H...B with B carrying two
    !* neighbours: as abhgfnff_eg2new plus an out-of-line damping towards a
    !* lone pair placed at lp_dist from B, opposite the neighbour sum vector.
    !* gdr and sigma as in abhgfnff_eg2new.
    !***********************************************************************
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    real(wp),intent(inout)        :: sigma(3,3)
    real(wp),intent(in)           :: mcf_ehb
    integer A,B,H,iTrA,iTrB,n,at(n)
    real(wp) xyz(3,n),energy,gdr(3,n)
    real(wp) q(n)
    real(wp) sqrab(n*(n+1)/2)   ! squared dist
    real(wp) srab(n*(n+1)/2)    ! dist

    real(wp) outl,dampl,damps,rdamp,damp
    real(wp) ddamp,rabdamp,rbhdamp
    real(wp) ratio1,ratio2,ratio2_lp,ratio2_nb(22),ratio3
    real(wp) rab,rah,rbh,rab2,rah2,rbh2
    real(wp) ranb(2),ranb2(2),rbnb(2),rbnb2(2)
    real(wp) drah(3),drbh(3),drab(3),dralp(3)
    real(wp) dranb(3,2),drbnb(3,2)
    real(wp) dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
    real(wp) ga(3),gb(3),gh(3),gnb(3,2),gnb_lp(3),glp(3)
    real(wp) qhoutl,radab
    real(wp) gi,gi_nb(2)
    real(wp) tmp1,tmp2(2),tmp3
    real(wp) rahprbh,ranbprbnb(2)
    real(wp) ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_lp,expo_nb(2)
    real(wp) aterm,dterm,nbterm,lpterm
    real(wp) qa,qb,qh
    real(wp) ca(2),cb(2)
    real(wp) shortcut
    real(wp) const
    real(wp) outl_nb(2),outl_nb_tot,outl_lp
    real(wp) vector(3),vnorm
    real(wp) gii(3,3)
    real(wp) unit_vec(3)
    real(wp) drnb(3,2)
    real(wp) lp(3)   !lonepair position
    real(wp) lp_dist !distance parameter between B and lonepair
    real(wp) ralp,ralp2,rblp,rblp2,ralpprblp
    logical mask_nb(2)
    logical lp_on_b

    !>-- proportion between the rbh and rab distance dependencies
    real(wp) :: p_bh
    real(wp) :: p_ab
    !>-- cutoff of the lone-pair out-of-line damping
    real(wp) hblpcut

    real(wp) :: vTrinb(3)
    integer i,nbb,inb,iTr

    p_bh = 1.d0+param%hbabmix
    p_ab = -param%hbabmix

    gdr = 0
    energy = 0
    vector = 0
    lp_dist = 0.50-0.018*param%repz(at(B))
    hblpcut = 56

    call hbonds(A,B,ca,cb,param,topo)

    nbb = 2 ! given through if condition before call
    do i = 1,nbb
      inb = 0; iTr = 0 ! jth_nb output
      call neigh%jth_nb(n,xyz,inb,i,B,iTr) ! inb is the i-th nb of B when shifted to iTr
      vTrinb = neigh%transVec(:,iTr)+neigh%transVec(:,iTrB)
      dranb(1:3,i) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,inb)+vTrinb)
      drbnb(1:3,i) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-(xyz(1:3,inb)+vTrinb)
      ranb2(i) = sum(dranb(1:3,i)**2)
      ranb(i) = sqrt(ranb2(i))
      rbnb2(i) = sum(drbnb(1:3,i)**2)
      rbnb(i) = sqrt(rbnb2(i))

      drnb(1:3,i) = (xyz(1:3,inb)+vTrinb)-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))
      vector = vector+drnb(1:3,i)
    end do

    vnorm = norm2(vector)
    !>-- lone pair position; when the neighbour vectors cancel it has no direction
    !>   and sits on B, the neighbour factors are kept
    lp_on_b = vnorm .le. 1.d-10
    if (.not.lp_on_b) then
      lp = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-lp_dist*(vector/vnorm)
    else
      lp = xyz(1:3,B)+neigh%transVec(1:3,iTrB)
    end if

    rab = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,B)+neigh%transVec(:,iTrB)))
    rab2 = rab**2
    rah = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,H)))
    rah2 = rah**2
    rbh = NORM2((xyz(:,B)+neigh%transVec(:,iTrB))-xyz(:,H))
    rbh2 = rbh**2

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))
    !>-- out-of-line damping: A-H...B
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

    !>-- out-of-line damping: A...LP-B
    rblp2 = sum(((xyz(1:3,B)+neigh%transVec(1:3,iTrB))-lp(1:3))**2)
    rblp = sqrt(rblp2)
    ralp2 = sum(((xyz(1:3,A)+neigh%transVec(1:3,iTrA))-lp(1:3))**2)
    ralp = sqrt(ralp2)
    ralpprblp = ralp+rblp+1.d-12
    expo_lp = (hblpcut/radab)*(ralpprblp/rab-1.d0)
    ratio2_lp = exp(expo_lp)
    outl_lp = 2.d0/(1.d0+ratio2_lp)

    !>-- out-of-line damping: A...nb(B)-B
    do i = 1,nbb
      ranbprbnb(i) = ranb(i)+rbnb(i)+1.d-12
      expo_nb(i) = (param%hbnbcut/radab)*(ranbprbnb(i)/rab-1.d0)
      ratio2_nb(i) = exp(-expo_nb(i))**(1.0)
      outl_nb(i) = (2.d0/(1.d0+ratio2_nb(i)))-1.0d0
    end do
    outl_nb_tot = product(outl_nb)

    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
    rbhdamp = damp*((p_bh/rbh2/rbh))
    rabdamp = damp*((p_ab/rab2/rab))
    rdamp = rbhdamp+rabdamp

    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

    qhoutl = qh*outl*outl_nb_tot*outl_lp

    const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

    energy = -rdamp*qhoutl*const
    drah(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-xyz(1:3,H)
    drbh(1:3) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-xyz(1:3,H)
    drab(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))
    dralp(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-lp(1:3)

    aterm = -rdamp*qh*outl_nb_tot*outl_lp*const
    nbterm = -rdamp*qh*outl*outl_lp*const
    lpterm = -rdamp*qh*outl*outl_nb_tot*const
    dterm = -qhoutl*const

    !>-- damping part: rab
    gi = ((rabdamp+rbhdamp)*ddamp-3.d0*rabdamp)/rab2
    gi = gi*dterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = dg(1:3)
    gb(1:3) = -dg(1:3)

    !>-- damping part: rbh
    gi = -3.d0*rbhdamp/rbh2
    gi = gi*dterm
    dg(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dg(1:3)
    gh(1:3) = -dg(1:3)

    !>-- angular A-H...B term: rab, then rah and rbh
    tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    gi = -tmp1*rahprbh/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

    gi = tmp1/rah
    dga(1:3) = gi*drah(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp1/rbh
    dgb(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

    !>-- angular A...LP-B term: rab, then ralp and rblp
    tmp3 = -2.d0*lpterm*ratio2_lp*expo_lp/(1+ratio2_lp)**2/(ralpprblp-rab)
    gi = -tmp3*ralpprblp/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

    gi = tmp3/ralp
    dga(1:3) = gi*dralp(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gb(1:3) = gb(1:3)-dga(1:3)
    glp(1:3) = -dga(1:3)!-dgb(1:3)

    !>-- LP gradient passed on to B and its neighbours; none if the LP sits on B
    gnb_lp = 0.0_wp
    if (.not.lp_on_b) then
      unit_vec = 0
      do i = 1,3
        unit_vec(i) = -1
        gii(1:3,i) = -lp_dist*dble(nbb)*(unit_vec/vnorm+(vector*vector(i)/sum(vector**2)**(1.5d0)))
        unit_vec = 0
      end do
      gnb_lp = matmul(gii,glp)
    end if

    !>-- angular A...nb(B)-B term: rab, then ranb and rbnb
    mask_nb = .true.
    do i = 1,nbb
      mask_nb(i) = .false.
      tmp2(i) = 2.d0*nbterm*product(outl_nb,mask_nb)*ratio2_nb(i)*expo_nb(i)/&
               & (1+ratio2_nb(i))**2/(ranbprbnb(i)-rab)
      gi_nb(i) = -tmp2(i)*ranbprbnb(i)/rab2
      dg(1:3) = gi_nb(i)*drab(1:3)
      ga(1:3) = ga(1:3)+dg(1:3)
      gb(1:3) = gb(1:3)-dg(1:3)
      mask_nb = .true.
    end do

    do i = 1,nbb
      gi_nb(i) = tmp2(i)/ranb(i)
      dga(1:3) = gi_nb(i)*dranb(1:3,i)
      ga(1:3) = ga(1:3)+dga(1:3)
      gi_nb(i) = tmp2(i)/rbnb(i)
      dgb(1:3) = gi_nb(i)*drbnb(1:3,i)
      gb(1:3) = gb(1:3)+dgb(1:3)
      dgnb(1:3) = -dga(1:3)-dgb(1:3)
      gnb(1:3,i) = dgnb(1:3)
    end do

    gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
    gdr(1:3,B) = gdr(1:3,B)+gb(1:3)+gnb_lp(1:3)
    gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
    do i = 1,nbb
      inb = 0; iTr = 0 ! jth_nb output
      call neigh%jth_nb(n,xyz,inb,i,B,iTr) ! inb is the i-th nb of B when shifted to iTr
      gdr(1:3,inb) = gdr(1:3,inb)+gnb(1:3,i)-gnb_lp(1:3)/dble(nbb)
    end do

    sigma(:,1) = sigma(:,1)+mcf_ehb*ga(1)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,2) = sigma(:,2)+mcf_ehb*ga(2)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,3) = sigma(:,3)+mcf_ehb*ga(3)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,1) = sigma(:,1)+mcf_ehb*gb(1)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,2) = sigma(:,2)+mcf_ehb*gb(2)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,3) = sigma(:,3)+mcf_ehb*gb(3)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,1) = sigma(:,1)+mcf_ehb*gh(1)*xyz(:,H)
    sigma(:,2) = sigma(:,2)+mcf_ehb*gh(2)*xyz(:,H)
    sigma(:,3) = sigma(:,3)+mcf_ehb*gh(3)*xyz(:,H)
    sigma(:,1) = sigma(:,1)+mcf_ehb*gnb_lp(1)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,2) = sigma(:,2)+mcf_ehb*gnb_lp(2)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,3) = sigma(:,3)+mcf_ehb*gnb_lp(3)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    gnb_lp = gnb_lp/dble(nbb)
    do i = 1,nbb
      inb = 0; iTr = 0 ! jth_nb output
      call neigh%jth_nb(n,xyz,inb,i,B,iTr) ! inb is the i-th nb of B when shifted to iTr
      vTrinb = neigh%transVec(:,iTr)+neigh%transVec(:,iTrB)
      sigma(:,1) = sigma(:,1)+mcf_ehb*(gnb(1,i)-gnb_lp(1))*(xyz(:,inb)+vTrinb)
      sigma(:,2) = sigma(:,2)+mcf_ehb*(gnb(2,i)-gnb_lp(2))*(xyz(:,inb)+vTrinb)
      sigma(:,3) = sigma(:,3)+mcf_ehb*(gnb(3,i)-gnb_lp(3))*(xyz(:,inb)+vTrinb)
    end do

  end subroutine abhgfnff_eg2_rnr

  subroutine abhgfnff_eg3(n,A,B,H,iTrA,iTrB,C,iTrC,at,xyz,q,sqrab,srab,energy,&
                  & gdr,param,topo,neigh,sigma,mcf_ehb)
    !***********************************************************************
    !* HB case 3, A-H...B with B the O of a carbonyl or nitro group and C
    !* its only neighbour (in cell iTrC); accounts for two in-plane LPs at B.
    !* Multiplicative form: an abhgfnff_eg2new-type energy times a torsion
    !* factor etors (R-C=O...H) and a bend factor eangl (C=O...H, 120 deg).
    !* gdr and sigma as in abhgfnff_eg2new.
    !***********************************************************************
    implicit none
    type(TGFFData),intent(in)     :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout)    :: neigh
    real(wp),intent(inout)        :: sigma(3,3)
    real(wp),intent(in)           :: mcf_ehb
    integer  :: A,B,H,iTrA,iTrB,n,at(n),C,iTrC
    real(wp) :: xyz(3,n),energy,gdr(3,n)
    real(wp) :: q(n)
    real(wp) :: sqrab(n*(n+1)/2)   ! squared dist
    real(wp) :: srab(n*(n+1)/2)    ! dist

    real(wp) :: outl,dampl,damps,rdamp,damp
    real(wp) :: ddamp,rabdamp,rbhdamp
    real(wp) :: ratio1,ratio2,ratio2_nb,ratio3
    real(wp) :: rab,rah,rbh,rab2,rah2,rbh2
    real(wp) :: vTrR(3),vTrB(3),vTrC(3)
    real(wp) :: ranb,ranb2,rbnb,rbnb2
    real(wp) :: drah(3),drbh(3),drab(3)
    real(wp) :: dranb(3),drbnb(3)
    real(wp) :: dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
    real(wp) :: ga(3),gb(3),gh(3),gnb(3)
    real(wp) :: phi,phi0,r0,fc,tshift,bshift
    real(wp) :: eangl,etors,gangl(3,n),gtors(3,n)
    real(wp) :: etmp(20),g3tmp(3,3),g4tmp(3,4,20)
    real(wp) :: qhoutl,radab
    real(wp) :: gi,gi_nb
    real(wp) :: tmp1,tmp2
    real(wp) :: rahprbh,ranbprbnb
    real(wp) :: ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb
    real(wp) :: aterm,dterm,nbterm,bterm,tterm
    real(wp) :: qa,qb,qh
    real(wp) :: ca(2),cb(2)
    real(wp) :: shortcut
    integer :: tlist(6,sum(neigh%nb(neigh%numnb,C,:)))
    real(wp) :: vtors(2,sum(neigh%nb(neigh%numnb,C,:)))
    real(wp) :: const
    real(wp) :: outl_nb_tot
    logical :: t_mask(20)

    !>-- proportion between the rbh and rab distance dependencies
    real(wp) :: p_bh
    real(wp) :: p_ab

    integer :: i,j,ii,jj,kk,ll,iTr,iTrR
    integer :: nbb,nbc
    integer :: ntors,rn

    p_bh = 1.d0+param%hbabmix
    p_ab = -param%hbabmix

    gdr = 0
    energy = 0
    etors = 0
    gtors = 0
    eangl = 0
    gangl = 0
    call hbonds(A,B,ca,cb,param,topo)

    nbb = 1 ! routine only called for nbb.eq.1
    nbc = sum(neigh%nb(neigh%numnb,C,:))
    ntors = nbc-nbb

    dranb(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,C)+neigh%transVec(1:3,iTrC))
    drbnb(1:3) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-(xyz(1:3,C)+neigh%transVec(1:3,iTrC))
    ranb2 = sum(dranb(1:3)**2)
    ranb = sqrt(ranb2)
    rbnb2 = sum(drbnb(1:3)**2)
    rbnb = sqrt(rbnb2)
    rab = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,B)+neigh%transVec(:,iTrB)))
    rab2 = rab**2
    rah = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,H)))
    rah2 = rah**2
    rbh = NORM2((xyz(:,B)+neigh%transVec(:,iTrB))-xyz(:,H))
    rbh2 = rbh**2

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))
    !>-- out-of-line damping: A-H...B
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

    !>-- out-of-line damping: A...nb(B)-B
    ranbprbnb = ranb+rbnb+1.d-12
    expo_nb = (param%hbnbcut/radab)*(ranbprbnb/rab-1.d0)
    ratio2_nb = exp(-expo_nb)
    outl_nb_tot = (2.d0/(1.d0+ratio2_nb))-1.0d0

    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
    rbhdamp = damp*((p_bh/rbh2/rbh))
    rabdamp = damp*((p_ab/rab2/rab))
    rdamp = rbhdamp+rabdamp
    !>-- torsion list: one entry per neighbour R of C other than B, atoms
    !>   ordered R, B(=O), C, H as ii, jj, kk, ll; tlist(6,:) is the cell of R
    j = 0
    do iTr = 1,neigh%numctr
      do i = 1,neigh%nb(neigh%numnb,C,iTr)
        if (neigh%nb(i,C,iTr) .eq. B) cycle
        j = j+1
        tlist(1,j) = neigh%nb(i,C,iTr) ! R
        tlist(2,j) = B
        tlist(3,j) = C
        tlist(4,j) = H
        tlist(5,j) = 2
        tlist(6,j) = neigh%fTrSum(iTr,iTrC)
        !>-- R outside the cutoff: vtors stays unset, the loops below skip it too
        if (tlist(6,j) .le. 0.or.tlist(6,j) .gt. neigh%nTrans) cycle
        vtors(1,j) = pi/2.0
        vtors(2,j) = param%tors_hb
      end do
    end do
    vTrB = neigh%transVec(:,iTrB)
    vTrC = neigh%transVec(:,iTrC)
    do i = 1,ntors
      ii = tlist(1,i) ! R
      jj = tlist(2,i) ! B
      kk = tlist(3,i) ! C
      ll = tlist(4,i) ! H
      rn = tlist(5,i)
      iTrR = tlist(6,i)
      if (iTrR .le. 0.or.iTrR .gt. neigh%nTrans) then
        g4tmp(:,:,i) = 0.0_wp
        etmp(i) = 0.0_wp
        cycle
      end if
      phi0 = vtors(1,i)
      tshift = vtors(2,i)
      vTrR = neigh%transVec(:,iTrR)

      phi = valijklffPBC(2,n,xyz,ii,jj,kk,ll,vTrR,vTrB,vTrC)
      call egtors_nci_mul(ii,jj,kk,ll,vTrR,vTrB,vTrC,rn,phi,phi0, &
                        & tshift,n,at,xyz,etmp(i),g4tmp(:,:,i))
    end do
    etors = product(etmp(1:ntors))
    !>-- product rule: gradient of factor i times all other torsion factors
    t_mask = .true.
    do i = 1,ntors
      t_mask(i) = .false.
      ii = tlist(1,i)
      jj = tlist(2,i)
      kk = tlist(3,i)
      ll = tlist(4,i)
      gtors(1:3,ii) = gtors(1:3,ii)+g4tmp(1:3,1,i)*product(etmp(1:ntors),t_mask(1:ntors))
      gtors(1:3,jj) = gtors(1:3,jj)+g4tmp(1:3,2,i)*product(etmp(1:ntors),t_mask(1:ntors))
      gtors(1:3,kk) = gtors(1:3,kk)+g4tmp(1:3,3,i)*product(etmp(1:ntors),t_mask(1:ntors))
      gtors(1:3,ll) = gtors(1:3,ll)+g4tmp(1:3,4,i)*product(etmp(1:ntors),t_mask(1:ntors))
      t_mask = .true.
    end do

    r0 = 120
    phi0 = r0*pi/180.
    bshift = param%bend_hb
    fc = 1.0d0-bshift
    call egbend_nci_mul(jj,kk,ll,vTrB,vTrC,phi0,fc,n,at,xyz,eangl,g3tmp)
    gangl(1:3,jj) = gangl(1:3,jj)+g3tmp(1:3,1)
    gangl(1:3,kk) = gangl(1:3,kk)+g3tmp(1:3,2)
    gangl(1:3,ll) = gangl(1:3,ll)+g3tmp(1:3,3)

    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

    qhoutl = qh*outl*outl_nb_tot

    const = ca(2)*qa*cb(1)*qb*param%xhaci_coh
    energy = -rdamp*qhoutl*eangl*etors*const
    drah(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-xyz(1:3,H)
    drbh(1:3) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-xyz(1:3,H)
    drab(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))

    aterm = -rdamp*qh*outl_nb_tot*eangl*etors*const
    nbterm = -rdamp*qh*outl*eangl*etors*const
    dterm = -qhoutl*eangl*etors*const
    tterm = -rdamp*qhoutl*eangl*const
    bterm = -rdamp*qhoutl*etors*const

    !>-- damping part: rab
    gi = ((rabdamp+rbhdamp)*ddamp-3.d0*rabdamp)/rab2
    gi = gi*dterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = dg(1:3)
    gb(1:3) = -dg(1:3)

    !>-- damping part: rbh
    gi = -3.d0*rbhdamp/rbh2
    gi = gi*dterm
    dg(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dg(1:3)
    gh(1:3) = -dg(1:3)

    !>-- angular A-H...B term: rab, then rah and rbh
    tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    gi = -tmp1*rahprbh/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

    gi = tmp1/rah
    dga(1:3) = gi*drah(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp1/rbh
    dgb(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

    !>-- angular A...nb(B)-B term: rab, then ranb and rbnb
    tmp2 = 2.d0*nbterm*ratio2_nb*expo_nb/(1+ratio2_nb)**2/(ranbprbnb-rab)
    gi_nb = -tmp2*ranbprbnb/rab2
    dg(1:3) = gi_nb*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

    gi_nb = tmp2/ranb
    dga(1:3) = gi_nb*dranb(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi_nb = tmp2/rbnb
    dgb(1:3) = gi_nb*drbnb(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgnb(1:3) = -dga(1:3)-dgb(1:3)
    gnb(1:3) = dgnb(1:3)

    do i = 1,ntors
      ii = tlist(1,i)
      gdr(1:3,ii) = gdr(1:3,ii)+gtors(1:3,ii)*tterm
    end do
    gdr(1:3,jj) = gdr(1:3,jj)+gtors(1:3,jj)*tterm
    gdr(1:3,kk) = gdr(1:3,kk)+gtors(1:3,kk)*tterm
    gdr(1:3,ll) = gdr(1:3,ll)+gtors(1:3,ll)*tterm

    gdr(1:3,jj) = gdr(1:3,jj)+gangl(1:3,jj)*bterm
    gdr(1:3,kk) = gdr(1:3,kk)+gangl(1:3,kk)*bterm
    gdr(1:3,ll) = gdr(1:3,ll)+gangl(1:3,ll)*bterm
    sigma(:,1) = sigma(:,1)+mcf_ehb*ga(1)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,2) = sigma(:,2)+mcf_ehb*ga(2)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,3) = sigma(:,3)+mcf_ehb*ga(3)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,1) = sigma(:,1)+mcf_ehb*gb(1)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,2) = sigma(:,2)+mcf_ehb*gb(2)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,3) = sigma(:,3)+mcf_ehb*gb(3)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,1) = sigma(:,1)+mcf_ehb*gh(1)*xyz(:,H)
    sigma(:,2) = sigma(:,2)+mcf_ehb*gh(2)*xyz(:,H)
    sigma(:,3) = sigma(:,3)+mcf_ehb*gh(3)*xyz(:,H)
    sigma(:,1) = sigma(:,1)+mcf_ehb*gnb(1)*(xyz(:,C)+neigh%transVec(1:3,iTrC))
    sigma(:,2) = sigma(:,2)+mcf_ehb*gnb(2)*(xyz(:,C)+neigh%transVec(1:3,iTrC))
    sigma(:,3) = sigma(:,3)+mcf_ehb*gnb(3)*(xyz(:,C)+neigh%transVec(1:3,iTrC))
    do i = 1,ntors
      ii = tlist(1,i)
      jj = tlist(2,i)
      kk = tlist(3,i)
      ll = tlist(4,i)
      iTrR = tlist(6,i)
      if (iTrR .le. 0.or.iTrR .gt. neigh%nTrans) then
        cycle
      end if
      sigma(:,1) = sigma(:,1)+mcf_ehb*tterm*gtors(1,ii)* &        ! R
                & (xyz(1:3,ii)+neigh%transVec(1:3,iTrR))
      sigma(:,2) = sigma(:,2)+mcf_ehb*tterm*gtors(2,ii)* &        ! R
                & (xyz(1:3,ii)+neigh%transVec(1:3,iTrR))
      sigma(:,3) = sigma(:,3)+mcf_ehb*tterm*gtors(3,ii)* &        ! R
                & (xyz(1:3,ii)+neigh%transVec(1:3,iTrR))
    end do
    !>-- jj, kk and ll are the same for every i above, only ii (R) changes
    sigma(:,1) = sigma(:,1)+mcf_ehb*tterm*gtors(1,jj)* &           ! B
              & (xyz(:,jj)+neigh%transVec(:,iTrB))
    sigma(:,2) = sigma(:,2)+mcf_ehb*tterm*gtors(2,jj)* &           ! B
              & (xyz(:,jj)+neigh%transVec(:,iTrB))
    sigma(:,3) = sigma(:,3)+mcf_ehb*tterm*gtors(3,jj)* &           ! B
              & (xyz(:,jj)+neigh%transVec(:,iTrB))
    sigma(:,1) = sigma(:,1)+mcf_ehb*tterm*gtors(1,kk)* &           ! C
              & (xyz(:,kk)+neigh%transVec(:,iTrC))
    sigma(:,2) = sigma(:,2)+mcf_ehb*tterm*gtors(2,kk)* &           ! C
              & (xyz(:,kk)+neigh%transVec(:,iTrC))
    sigma(:,3) = sigma(:,3)+mcf_ehb*tterm*gtors(3,kk)* &           ! C
              & (xyz(:,kk)+neigh%transVec(:,iTrC))
    sigma(:,1) = sigma(:,1)+mcf_ehb*tterm*gtors(1,ll)* &           ! H
              & xyz(:,ll)
    sigma(:,2) = sigma(:,2)+mcf_ehb*tterm*gtors(2,ll)* &           ! H
              & xyz(:,ll)
    sigma(:,3) = sigma(:,3)+mcf_ehb*tterm*gtors(3,ll)* &           ! H
              & xyz(:,ll)
    sigma(:,1) = sigma(:,1)+mcf_ehb*bterm*gangl(1,jj)* &           ! B
              & (xyz(:,jj)+neigh%transVec(:,iTrB))
    sigma(:,2) = sigma(:,2)+mcf_ehb*bterm*gangl(2,jj)* &           ! B
              & (xyz(:,jj)+neigh%transVec(:,iTrB))
    sigma(:,3) = sigma(:,3)+mcf_ehb*bterm*gangl(3,jj)* &           ! B
              & (xyz(:,jj)+neigh%transVec(:,iTrB))
    sigma(:,1) = sigma(:,1)+mcf_ehb*bterm*gangl(1,kk)* &           ! C
              & (xyz(:,kk)+neigh%transVec(:,iTrC))
    sigma(:,2) = sigma(:,2)+mcf_ehb*bterm*gangl(2,kk)* &           ! C
              & (xyz(:,kk)+neigh%transVec(:,iTrC))
    sigma(:,3) = sigma(:,3)+mcf_ehb*bterm*gangl(3,kk)* &           ! C
              & (xyz(:,kk)+neigh%transVec(:,iTrC))
    sigma(:,1) = sigma(:,1)+mcf_ehb*bterm*gangl(1,ll)* &           ! H
              & xyz(:,ll)
    sigma(:,2) = sigma(:,2)+mcf_ehb*bterm*gangl(2,ll)* &           ! H
              & xyz(:,ll)
    sigma(:,3) = sigma(:,3)+mcf_ehb*bterm*gangl(3,ll)* &           ! H
              & xyz(:,ll)

    gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
    gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
    gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
    gdr(1:3,C) = gdr(1:3,C)+gnb(1:3)

  end subroutine abhgfnff_eg3

  subroutine rbxgfnff_eg(n,A,B,X,iTrB,iTrX,at,xyz,q,energy,gdr,param,neigh,sigma)
    !***********************************************************************
    !* Halogen bond A...X-B: energy and gradient of a single bond.
    !* A is in the central cell, B and X in cells iTrB and iTrX.
    !* gdr(:,1:3) returns the A, B, X gradients; sigma is incremented.
    !***********************************************************************

    implicit none
    type(TGFFData),intent(in) :: param
    integer               :: A,B,X,iTrB,iTrX,n,at(n)
    real(wp)                :: xyz(3,n)
    real(wp),intent(inout)  :: energy,gdr(3,3)
    real(wp)                :: q(n)
    type(TNeigh),intent(inout) :: neigh
    real(wp),intent(inout)        :: sigma(3,3)

    real(wp) ::  outl,dampl,damps,rdamp,damp
    real(wp) ::  ratio1,ratio2,ratio3
    real(wp) :: rab,rax,rbx,rab2,rax2,rbx2
    real(wp) :: drax(3),drbx(3),drab(3)
    real(wp) ::  dg(3),dga(3),dgb(3),dgx(3)
    real(wp) ::  gi,ga(3),gb(3),gx(3)
    real(wp) :: ex1_b,ex2_b,ex1_x,ex2_x,expo
    real(wp) ::  aterm,dterm
    real(wp) :: qb,qx
    real(wp) ::  cx,cb
    real(wp) ::  shortcut,const


    gdr = 0
    energy = 0

    cb = 1. ! param%xhbas(at(B)) !
    cx = param%xbaci(at(X))

    drax(1:3) = xyz(1:3,A)-(xyz(1:3,X)+neigh%transVec(1:3,iTrX))
    drbx(1:3) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-(xyz(1:3,X)+neigh%transVec(1:3,iTrX))
    drab(1:3) = xyz(1:3,A)-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))

    rab2 = sum(drab**2)
    rab = sqrt(rab2)

    rax2 = sum(drax**2)
    rax = sqrt(rax2)+1.d-12

    rbx2 = sum(drbx**2)
    rbx = sqrt(rbx2)+1.d-12

    !>-- out-of-line damping
    expo = param%xbacut*((rax+rbx)/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow !
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

    ratio1 = (rbx2/param%hblongcut_xb)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

    shortcut = param%xbscut*(param%rad(at(A))+param%rad(at(B)))
    ratio3 = (shortcut/rbx2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    rdamp = damp/rbx2/rbx ! **2

    ex1_x = exp(param%xbst*q(X))
    ex2_x = ex1_x+param%xbsf
    qx = ex1_x/ex2_x

    ex1_b = exp(-param%xbst*q(B))
    ex2_b = ex1_b+param%xbsf
    qb = ex1_b/ex2_b

    const = cb*qb*cx*qx

    !>-- r^3 only slightly better than r^4
    aterm = -rdamp*const
    dterm = -outl*const
    energy = -rdamp*outl*const

    !>-- damping part: rbx
    gi = rdamp*(-(2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3&
    &     /(1.d0+ratio3))-3.d0)/rbx2   ! 4,5,6 instead of 3 !
    gi = gi*dterm
    dg(1:3) = gi*drbx(1:3)
    gb(1:3) = dg(1:3)
    gx(1:3) = -dg(1:3)

    !>-- out-of-line term: rab, then rax and rbx
    gi = 2.d0*ratio2*expo*(rax+rbx)/(1.d0+ratio2)**2/(rax+rbx-rab)/rab2
    gi = gi*aterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = +dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

    gi = -2.d0*ratio2*expo/(1.d0+ratio2)**2/(rax+rbx-rab)/rax
    gi = gi*aterm
    dga(1:3) = gi*drax(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = -2.d0*ratio2*expo/(1.d0+ratio2)**2/(rax+rbx-rab)/rbx
    gi = gi*aterm
    dgb(1:3) = gi*drbx(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgx(1:3) = -dga(1:3)-dgb(1:3)
    gx(1:3) = gx(1:3)+dgx(1:3)

    gdr(1:3,1) = ga(1:3)
    gdr(1:3,2) = gb(1:3)
    gdr(1:3,3) = gx(1:3)
    sigma(:,1) = sigma(:,1)+ga(1)*xyz(:,A)
    sigma(:,2) = sigma(:,2)+ga(2)*xyz(:,A)
    sigma(:,3) = sigma(:,3)+ga(3)*xyz(:,A)
    sigma(:,1) = sigma(:,1)+gb(1)*(xyz(:,B)+neigh%transVec(:,iTrB))
    sigma(:,2) = sigma(:,2)+gb(2)*(xyz(:,B)+neigh%transVec(:,iTrB))
    sigma(:,3) = sigma(:,3)+gb(3)*(xyz(:,B)+neigh%transVec(:,iTrB))
    sigma(:,1) = sigma(:,1)+gx(1)*(xyz(:,X)+neigh%transVec(:,iTrX))
    sigma(:,2) = sigma(:,2)+gx(2)*(xyz(:,X)+neigh%transVec(:,iTrX))
    sigma(:,3) = sigma(:,3)+gx(3)*(xyz(:,X)+neigh%transVec(:,iTrX))

    return

  end subroutine rbxgfnff_eg

  subroutine egbend_nci_mul(j,i,k,vTrB,vTrC,c0,fc,n,at,xyz,e,g)
    !***********************************************************************
    !* Multiplicative bend factor e = 1-ea for the angle C(i)-B(j)...H(k) at
    !* vertex j, without distance damping; kijk is set so that the cosine
    !* form gives ea = fc at theta = 0.
    !* g(:,1:3) holds the gradients of j, i, k in that order.
    !***********************************************************************
    implicit none

    integer n,at(n)
    integer j,i,k ! B,C,H
    real(wp) vTrB(3),vTrC(3)
    real(wp) c0,fc
    real(wp) xyz(3,n),g(3,3),e

    real(wp) kijk,va(3),vb(3),vc(3),cosa
    real(wp) dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) rab2,vab(3),vcb(3),rp
    real(wp) rcb2
    real(wp) theta,deda(3),vp(3)

    kijk = fc/(cos(0.0d0)-cos(c0))**2
    va(1:3) = xyz(1:3,i)+vTrC
    vb(1:3) = xyz(1:3,j)+vTrB
    vc(1:3) = xyz(1:3,k)
    call vsub(va,vb,vab,3)
    call vsub(vc,vb,vcb,3)
    rab2 = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
    rcb2 = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
    call crprod(vcb,vab,vp)
    rp = vlen(vp)+1.d-14
    call impsc(vab,vcb,cosa)
    cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
    theta = dacos(cosa)

    if (pi-c0 .lt. 1.d-6) then     ! linear
      dt = theta-c0
      ea = kijk*dt**2
      deddt = 2.d0*kijk*dt
    else
      ea = kijk*(cosa-cos(c0))**2  ! not linear
      deddt = 2.0d0*kijk*sin(theta)*(cos(c0)-cosa)
    end if

    e = (1.0d0-ea)
    call crprod(vab,vp,deda)
    rmul1 = -deddt/(rab2*rp)
    deda = deda*rmul1
    call crprod(vcb,vp,dedc)
    rmul2 = deddt/(rcb2*rp)
    dedc = dedc*rmul2
    dedb = deda+dedc
    g(1:3,1) = dedb(1:3)
    g(1:3,2) = -deda(1:3)
    g(1:3,3) = -dedc(1:3)

  end subroutine egbend_nci_mul

  subroutine egtors_nci_mul(i,j,k,l,vTrR,vTrB,vTrC,rn,phi,phi0,tshift,n,at,xyz,e,g)
    !***********************************************************************
    !* Multiplicative torsion factor for i-j-k-l (called with R, B, C, H),
    !* without distance damping, which is inherent in the HB term.
    !* e runs from tshift to 1. g(:,1:4) holds the gradients of i, j, k, l.
    !***********************************************************************
    implicit none
    integer n,at(n)
    integer i,j,k,l
    integer rn
    real(wp) :: vTrR(3),vTrB(3),vTrC(3)
    real(wp) :: phi,phi0,tshift
    real(wp) :: xyz(3,n),g(3,4),e
    real(wp) :: fc
    real(wp) :: vab(3),vcb(3)
    real(wp) :: et,dij,c1
    real(wp) :: x1sin,x1cos,dphi1,vdc(3)
    real(wp) :: ddd(3),ddc(3),ddb(3),dda(3)
    real(wp) :: rij,rkl,rjk

    fc = (1.0d0-tshift)/2.0d0
    vab(1:3) = xyz(1:3,i)+vTrR-xyz(1:3,j)-vTrB
    vcb(1:3) = xyz(1:3,j)+vTrB-xyz(1:3,k)-vTrC
    vdc(1:3) = xyz(1:3,k)+vTrC-xyz(1:3,l)
    rij = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
    rjk = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
    rkl = vdc(1)*vdc(1)+vdc(2)*vdc(2)+vdc(3)*vdc(3)

    call dphidrPBC(1,n,xyz,i,j,k,l,vTrR,vTrB,vTrC,phi,dda,ddb,ddc,ddd)

    dphi1 = phi-phi0
    c1 = rn*dphi1+pi
    x1cos = cos(c1)
    x1sin = sin(c1)
    et = (1.+x1cos)*fc+tshift
    dij = -rn*x1sin*fc
    g(1:3,1) = dij*dda(1:3)
    g(1:3,2) = dij*ddb(1:3)
    g(1:3,3) = dij*ddc(1:3)
    g(1:3,4) = dij*ddd(1:3)
    e = et !*damp
  end subroutine egtors_nci_mul

end module gfnff_eg_hb
