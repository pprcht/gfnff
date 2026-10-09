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
!> Bonded energy terms of GFN-FF with gradient and strain: bonds (plain and HB
!> corrected), angles, torsions, the triple-bond torsion, the bonded ATM term.
module gfnff_eg_bonded

  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_data_types,only:TGFFData,TGFFTopology
  use gfnff_neighbor,only:TNeigh
  use gfnff_geometry,only:crprod,domegadrPBC,dphidr,impsc,omegaPBC,torsPBC,valijklff, &
    &                     vlen,vsub
  use gfnff_math_wrapper,only:gemv
  use gfnff_rab,only:gfnffdrab
  use gfnff_bond_potential,only:bond_potential,bond_potential_hb
  use gfnff_eg_terms,only:gfnffdampa,gfnffdampt
  implicit none
  private

  !> one driver per term; single-coordinate routines (egbond, ...) stay private
  public :: eg_bonds,eg_bonds_harmonic
  public :: eg_angles
  public :: eg_torsions,eg_storsions
  public :: eg_batm

  real(wp),private,parameter :: pi = 3.1415926535897932385_wp

contains  !> MODULE PROCEDURES START HERE

  subroutine eg_bonds(n,at,xyz,cn,dcn,dcndL,hb_cn,hb_dcn,dhbcndL, &
        & param,topo,neigh,version,ebond,g,sigma)
    !***********************************************************************
    !* Bond stretching energy over the whole bond list. Bonds that carry a
    !* hydrogen bond go to egbond_hb, all others to egbond; the reference
    !* distances rab0 and their CN derivatives are set up here for both.
    !* dcn/dcndL and hb_dcn/dhbcndL are the Cartesian and strain derivatives
    !* of cn and of hb_cn, the ERF coordination number of the HB correction.
    !* version selects the well shape. ebond, g and sigma are incremented.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: cn(n),dcn(3,n,n),dcndL(3,3,n)
    real(wp),intent(in) :: hb_cn(n),hb_dcn(3,n,n),dhbcndL(3,3,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    integer,intent(in) :: version
    real(wp),intent(inout) :: ebond
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)

    integer :: i,iat,jat,iTr
    real(wp) :: rab,rij,drijdcn(2),gf(3)
    real(wp),allocatable :: gfac(:,:),rab0(:),rabdcn(:,:),dEdcn(:)
    logical,allocatable :: considered_ABH(:,:,:)

    if (neigh%nbond .le. 0) return

    !>-- gfac is (3,nbond), not (3,n,nbond): drab/dr is rebuilt from dcn inside
    !>   egbond/egbond_hb rather than materialised for every bond
    allocate (gfac(3,neigh%nbond),rab0(neigh%nbond),rabdcn(2,neigh%nbond))
    rab0(:) = neigh%vbond(1,:)
    call gfnffdrab(n,at,cn,neigh%nbond,neigh%blist,rab0,gfac,rabdcn)
    allocate (dEdcn(n),source=0.0_wp)
    allocate (considered_ABH(topo%hb_mapNAB,topo%hb_mapNAB,topo%hb_mapNH),source=.false.)

    !$omp parallel do default(none) reduction(+:g, ebond, sigma, dEdcn) &
    !$omp shared(gfac, dcn, topo, neigh, param, rab0, rabdcn, xyz, at, hb_cn, hb_dcn, n, dhbcndL, considered_ABH, version) &
    !$omp private(i, iat, jat, rab, rij, gf, drijdcn, iTr)
    !> Atom indices come out of a list, so two iterations can hit the same
    !> entry of the gradient. Vectorising that scatter needs conflict
    !> detection; ifx omits it above -O1 and silently drops contributions.
    !> A comment to compilers that do not know the directive.
!DIR$ NOVECTOR
    do i = 1,neigh%nbond
      jat = neigh%blist(1,i)
      iat = neigh%blist(2,i)
      iTr = neigh%blist(3,i)
      if (iTr .gt. neigh%nTrans) cycle
      rab = NORM2(xyz(:,jat)-xyz(:,iat)+neigh%transVec(:,iTr))
      rij = rab0(i)
      gf = gfac(:,i)
      drijdcn = rabdcn(:,i)
      if (neigh%nr_hb(i) .ge. 1) then
        call egbond_hb(i,iat,jat,iTr,rab,rij,gf,dcn,drijdcn,hb_cn,hb_dcn,n,at,xyz,&
             &ebond,g,sigma,param,topo,neigh,version,dEdcn,dhbcndL,considered_ABH)
      else
        call egbond(i,iat,jat,iTr,rab,rij,gf,dcn,drijdcn,n,at,xyz,ebond,g,sigma,neigh,version,dEdcn)
      end if
    end do
    !$omp end parallel do

    !>-- strain part of the CN chain rule
    call gemv(dcndL,dEdcn,sigma,alpha=1.0_wp,beta=1.0_wp)

  end subroutine eg_bonds

  subroutine eg_angles(n,at,xyz,param,topo,neigh,eangl,g,sigma)
    !***********************************************************************
    !* Angle bending energy over the whole angle list topo%alist.
    !* eangl and g are incremented, sigma in the periodic case only.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(inout) :: eangl
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)

    integer :: m,i,j,k
    real(wp) :: etmp,g3tmp(3,3),ds(3,3)

    if (topo%nangl .le. 0) return

    !$omp parallel do default(none) reduction (+:eangl, g, sigma) &
    !$omp shared(n, at, xyz, topo, neigh, param) &
    !$omp private(m, j, i, k, etmp, g3tmp,ds)
    !>-- gradient scatter must not be vectorised, see eg_bonds
!DIR$ NOVECTOR
    do m = 1,topo%nangl
      i = topo%alist(1,m)
      j = topo%alist(2,m)
      k = topo%alist(3,m)
      call egbend(m,j,i,k,n,at,xyz,etmp,g3tmp,ds,param,topo,neigh)
      g(1:3,i) = g(1:3,i)+g3tmp(1:3,1) ! alist has swapped i and j
      g(1:3,j) = g(1:3,j)+g3tmp(1:3,2) ! compared to orig gfnff
      g(1:3,k) = g(1:3,k)+g3tmp(1:3,3) ! therefore swapped i and j here too
      if (neigh%nTrans .ne. 1) sigma = sigma+ds
      eangl = eangl+etmp
    end do
    !$omp end parallel do

  end subroutine eg_angles

  subroutine eg_torsions(n,at,xyz,param,topo,neigh,etors,g,sigma)
    !***********************************************************************
    !* Torsion energy over topo%tlist; entries whose translation indices exceed
    !* neigh%nTrans are skipped. etors and g are incremented, sigma if periodic.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(inout) :: etors
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)

    integer :: m,i,j,k,l,iTrl,iTrj,iTrk
    real(wp) :: etmp,g4tmp(3,4),ds(3,3)

    if (topo%ntors .le. 0) return

    !$omp parallel do default(none) reduction(+:etors, g, sigma) &
    !$omp shared(param, topo, neigh, n, at, xyz) &
    !$omp private(m, i, j, k, l,iTrl,iTrj,iTrk, etmp, g4tmp,ds)
    !>-- gradient scatter must not be vectorised, see eg_bonds
!DIR$ NOVECTOR
    do m = 1,topo%ntors
      i = topo%tlist(1,m)  ! is actually l  ! for out-of-plane it is correct
      j = topo%tlist(2,m)  ! is actually i  ! for out-of-plane it is correct
      k = topo%tlist(3,m)  ! is actually j  ! for out-of-plane it is correct
      l = topo%tlist(4,m)  ! is actually k  ! for out-of-plane it is correct
      iTrl = topo%tlist(6,m)
      iTrj = topo%tlist(7,m)
      iTrk = topo%tlist(8,m)
      if (iTrj .gt. neigh%nTrans.or.iTrk .gt. neigh%nTrans.or.iTrl .gt. neigh%nTrans) cycle
      call egtors(m,i,j,k,l,iTrl,iTrj,iTrk,n,at,xyz,etmp,g4tmp,ds,param,topo,neigh)
      g(1:3,i) = g(1:3,i)+g4tmp(1:3,1)
      g(1:3,j) = g(1:3,j)+g4tmp(1:3,2)
      g(1:3,k) = g(1:3,k)+g4tmp(1:3,3)
      g(1:3,l) = g(1:3,l)+g4tmp(1:3,4)
      if (neigh%nTrans .ne. 1) sigma = sigma+ds
      etors = etors+etmp
    end do
    !$omp end parallel do

  end subroutine eg_torsions

  subroutine eg_storsions(n,xyz,topo,etors,g)
    !***********************************************************************
    !* Special torsion potential for rotation around triple bonded carbon,
    !* over topo%sTorsl. The regular torsion energy etors and g are incremented.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFTopology),intent(in) :: topo
    real(wp),intent(inout) :: etors
    real(wp),intent(inout) :: g(3,n)

    integer :: i,m
    real(wp) :: etmp
    real(wp),allocatable :: g5tmp(:,:)

    if (.not.allocated(topo%sTorsl)) return
    m = size(topo%sTorsl(1,:))
    if (m .eq. 0) return

    allocate (g5tmp(3,n))
    do i = 1,m
      call sTors_eg(i,n,xyz,topo,etmp,g5tmp)
      etors = etors+etmp
      g = g+g5tmp
    end do

  end subroutine eg_storsions

  subroutine eg_batm(n,at,xyz,param,topo,neigh,ebatm,g,sigma)
    !***********************************************************************
    !* Bonded three-body (ATM) term over the triples in topo%b3list.
    !* ebatm and g are incremented, sigma in the periodic case only.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(inout) :: ebatm
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)

    integer :: i,j,k,l,iTrk,iTrl
    real(wp) :: etmp,g3tmp(3,3),ds(3,3)

    if (topo%nbatm .le. 0) return

    !$omp parallel do default(none) reduction(+:ebatm, g, sigma) &
    !$omp shared(n, at, xyz, topo, neigh, param) &
    !$omp private(i, j, k, l, iTrk, iTrl, etmp, g3tmp, ds)
    !>-- gradient scatter must not be vectorised, see eg_bonds
!DIR$ NOVECTOR
    do i = 1,topo%nbatm
      j = topo%b3list(1,i)
      k = topo%b3list(2,i)
      l = topo%b3list(3,i)
      iTrk = topo%b3list(4,i)
      iTrl = topo%b3list(5,i)
      call batmgfnff_eg(n,j,k,l,iTrk,iTrl,at,xyz,topo%qa,etmp,g3tmp,ds,param,neigh)
      g(1:3,j) = g(1:3,j)+g3tmp(1:3,1)
      g(1:3,k) = g(1:3,k)+g3tmp(1:3,2)
      g(1:3,l) = g(1:3,l)+g3tmp(1:3,3)
      if (neigh%nTrans .ne. 1) sigma = sigma+ds
      ebatm = ebatm+etmp
    end do
    !$omp end parallel do

  end subroutine eg_batm

  subroutine eg_bonds_harmonic(n,at,xyz,param,neigh,ebond,g)
    !***********************************************************************
    !* Extremely crude harmonic bond potential used only by the harmonic2020
    !* version for 2D-to-3D structure conversion: equilibrium distance from
    !* covalent radii, force constant fixed at 0.1. The bond list is
    !* neigh%blist; topo%blist is never filled and would silently reduce
    !* harmonic2020 to its repulsion term. ebond is overwritten, g incremented.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(out) :: ebond
    real(wp),intent(inout) :: g(3,n)

    integer :: i,iat,jat,iTr
    real(wp) :: rab,r2,r3(3),rn,dum

    ebond = 0
    !$omp parallel do default(none) reduction(+:ebond, g) &
    !$omp shared(neigh, param, xyz, at) private(i, iat, jat, iTr, rab, r2, r3, rn, dum)
    !>-- gradient scatter must not be vectorised, see eg_bonds
!DIR$ NOVECTOR
    do i = 1,neigh%nbond
      iat = neigh%blist(2,i)
      jat = neigh%blist(1,i)
      iTr = neigh%blist(3,i)
      if (iTr .gt. neigh%nTrans) cycle
      r3 = xyz(:,iat)-xyz(:,jat)-neigh%transVec(:,iTr)
      rab = sqrt(sum(r3*r3))
      rn = 0.7*(param%rcov(at(iat))+param%rcov(at(jat)))
      r2 = rn-rab
      ebond = ebond+0.1d0*r2**2  ! fixfc = 0.1
      dum = 0.1d0*2.0d0*r2/rab
      g(:,jat) = g(:,jat)+dum*r3
      g(:,iat) = g(:,iat)-dum*r3
    end do
    !$omp end parallel do

  end subroutine eg_bonds_harmonic

  subroutine egbond(i,iat,jat,iTr,rab,rij,gf,dcn,drijdcn,n,at,xyz,e,g,sigma,neigh,version,dEdcn)
    !***********************************************************************
    !* Energy and gradient of a single bond. The well itself lives in
    !* bond_potential; the Cartesian projection, the strain contribution and
    !* the CN chain rule (the three-body part) are independent of its shape.
    !***********************************************************************
    implicit none
    type(TNeigh),intent(in) :: neigh
    integer,intent(in)   :: i
    integer,intent(in)   :: n
    integer,intent(in)   :: iat
    integer,intent(in)   :: jat
    integer,intent(in)   :: at(n)
    integer,intent(in)   :: iTr
    real(wp),intent(in)    :: rab
    real(wp),intent(in)    :: rij
    real(wp),intent(in)    :: gf(3)      !> (scaleF*ff, cnfak(ati), cnfak(atj))
    real(wp),intent(in)    :: dcn(3,n,n)
    real(wp),intent(in)    :: drijdcn(2)
    real(wp),intent(in)    :: xyz(3,n)
    integer,intent(in)   :: version
    real(wp),intent(inout) :: dEdcn(n)
    real(wp),intent(inout) :: e
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)
    integer k
    real(wp) dr,dum,vrab(3)
    real(wp) dx,dy,dz,dg(3)
    real(wp) yy,ed
    real(wp) t4,t5,t6,t8

    t8 = neigh%vbond(2,i)
    dr = rab-rij
    call bond_potential(version,t8,neigh%vbond(3,i),dr,dum,ed=ed)
    e = e+dum
    !>-- yy is dE/drij, the reference length being what the CN chain rule sees
    yy = -ed
    dx = -xyz(1,jat)+xyz(1,iat)-neigh%transVec(1,iTr)
    dy = -xyz(2,jat)+xyz(2,iat)-neigh%transVec(2,iTr)
    dz = -xyz(3,jat)+xyz(3,iat)-neigh%transVec(3,iTr)
    vrab(1) = dx
    vrab(2) = dy
    vrab(3) = dz
    t4 = -yy*dx/rab
    t5 = -yy*dy/rab
    t6 = -yy*dz/rab
    dg(1) = t4
    dg(2) = t5
    dg(3) = t6
    g(1,iat) = g(1,iat)+t4 ! to avoid if in loop below
    g(2,iat) = g(2,iat)+t5
    g(3,iat) = g(3,iat)+t6
    dEdcn(iat) = dEdcn(iat)+yy*drijdcn(1)
    if (neigh%nTrans .ne. 1) then   !> stress only wanted for PBC
      sigma(:,1) = sigma(:,1)+dg(1)*vrab
      sigma(:,2) = sigma(:,2)+dg(2)*vrab
      sigma(:,3) = sigma(:,3)+dg(3)*vrab
    end if

    t4 = yy*(dx/rab)
    t5 = yy*(dy/rab)
    t6 = yy*(dz/rab)
    dg(1) = t4
    dg(2) = t5
    dg(3) = t6
    g(1,jat) = g(1,jat)+t4 ! to avoid if in loop below
    g(2,jat) = g(2,jat)+t5
    g(3,jat) = g(3,jat)+t6
    !>-- rij depends on cn; dEdcn feeds the strain chain rule in eg_bonds
    dEdcn(jat) = dEdcn(jat)+yy*drijdcn(2)

    if (neigh%nTrans .eq. 1) then   !> molecular: no stress tensor wanted
      do k = 1,n !3B gradient
        g(:,k) = g(:,k)+(gf(1)*(gf(2)*dcn(:,k,iat)+gf(3)*dcn(:,k,jat)))*yy
      end do
    else
      do k = 1,n !3B gradient
        dg = (gf(1)*(gf(2)*dcn(:,k,iat)+gf(3)*dcn(:,k,jat)))*yy
        g(:,k) = g(:,k)+dg
        sigma(:,1) = sigma(:,1)+dg(1)*vrab
        sigma(:,2) = sigma(:,2)+dg(2)*vrab
        sigma(:,3) = sigma(:,3)+dg(3)*vrab
      end do
    end if

  end subroutine egbond

  subroutine egbond_hb(i,iat,jat,iTr,rab,rij,gf,dcn,drijdcn,hb_cn,hb_dcn,n,at,xyz,e,&
                  &g,sigma,param,topo,neigh,version,dEdcn,dhbcndL,considered_ABH)
    !***********************************************************************
    !* As egbond, for a bond that takes part in a hydrogen bridge. The
    !* steepness is also scaled by the hydrogen-bond coordination number of
    !* the bridging hydrogen, which adds a second chain rule.
    !***********************************************************************
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    integer,intent(in)   :: i
    integer,intent(in)   :: n
    integer,intent(in)   :: iat
    integer,intent(in)   :: jat
    integer,intent(in)   :: iTr ! transVec index for jat
    integer,intent(in)   :: at(n)
    real(wp),intent(in)    :: rab
    real(wp),intent(in)    :: rij
    real(wp),intent(in)    :: gf(3)      !> (scaleF*ff, cnfak(ati), cnfak(atj))
    real(wp),intent(in)    :: dcn(3,n,n)
    real(wp),intent(in)    :: drijdcn(2)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(in)    :: hb_cn(n)
    real(wp),intent(in)    :: hb_dcn(3,n,n)
    real(wp),intent(in) :: dhbcndL(3,3,n)
    integer,intent(in)   :: version
    logical,intent(inout)  :: considered_ABH(topo%hb_mapNAB,topo%hb_mapNAB,topo%hb_mapNH)! only consider ABH triplets once; indep of iTr
    real(wp),intent(inout) :: e
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)
    real(wp),intent(inout) :: dEdcn(n)
    integer j,k
    integer jA,jH,iTrA,iTrH,iTrB
    integer hbH,hbB,hbA
    integer mapA,mapB,mapH
    real(wp) dr,dum
    real(wp) dx,dy,dz,vrab(3),dg(3)
    real(wp) yy,zz,ed
    real(wp) t1,t4,t5,t6

    if (at(iat) .eq. 1) then
      hbH = iat
      hbA = jat
    else if (at(jat) .eq. 1) then
      hbH = jat
      hbA = iat
    else
      write (stdout,'(10x,"No H-atom found in this bond ",i0,1x,i0)') iat,jat
      return
    end if

    t1 = 1.0-param%vbond_scale
    dr = rab-rij
    call bond_potential_hb(version,neigh%vbond(2,i),t1,hb_cn(hbH), &
       & neigh%vbond(3,i),dr,dum,ed,zz)
    e = e+dum
    !>-- yy is dE/drij as in egbond, zz is dE/dhbcn for the second chain rule
    yy = -ed
    dx = -xyz(1,jat)+xyz(1,iat)-neigh%transVec(1,iTr)
    dy = -xyz(2,jat)+xyz(2,iat)-neigh%transVec(2,iTr)
    dz = -xyz(3,jat)+xyz(3,iat)-neigh%transVec(3,iTr)
    vrab(1) = dx
    vrab(2) = dy
    vrab(3) = dz
    t4 = -yy*dx/rab
    t5 = -yy*dy/rab
    t6 = -yy*dz/rab
    dg(1) = t4
    dg(2) = t5
    dg(3) = t6
    g(1,iat) = g(1,iat)+t4 ! to avoid if in loop below
    g(2,iat) = g(2,iat)+t5
    g(3,iat) = g(3,iat)+t6
    dEdcn(iat) = dEdcn(iat)+yy*drijdcn(1)
    if (neigh%nTrans .ne. 1) then   !> stress only wanted for PBC
      sigma(:,1) = sigma(:,1)+dg(1)*vrab
      sigma(:,2) = sigma(:,2)+dg(2)*vrab
      sigma(:,3) = sigma(:,3)+dg(3)*vrab
    end if

    t4 = yy*(dx/rab)
    t5 = yy*(dy/rab)
    t6 = yy*(dz/rab)
    dg(1) = t4
    dg(2) = t5
    dg(3) = t6
    g(1,jat) = g(1,jat)+t4 ! to avoid if in loop below
    g(2,jat) = g(2,jat)+t5
    g(3,jat) = g(3,jat)+t6
    dEdcn(jat) = dEdcn(jat)+yy*drijdcn(2)

    if (neigh%nTrans .eq. 1) then   !> molecular: no stress tensor wanted
      do k = 1,n !3B gradient
        g(:,k) = g(:,k)+(gf(1)*(gf(2)*dcn(:,k,iat)+gf(3)*dcn(:,k,jat)))*yy
      end do
    else
      do k = 1,n !3B gradient
        dg = (gf(1)*(gf(2)*dcn(:,k,iat)+gf(3)*dcn(:,k,jat)))*yy
        g(:,k) = g(:,k)+dg
        sigma(:,1) = sigma(:,1)+dg(1)*vrab
        sigma(:,2) = sigma(:,2)+dg(2)*vrab
        sigma(:,3) = sigma(:,3)+dg(3)*vrab
      end do
    end if

    do j = 1,topo%bond_hb_nr !CN gradient
      jA = topo%bond_hb_AH(1,j)
      jH = topo%bond_hb_AH(2,j)
      iTrA = topo%bond_hb_AH(3,j)
      iTrH = topo%bond_hb_AH(4,j)
      if ((jH .eq. hbH.and.jA .eq. hbA.and.iTrA .eq. 1.and.iTrH .eq. iTr).or.&
         &(jH .eq. hbH.and.jA .eq. hbA.and.iTrH .eq. 1.and.iTrA .eq. iTr)) then
        dg = hb_dcn(:,hbH,hbH)*zz
        g(:,hbH) = g(:,hbH)+dg
        if (hbH .eq. iat) then  ! only hbH cf. energy -> dum above
          sigma = sigma+zz*dhbcndL(:,:,iat)
        end if
        if (hbH .eq. jat) then ! only hbH cf. energy -> dum above
          sigma = sigma+zz*dhbcndL(:,:,jat)
        end if
        do k = 1,topo%bond_hb_Bn(j)
          hbB = topo%bond_hb_B(1,k,j)
          iTrB = topo%bond_hb_B(2,k,j)
          mapA = topo%hb_mapABH(hbA)
          mapB = topo%hb_mapABH(hbB)
          mapH = topo%hb_mapABH(hbH)
          if (.not.considered_ABH(mapA,mapB,mapH)) then
            considered_ABH(mapA,mapB,mapH) = .true.
            dg = hb_dcn(:,hbB,hbH)*zz
            g(:,hbB) = g(:,hbB)-dg
          end if
        end do
      end if
    end do
  end subroutine egbond_hb

  subroutine egbend(m,j,i,k,n,at,xyz,e,g,ds,param,topo,neigh)
    !***********************************************************************
    !* Energy and gradient of bending angle m (j-i-k, i central), damped along
    !* both bonds; harmonic in theta for a linear reference, else in cos(theta).
    !* g holds the columns for i, j, k; ds stays zero in the molecular case.
    !***********************************************************************
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    integer m,n,at(n)
    integer i,j,k,dim1,dim2
    real(wp) xyz(3,n),g(3,3),e,ds(3,3)

    real(wp) c0,kijk,va(3),vb(3),vc(3),cosa
    real(wp) dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) theta,deda(3),vp(3)
    integer :: iTrj,iTrk
    real(wp)  :: vTrj(3),vTrk(3)

    ds = 0.0_wp
    c0 = topo%vangl(1,m)
    kijk = topo%vangl(2,m)
    iTrj = topo%alist(4,m)
    vTrj = neigh%transVec(:,iTrj)
    iTrk = topo%alist(5,m)
    vTrk = neigh%transVec(:,iTrk)
    va(1:3) = xyz(1:3,j)+vTrj
    vb(1:3) = xyz(1:3,i)
    vc(1:3) = xyz(1:3,k)+vTrk
    call vsub(va,vb,vab,3)
    call vsub(vc,vb,vcb,3)
    rab2 = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
    rcb2 = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
    call crprod(vcb,vab,vp)
    rp = vlen(vp)+1.d-14
    call impsc(vab,vcb,cosa)
    cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
    theta = dacos(cosa)  ! angle for bond j-i-k  => va-vb-vc  (vb is in the middle)

    call gfnffdampa(at(j),at(i),rab2,dampij,damp2ij,param)
    call gfnffdampa(at(k),at(i),rcb2,dampjk,damp2jk,param)
    damp = dampij*dampjk

    if (pi-c0 .lt. 1.d-6) then ! linear
      dt = theta-c0
      ea = kijk*dt**2
      deddt = 2.d0*kijk*dt
    else
      ea = kijk*(cosa-cos(c0))**2
      deddt = 2.0d0*kijk*sin(theta)*(cos(c0)-cosa)
    end if

    e = ea*damp
    call crprod(vab,vp,deda)
    rmul1 = -deddt/(rab2*rp)
    deda = deda*rmul1
    call crprod(vcb,vp,dedc)
    rmul2 = deddt/(rcb2*rp)
    dedc = dedc*rmul2
    dedb = deda+dedc
    term1(1:3) = ea*damp2ij*dampjk*vab(1:3)
    term2(1:3) = ea*damp2jk*dampij*vcb(1:3)
    g(1:3,1) = -dedb(1:3)*damp-term1(1:3)-term2(1:3)
    g(1:3,2) = deda(1:3)*damp+term1(1:3)
    g(1:3,3) = dedc(1:3)*damp+term2(1:3)
    if (neigh%nTrans .ne. 1) then
      do dim1 = 1,3
        do dim2 = dim1,3
          ds(dim1,dim2) = g(dim2,1)*vb(dim1)+g(dim2,2)*va(dim1)+g(dim2,3)*vc(dim1)
        end do
      end do
      do dim1 = 1,3
        do dim2 = 1,dim1-1
          ds(dim1,dim2) = ds(dim2,dim1)
        end do
      end do
    end if

  end subroutine egbend

  subroutine egtors(m,i,j,k,l,iTrl,iTrj,iTrk,n,at,xyz,e,g,ds,param,topo,neigh)
    !***********************************************************************
    !* Energy and gradient of torsion m. tlist(5,m) > 0 is the periodicity
    !* of a proper torsion; otherwise the entry is an out-of-plane (improper)
    !* term around atom j. Both are damped along their three bonds. g holds
    !* the columns for i, j, k, l; ds is set in the periodic case only.
    !***********************************************************************
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    integer :: m,n,at(n)
    integer :: i,j,k,l,iTrl,iTrj,iTrk,dim1,dim2
    real(wp) :: xyz(3,n),g(3,4),e,ds(3,3)
    real(wp) :: vTrl(3),vTrj(3),vTrk(3)

    real(wp) :: va(3),vb(3),vc(3),vd(3)
    real(wp) :: term1(3),term2(3),vab(3),vcb(3)
    real(wp) :: damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) :: et,dij,c1
    real(wp) :: term3(3),x1sin,x1cos,dphi1,vdc(3)
    real(wp) :: ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) :: rij,phi0,rkl,rjk,dampkl,damp2kl
    real(wp) :: dampjl,damp2jl,rn

    rn = dble(topo%tlist(5,m))
    phi0 = topo%vtors(1,m)

    if (topo%tlist(5,m) .gt. 0) then
      vTrl = neigh%transVec(:,iTrl)
      vTrj = neigh%transVec(:,iTrj)
      vTrk = neigh%transVec(:,iTrk)
      va = xyz(1:3,i)+vTrl
      vb = xyz(1:3,j)
      vc = xyz(1:3,k)+vTrj
      vd = xyz(1:3,l)+vTrk
      vab(1:3) = va-vb
      vcb(1:3) = vb-vc
      vdc(1:3) = vc-vd
      rij = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
      rjk = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
      rkl = vdc(1)*vdc(1)+vdc(2)*vdc(2)+vdc(3)*vdc(3)

      call gfnffdampt(at(i),at(j),rij,dampij,damp2ij,param)
      call gfnffdampt(at(k),at(j),rjk,dampjk,damp2jk,param)
      call gfnffdampt(at(k),at(l),rkl,dampkl,damp2kl,param)
      damp = dampjk*dampij*dampkl

      !>-- angle and derivatives in one pass, sharing bond vectors and normals
      call torsPBC(1,n,xyz,i,j,k,l,vTrj,vTrk,vTrl,phi,dda,ddb,ddc,ddd)
      dphi1 = phi-phi0
      c1 = rn*dphi1+pi
      x1cos = cos(c1)
      x1sin = sin(c1)
      et = (1.+x1cos)*topo%vtors(2,m)
      dij = -rn*x1sin*topo%vtors(2,m)*damp
      term1(1:3) = et*damp2ij*dampjk*dampkl*vab(1:3)
      term2(1:3) = et*damp2jk*dampij*dampkl*vcb(1:3)
      term3(1:3) = et*damp2kl*dampij*dampjk*vdc(1:3)
      g(1:3,1) = dij*dda(1:3)+term1
      g(1:3,2) = dij*ddb(1:3)-term1+term2
      g(1:3,3) = dij*ddc(1:3)+term3-term2
      g(1:3,4) = dij*ddd(1:3)-term3
      if (neigh%nTrans .ne. 1) then
        do dim1 = 1,3
          do dim2 = dim1,3
            ds(dim1,dim2) = g(dim2,1)*va(dim1) &
                    &      +g(dim2,2)*vb(dim1) &
                    &      +g(dim2,3)*vc(dim1) &
                    &      +g(dim2,4)*vd(dim1)
            ds(dim2,dim1) = ds(dim1,dim2)
          end do
        end do
        do dim1 = 1,3
          do dim2 = 1,dim1-1
            ds(dim1,dim2) = ds(dim2,dim1)
          end do
        end do
      end if
      e = et*damp
    else  ! out-of-plane, improper
      vTrl = neigh%transVec(:,iTrl)
      vTrj = neigh%transVec(:,iTrj)
      vTrk = neigh%transVec(:,iTrk)
      va = xyz(1:3,i)
      vb = xyz(1:3,j)+vTrj
      vc = xyz(1:3,k)+vTrk
      vd = xyz(1:3,l)+vTrl
      vab(1:3) = vb-va
      vcb(1:3) = vb-vc
      vdc(1:3) = vb-vd
      rij = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
      rjk = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
      rjl = vdc(1)*vdc(1)+vdc(2)*vdc(2)+vdc(3)*vdc(3)

      call gfnffdampt(at(i),at(j),rij,dampij,damp2ij,param)
      call gfnffdampt(at(k),at(j),rjk,dampjk,damp2jk,param)
      call gfnffdampt(at(j),at(l),rjl,dampjl,damp2jl,param)
      damp = dampjk*dampij*dampjl

      phi = omegaPBC(n,xyz,i,j,k,l,vTrl,vTrj,vTrk)
      call domegadrPBC(n,xyz,i,j,k,l,vTrl,vTrj,vTrk,&
                     & phi,dda,ddb,ddc,ddd)

      if (topo%tlist(5,m) .eq. 0) then  ! phi0=0 case
        dphi1 = phi-phi0
        c1 = dphi1+pi
        x1cos = cos(c1)
        x1sin = sin(c1)
        et = (1.+x1cos)*topo%vtors(2,m)
        dij = -x1sin*topo%vtors(2,m)*damp
      else                     ! double min at phi0,-phi0
        et = topo%vtors(2,m)*(cos(phi)-cos(phi0))**2
        dij = 2.*topo%vtors(2,m)*sin(phi)*(cos(phi0)-cos(phi))*damp
      end if
      term1(1:3) = et*damp2ij*dampjk*dampjl*vab(1:3)
      term2(1:3) = et*damp2jk*dampij*dampjl*vcb(1:3)
      term3(1:3) = et*damp2jl*dampij*dampjk*vdc(1:3)
      g(1:3,1) = dij*dda(1:3)-term1
      g(1:3,2) = dij*ddb(1:3)+term1+term2+term3
      g(1:3,3) = dij*ddc(1:3)-term2
      g(1:3,4) = dij*ddd(1:3)-term3
      if (neigh%nTrans .ne. 1) then
        do dim1 = 1,3
          do dim2 = dim1,3
            ds(dim1,dim2) = g(dim2,1)*va(dim1) &
                    &      +g(dim2,2)*vb(dim1) &
                    &      +g(dim2,3)*vc(dim1) &
                    &      +g(dim2,4)*vd(dim1)
          end do
        end do
        do dim1 = 1,3
          do dim2 = 1,dim1-1
            ds(dim1,dim2) = ds(dim2,dim1)
          end do
        end do
      end if
      e = et*damp
    end if

  end subroutine egtors

  subroutine batmgfnff_eg(n,iat,jat,kat,iTrj,iTrk,at,xyz,q,energy,g,ds,param,neigh)
    !***********************************************************************
    !* Bonded three-body ATM term of the triple iat, jat, kat, scaled by a
    !* charge-dependent factor. Taken from the D3 ATM code. g holds the
    !* columns for iat, jat, kat; ds stays zero in the molecular case.
    !***********************************************************************
    implicit none
    type(TGFFData),intent(in) :: param
    type(TNeigh),intent(in) :: neigh
    integer,intent(in) :: iat,jat,kat,n,at(n),iTrj,iTrk
    real(wp),intent(in) :: xyz(3,n),q(n)
    real(wp),intent(out) :: energy,g(3,3),ds(3,3)

    real(wp) :: r2ij,r2jk,r2ik,sr2ij,sr2jk,sr2ik,invsr2ij,invsr2jk,invsr2ik
    real(wp) :: c9,mijk,imjk,ijmk,rijk3,ang,angr9,rav3
    real(wp) :: rij(3),rik(3),rjk(3),ri(3),rj(3),rk(3),drij,drik,drjk,dang,ff,fi,fj,fk
    real(wp),parameter :: fqq = 3.0_wp
    integer :: iTrDum,dm1,dm2

    energy = 0.0_wp
    g = 0.0_wp
    ds = 0.0_wp

    fi = (1.0_wp-fqq*q(iat))
    fi = min(max(fi,-4.0_wp),4.0_wp)
    fj = (1.0_wp-fqq*q(jat))
    fj = min(max(fj,-4.0_wp),4.0_wp)
    fk = (1.0_wp-fqq*q(kat))
    fk = min(max(fk,-4.0_wp),4.0_wp)
    ff = fi*fj*fk ! charge term
    c9 = ff*param%zb3atm(at(iat))*param%zb3atm(at(jat))*param%zb3atm(at(kat)) ! strength of interaction
    r2ij = sum((xyz(1:3,iat)-(xyz(1:3,jat)+neigh%transVec(1:3,iTrj)))**2)
    r2ik = sum((xyz(1:3,iat)-(xyz(1:3,kat)+neigh%transVec(1:3,iTrk)))**2)
    !>-- iTrDum indexes transVec(iTrk)-transVec(iTrj) if that is tabulated
    iTrDum = neigh%fTrSum(neigh%iTrNeg(iTrj),iTrk)
    if (iTrDum <= 0.or.iTrDum > neigh%numctr) then
      r2jk = sum(((xyz(:,kat)+neigh%transVec(:,iTrk)) &
                  -(xyz(:,jat)+neigh%transVec(:,iTrj)))**2)
    else
      r2jk = sum((xyz(:,jat)-(xyz(:,kat)+neigh%transVec(:,iTrDum)))**2)
    end if
    sr2ij = sqrt(r2ij)
    sr2ik = sqrt(r2ik)
    sr2jk = sqrt(r2jk)
    invsr2ij = 1._wp/sr2ij
    invsr2ik = 1._wp/sr2ik
    invsr2jk = 1._wp/sr2jk
    mijk = -r2ij+r2jk+r2ik
    imjk = r2ij-r2jk+r2ik
    ijmk = r2ij+r2jk-r2ik
    rijk3 = r2ij*r2jk*r2ik
    rav3 = rijk3*sr2ij*sr2jk*sr2ik ! R^9
    ang = 0.375_wp*ijmk*imjk*mijk/rijk3
    angr9 = (ang+1.0_wp)/rav3
    energy = c9*angr9

    dang = -0.375_wp*(r2ij**3+r2ij**2*(r2jk+r2ik) &
        & +r2ij*(3.0_wp*r2jk**2+2.0_wp*r2jk*r2ik+3.0_wp*r2ik**2) &
        & -5.0_wp*(r2jk-r2ik)**2*(r2jk+r2ik)) &
        & /(sr2ij*rijk3*rav3)
    drij = -dang*c9
    dang = -0.375_wp*(r2jk**3+r2jk**2*(r2ik+r2ij) &
        & +r2jk*(3.0_wp*r2ik**2+2.0_wp*r2ik*r2ij+3.0_wp*r2ij**2) &
        & -5.0_wp*(r2ik-r2ij)**2*(r2ik+r2ij)) &
        & /(sr2jk*rijk3*rav3)
    drjk = -dang*c9
    dang = -0.375_wp*(r2ik**3+r2ik**2*(r2jk+r2ij) &
        & +r2ik*(3.0_wp*r2jk**2+2.0_wp*r2jk*r2ij+3.0_wp*r2ij**2) &
        & -5.0_wp*(r2jk-r2ij)**2*(r2jk+r2ij)) &
        & /(sr2ik*rijk3*rav3)
    drik = -dang*c9

    rij = xyz(:,jat)-xyz(:,iat)+neigh%transVec(:,iTrj)
    rik = xyz(:,kat)-xyz(:,iat)+neigh%transVec(:,iTrk)
    if (iTrDum <= 0.or.iTrDum > neigh%numctr) then
      rjk = (xyz(:,kat)+neigh%transVec(:,iTrk))-(xyz(:,jat)+neigh%transVec(:,iTrj))
    else
      rjk = xyz(:,kat)-xyz(:,jat)+neigh%transVec(:,iTrDum)
    end if
    g(:,1) = drij*rij*invsr2ij
    g(:,1) = g(:,1)+drik*rik*invsr2ik
    g(:,2) = drjk*rjk*invsr2jk
    g(:,2) = g(:,2)-drij*rij*invsr2ij
    g(:,3) = -drik*rik*invsr2ik
    g(:,3) = g(:,3)-drjk*rjk*invsr2jk

    if (neigh%nTrans /= 1) then
      ri = xyz(:,iat)
      rj = xyz(:,jat)+neigh%transVec(:,iTrj)
      rk = xyz(:,kat)+neigh%transVec(:,iTrk)
      do dm1 = 1,3
        do dm2 = dm1,3
          ds(dm1,dm2) = (drij*rij(dm2)*invsr2ij)*ri(dm1) & ! i derivatives
              & +(drik*rik(dm2)*invsr2ik)*ri(dm1) &
              & +(drjk*rjk(dm2)*invsr2jk)*rj(dm1) & ! j derivatives
              & -(drij*rij(dm2)*invsr2ij)*rj(dm1) &
              & -(drik*rik(dm2)*invsr2ik)*rk(dm1) & ! k derivatives
              & -(drjk*rjk(dm2)*invsr2jk)*rk(dm1)
          ds(dm2,dm1) = ds(dm1,dm2)
        end do
      end do
    end if

  end subroutine batmgfnff_eg

  subroutine sTors_eg(m,n,xyz,topo,energy,dg)
    !***********************************************************************
    !* Torsion term m around triple bonded carbon on atoms sTorsl((/1,2,5,6/),m).
    !* energy and dg are overwritten, zero if any entry of sTorsl(:,m) is zero.
    !***********************************************************************
    integer,intent(in) :: m
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFTopology),intent(in) :: topo
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: dg(3,n)
    integer :: c1,c2,c3,c4
    integer :: i

    real(wp) :: phi
    real(wp) :: erefhalf
    real(wp) :: dp1(3),dp2(3),dp3(3),dp4(3)

    energy = 0.0_wp
    dg(:,:) = 0.0_wp

    if (.not.any(topo%sTorsl(:,m) .eq. 0)) then

      c1 = topo%sTorsl(1,m)
      c2 = topo%sTorsl(2,m)
      c3 = topo%sTorsl(5,m)
      c4 = topo%sTorsl(6,m)

      phi = valijklff(n,xyz,c1,c2,c3,c4)
      call dphidr(n,xyz,c1,c2,c3,c4,phi,dp1,dp2,dp3,dp4)

      !>-- reference energy for torsion of 90 deg,
      !>   calculated with DLPNO-CCSD(T) CBS on diphenylacetylene
      erefhalf = 3.75_wp*1.0e-4_wp  ! approx 1.97 kJ/mol !
      energy = -erefhalf*cos(2.0_wp*phi)+erefhalf
      do i = 1,3
        dg(i,c1) = dg(i,c1)+erefhalf*2.0_wp*sin(2.0_wp*phi)*dp1(i)
        dg(i,c2) = dg(i,c2)+erefhalf*2.0_wp*sin(2.0_wp*phi)*dp2(i)
        dg(i,c3) = dg(i,c3)+erefhalf*2.0_wp*sin(2.0_wp*phi)*dp3(i)
        dg(i,c4) = dg(i,c4)+erefhalf*2.0_wp*sin(2.0_wp*phi)*dp4(i)
      end do
    end if
  end subroutine sTors_eg

end module gfnff_eg_bonded
