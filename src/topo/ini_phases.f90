! ──────────────────────────────────────────────────────────────────────────────
! This file is part of gfnff.
!
! Copyright (C) 2023-2026 Philipp Pracht
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
! along with gfnff. If not, see <https://www.gnu.org/licenses/>.
! ──────────────────────────────────────────────────────────────────────────────
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert, Sebastian Spicher, Stefan Grimme
!> at https://github.com/grimme-lab/xtb
! ──────────────────────────────────────────────────────────────────────────────

!> Phases lifted out of gfnff_ini.
!>
!> Each routine here was a block inside that one very long subroutine, marked
!> off by a banner comment and communicating with the rest only through a
!> handful of variables. They are gathered in one module because that is what
!> they are -- consecutive steps of one setup -- and separated into routines
!> because their inputs and outputs are now stated instead of implied.
!>
!> ── hydrogen and halogen bond participant lists ─────────────────────────────
!>
!> Three lists come out of here and all three are consumed by the HB/XB energy
!> terms rather than by the rest of the topology:
!>
!>   hbbas/hbaci  per-atom basicity and acidity, element parameters adjusted
!>                for the local environment (carbene, carbonyl, nitro, amide)
!>   hbatHl       the hydrogens that can donate, i.e. bonded to N/O/C and
!>                carrying enough positive charge
!>   hbatABl      the acceptor/donor heavy atom pairs
!>   xbatABl      the A-X...B triples for halogen bonds, which unlike the
!>                hydrogen bond lists also carry the two cell translations
!>
!> The charge thresholds are applied against the topological charges qa, not
!> the geometry-dependent ones, so the lists are a property of the topology
!> and do not have to be rebuilt when the geometry moves a little.
module gfnff_topo_iniphases
  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData,TGFFTopology,TGFFGenerator,TCell
  use gfnff_neighbor,only:TNeigh
  use gfnff_topo_hbset,only:hbonds
  use gfnff_topo_predicates,only:xatom,amideH
  use gfnff_topo_rings,only:getring36
  use gfnff_geometry,only:lin,valijklff,omegaPBC,banglPBC
  use gfnff_rab,only:gfnffrab,itabrow6
  use gfnff_param,only:pse
  use gfnff_qm,only:gfnffqmsolve
  use gfnff_topo_hbset,only:bond_hbset,bond_hbset0, &
    &                       bond_hb_AHB_set,bond_hb_AHB_set1,bond_hb_AHB_set0
  use gfnff_topo_rings,only:ringsbond,ringsbend,ringstors,ringstorl,chktors,ringsatom, &
    &                       ssort
  use gfnff_topo_predicates,only:ctype,alphaCO,amide
!$ use omp_lib
  implicit none(type,external)
  private

  public :: set_hb_xb_lists
  public :: perceive_rings,set_bonded_triples,set_pair_exponents
  public :: set_bonded_parameters
  public :: hueckel_solve

!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!

  subroutine set_hb_xb_lists(nat,at,xyz,cell,param,gen,itag,piadr2,hbthr2, &
        & printlevel,printunit,neigh,topo)
    !***********************************************************************
    !* Build the hydrogen and halogen bond participant lists.
    !* Input:
    !*   nat/at/xyz - system definition
    !*   cell       - lattice, needed for the translation vector refresh
    !*   param/gen  - GFN-FF parameters and generator thresholds
    !*   itag       - carbene tag from the neighbour list setup
    !*   piadr2     - pi atom assignment, gates carbon acceptors
    !*   hbthr2     - squared long-range HB cutoff
    !*   printlevel - verbosity
    !*   printunit  - output unit
    !* In/out:
    !*   neigh      - neighbour data; the translation vectors are regenerated
    !*                for the HB cutoff at the end
    !*   topo       - topology; hbbas, hbaci, hbatHl, hbatABl and xbatABl are
    !*                allocated and filled, with their counts
    !***********************************************************************
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TCell),intent(in) :: cell
    type(TGFFData),intent(in) :: param
    type(TGFFGenerator),intent(in) :: gen
    integer,intent(in) :: itag(nat)
    integer,intent(in) :: piadr2(nat)
    real(wp),intent(in) :: hbthr2
    integer,intent(in) :: printlevel,printunit
    type(TNeigh),intent(inout) :: neigh
    type(TGFFTopology),intent(inout) :: topo

    integer :: i,j,m,ia,ix,ati,nn,iTr,iTrj,iTrDum
    integer :: myunit
    real(wp) :: ff,hbpi(2),hbpj(2)
    integer,allocatable :: locarr(:,:)

    myunit = printunit

    allocate (topo%hbbas(nat),source=1.0d0)
    do i = 1,nat
      nn = sum(neigh%nb(neigh%numnb,i,:))
      ati = at(i)
      topo%hbbas(i) = param%xhbas(at(i))
      ! Carbene:
      if (ati .eq. 6.and.nn .eq. 2.and.itag(i) .eq. 1) topo%hbbas(i) = 1.46
      iTr = 0
      if (ati .eq. 8.and.nn .eq. 1) then
        call neigh%nbLoc(nat,neigh%nb,i,locarr)
        iTr = locarr(neigh%numnb,1)
        deallocate (locarr)
      end if
      if (iTr .eq. 0) cycle
      ! Carbonyl R-C=O
      if (ati .eq. 8.and.nn .eq. 1.and.at(neigh%nb(nn,i,iTr)) .eq. 6) topo%hbbas(i) = 0.68
      ! Nitro R-N=O
      if (ati .eq. 8.and.nn .eq. 1.and.at(neigh%nb(nn,i,iTr)) .eq. 7) topo%hbbas(i) = 0.47
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make list of HB donor acidity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !atom specific (not element) acidity parameters
    allocate (topo%hbaci(nat),source=1.0d0)
    do i = 1,nat
      topo%hbaci(i) = param%xhaci(at(i))
    end do
    do i = 1,nat
      ! get first nb, iTr expected to be irrelevant since same in each cell
      call neigh%jth_nb(nat,xyz,nn,1,i,iTr)          !  R - N/C/O - H  ! terminal  H
      if (nn .le. 0) cycle  ! cycle if there is no nb
      topo%hbaci(i) = param%xhaci(at(i))   ! only first neighbor of interest
      if (amideH(nat,at,topo%hyb,neigh%numnb,neigh%numctr,neigh%nb,piadr2,i,neigh)) &
              &topo%hbaci(nn) = topo%hbaci(nn)*0.80
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make list of ABs for HAB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate (topo%hbatHl(2,nat),topo%hbatABl(2,nat*(nat+1)/2),source=0)

    topo%nathbH = 0
    do i = 1,nat
      if (at(i) .ne. 1) cycle
      if (topo%hyb(i) .eq. 1) cycle      ! exclude bridging hydrogens from HB correction
      ff = gen%hqabthr
      call neigh%jth_nb(nat,xyz,j,1,i,iTr)  ! get j, first neighbor of H when j is shifted to iTr
      if (j .le. 0) cycle
      if (at(j) .gt. 10) ff = ff-0.20                ! H on heavy atoms may be negatively charged
      if (at(j) .eq. 6.and.topo%hyb(j) .eq. 3) ff = ff+0.05 ! H on sp3 C must be really positive 0.05
      if (topo%qa(i) .gt. ff) then                       ! make list of HB H atoms but only if they have a positive charge
        topo%nathbH = topo%nathbH+1
        topo%hbatHl(1,topo%nathbH) = i
        topo%hbatHl(2,topo%nathbH) = iTr
      end if
    end do
    if (printlevel >= 2) write (myunit,'(10x,"# H in HB",3x,i0)') topo%nathbH

    topo%nathbAB = 0
    do i = 1,nat
      if (at(i) .eq. 6.and.piadr2(i) .eq. 0) cycle ! C sp or sp2 pi
      ff = gen%qabthr
      if (at(i) .gt. 10) ff = ff+0.2   ! heavy atoms may be positively charged
      if (topo%qa(i) .gt. ff) cycle
      do j = 1,i-1
        ff = gen%qabthr
        if (at(j) .gt. 10) ff = ff+0.2  ! heavy atoms may be positively charged
        if (topo%qa(j) .gt. ff) cycle
        call hbonds(i,j,hbpi,hbpj,param,topo)
        if (hbpi(1)*hbpj(2) .lt. 1.d-6.and.hbpi(2)*hbpj(1) .lt. 1.d-6) cycle
        if (at(j) .eq. 6.and.piadr2(j) .eq. 0) cycle ! C sp or sp2 pi
        topo%nathbAB = topo%nathbAB+1
        topo%hbatABl(1,topo%nathbAB) = i
        topo%hbatABl(2,topo%nathbAB) = j
      end do
    end do

! make ABX list
    m = 0
    do i = 1,nat  ! potential A
      do iTr = 1,neigh%numctr
        do ia = 1,neigh%nb(neigh%numnb,i,iTr)
          ix = neigh%nb(ia,i,iTr) ! potential X (P,S,Cl,As,Se,Br,Sb,Te and I)
          if (xatom(at(ix))) then
            if (at(ix) .eq. 16.and.sum(neigh%nb(neigh%numnb,ix,:)) .gt. 2) cycle ! no sulphoxide etc S
            do iTrj = 1,neigh%numctr
              do j = 1,nat  ! loop over B:  from group 15 - group 17| in code group 5,6 or 7
                if ((i .eq. j.and.iTrj .eq. 1).or.(j .eq. ix.and.iTrj .eq. iTr)) cycle
                iTrDum = neigh%fTrSum(neigh%iTrNeg(iTrj),iTr) ! j and ix shifted -> need dummy
                if (iTrDum .gt. neigh%numctr.or.iTrDum .eq. -1) cycle  ! bpair cant handle this
                if (neigh%bpair(ix,j,iTrDum) .le. 3) cycle ! must be A...B and not X-B i.e. A-X...B
                if (param%xhbas(at(j)) .lt. 1.d-6) cycle   ! B must be O,N,...
                if (param%group(at(j)) .eq. 4) then
                  if (piadr2(j) .eq. 0.or.topo%qa(j) .gt. 0.05) cycle   ! must be a (pi)base
                end if
                m = m+1
              end do
            end do
          end if
        end do
      end do
    end do
    topo%natxbAB = m
    allocate (topo%xbatABl(5,topo%natxbAB),source=0)
    m = 0
    do i = 1,nat
      do iTr = 1,neigh%numctr
        do ia = 1,neigh%nb(neigh%numnb,i,iTr)
          ix = neigh%nb(ia,i,iTr)
          if (xatom(at(ix))) then
            if (at(ix) .eq. 16.and.sum(neigh%nb(neigh%numnb,ix,:)) .gt. 2) cycle ! no sulphoxide etc S
            do iTrj = 1,neigh%numctr
              do j = 1,nat
                if (i .eq. j.and.iTrj .eq. 1.or.j .eq. ix.and.iTrj .eq. iTr) cycle
                iTrDum = neigh%fTrSum(neigh%iTrNeg(iTrj),iTr) ! j and ix shifted -> need dummy
                if (iTrDum .gt. neigh%numctr.or.iTrDum .le. 0) cycle  ! bpair cant handle this
                if (neigh%bpair(ix,j,iTrDum) .le. 3) cycle  ! must be A...B and not X-B i.e. A-X...B
                if (param%xhbas(at(j)) .lt. 1.d-6) cycle  ! B must be O,N,...
                if (param%group(at(j)) .eq. 4) then
                  if (piadr2(j) .eq. 0.or.topo%qa(j) .gt. 0.05) cycle   ! must be a (pi)base
                end if
                m = m+1
                topo%xbatABl(1,m) = i    ! A
                topo%xbatABl(2,m) = j    ! B
                topo%xbatABl(3,m) = ix   ! X
                topo%xbatABl(4,m) = iTrj ! iTrB
                topo%xbatABl(5,m) = iTr  ! iTrX
              end do
            end do
          end if
        end do
      end do
    end do
    call neigh%getTransVec(nat,at,xyz,cell,sqrt(hbthr2))

  end subroutine set_hb_xb_lists

!========================================================================================!

  subroutine perceive_rings(nat,at,xyz,cell,printlevel,printunit,neigh,sring,cring)
    !***********************************************************************
    !* Smallest ring through every atom.
    !* getring36 walks the neighbour list, so for a periodic system the
    !* neighbours living in other cells have to be folded into one list
    !* first: nbrngs is neigh%nbm with the images of each atom appended to
    !* its central cell entry. A ring that closes through a cell boundary is
    !* invisible otherwise.
    !* Input:
    !*   nat/at/xyz - system definition
    !*   cell       - lattice; only npbc is read, to decide whether to fold
    !*   printlevel/printunit - verbosity and output unit
    !* In/out:
    !*   neigh      - neighbour data; nbm is consumed and deallocated here
    !* Output:
    !*   sring      - ring size per atom and ring index
    !*   cring      - the member atoms of each of those rings
    !***********************************************************************
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TCell),intent(in) :: cell
    integer,intent(in) :: printlevel,printunit
    type(TNeigh),intent(inout) :: neigh
    integer,intent(inout) :: sring(:,:)
    integer,intent(inout) :: cring(:,:,:)

    integer :: i,j,nni,iTr
    integer :: cr(10,20),sr(20)
    integer :: myunit
    integer,allocatable :: nbrngs(:,:)

    myunit = printunit

    allocate (nbrngs(neigh%numnb,nat),source=0)
    nbrngs = neigh%nbm(:,:,1)
    if (cell%npbc .ne. 0) then
      do i = 1,nat
        nni = neigh%nbm(neigh%numnb,i,1)
        do iTr = 2,neigh%numctr
          do j = 1,neigh%nbm(neigh%numnb,i,iTr)
            ! append neighbors from other cells to cell one
            nbrngs(nni+j,i) = neigh%nbm(j,i,iTr)
            ! adjust number of nb
            nbrngs(neigh%numnb,i) = nbrngs(neigh%numnb,i)+1
          end do
          nni = nni+neigh%nbm(neigh%numnb,i,iTr)
        end do
      end do
    end if
    if (printlevel >= 2) write (myunit,'(10x,"rings ...")')
!$omp parallel default(none) private(i,cr,sr) shared(nat,at,xyz,neigh,nbrngs,cring,sring)
!$omp do
    do i = 1,nat
      call getring36(nat,at,neigh%numnb,neigh%numctr,nbrngs,i,cr,sr)
      cring(1:10,1:20,i) = cr(1:10,1:20)
      sring(1:20,i) = sr(1:20)
    end do
!$omp end do
!$omp end parallel
    deallocate (neigh%nbm,nbrngs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine perceive_rings

!========================================================================================!

  subroutine set_bonded_triples(nat,printlevel,printunit,neigh,topo,io)
    !***********************************************************************
    !* Bonded atom triples for the three-body (ATM) term.
    !* Collects i-j-k where i and j are exactly three bonds apart and k is a
    !* neighbour of either, i.e. the triples that the bend and torsion terms
    !* do not already cover.
    !* Input:
    !*   nat        - number of atoms
    !*   printlevel/printunit - verbosity and output unit
    !*   neigh      - neighbour data, provides bpair and the neighbour lists
    !* In/out:
    !*   topo       - topology; b3list is allocated and filled, nbatm set
    !* Output:
    !*   io         - non-zero if the triple count overran the estimate
    !***********************************************************************
    integer,intent(in) :: nat
    integer,intent(in) :: printlevel,printunit
    type(TNeigh),intent(in) :: neigh
    type(TGFFTopology),intent(inout) :: topo
    integer,intent(out) :: io

    character(len=*),parameter :: source = 'set_bonded_triples'
    integer :: i,j,k,m,idum,iTr,iTr2
    integer :: myunit

    myunit = printunit
    io = 0

    idum = 1000*nat
    allocate (topo%b3list(5,idum),source=0)
    topo%nbatm = 0
    do i = 1,nat
      do j = 1,i-1 !
        do iTr = 1,neigh%numctr
          if (neigh%bpair(j,i,iTr) .eq. 3) then  ! 1,4 exclusion of back-pair makes it worse, 1,3 makes little effect
            do iTr2 = 1,neigh%numctr
              do m = 1,neigh%nb(neigh%numnb,j,iTr2)
                k = neigh%nb(m,j,iTr2)
                if ((i == k.or.j == k).and.iTr == 1.and.iTr2 == 1) cycle
                topo%nbatm = topo%nbatm+1
                topo%b3list(1,topo%nbatm) = i                      !in central cell
                topo%b3list(2,topo%nbatm) = j
                topo%b3list(3,topo%nbatm) = k
                topo%b3list(4,topo%nbatm) = iTr                    !iTrj
                topo%b3list(5,topo%nbatm) = neigh%fTrSum(iTr,iTr2) !iTrk
              end do
              do m = 1,neigh%nb(neigh%numnb,i,iTr2)
                k = neigh%nb(m,i,iTr2)
                if ((i == k.or.j == k).and.iTr == 1.and.iTr2 == 1) cycle
                topo%nbatm = topo%nbatm+1
                topo%b3list(1,topo%nbatm) = i
                topo%b3list(2,topo%nbatm) = j
                topo%b3list(3,topo%nbatm) = k
                topo%b3list(4,topo%nbatm) = iTr
                topo%b3list(5,topo%nbatm) = iTr2
              end do
            end do
          end if
        end do
      end do
    end do
    if (topo%nbatm .gt. idum) then
      if (printlevel >= 1) then
        write (myunit,*) idum,topo%nbatm
        write (myunit,'("**ERROR** ",a,1x,a)') 'overflow in ini',source
      end if
      return
    end if
    if (printlevel >= 2) write (myunit,'(10x,"# BATM",3x,i0)') topo%nbatm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine set_bonded_triples

!========================================================================================!

  subroutine set_pair_exponents(nat,at,param,gen,neigh,topo)
    !***********************************************************************
    !* Screening exponents for the non-bonded repulsion, one per atom pair,
    !* and the D4 zeta charge scaling for the dispersion C6.
    !*
    !* The exponent carries no cell dependence. Only the H...H case ever
    !* varied with the cell, and only through the number of bonds between the
    !* two atoms; that factor is applied where the exponent is used and lives
    !* in topo%hhrep, so what is stored here is the bare product.
    !* Input:
    !*   nat/at     - system definition
    !*   param/gen  - GFN-FF parameters and generator scaling factors
    !*   neigh      - neighbour data, for the coordination-dependent decrease
    !* In/out:
    !*   topo       - topology; alphanb, zetac6 and hhrep are filled. qa must
    !*                already hold the topological charges.
    !***********************************************************************
    integer,intent(in) :: nat,at(nat)
    type(TGFFData),intent(in) :: param
    type(TGFFGenerator),intent(in) :: gen
    type(TNeigh),intent(in) :: neigh
    type(TGFFTopology),intent(inout) :: topo

    integer :: i,j,ij,ati,atj
    real(wp) :: f1,f2,ff,fn,dum1,dum2

    ! H...H repulsion scaling, looked up by the number of bonds between the
    ! pair. Only the 1,3 and 1,4 cases are scaled beyond the plain H...H factor
    topo%hhrep(:) = gen%hhfac
    topo%hhrep(2) = gen%hhfac*gen%hh13rep   ! 1,3 case
    topo%hhrep(3) = gen%hhfac*gen%hh14rep   ! 1,4 case, important for torsions

    do i = 1,nat
      ati = at(i)
      fn = 1.0d0+gen%nrepscal/(1.0d0+dble(sum(neigh%nb(neigh%numnb,i,:)))**2)
      dum1 = param%repan(ati)*(1.d0+topo%qa(i)*gen%qrepscal)*fn ! a small but physically correct decrease of repulsion with q
      f1 = zeta(ati,topo%qa(i))
      do j = 1,i
        atj = at(j)
        fn = 1.0d0+gen%nrepscal/(1.0d0+dble(sum(neigh%nb(neigh%numnb,j,:))**2))
        dum2 = param%repan(atj)*(1.d0+topo%qa(j)*gen%qrepscal)*fn
        f2 = zeta(atj,topo%qa(j))
        ij = lin(j,i) ! for zetac6
        ! ── one exponent per pair, no cell dependence ───────────────────────
        ! Only the H...H case ever varied with the cell, and only through the
        ! bond count between the two atoms. That factor is applied where the
        ! exponent is used, so for H...H the bare product is stored here and
        ! the whole H...H scaling, gen%hhfac included, sits in topo%hhrep.
        if (ati .eq. 1.and.atj .eq. 1) then
          ff = 1.0d0                                 ! see topo%hhrep
        else
          ff = 1.0d0
if ((ati .eq. 1.and.param%metal(atj) .gt. 0).or.(atj .eq. 1.and.param%metal(ati) .gt. 0)) ff = 0.85 ! M...H
          if ((ati .eq. 1.and.atj .eq. 6).or.(atj .eq. 1.and.ati .eq. 6)) ff = 0.91 ! C...H, good effect
          if ((ati .eq. 1.and.atj .eq. 8).or.(atj .eq. 1.and.ati .eq. 8)) ff = 1.04 ! O...H, good effect
        end if
        topo%alphanb(ij) = sqrt(dum1*dum2)*ff
        topo%zetac6(ij) = f1*f2  ! D4 zeta scaling using qref=0
      end do
    end do

  end subroutine set_pair_exponents

!========================================================================================!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)

!> @brief charge scaling function
  pure elemental function zeta(at,q)
    implicit none
    integer,intent(in) :: at
    real(wp),intent(in) :: q

    real(wp)           :: zeta,qmod
    real(wp),parameter :: zeff(103) = (/ &
    &   1,2,  & ! H-He
    &   3,4,5,6,7,8,9,10,  & ! Li-Ne
    &  11,12,13,14,15,16,17,18,  & ! Na-Ar
    &  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,  & ! K-Kr
    &   9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,  & ! Rb-Xe
    &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Cs-Lu
    &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26, &  ! Hf-Rn
    &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43  & ! Fr-Lr
    &/)
!! Semiempirical Evaluation of the GlobalHardness of the Atoms of 103
!! Elements of the Periodic Table Using the Most Probable Radii as
!! their Size Descriptors DULAL C. GHOSH, NAZMUL ISLAM 2009 in
!! Wiley InterScience (www.inte"rscience.wiley.com).
!! DOI 10.1002/qua.22202
!! values in the paper multiplied by two because
!! (ii:ii)=(IP-EA)=d^2 E/dN^2 but the hardness
!! definition they use is 1/2d^2 E/dN^2 (in Eh)
    real(wp),parameter :: c(1:103) = (/ &
   &0.47259288_wp,0.92203391_wp,0.17452888_wp,0.25700733_wp,0.33949086_wp,0.42195412_wp, & ! H-C
   &0.50438193_wp,0.58691863_wp,0.66931351_wp,0.75191607_wp,0.17964105_wp,0.22157276_wp, & ! N-Mg
   &0.26348578_wp,0.30539645_wp,0.34734014_wp,0.38924725_wp,0.43115670_wp,0.47308269_wp, & ! Al-Ar
   &0.17105469_wp,0.20276244_wp,0.21007322_wp,0.21739647_wp,0.22471039_wp,0.23201501_wp, & ! Ca-Cr
   &0.23933969_wp,0.24665638_wp,0.25398255_wp,0.26128863_wp,0.26859476_wp,0.27592565_wp, & ! Mn-Zn
   &0.30762999_wp,0.33931580_wp,0.37235985_wp,0.40273549_wp,0.43445776_wp,0.46611708_wp, & ! Ga-Kr
   &0.15585079_wp,0.18649324_wp,0.19356210_wp,0.20063311_wp,0.20770522_wp,0.21477254_wp, & ! Rb-Mo
   &0.22184614_wp,0.22891872_wp,0.23598621_wp,0.24305612_wp,0.25013018_wp,0.25719937_wp, & ! Tc-Cd
   &0.28784780_wp,0.31848673_wp,0.34912431_wp,0.37976593_wp,0.41040808_wp,0.44105777_wp, & ! In-Xe
   &0.05019332_wp,0.06762570_wp,0.08504445_wp,0.10247736_wp,0.11991105_wp,0.13732772_wp, & ! Cs-Nd
   &0.15476297_wp,0.17218265_wp,0.18961288_wp,0.20704760_wp,0.22446752_wp,0.24189645_wp, & ! Pm-Dy
   &0.25932503_wp,0.27676094_wp,0.29418231_wp,0.31159587_wp,0.32902274_wp,0.34592298_wp, & ! Ho-Hf
   &0.36388048_wp,0.38130586_wp,0.39877476_wp,0.41614298_wp,0.43364510_wp,0.45104014_wp, & ! Ta-Pt
   &0.46848986_wp,0.48584550_wp,0.12526730_wp,0.14268677_wp,0.16011615_wp,0.17755889_wp, & ! Au-Po
   &0.19497557_wp,0.21240778_wp,& ! At, Rn
   &0.07263133_wp,0.09421788_wp,0.09920108_wp,0.10418429_wp,0.14235212_wp,0.16393866_wp, & ! Fr-U
   &0.18675998_wp,0.22370039_wp,0.25113742_wp,0.25026279_wp,0.28843797_wp,0.31002451_wp, & ! Np-Cf
   &0.33159636_wp,0.35316820_wp,0.36822807_wp,0.39634864_wp,0.40135389_wp & ! Es-Lr
   &/)

    intrinsic :: exp

    qmod = zeff(at)+q
    if (qmod .lt. 0._wp) then
      zeta = exp(3.0d0)
    else
      zeta = exp(3.0d0*(1._wp-exp(c(at)*(1._wp-zeff(at)/qmod))))
    end if

  end function zeta

!========================================================================================!

  subroutine set_bonded_parameters(nat,at,xyz,cell,param,gen,cn,rab,rtmp,mchar, &
        & pbo,pibo,piadr,imetal,itag,btyp,sring,cring,hbthr1,hbthr2, &
        & printlevel,printunit,neigh,topo,io)
    !***********************************************************************
    !* Assign the parameters of every bonded term: bond stretch, bend,
    !* torsion and out-of-plane, plus the special torsion around carbon
    !* triple bonds and the hydrogen-bridge force constant scaling.
    !*
    !* This is one routine and not five because the five steps share their
    !* scratch. They were consecutive blocks of gfnff_ini reusing the same
    !* loop counters and the same temporaries, and splitting them apart
    !* before that scratch is scoped would mean threading a dozen meaningless
    !* variables through the interface. Scoped here first, they can be split
    !* further without paying that cost.
    !*
    !* What actually crosses the boundary is listed below; everything else
    !* these blocks touch is local to this routine.
    !* Input:
    !*   nat/at/xyz - system definition
    !*   cell       - lattice, npbc gates the periodic hydrogen bridge setup
    !*   param/gen  - GFN-FF parameters and generator thresholds
    !*   cn         - coordination numbers
    !*   rab        - packed interatomic distances
    !*   mchar      - metal character per atom
    !*   pbo/pibo   - pi bond orders from the Hueckel treatment
    !*   piadr      - pi atom assignment
    !*   imetal     - metal classification per atom
    !*   itag       - carbene/eta tag from the neighbour setup
    !*   sring/cring- smallest ring size and members per atom
    !*   hbthr1/2   - squared hydrogen bond cutoffs
    !*   printlevel/printunit - verbosity and output unit
    !* In/out:
    !*   rtmp       - scratch for the RAB guess, sized nat*(nat+1)/2
    !*   btyp       - bond type per bond; allocated by the caller, filled here
    !*   neigh      - neighbour data; vbond and the angle/torsion lists grow
    !*   topo       - topology; alist, vangl, tlist, vtors, sTorsl and the
    !*                hydrogen bridge maps are filled
    !* Output:
    !*   io         - non-zero if a list overran its allocation
    !***********************************************************************
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TCell),intent(in) :: cell
    type(TGFFData),intent(in) :: param
    type(TGFFGenerator),intent(in) :: gen
    real(wp),intent(in) :: cn(nat)
    real(wp),intent(in) :: rab(:)
    real(wp),intent(inout) :: rtmp(:)
    real(wp),intent(in) :: mchar(nat)
    real(wp),intent(in) :: pbo(:),pibo(:)
    integer,intent(in) :: piadr(nat)
    integer,intent(in) :: imetal(nat),itag(nat)
    integer,intent(inout) :: btyp(:)
    integer,intent(in) :: sring(:,:),cring(:,:,:)
    real(wp),intent(in) :: hbthr1,hbthr2
    integer,intent(in) :: printlevel,printunit
    type(TNeigh),intent(inout) :: neigh
    type(TGFFTopology),intent(inout) :: topo
    integer,intent(out) :: io

    character(len=*),parameter :: source = 'set_bonded_parameters'
    real(wp),parameter :: pi = 3.1415926535897932385_wp

    integer :: ati,atj,atk,i,j,k,nn,ii,jj,kk,ll,m,rings,ia,ja,ij,nnn,idum,no,nbi
    integer :: nni,nnj
    integer :: ineig,jneig,nrot,bbtyp,hybi,hybj,nh,nc
    integer :: ringsi,ringsj,ringsk,ringl,npi,maxtors,rings4,nheav
    integer :: ncarbo,mtyp1,mtyp2
    integer :: nf,nsi,nmet,nhi,nhj
    integer :: AHB_nr,bond_hbn
    integer :: iTr,iTr2,iTrj,iTrk,iTrlDum,iTrl
    integer :: myunit
    real(wp) :: vTrl(3),vTrj(3),vTrk(3)
    real(wp) :: r0,ff,f1,f2,phi,ringf,fcn
    real(wp) :: shift,dum,qafac,fqq,feta
    real(wp) :: sumppi,fpi,fxh,fijk,fsrb2
    real(wp) :: fheavy,fn,fctot,fij
    real(wp) :: bstrength
    real(wp) :: fkl,fbsmall
    logical :: lring,picon,notpicon,bridge,sp3ij,ccij
    logical :: triple,sp3kl
    integer :: cDbl(4,42),cd,cdi
    integer :: ind3(3)
    real(wp) :: sdum3(3)
    logical :: tDbl
    integer,allocatable :: locarr(:,:)
    integer,allocatable :: lin_AHB(:,:)
    integer,allocatable :: bond_hbl(:,:)

    myunit = printunit
    io = 0

    call gfnffrab(nat,at,cn,rtmp)           ! guess RAB for output

    topo%nbond_vbond = neigh%nbond
    allocate (neigh%vbond(3,neigh%nbond),source=0.0d0)

    if (printlevel >= 2) then
      write (myunit,*)
      write (myunit,'(10x,"#atoms :",3x,i0)') nat
      write (myunit,'(10x,"#bonds :",3x,i0)') neigh%nbond
    end if
    if (printlevel >= 3) then
      write (myunit,*)
     write (myunit,*) 'bond atoms        type  in ring    R      R0    piBO    fqq  kbond(tot)  alp'
    end if

    do i = 1,neigh%nbond
      jj = neigh%blist(1,i)
      ii = neigh%blist(2,i)
      nni = sum(neigh%nb(neigh%numnb,ii,:))
      nnj = sum(neigh%nb(neigh%numnb,jj,:))
      ij = lin(ii,jj)
      ia = at(ii)
      ja = at(jj)
      call ringsbond(nat,ii,jj,cring,sring,rings)
      shift = 0.d0
      fxh = 1.d0
      ringf = 1.d0
      fqq = 1.d0
      fpi = 1.d0
      fheavy = 1.d0
      fheavy = 1.d0
      fcn = 1.d0
      fsrb2 = gen%srb2
      bridge = .false.
      shift = 0.d0
! assign bond type
      btyp(i) = 1 ! single
      if (topo%hyb(ii) .eq. 2.and.topo%hyb(jj) .eq. 2) btyp(i) = 2 ! sp2-sp2 = pi
      if (topo%hyb(ii) .eq. 3.and.topo%hyb(jj) .eq. 2.and.ia .eq. 7) btyp(i) = 2 ! N-sp2
      if (topo%hyb(jj) .eq. 3.and.topo%hyb(ii) .eq. 2.and.ja .eq. 7) btyp(i) = 2 ! N-sp2
      if (topo%hyb(ii) .eq. 1.or.topo%hyb(jj) .eq. 1) btyp(i) = 3 ! sp-X i.e. no torsion
      if ((param%group(ia) .eq. 7.or.ia .eq. 1).and.topo%hyb(ii) .eq. 1) then
        btyp(i) = 3 ! linear halogen i.e. no torsion
        bridge = .true.
      end if
      if ((param%group(ja) .eq. 7.or.ja .eq. 1).and.topo%hyb(jj) .eq. 1) then
        btyp(i) = 3 ! linear halogen i.e. no torsion
        bridge = .true.
      end if
      if (topo%hyb(ii) .eq. 5.or.topo%hyb(jj) .eq. 5) btyp(i) = 4 ! hypervalent
      if (imetal(ii) .gt. 0.or.imetal(jj) .gt. 0) btyp(i) = 5 ! metal
      if (imetal(ii) .eq. 2.and.imetal(jj) .eq. 2) btyp(i) = 7 ! TM metal-metal
      if (imetal(jj) .eq. 2.and.itag(ii) .eq. -1.and.piadr(ii) .gt. 0) btyp(i) = 6 ! eta
      if (imetal(ii) .eq. 2.and.itag(jj) .eq. -1.and.piadr(jj) .gt. 0) btyp(i) = 6 ! eta
      bbtyp = btyp(i)
! normal bond
      if (bbtyp .lt. 5) then
        hybi = max(topo%hyb(ii),topo%hyb(jj))
        hybj = min(topo%hyb(ii),topo%hyb(jj))
        if (hybi .eq. 5.or.hybj .eq. 5) then
          bstrength = gen%bstren(4)                                       ! base value hypervalent
        else
          bstrength = gen%bsmat(hybi,hybj)                                ! base value normal hyb
        end if
        if (hybi .eq. 3.and.hybj .eq. 2.and.(ia .eq. 7.or.ja .eq. 7)) &
 &                                      bstrength = gen%bstren(2)*1.04   ! N-sp2

        if (bridge) then
          if (param%group(ia) .eq. 7) bstrength = gen%bstren(1)*0.50d0 ! bridging X
          if (param%group(ja) .eq. 7) bstrength = gen%bstren(1)*0.50d0 ! bridging X
          if (ia .eq. 1.or.ia .eq. 9) bstrength = gen%bstren(1)*0.30d0 ! bridging H/F
          if (ja .eq. 1.or.ja .eq. 9) bstrength = gen%bstren(1)*0.30d0 ! bridging H/F
        end if
        if (bbtyp .eq. 4) shift = gen%hyper_shift          ! hypervalent
        if (ia .eq. 1.or.ja .eq. 1) shift = gen%rabshifth            ! XH
        if (ia .eq. 9.and.ja .eq. 9) shift = 0.22                 ! f2
        if (topo%hyb(ii) .eq. 3.and.topo%hyb(jj) .eq. 0) shift = shift-0.022         ! X-sp3
        if (topo%hyb(ii) .eq. 0.and.topo%hyb(jj) .eq. 3) shift = shift-0.022         ! X-sp3
        if (topo%hyb(ii) .eq. 1.and.topo%hyb(jj) .eq. 0) shift = shift+0.14          ! X-sp
        if (topo%hyb(ii) .eq. 0.and.topo%hyb(jj) .eq. 1) shift = shift+0.14          ! X-sp
        if ((ia .eq. 1.and.ja .eq. 6)) then
          call ringsatom(nat,jj,cring,sring,ringsj)
          if (ringsj .eq. 3) fxh = 1.05    ! 3-ring CH
          if (ctype(nat,at,neigh%numnb,neigh%numctr,neigh%nb,piadr,jj) .eq. 1) fxh = 0.95    ! aldehyd CH
        end if
        if ((ia .eq. 6.and.ja .eq. 1)) then
          call ringsatom(nat,ii,cring,sring,ringsi)
          if (ringsi .eq. 3) fxh = 1.05    ! 3-ring CH
          if (ctype(nat,at,neigh%numnb,neigh%numctr,neigh%nb,piadr,ii) .eq. 1) fxh = 0.95    ! aldehyd CH
        end if
        if ((ia .eq. 1.and.ja .eq. 5)) fxh = 1.10    ! BH
        if ((ja .eq. 1.and.ia .eq. 5)) fxh = 1.10    !
        if ((ia .eq. 1.and.ja .eq. 7)) fxh = 1.06    ! NH
        if ((ja .eq. 1.and.ia .eq. 7)) fxh = 1.06    !
        if ((ia .eq. 1.and.ja .eq. 8)) fxh = 0.93    ! OH
        if ((ja .eq. 1.and.ia .eq. 8)) fxh = 0.93    !
        if (bbtyp .eq. 3.and.ia .eq. 6.and.ja .eq. 8) bstrength = gen%bstren(3)*0.90d0 ! makes CO right and M-CO reasonable
        if (bbtyp .eq. 3.and.ia .eq. 8.and.ja .eq. 6) bstrength = gen%bstren(3)*0.90d0 !
!           modify locally for triple bonds
        if (bbtyp .eq. 3.and.(topo%hyb(ii) .eq. 0.or.topo%hyb(jj) .eq. 0)) bbtyp = 1 ! sp-sp3
        if (bbtyp .eq. 3.and.(topo%hyb(ii) .eq. 3.or.topo%hyb(jj) .eq. 3)) bbtyp = 1 ! sp-sp3
        if (bbtyp .eq. 3.and.(topo%hyb(ii) .eq. 2.or.topo%hyb(jj) .eq. 2)) bbtyp = 2 ! sp-sp2
!           Pi stuff
        if (pibo(i) .gt. 0) then
          shift = gen%hueckelp*(gen%bzref-pibo(i)) ! ref value = no correction is benzene, P=2/3
          if (bbtyp .ne. 3.and.pibo(i) .gt. 0.1) then
            btyp(i) = 2
            bbtyp = 2
          end if
          fpi = 1.0d0-gen%hueckelp2*(gen%bzref2-pibo(i)) ! deepness
        end if
        if (ia .gt. 10.and.ja .gt. 10) then
          fcn = fcn/(1.0d0+0.007*dble(nni)**2)
          fcn = fcn/(1.0d0+0.007*dble(nnj)**2)
        end if
        qafac = topo%qa(ii)*topo%qa(jj)*70.0d0
        fqq = 1.0_wp+gen%qfacbm0/(1.0_wp+exp(15.0_wp*qafac))
! metal involed
      else
        shift = 0
        bstrength = gen%bstren(bbtyp)
        if (bbtyp .eq. 7) then ! TM-TM
          if (itabrow6(ia) .gt. 4.and.itabrow6(ja) .gt. 4) bstrength = gen%bstren(8) ! 4/5d-4/5d
          if (itabrow6(ia) .eq. 4.and.itabrow6(ja) .gt. 4) bstrength = gen%bstren(9) ! 3d-4/5d
          if (itabrow6(ja) .eq. 4.and.itabrow6(ia) .gt. 4) bstrength = gen%bstren(9) ! 3d-4/5d
          dum = 2.0d0*mchar(ii)+2.0d0*mchar(jj)
          dum = min(dum,0.5d0)  ! limit the "metallic" correction
          bstrength = bstrength*(1.0d0-dum)
        end if
        mtyp1 = 0  ! no metal
        mtyp2 = 0
        if (param%group(ia) .eq. 1) mtyp1 = 1  ! Li...
        if (param%group(ia) .eq. 2) mtyp1 = 2  ! Be...
        if (param%group(ia) .gt. 2.and.imetal(ii) .eq. 1) mtyp1 = 3  ! main group
        if (imetal(ii) .eq. 2) mtyp1 = 4  ! TM
        if (param%group(ja) .eq. 1) mtyp2 = 1  ! Li...
        if (param%group(ja) .eq. 2) mtyp2 = 2  ! Be...
        if (param%group(ja) .gt. 2.and.imetal(jj) .eq. 1) mtyp2 = 3  ! main group
        if (imetal(jj) .eq. 2) mtyp2 = 4  ! TM
        qafac = topo%qa(ii)*topo%qa(jj)*25.0d0
        dum = 1.0_wp/(1.0_wp+exp(15.0_wp*qafac))
        fqq = 1.0d0+dum*(gen%qfacbm(mtyp1)+gen%qfacbm(mtyp2))*0.5   ! metal charge corr.
        if (imetal(ii) .eq. 2.and.ja .gt. 10) fheavy = 0.65d0 ! heavy gen. ligand
        if (imetal(jj) .eq. 2.and.ia .gt. 10) fheavy = 0.65d0
        if (imetal(ii) .eq. 2.and.ja .eq. 15) fheavy = 1.60d0 ! P ligand
        if (imetal(jj) .eq. 2.and.ia .eq. 15) fheavy = 1.60d0
        if (imetal(ii) .eq. 2.and.param%group(ja) .eq. 6) fheavy = 0.85d0 ! chalcogen ligand
        if (imetal(jj) .eq. 2.and.param%group(ia) .eq. 6) fheavy = 0.85d0
        if (imetal(ii) .eq. 2.and.param%group(ja) .eq. 7) fheavy = 1.30d0 ! halogen ligand
        if (imetal(jj) .eq. 2.and.param%group(ia) .eq. 7) fheavy = 1.30d0
        if (imetal(ii) .eq. 2.and.ja .eq. 1.and.itabrow6(ia) .le. 5) fxh = 0.80d0 ! hydrogen 3d/4d
        if (imetal(jj) .eq. 2.and.ia .eq. 1.and.itabrow6(ja) .le. 5) fxh = 0.80d0 ! hydrogen 3d/4d
        if (imetal(ii) .eq. 2.and.ja .eq. 1.and.itabrow6(ia) .gt. 5) fxh = 1.00d0 ! hydrogen 5d
        if (imetal(jj) .eq. 2.and.ia .eq. 1.and.itabrow6(ja) .gt. 5) fxh = 1.00d0 ! hydrogen 5d
        if (imetal(ii) .eq. 1.and.ja .eq. 1) fxh = 1.20d0
        if (imetal(jj) .eq. 1.and.ia .eq. 1) fxh = 1.20d0
        if (imetal(jj) .eq. 2.and.topo%hyb(ii) .eq. 1) then !CO/CN/NC...
          if (ia .eq. 6) then
            fpi = 1.5d0
            shift = -0.45d0
          end if
          if (ia .eq. 7.and.nni .ne. 1) then
            fpi = 0.4d0
            shift = 0.47d0
          end if
        end if
        if (imetal(ii) .eq. 2.and.topo%hyb(jj) .eq. 1) then !CO/CN/NC...
          if (ja .eq. 6) then
            fpi = 1.5d0
            shift = -0.45d0
          end if
          if (ja .eq. 7.and.nnj .ne. 1) then
            fpi = 0.4d0
            shift = 0.47d0
          end if
        end if
        if (imetal(ii) .eq. 2) shift = shift+gen%metal2_shift   ! metal shift TM
        if (imetal(jj) .eq. 2) shift = shift+gen%metal2_shift   !
        if (imetal(ii) .eq. 1.and.param%group(ia) .le. 2) shift = shift+gen%metal1_shift   ! metal shift group 1+2
        if (imetal(jj) .eq. 1.and.param%group(ja) .le. 2) shift = shift+gen%metal1_shift   !
        if (mtyp1 .eq. 3) shift = shift+gen%metal3_shift   ! metal shift MG
        if (mtyp2 .eq. 3) shift = shift+gen%metal3_shift   !
        if (bbtyp .eq. 6.and.param%metal(ia) .eq. 2) shift = shift+gen%eta_shift*nni! eta coordinated
        if (bbtyp .eq. 6.and.param%metal(ja) .eq. 2) shift = shift+gen%eta_shift*nnj! eta coordinated
        if (mtyp1 .gt. 0.and.mtyp1 .lt. 3) fcn = fcn/(1.0d0+0.100*dble(nni)**2)
        if (mtyp2 .gt. 0.and.mtyp2 .lt. 3) fcn = fcn/(1.0d0+0.100*dble(nnj)**2)
        if (mtyp1 .eq. 3) fcn = fcn/(1.0d0+0.030*dble(nni)**2)
        if (mtyp2 .eq. 3) fcn = fcn/(1.0d0+0.030*dble(nnj)**2)
        if (mtyp1 .eq. 4) fcn = fcn/(1.0d0+0.036*dble(nni)**2)
        if (mtyp2 .eq. 4) fcn = fcn/(1.0d0+0.036*dble(nnj)**2)
        if (mtyp1 .eq. 4.or.mtyp2 .eq. 4) then
          fsrb2 = -gen%srb2*0.22! weaker, inverse EN dep. for TM metals
        else
          fsrb2 = gen%srb2*0.28! "normal" for other metals
        end if
      end if

      if (ia .gt. 10.and.ja .gt. 10) then  ! both atoms are heavy
        shift = shift+gen%hshift3
        if (ia .gt. 18) shift = shift+gen%hshift4
        if (ja .gt. 18) shift = shift+gen%hshift4
        if (ia .gt. 36) shift = shift+gen%hshift5
        if (ja .gt. 36) shift = shift+gen%hshift5
      end if

! shift
      neigh%vbond(1,i) = gen%rabshift+shift   ! value for all bonds + special part

! RINGS prefactor
      if (rings .gt. 0) ringf = 1.0d0+gen%fringbo*(6.0d0-dble(rings))**2  ! max ring size is 6

! steepness
      neigh%vbond(2,i) = gen%srb1*(1.0d0+fsrb2*(param%en(ia)-param%en(ja))**2+gen%srb3*bstrength)

! tot prefactor        atoms              spec     typ       qterm    heavy-M  pi   XH(3ring,OH...) CN for M
      neigh%vbond(3,i) = -param%bond(ia)*param%bond(ja)*ringf*bstrength*fqq*fheavy*fpi*fxh*fcn
!        write(myunit,*) bond(ia),bond(ja),ringf,bstrength,fqq,fheavy,fpi,fxh

! output
      r0 = (rtmp(ij)+neigh%vbond(1,i))*0.529167
      if (printlevel >= 3) write (myunit,'(2a3,2i5,2x,2i5,2x,6f8.3)') &
  &   pse(at(ii)),pse(at(jj)),ii,jj,bbtyp,rings,0.529167*rab(ij),r0,pibo(i),fqq,neigh%vbond(3,i),neigh%vbond(2,i)
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     scale FC if bond is part of hydrogen bridge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     scale FC if bond is part of hydrogen bridge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Set up fix hblist just like for the HB term
    allocate (topo%isABH(nat),source=.false.)
    call bond_hbset0(nat,at,xyz,cell%npbc,bond_hbn,topo,neigh,hbthr1,hbthr2)
    allocate (bond_hbl(6,bond_hbn))
    allocate (neigh%nr_hb(neigh%nbond),source=0)
    call bond_hbset(nat,at,xyz,cell%npbc,bond_hbn,bond_hbl,&
         &           topo,neigh,hbthr1,hbthr2)
    !Set up AH, B and nr. of B list
    call bond_hb_AHB_set0(nat,at,neigh%nbond,bond_hbn,bond_hbl,AHB_nr,neigh)
    allocate (lin_AHB(4,0:AHB_nr),source=0)
    call bond_hb_AHB_set1(nat,at,neigh%nbond,bond_hbn,bond_hbl,AHB_nr,lin_AHB,topo%bond_hb_nr,topo%b_max,topo,neigh)
    allocate (topo%bond_hb_AH(4,topo%bond_hb_nr),source=0)
    allocate (topo%bond_hb_B(2,topo%b_max,topo%bond_hb_nr),source=0)
    allocate (topo%bond_hb_Bn(topo%bond_hb_nr),source=0)
    call bond_hb_AHB_set(nat,at,neigh%nbond,bond_hbn,bond_hbl,AHB_nr,lin_AHB,topo,neigh)

    ! create mapping from atom index to hb index, for AB and H seperately
    allocate (topo%hb_mapABH(nat),source=0)
    j = 0 ! H counter
    k = 0 ! AB counter
    do i = 1,nat
      ! check if atom i is A,B, or H
      if (topo%isABH(i)) then
        ! check if it is H
        if (at(i) .eq. 1) then
          j = j+1
          topo%hb_mapABH(i) = j
          ! then it is A or B
        else
          k = k+1
          topo%hb_mapABH(i) = k
        end if
      end if
    end do
    topo%hb_mapNAB = k
    topo%hb_mapNH = j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!               bend
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    topo%nangl = 0
    do i = 1,nat
      nn = sum(neigh%nb(neigh%numnb,i,:))                  ! take full set to include M-X-Y
      if (nn .le. 1) cycle                                  !
      if (nn .gt. 6) cycle     ! no highly coordinated atom
      cDbl = 0 ! cDbl stores j,k,iTr,iTr2 that were already considered
      cdi = 0
      ati = at(i)
      do iTr = 1,neigh%numctr
        do j = 1,neigh%nb(neigh%numnb,i,iTr)
          do iTr2 = 1,neigh%numctr !
            do k = 1,j !
              if (iTr .eq. iTr2.and.k .eq. j) cycle  !dont use same atom as both neighbors
              jj = neigh%nb(j,i,iTr)
              kk = neigh%nb(k,i,iTr2)
              if (kk .eq. 0) cycle  !only j goes over nb so k or kk might be "out of bounds"
              atj = at(jj)
              atk = at(kk)
              fijk = param%angl(ati)*param%angl2(atj)*param%angl2(atk)
              if (fijk .lt. gen%fcthr) cycle     ! too small
              ! check for double counting
              tDbl = .false.
              do cd = 1,42
                if (cDbl(1,cd) .eq. 0) exit
if (cDbl(1,cd) .eq. kk.and.cDbl(2,cd) .eq. jj.and.cDbl(3,cd) .eq. iTr2.and.cDbl(4,cd) .eq. iTr) then
                  tDbl = .true.
                  exit
                end if
              end do
              if (tDbl) cycle
              cdi = cdi+1
              cDbl(1,cdi) = jj
              cDbl(2,cdi) = kk
              cDbl(3,cdi) = iTr
              cDbl(4,cdi) = iTr2
              topo%nangl = topo%nangl+1
            end do
          end do
        end do
      end do
    end do

    if (printlevel >= 2) write (myunit,'(10x,"#angl  :",3x,i0)') topo%nangl
    if (printlevel >= 3) then
      write (myunit,*)
      write (myunit,*) 'angle atoms        phi0    phi      FC  pi rings'
    end if

    topo%nangl_alloc = topo%nangl
    allocate (topo%alist(5,topo%nangl),source=0)
    allocate (topo%vangl(2,topo%nangl),source=0.0d0)
    topo%nangl = 0
    do i = 1,nat  ! start angl_loop
      nn = sum(neigh%nb(neigh%numnb,i,:))
      if (nn .le. 1) cycle  ! no angle with only one neighbor
      if (nn .gt. 6) cycle  ! no highly coordinated systems
      cDbl = 0 ! cDbl(j,k,iTr,iTr2) with iTr->j and iTr2->k
      cdi = 0
      ii = i
      ati = at(i)
      do iTr = 1,neigh%numctr
        do j = 1,nn
          do iTr2 = 1,neigh%numctr
            do k = 1,j
              if (iTr .eq. iTr2.and.k .eq. j) cycle !dont use same atom as both neighbors
              jj = neigh%nb(j,i,iTr)
              kk = neigh%nb(k,i,iTr2)
              if (kk .eq. 0.or.jj .eq. 0) cycle !only j goes over nb so k or kk might be "out of bounds"
              atj = at(jj)
              atk = at(kk)
              fijk = param%angl(ati)*param%angl2(atj)*param%angl2(atk)
              if (fijk .lt. gen%fcthr) cycle     ! too small
              ! check for double counting
              tDbl = .false.
              do cd = 1,42
                if (cDbl(1,cd) .eq. 0) exit
if (cDbl(1,cd) .eq. kk.and.cDbl(2,cd) .eq. jj.and.cDbl(3,cd) .eq. iTr2.and.cDbl(4,cd) .eq. iTr) then
                  tDbl = .true.
                  exit
                end if
              end do
              if (tDbl) cycle
              cdi = cdi+1
              cDbl(1,cdi) = jj
              cDbl(2,cdi) = kk
              cDbl(3,cdi) = iTr
              cDbl(4,cdi) = iTr2
              call banglPBC(1,xyz,jj,i,kk,iTr,iTr2,neigh%transVec,phi)
              if (param%metal(ati) .gt. 0.and.phi*180./pi .lt. 60.) cycle ! skip eta cases even if CN < 6 (e.g. CaCp+)
              feta = 1.0d0
              if (imetal(ii) .eq. 2.and.itag(jj) .eq. -1.and.piadr(jj) .gt. 0) feta = 0.3d0       ! eta coord.
              if (imetal(ii) .eq. 2.and.itag(kk) .eq. -1.and.piadr(kk) .gt. 0) feta = feta*0.3d0  !
              nh = 0
              if (atj .eq. 1) nh = nh+1
              if (atk .eq. 1) nh = nh+1
              nnn = 0
              if (atj .eq. 7) nnn = nnn+1
              if (atk .eq. 7) nnn = nnn+1
              no = 0
              if (atj .eq. 8) no = no+1
              if (atk .eq. 8) no = no+1
              nheav = 0
              if (atj .gt. 14) nheav = nheav+1
              if (atk .gt. 14) nheav = nheav+1
              nsi = 0
              if (atj .eq. 14) nsi = nsi+1
              if (atk .eq. 14) nsi = nsi+1
              nc = 0
              if (atj .eq. 6) nc = nc+1
              if (atk .eq. 6) nc = nc+1
              nmet = 0
              if (param%metal(atj) .ne. 0) nmet = nmet+1
              if (param%metal(atk) .ne. 0) nmet = nmet+1
              npi = 0
              if (piadr(jj) .ne. 0) npi = npi+1
              if (piadr(kk) .ne. 0) npi = npi+1
              topo%nangl = topo%nangl+1
              topo%alist(1,topo%nangl) = ii
              topo%alist(2,topo%nangl) = jj  ! jj is shifted to iTr
              topo%alist(3,topo%nangl) = kk  ! kk is shifted to iTr2
              topo%alist(4,topo%nangl) = iTr
              topo%alist(5,topo%nangl) = iTr2
              call ringsbend(nat,ii,jj,kk,cring,sring,rings)
              triple = (topo%hyb(ii) .eq. 1.or.topo%hyb(jj) .eq. 1).or. &
    &                (topo%hyb(ii) .eq. 1.or.topo%hyb(kk) .eq. 1)
              if (imetal(ii) .eq. 0.and.imetal(jj) .eq. 0.and.imetal(kk) .eq. 0) then
                fqq = 1.0d0-(topo%qa(ii)*topo%qa(jj)+topo%qa(ii)*topo%qa(kk))*gen%qfacBEN      ! weaken it
              else
                fqq = 1.0d0-(topo%qa(ii)*topo%qa(jj)+topo%qa(ii)*topo%qa(kk))*gen%qfacBEN*2.5
              end if
              f2 = 1.0d0
              fn = 1.0d0

!-------------------------
! definitions come here
!-------------------------

!!!!!!!!!!
! DEFAULT
!!!!!!!!!!
              r0 = 100.0

              if (topo%hyb(i) .eq. 1) r0 = 180.
              if (topo%hyb(i) .eq. 2) r0 = 120.
              if (topo%hyb(i) .eq. 3) r0 = 109.5
              if (topo%hyb(i) .eq. 3.and.at(i) .gt. 10) then
                if (nn .le. 3) r0 = gen%aheavy3    ! heavy maingroup three coordinated
                if (nn .ge. 4) r0 = gen%aheavy4    ! heavy maingroup four  coordinated
                if (nn .eq. 4.and.param%group(ati) .eq. 5) r0 = 109.5      ! four coordinated group 5
                if (nn .eq. 4.and.param%group(ati) .eq. 4.and.ati .gt. 49) r0 = 109.5      ! four coordinated Sn, Pb
                if (param%group(ati) .eq. 4) r0 = r0-nh*5.   ! smaller angles for XHn Si...
                if (param%group(ati) .eq. 5) r0 = r0-nh*5.   ! smaller angles for XHn P..
                if (param%group(ati) .eq. 6) r0 = r0-nh*5.   ! smaller angles for XHn S..
              end if
              if (topo%hyb(i) .eq. 5) then
                r0 = 90.
                f2 = 0.11       ! not very important
                if (phi*180./pi .gt. gen%linthr) r0 = 180.       ! hypervalent coordination can be linear GEODEP
              end if
!!!!!!!!!!
! B
!!!!!!!!!!
              if (ati .eq. 5) then
                if (topo%hyb(i) .eq. 3) r0 = 115.
                if (topo%hyb(i) .eq. 2) r0 = 115.
              end if
!!!!!!!!!!
! C cases
!!!!!!!!!!
              if (ati .eq. 6) then
                if (topo%hyb(i) .eq. 3.and.nh .eq. 2) r0 = 108.6  ! CHH
                if (topo%hyb(i) .eq. 3.and.no .eq. 1) r0 = 108.5  ! COR
                if (topo%hyb(i) .eq. 2.and.no .eq. 2) r0 = 122.   ! COO
                if (topo%hyb(i) .eq. 2.and.no .eq. 1) f2 = 0.7    ! C=O
                if (topo%hyb(i) .eq. 1.and.no .eq. 2) then
                  triple = .false.   ! CO2
                  f2 = 2.0
                end if
                if (topo%hyb(i) .eq. 3.and.nn .gt. 4) then
                  if (phi*180./pi .gt. gen%linthr) r0 = 180.       ! hypervalent coordination can be linear GEODEP
                end if
              end if
!!!!!!!!!!
! O cases
!!!!!!!!!!
              if (ati .eq. 8.and.nn .eq. 2) then
                r0 = 104.5
!                   H2O
                if (nh .eq. 2) then
                  r0 = 100. ! compensate ES of the Hs
                  f2 = 1.20 ! H2O is better with 1.2-1.3 but H2O in fit behaves differently
                end if
                r0 = r0+7.*nsi   ! O angles widen with Si attached
                r0 = r0+14.*nmet  ! O angles widen with M attached
                if (npi .eq. 2) then
                  r0 = 109. ! e.g. Ph-O-Ph
                end if
                if (nmet .gt. 0.and.phi*180./pi .gt. gen%linthr) then
                  r0 = 180. ! metal coordination can be linear GEODEP
                  f2 = 0.3
                end if
              end if
!!!!!!!!!!
! N cases
!!!!!!!!!!
              if (ati .eq. 7.and.nn .eq. 2) then
                f2 = 1.4
                r0 = 115.
                if (rings .ne. 0) r0 = 105.
                if (at(kk) .eq. 8.or.at(jj) .eq. 8) r0 = 103.
                if (at(kk) .eq. 9.or.at(jj) .eq. 9) r0 = 102.
                if (topo%hyb(i) .eq. 1) r0 = 180.   ! NC or NNN
                if (imetal(jj) .eq. 2.and.topo%hyb(i) .eq. 1.and.at(kk) .eq. 7) r0 = 135.   ! NN on M
                if (imetal(kk) .eq. 2.and.topo%hyb(i) .eq. 1.and.at(jj) .eq. 7) r0 = 135.   ! NN on M
              end if
! NR3
              if (ati .eq. 7.and.topo%hyb(i) .eq. 3) then
!                 in pi system
                if (npi .gt. 0) then
                  if (amide(nat,at,topo%hyb,neigh%numnb,neigh%numctr,neigh%nb,piadr,i)) then
                    r0 = 115.
                    f2 = 1.2d0
                  else
                    sumppi = pbo(lin(ii,jj))+pbo(lin(ii,kk))
                    r0 = 113.
                    f2 = 1.d0-sumppi*0.7d0 ! must be -!
                  end if
                else
                  r0 = 104. ! sat. pyr. N, steep around 106
                  f2 = 0.40 ! 1.0 is better for NH3
                  f2 = f2+nh*0.19
                  f2 = f2+no*0.25
                  f2 = f2+nc*0.01
                end if
              end if
!!!!!!!!!!
! RING < 5
!!!!!!!!!!
              if (rings .eq. 3) r0 = 82. ! 60 gives too little strain
              if (rings .eq. 4) r0 = 96.
              if (rings .eq. 5.and.ati .eq. 6) r0 = 109.
!!!!!!!!!!
! specials
!!!!!!!!!!
! R-X in 3-rings e.g. cyclopropene
              if (rings .eq. 0) then
                call ringsatom(nat,i,cring,sring,idum)
                if (idum .eq. 3) then
                  call ringsatom(nat,jj,cring,sring,ringsj)
                  call ringsatom(nat,kk,cring,sring,ringsk)
                  if (ringsj+ringsk .eq. 102) r0 = r0+4.d0
                end if
              end if

! triple bonds
              if (triple) then
                f2 = 0.60d0  ! complex 7 in S30L makes artificial torsions if this is 0.4 which is
                ! slightly better for the phenylmethylethyne bending pot.
                if (atj .eq. 7.or.atk .eq. 7) f2 = 1.00d0
                if ((imetal(jj) .eq. 2.or.imetal(kk) .eq. 2).and.phi*180./pi .gt. gen%linthr) then
                  if (ati .eq. 6.and.atj .eq. 6) f2 = 3.   ! M-CC
                  if (ati .eq. 6.and.atk .eq. 6) f2 = 3.   ! M-CC
                  if (ati .eq. 6.and.atj .eq. 7) f2 = 3.   ! M-CN
                  if (ati .eq. 6.and.atk .eq. 7) f2 = 3.   ! M-CN
                  if (ati .eq. 6.and.param%group(atj) .eq. 6) f2 = 14.  ! M-CO or CS
                  if (ati .eq. 6.and.param%group(atk) .eq. 6) f2 = 14.  ! M-CO or CS
                  if (ati .eq. 7.and.atj .eq. 7) f2 = 10.  ! M-NN
                  if (ati .eq. 7.and.atj .eq. 6) f2 = 10.  ! M-NC
                  if (ati .eq. 7.and.atk .eq. 6) f2 = 10.  ! M-NC
                  if (ati .eq. 7.and.atj .eq. 8) then; r0 = 180.; f2 = 12.; end if  ! M-NO
                  if (ati .eq. 7.and.atk .eq. 8) then; r0 = 180.; f2 = 12.; end if  ! M-NO
                end if
              end if
! carbene analogous
              if (param%group(ati) .eq. 4.and.nn .eq. 2.and.itag(i) .eq. 1) then
                if (ati .eq. 6) r0 = 145.
                if (ati .gt. 6) r0 = 90.
              end if
! SO3X
              if (param%group(ati) .eq. 6.and.nn .eq. 4.and.no .ge. 1) r0 = 115.
! halogens CN=2
              if (param%group(ati) .eq. 7.and.topo%hyb(i) .eq. 1) then
                if (ati .eq. 9) r0 = 90.
                if (ati .eq. 17) r0 = 90.
                if (ati .eq. 35) r0 = 90.
                if (ati .eq. 53) r0 = 90.
                if (ati .gt. 9.and.phi*180./pi .gt. gen%linthr) r0 = 180. ! change to linear if linear coordinated, GEODEP
                f2 = 0.6/dble(ati)**0.15
              end if
! PB or Sn can be pyramidal
    if (topo%hyb(i) .eq. 3.and.param%group(ati) .eq. 4.and.ati .gt. 32.and.topo%qa(i) .gt. 0.4) then
                if (phi*180./pi .gt. 140.) then
                  r0 = 180. ! change to linear
                end if
                if (phi*180./pi .lt. 100.) then
                  r0 = 90.
                end if
                f2 = 1.0
              end if
! METAL
              if (imetal(ii) .gt. 0) then
                if (topo%hyb(i) .eq. 0) then
                  r0 = 90.
                  f2 = 1.35  ! important difference to other bends, big effect 1.15,1.25,1.35
                end if
                if (topo%hyb(i) .eq. 1) r0 = 180.
                if (topo%hyb(i) .eq. 2) r0 = 120.
                if (topo%hyb(i) .eq. 3) r0 = 109.5
                if (phi*180./pi .gt. gen%linthr) r0 = 180. ! change to linear
              end if

              fn = 1.0d0-2.36d0/dble(nn)**2

!----------------------
! end of definitions
!----------------------
              topo%vangl(1,topo%nangl) = r0*pi/180.
              fbsmall = (1.0d0-gen%fbs1*exp(-0.64*(topo%vangl(1,topo%nangl)-pi)**2))

!              central*neigbor charge spec. met.  small angle corr.
              topo%vangl(2,topo%nangl) = fijk*fqq*f2*fn*fbsmall*feta
              if (printlevel >= 3) write (myunit,'(3i5,2x,3f8.3,l2,i4)') ii,jj,kk,r0,phi*180./pi,topo%vangl(2,topo%nangl),picon,rings
            end do
          end do
        end do
      end do
    end do  ! end angl_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!              torsion
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    topo%ntors = sum(piadr)+nat
    do m = 1,neigh%nbond
      ii = neigh%blist(1,m)
      jj = neigh%blist(2,m)
      nni = sum(neigh%nb(neigh%numnb,ii,:))
      nnj = sum(neigh%nb(neigh%numnb,jj,:))
      if (btyp(m) .eq. 3.or.btyp(m) .eq. 6) cycle ! no sp-sp or metal eta
      if (param%tors(at(ii)) .lt. 0.or.param%tors(at(jj)) .lt. 0) cycle ! no negative values
      if (param%tors(at(ii))*param%tors(at(jj)) .lt. 1.d-3) cycle ! no small values
      if (param%metal(at(ii)) .gt. 1.and.nni .gt. 4) cycle ! no HC metals
      if (param%metal(at(jj)) .gt. 1.and.nnj .gt. 4) cycle !
      topo%ntors = topo%ntors+nni*nnj*2 ! upper limit
    end do
    maxtors = topo%ntors
    if (printlevel >= 3) write (myunit,*) 'torsion atoms        nrot   rings    phi0    phi      FC'

    topo%ntors_alloc = topo%ntors
    allocate (topo%tlist(8,topo%ntors),source=0)
    allocate (topo%vtors(2,topo%ntors),source=0.0d0)
    topo%ntors = 0
    do m = 1,neigh%nbond
      jj = neigh%blist(1,m)
      ii = neigh%blist(2,m)
      iTrj = neigh%blist(3,m)
      nni = sum(neigh%nb(neigh%numnb,ii,:))
      nnj = sum(neigh%nb(neigh%numnb,jj,:))
      if (btyp(m) .eq. 3.or.btyp(m) .eq. 6) cycle    ! metal eta or triple
      fij = param%tors(at(ii))*param%tors(at(jj))             ! atom contribution, central bond
      if (fij .lt. gen%fcthr) cycle
      if (param%tors(at(ii)) .lt. 0.or.param%tors(at(jj)) .lt. 0) cycle ! no negative values
      if (param%metal(at(ii)) .gt. 1.and.nni .gt. 4) cycle ! no HC metals
      if (param%metal(at(jj)) .gt. 1.and.nnj .gt. 4) cycle !
      fqq = 1.0d0+abs(topo%qa(ii)*topo%qa(jj))*gen%qfacTOR      ! weaken it for e.g. CF-CF and similar
      call ringsbond(nat,ii,jj,cring,sring,rings) ! i and j in same ring
      lring = .false.
      ccij = .false.
      if (rings .gt. 0) lring = .true.
      sp3ij = topo%hyb(ii) .eq. 3.and.topo%hyb(jj) .eq. 3
      if (at(ii) .eq. 6.and.at(jj) .eq. 6) ccij = .true.
      nhi = 1
      nhj = 1
      do iTr = 1,neigh%numctr
        do ineig = 1,neigh%nb(neigh%numnb,ii,iTr)
          if (at(neigh%nb(ineig,ii,iTr)) .eq. 1) nhi = nhi+1
        end do
      end do
      do iTr = 1,neigh%numctr
        do jneig = 1,neigh%nb(neigh%numnb,jj,iTr)
          if (at(neigh%nb(jneig,jj,iTr)) .eq. 1) nhj = nhj+1
        end do
      end do
      fij = fij*(dble(nhi)*dble(nhj))**0.07 ! n H term
      ! amides and alpha carbons in peptides/proteins
      if (alphaCO(nat,at,topo%hyb,neigh%numnb,neigh%numctr,neigh%nb,piadr,ii,jj)) fij = fij*1.3d0
      if (amide(nat,at,topo%hyb,neigh%numnb,neigh%numctr,neigh%nb,piadr,ii).and.topo%hyb(jj) .eq. 3.and.at(jj) .eq. 6) fij = fij*1.3d0
      if (amide(nat,at,topo%hyb,neigh%numnb,neigh%numctr,neigh%nb,piadr,jj).and.topo%hyb(ii) .eq. 3.and.at(ii) .eq. 6) fij = fij*1.3d0
      ! hypervalent
      if (btyp(m) .eq. 4) fij = fij*0.2d0
!     loop over neighbors of ij
      do iTrk = 1,neigh%numctr
        do ineig = 1,neigh%nb(neigh%numnb,ii,iTrk)
          kk = neigh%nb(ineig,ii,iTrk)  !neighbors of ii that are not jj
          if (kk .eq. jj.and.iTrk .eq. iTrj) cycle
          do iTrlDum = 1,neigh%numctr
            do jneig = 1,neigh%nb(neigh%numnb,jj,iTrlDum)
              ll = neigh%nb(jneig,jj,iTrlDum)   !neighbors of jj that are neither ii or kk
              iTrl = neigh%fTrSum(iTrlDum,iTrj) ! ll has to be shifted if jj is shifted
              if (iTrl .eq. -1.or.iTrl .gt. neigh%numctr) cycle
              if (ll .eq. ii.and.iTrl .eq. 1) cycle
              if (ll .eq. kk.and.iTrl .eq. iTrk) cycle
              if (chktors(nat,xyz,ii,jj,kk,ll,iTrj,iTrk,iTrl,neigh)) cycle  ! near 180
              fkl = param%tors2(at(kk))*param%tors2(at(ll))       ! outer kl term
           if (at(kk) .eq. 7.and.piadr(kk) .eq. 0) fkl = param%tors2(at(kk))*param%tors2(at(ll))*0.5
           if (at(ll) .eq. 7.and.piadr(ll) .eq. 0) fkl = param%tors2(at(kk))*param%tors2(at(ll))*0.5
              if (fkl .lt. gen%fcthr) cycle
              if (param%tors(at(kk)) .lt. 0.or.param%tors(at(ll)) .lt. 0) cycle ! no negative values
              f1 = gen%torsf(1)
              f2 = 0.0d0
fkl = fkl*(dble(sum(neigh%nb(neigh%numnb,kk,:)))*dble(sum(neigh%nb(neigh%numnb,ll,:))))**(-0.14)  ! CN term

!-----------------------
! definitions come here
!-----------------------
              if (lring) then
                if (rings .gt. 3) then
                  call ringstors(nat,ii,jj,kk,ll,cring,sring,rings4) ! smallest ring in which i,j,k,l are
                else
                  rings4 = 3 ! the 3-ring is special
                end if
! RING CASE
                nrot = 1
                if (btyp(m) .eq. 2) nrot = 2 ! max at 90 for pi and symmetric at 0,-180,180
                phi = 0  ! cis
                if (btyp(m) .eq. 1.and.rings4 .gt. 0) then
                  call ringstorl(nat,ii,jj,kk,ll,cring,sring,ringl)  ! largest ring in which i,j,k,l are
                  notpicon = piadr(kk) .eq. 0.and.piadr(ll) .eq. 0                      ! do it only for sat. rings
                  if (rings4 .eq. 3.and.notpicon) then; nrot = 1; phi = 0.d0; f1 = gen%fr3; end if
                  if (rings4 .eq. 4.and.ringl .eq. rings4.and.notpicon) then; nrot = 6; phi = 30.d0; f1 = gen%fr4; end if
                  if (rings4 .eq. 5.and.ringl .eq. rings4.and.notpicon) then; nrot = 6; phi = 30.d0; f1 = gen%fr5; end if
                  if (rings4 .eq. 6.and.ringl .eq. rings4.and.notpicon) then; nrot = 3; phi = 60.d0; f1 = gen%fr6; end if
                end if
                if (rings4 .eq. 0.and.btyp(m) .eq. 1.and.sum(neigh%nb(neigh%numnb,kk,:)) .eq. 1&
                        &.and.sum(neigh%nb(neigh%numnb,ll,:)) .eq. 1) then; nrot = 6; phi = 30.d0; f1 = 0.30; end if
                if (btyp(m) .eq. 2.and.rings .eq. 5.and.at(ii)*at(jj) .eq. 42) then
                  if (amide(nat,at,topo%hyb,neigh%numnb,neigh%numctr,neigh%nb,piadr,ii).or.&
        &amide(nat,at,topo%hyb,neigh%numnb,neigh%numctr,neigh%nb,piadr,jj)) f1 = 5.  ! improving CB7
                end if
              else
! ACYCLIC
                phi = 180.d0 ! trans
                nrot = 1
                if (topo%hyb(ii) .eq. 3.and.topo%hyb(jj) .eq. 3) nrot = 3 ! Me case
                if (btyp(m) .eq. 2) nrot = 2 ! max at 90 for pi and symmetric at 0,-180,180
                if (piadr(ii) .gt. 0.and.(piadr(jj) .eq. 0.and.topo%hyb(jj) .eq. 3)) then  ! pi-sp3
                  f1 = 0.5d0
                  if (at(ii) .eq. 7) f1 = 0.2d0 ! important for CB7 conf.
                  phi = 180.d0
                  nrot = 3
                end if
                if (piadr(jj) .gt. 0.and.(piadr(ii) .eq. 0.and.topo%hyb(ii) .eq. 3)) then
                  f1 = 0.5d0
                  if (at(jj) .eq. 7) f1 = 0.2d0 ! important for CB7 conf.
                  phi = 180.d0
                  nrot = 3
                end if
              end if
! SP3 specials
              if (topo%hyb(ii) .eq. 3.and.topo%hyb(jj) .eq. 3) then
! N-N, P-P ...
                if (param%group(at(ii)) .eq. 5.and.param%group(at(jj)) .eq. 5) then
                  nrot = 3
                  phi = 60.d0
                  f1 = 3.0d0
                end if
! 5-6
                if ((param%group(at(ii)) .eq. 5.and.param%group(at(jj)) .eq. 6).or. &
      &            (param%group(at(ii)) .eq. 6.and.param%group(at(jj)) .eq. 5)) then
                  nrot = 2
                  phi = 90.d0
                  f1 = 1.0d0
                  if (at(ii) .ge. 15.and.at(jj) .ge. 15) f1 = 20.0d0
                end if
! O-O, S-S ...
                if (param%group(at(ii)) .eq. 6.and.param%group(at(jj)) .eq. 6) then
                  nrot = 2
                  phi = 90.d0
                  f1 = 5.0d0
                  if (at(ii) .ge. 16.and.at(jj) .ge. 16) f1 = 25.0d0 ! better for h2s2
                end if
              end if
! pi system
              if (pibo(m) .gt. 0) then
                f2 = pibo(m)*exp(-2.5d0*(1.24d0-pibo(m))**14)  ! decrease to very small values for P < 0.3
                ! values of 2.5 instead of 2.4 give larger tangles
                ! the parameter 1.24 is very sensitive ie 1.25 yield 5 deg more in 1,3cB
                if (piadr(kk) .eq. 0.and.at(kk) .gt. 10) f2 = f2*1.3! the pi BO becomes more significant if heavies are attached
                if (piadr(ll) .eq. 0.and.at(ll) .gt. 10) f2 = f2*1.3
                f1 = f1*0.55
              end if

              if (topo%hyb(kk) .eq. 5.or.topo%hyb(ll) .eq. 5) fkl = fkl*1.5 ! hypervalent corr.
!--------------------
! end of definitions
!-------------------

! total FC            sigma       pi             charge central outer kl
              fctot = (f1+10.d0*gen%torsf(2)*f2)*fqq*fij*fkl

              if (fctot .gt. gen%fcthr) then ! avoid tiny potentials
                topo%ntors = topo%ntors+1
                if (topo%ntors .gt. maxtors) then
 if (printlevel >= 1) write (myunit,'("**ERROR** ",a,1x,a)') 'internal (torsion setup) error',source
                  return
                end if
                topo%tlist(1,topo%ntors) = ll
                topo%tlist(2,topo%ntors) = ii
                topo%tlist(3,topo%ntors) = jj
                topo%tlist(4,topo%ntors) = kk
                topo%tlist(5,topo%ntors) = nrot
                topo%tlist(6,topo%ntors) = iTrl
                topo%tlist(7,topo%ntors) = iTrj
                topo%tlist(8,topo%ntors) = iTrk
                topo%vtors(1,topo%ntors) = phi*pi/180.0d0
                topo%vtors(2,topo%ntors) = fctot
!                 printout
                phi = valijklff(nat,xyz,ll,ii,jj,kk)
                if (printlevel >= 3) write (myunit,'(4i5,2x,i2,5x,i2,4x,3f8.3)') &
   &            ii,jj,kk,ll,topo%tlist(5,topo%ntors),rings,topo%vtors(1,topo%ntors)*180./pi,phi*180./pi,topo%vtors(2,topo%ntors)
              end if

! extra rot=1 torsion potential for sp3-sp3 to get gauche conf energies well
              sp3kl = topo%hyb(kk) .eq. 3.and.topo%hyb(ll) .eq. 3
              if (sp3kl.and.sp3ij.and.(.not.lring).and.btyp(m) .lt. 5) then
                topo%ntors = topo%ntors+1
                if (topo%ntors .gt. maxtors) then
 if (printlevel >= 1) write (myunit,'("**ERROR** ",a,1x,a)') 'internal (torsion setup) error',source
                  return
                end if
                ff = gen%torsf(6)
                if (at(ii) .eq. 7.or.at(jj) .eq. 7) ff = gen%torsf(7)
                if (at(ii) .eq. 8.or.at(jj) .eq. 8) ff = gen%torsf(8)
                topo%tlist(1,topo%ntors) = ll
                topo%tlist(2,topo%ntors) = ii
                topo%tlist(3,topo%ntors) = jj
                topo%tlist(4,topo%ntors) = kk
                topo%tlist(6,topo%ntors) = iTrl
                topo%tlist(7,topo%ntors) = iTrj
                topo%tlist(8,topo%ntors) = iTrk
                topo%tlist(5,topo%ntors) = 1
                topo%vtors(1,topo%ntors) = pi
                topo%vtors(2,topo%ntors) = ff*fij*fkl*fqq
                if (printlevel >= 3) write (myunit,'(4i5,2x,i2,5x,i2,4x,3f8.3)') &
   &            ii,jj,kk,ll,topo%tlist(5,topo%ntors),rings,topo%vtors(1,topo%ntors)*180./pi,phi*180./pi,topo%vtors(2,topo%ntors)
              end if

            end do ! neighbors ij
          end do
        end do
      end do ! bond loop
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! out-of-plane, improper (three-fold coordinated central pi atom i or an N)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (printlevel >= 3) write (myunit,*) 'out-of-plane atoms          phi0    phi      FC'
    do i = 1,nat
      if (sum(neigh%nb(neigh%numnb,i,:)) .ne. 3) cycle
      if (piadr(i) .eq. 0) then
        if (at(i) .ne. 7) cycle
      end if
      topo%ntors = topo%ntors+1
      ! get those 3 neighbors
      call neigh%nbLoc(nat,neigh%nb,i,locarr)
      jj = 0
      kk = 0
      ll = 0
      iTrl = 0
      if (size(locarr,dim=2) .eq. 1) then     ! all 3 in same cell
        jj = locarr(1,1)
        kk = locarr(2,1)
        ll = locarr(3,1)
        iTrj = locarr(neigh%numnb,1)
        iTrk = locarr(neigh%numnb,1)
        iTrl = locarr(neigh%numnb,1)
      elseif (size(locarr,dim=2) .eq. 2) then ! in 2 cells
        jj = locarr(1,1)
        kk = locarr(1,2)
        ll = locarr(2,1)             !  guess its here
        iTrj = locarr(neigh%numnb,1)
        iTrk = locarr(neigh%numnb,2)
        iTrl = locarr(neigh%numnb,1)
        if (ll .eq. 0) then
          ll = locarr(2,2) ! if ll was zero its in this cell
          iTrl = locarr(neigh%numnb,2)
        end if
      else                                  ! in 3 cells
        jj = locarr(1,1)
        kk = locarr(1,2)
        ll = locarr(1,3)
        iTrj = locarr(neigh%numnb,1)
        iTrk = locarr(neigh%numnb,2)
        iTrl = locarr(neigh%numnb,3)
      end if
      deallocate (locarr)
!        sort atoms according to distance to central atom such that the same inversion angle def. always results
      sdum3(1) = NORM2(xyz(:,i)-(xyz(:,jj)+neigh%transVec(:,iTrj)))
      sdum3(2) = NORM2(xyz(:,i)-(xyz(:,kk)+neigh%transVec(:,iTrk)))
      sdum3(3) = NORM2(xyz(:,i)-(xyz(:,ll)+neigh%transVec(:,iTrl)))
      ind3(1) = jj
      ind3(2) = kk
      ind3(3) = ll
      call ssort(3,sdum3,ind3)
      sdum3(1) = NORM2(xyz(:,i)-(xyz(:,jj)+neigh%transVec(:,iTrj)))
      sdum3(2) = NORM2(xyz(:,i)-(xyz(:,kk)+neigh%transVec(:,iTrk)))
      sdum3(3) = NORM2(xyz(:,i)-(xyz(:,ll)+neigh%transVec(:,iTrl)))
      jj = ind3(1)  ! assign sorted indices
      kk = ind3(2)
      ll = ind3(3)
      ! now sort iTr's
      ind3(1) = iTrj
      ind3(2) = iTrk
      ind3(3) = iTrl
      call ssort(3,sdum3,ind3)
      iTrj = ind3(1)
      iTrk = ind3(2)
      iTrl = ind3(3)
      ! save to tlist
      topo%tlist(1,topo%ntors) = i
      topo%tlist(2,topo%ntors) = jj
      topo%tlist(3,topo%ntors) = kk
      topo%tlist(4,topo%ntors) = ll
      topo%tlist(6,topo%ntors) = iTrl
      topo%tlist(7,topo%ntors) = iTrj
      topo%tlist(8,topo%ntors) = iTrk
      if (piadr(i) .eq. 0.and.at(i) .eq. 7) then  ! sat N case
        r0 = 80.0d0
        ff = 0.60d0
        topo%tlist(5,topo%ntors) = -1
        topo%vtors(1,topo%ntors) = r0*pi/180. ! double min at +/- phi0
        topo%vtors(2,topo%ntors) = 0.0d0
        do iTr = 1,neigh%numctr
          do m = 1,neigh%nb(neigh%numnb,i,iTr)
            idum = neigh%nb(m,i,iTr)
            topo%vtors(2,topo%ntors) = topo%vtors(2,topo%ntors)+ff*sqrt(param%repz(at(idum)))  ! NX3 has higher inv barr. than NH3
          end do
        end do
      else
        ncarbo = 0
        nf = 0
        do iTr = 1,neigh%numctr
          do m = 1,neigh%nb(neigh%numnb,i,iTr)
            idum = neigh%nb(m,i,iTr)
            if (at(idum) .eq. 8.or.at(idum) .eq. 16) ncarbo = ncarbo+1
            if (param%group(at(idum)) .eq. 7) nf = nf+1
          end do
        end do
        fqq = 1.0d0+topo%qa(i)*5.0d0
        topo%tlist(5,topo%ntors) = 0         ! phi0=0 case (pi)
        topo%vtors(1,topo%ntors) = 0.0d0     !  "      "
        sumppi = pbo(lin(i,jj))+pbo(lin(i,kk))+pbo(lin(i,ll))
        f2 = 1.0d0-sumppi*gen%torsf(5)
!                         base val  piBO  charge term
        topo%vtors(2,topo%ntors) = gen%torsf(3)*f2*fqq
!          carbonyl corr.
        if (at(i) .eq. 5.and.ncarbo .gt. 0) topo%vtors(2,topo%ntors) = topo%vtors(2,topo%ntors)*38.
        if (at(i) .eq. 6.and.ncarbo .gt. 0) topo%vtors(2,topo%ntors) = topo%vtors(2,topo%ntors)*38.
        if (at(i) .eq. 6.and.nf .gt. 0.and.ncarbo .eq. 0) topo%vtors(2,topo%ntors) = topo%vtors(2,topo%ntors)*10.
if (at(i) .eq. 7.and.ncarbo .gt. 0) topo%vtors(2,topo%ntors) = topo%vtors(2,topo%ntors)*10./f2 ! no pi dep
      end if
!        printout
      vTrl = neigh%transVec(:,iTrl)
      vTrj = neigh%transVec(:,iTrj)
      vTrk = neigh%transVec(:,iTrk)
      phi = omegaPBC(nat,xyz,i,jj,kk,ll,vTrl,vTrj,vTrk)
      if (printlevel >= 3) write (myunit,'(4i5,7x,3f8.3)') i,jj,kk,ll,topo%vtors(1,topo%ntors)*180./pi,phi*180./pi,topo%vtors(2,topo%ntors)
    end do

    if (printlevel >= 2) then
      write (myunit,'(10x,"#tors  :",3x,i0)') topo%ntors
      write (myunit,'(10x,"#nmol  :",3x,i0)') topo%nfrag
    end if

! all done

    topo%maxsystem = 5000

    if (.false.) then
      !call fragmentize(nat,at,xyz,topo%maxsystem,500,rab,neigh%numnb,neigh%numctr,neigh%nb, &
      !   & topo%ispinsyst,topo%nspinsyst,topo%nsystem)
    else
      topo%nsystem = 1
    end if

    if (printlevel >= 2) write (myunit,'(10x,"#optfrag :",3x,i0)') topo%nfrag

    ! check if triple bonded carbon is present (for torsion term)
    nn = 0
    do i = 1,nat
      if (at(i) .eq. 6.and.sum(neigh%nb(neigh%numnb,i,:)) .eq. 2) then
        do j = 1,2
          nbi = 0
          call neigh%jth_nb(nat,xyz,nbi,j,i,iTr)  ! nbi is the jth nb of i in cell iTr
          if (nbi .eq. 0) cycle
          if (at(nbi) .eq. 6.and.sum(neigh%nb(neigh%numnb,nbi,:)) .eq. 2) then
            nn = nn+1
          end if
        end do
      end if
    end do
    nn = nn/2
    allocate (topo%sTorsl(6,nn),source=0)
    topo%nstors = nn
    if (nn .ne. 0) then
      ! fix double counting
      call specialTorsList(nn,nat,at,xyz,topo,neigh,topo%sTorsl)
    end if

  end subroutine set_bonded_parameters

! using C1=ii, C2=jj, C3=kk, C4=ll
  subroutine specialTorsList(nst,nat,at,xyz,topo,neigh,sTorsList)
    !***********************************************************************
    !* Collect the C1-C2-C(sp)#C(sp)-C3-C4 units that get the special
    !* torsion potential, i.e. a carbon triple bond whose two carbons each
    !* carry an sp2 carbon that is itself bonded to a further sp2 carbon.
    !* Input:
    !*   nst      - number of triple bonded carbon pairs found by the caller
    !*   nat/at/xyz - system definition
    !*   topo     - topology, provides the hybridisation
    !*   neigh    - neighbour data
    !* Output:
    !*   sTorsList - (C1,C2,Ci,Cnbi,C3,C4) per unit; entries that could not
    !*               be completed are left at zero and skipped by sTors_eg
    !***********************************************************************
    integer,intent(in) :: nst
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh ! main type for introducing PBC
    integer,intent(inout) :: sTorsList(6,nst)
    integer :: i,j,k,ii,jj,kk,ll,idx,iTr,nbi,nbk
    logical :: iiok,llok
    ! initialize variables
    idx = 0
    ii = -1
    jj = -1
    kk = -1
    ll = -1

    do i = 1,nat
      ! carbon with two neighbors bonded to other carbon* with two neighbors
      if (at(i) .eq. 6.and.sum(neigh%nb(neigh%numnb,i,:)) .eq. 2) then
        do j = 1,2
          call neigh%jth_nb(nat,xyz,nbi,j,i,iTr)  ! nbi is the jth nb of i in cell iTr
          if (nbi .eq. 0.or.iTr .eq. 0) cycle
          if (at(nbi) .eq. 6.and.sum(neigh%nb(neigh%numnb,nbi,:)) .eq. 2) then  ! *other carbon
            ! check carbon triple bond distance
            if (NORM2(xyz(1:3,i)-xyz(1:3,nbi)) .le. 2.37) then
              ! at this point we know that i and nbi are carbons bonded through triple bond
              ! check C2 and C3
              ! reset per candidate pair, otherwise a match found for an
              ! earlier atom i would still be in jj/kk here
              jj = -1
              kk = -1
              do k = 1,2  ! C2 is other nb of Ci
                nbk = 0
                call neigh%jth_nb(nat,xyz,nbk,k,i,iTr)  ! nbk is the kth nb of i in cell iTr .ne.nbi
                if (nbk .ne. nbi.and.nbk .ne. 0) then
                  if (at(nbk) .eq. 6) jj = nbk ! C2 index
                end if
              end do
              do k = 1,2  ! C3 is other nb of Cnbi
                nbk = 0
                call neigh%jth_nb(nat,xyz,nbk,k,nbi,iTr)  ! nbk is the kth nb of nbi in cell iTr .ne.i
                if (nbk .ne. i.and.nbk .ne. 0) then
                  if (at(nbk) .eq. 6) kk = nbk ! C3 index
                end if
              end do
              if (jj .eq. -1.or.kk .eq. -1) then
                exit ! next atom i
              end if
              ! check C1 through C4 are sp2 carbon
              if (topo%hyb(jj) .eq. 2.and.topo%hyb(kk) .eq. 2 &
              &   .and.at(jj) .eq. 6.and.at(kk) .eq. 6) then
                iiok = .false.
                llok = .false.
                ! which of the two valid neighbors is picked as C1 depends
                !  on atom sorting in input file !!! The last one in file.
                do k = 1,sum(neigh%nb(neigh%numnb,jj,:))
                  nbk = 0
                  call neigh%jth_nb(nat,xyz,nbk,k,jj,iTr)  ! nbk is the kth nb of C2
                  if (nbk .eq. 0.or.nbk .eq. i) cycle
                  if (topo%hyb(nbk) .eq. 2.and.at(nbk) .eq. 6 &
                     & .and.sum(neigh%nb(neigh%numnb,nbk,:)) .eq. 3) then
                    ii = nbk
                    iiok = .true.
                  end if
                end do
                ! which of the two valid neighbors is picked as C4 depends
                !  on atom sorting in input file !!! The last one in file.
                do k = 1,sum(neigh%nb(neigh%numnb,kk,:))
                  nbk = 0
                  call neigh%jth_nb(nat,xyz,nbk,k,kk,iTr)  ! nbk is the kth nb of C3
                  if (nbk .eq. 0.or.nbk .eq. nbi) cycle
                  if (topo%hyb(nbk) .eq. 2.and.at(nbk) .eq. 6 &
                     & .and.sum(neigh%nb(neigh%numnb,nbk,:)) .eq. 3) then
                    ll = nbk
                    llok = .true.
                  end if
                end do
                if (nbi .gt. i.and.iiok.and.llok) then ! to avoid double counting
                  idx = idx+1
                  sTorsList(1,idx) = ii  ! C1
                  sTorsList(2,idx) = jj  ! C2
                  sTorsList(3,idx) = i   ! Ci
                  sTorsList(4,idx) = nbi ! Cnbi
                  sTorsList(5,idx) = kk  ! C3
                  sTorsList(6,idx) = ll  ! C4
                end if
              end if ! C1-C4 are sp2 carbon
            end if  ! CC distance
          end if  ! other carbon
        end do
      end if ! is carbon with nnb=2
    end do
  end subroutine specialTorsList

!========================================================================================!

  subroutine hueckel_solve(nat,at,xyz,param,gen,rab,itag,npiall,picount,pimvec, &
        & ipis,piadr3,piadr4,itmp,piadr,pibo,pbo,printlevel,printunit,neigh,topo,io)
    !***********************************************************************
    !* Iterative Hueckel treatment of the pi subsystems.
    !*
    !* Each pi system found by the setup is solved on its own: build the
    !* Hueckel matrix over its atoms, diagonalise, occupy, and feed the
    !* resulting bond orders back into the off-diagonal elements. The
    !* iteration damps the off-diagonals by the density so that a system
    !* like cyclooctatetraene localises onto the right bonds instead of
    !* delocalising over the whole ring.
    !*
    !* The pi bond orders that come out drive the bond, torsion and
    !* out-of-plane parameters, and piadr is overwritten with the final
    !* pi-atom marker the rest of the setup reads.
    !* Input:
    !*   nat/at/xyz - system definition
    !*   param/gen  - GFN-FF parameters and generator thresholds
    !*   rab        - packed interatomic distances
    !*   itag       - carbene tag; a carbene carbon contributes no pi electron
    !*   npiall     - number of candidate pi atoms
    !*   picount    - number of separate pi subsystems
    !*   pimvec     - which subsystem each candidate belongs to; released here,
    !*                together with ipis, once the last subsystem is done
    !*   ipis        - net charge carried by each pi subsystem
    !*   printlevel/printunit - verbosity and output unit
    !* In/out:
    !*   piadr3/piadr4/itmp - setup scratch, allocated by the caller
    !*   piadr      - candidate pi atoms on entry, pi atom marker on exit
    !*   pibo       - pi bond order per bond
    !*   pbo        - pi bond order per atom pair
    !*   neigh      - neighbour data, provides the bond list
    !*   topo       - topology; hyb is read, nothing else is written here
    !* Output:
    !*   io         - non-zero if a Hueckel diagonalisation failed
    !***********************************************************************
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TGFFData),intent(in) :: param
    type(TGFFGenerator),intent(in) :: gen
    real(wp),intent(in) :: rab(:)
    integer,intent(in) :: itag(nat)
    integer,intent(in) :: npiall,picount
    integer,allocatable,intent(inout) :: pimvec(:)
    integer,allocatable,intent(inout) :: ipis(:)
    integer,intent(inout) :: piadr3(:),piadr4(:),itmp(:)
    integer,intent(inout) :: piadr(:)
    real(wp),intent(inout) :: pibo(:),pbo(:)
    integer,intent(in) :: printlevel,printunit
    type(TNeigh),intent(inout) :: neigh
    type(TGFFTopology),intent(inout) :: topo
    integer,intent(out) :: io

    character(len=*),parameter :: source = 'hueckel_solve'
    integer :: ati,i,k,nn,ii,jj,ia,ja
    integer :: hybi,pis,hcalc
    integer :: npi,nelpi
    integer :: iTr
    integer :: myunit
    real(wp) :: dum,dum2
    real(wp) :: eold
    integer,allocatable :: piel(:)
    real(wp),allocatable :: Api(:,:),S(:,:),Pold(:,:),occ(:),eps(:),apisave(:,:)
    real(wp),allocatable :: pispop(:),pisea(:),pisip(:)

    myunit = printunit
    io = 0

    if (picount .gt. 0) then
      if (printlevel >= 2) write (myunit,'(10x,"doing iterative Hueckel for ",i0," subsystem(s) ...")') picount
      allocate (pispop(picount),pisip(picount),pisea(picount),source=0.0d0)
      allocate (piel(nat),source=0)
      itmp = 0 ! save pi atom info
      hcalc = 0
      pisip = 0
      pisea = 0

      if (printlevel >= 2) then
        write (myunit,'(10x,"iterative Hueckel run to get P ...")')
      end if
      do pis = 1,picount ! loop over pi systems
        npi = 0
        nelpi = 0
        piadr3 = 0
        piadr4 = 0
        piel = 0

        do k = 1,npiall
          if (pimvec(k) .eq. pis) then
            npi = npi+1
            ati = at(piadr(k))
            hybi = topo%hyb(piadr(k))
            ii = nelpi
            if (ati .eq. 5.and.hybi .eq. 1) nelpi = nelpi+1  ! B in borine
            if (ati .eq. 6.and.itag(piadr(k)) .ne. 1) nelpi = nelpi+1  ! skip if its a carbene (tag itag=1)
            if (ati .eq. 7.and.hybi .eq. 2.and.itag(piadr(k)) .eq. 1) &
     &                                           nelpi = nelpi+1  ! the itag=1 avoids an odd el number for the nitro group (its 4)
            if (ati .eq. 7.and.hybi .le. 2) nelpi = nelpi+1
            if (ati .eq. 7.and.hybi .eq. 3) nelpi = nelpi+2
            if (ati .eq. 8.and.hybi .eq. 1) nelpi = nelpi+1
            if (ati .eq. 8.and.hybi .eq. 2) nelpi = nelpi+1
            if (ati .eq. 8.and.hybi .eq. 3) nelpi = nelpi+2
            if (ati .eq. 9.and.hybi .ne. 1) nelpi = nelpi+2
            if (ati .eq. 9.and.hybi .eq. 1) nelpi = nelpi+3 !??? otherwise fluor-furan+ is wrong
            if (ati .eq. 16.and.hybi .eq. 1) nelpi = nelpi+1
            if (ati .eq. 16.and.hybi .eq. 2) nelpi = nelpi+1
            if (ati .eq. 16.and.hybi .eq. 3) nelpi = nelpi+2
            if (ati .eq. 17.and.hybi .eq. 0) nelpi = nelpi+2
            if (ati .eq. 17.and.hybi .eq. 1) nelpi = nelpi+3
            piadr3(npi) = piadr(k) ! map to original, full atom set
            piadr4(piadr(k)) = npi
            piel(piadr(k)) = nelpi-ii
            if (piel(piadr(k)) .gt. 2) piel(piadr(k)) = 2
          end if
        end do
        nelpi = nelpi-ipis(pis)
        if (npi .lt. 2.or.nelpi .lt. 1) cycle
        allocate (Api(npi,npi),apisave(npi,npi),Pold(npi,npi),S(npi,npi),occ(npi),eps(npi)) ! S is just scratch here

        eold = 0
        Pold = 2.d0/3.d0
! iterative Hueckel loop, off-diag terms are reduced depending on P to avoid overdelocalization
        do nn = 1,nint(gen%maxhiter)      ! just some iterations
          Api = 0
          do i = 1,npi
            ii = piadr3(i)
            Api(i,i) = gen%hdiag(at(ii))+topo%qa(ii)*gen%hueckelp3-dble(piel(ii)-1)*gen%pilpf
          end do
!     loop over bonds for pair interactions
          do i = 1,neigh%nbond
            jj = neigh%blist(1,i)
            ii = neigh%blist(2,i)
            iTr = neigh%blist(3,i)
            ia = piadr4(ii)
            ja = piadr4(jj)
            if (ia .gt. 0.and.ja .gt. 0) then
              !dum=1.d-9*rab(lin(ii,jj))                                 ! distort so that Huckel for e.g. COT localizes to right bonds
              dum = 1.d-9*NORM2(xyz(:,ii)-(xyz(:,jj)+neigh%transVec(:,iTr)))  ! distort so that Huckel for e.g. COT localizes to right bonds
              dum = sqrt(gen%hoffdiag(at(ii))*gen%hoffdiag(at(jj)))-dum           ! better than arithmetic
              dum2 = gen%hiter
              if (topo%hyb(ii) .eq. 1) dum2 = dum2*gen%htriple        ! triple bond is different
              if (topo%hyb(jj) .eq. 1) dum2 = dum2*gen%htriple        ! triple bond is different
              Api(ja,ia) = -dum*(1.0d0-dum2*(2.0d0/3.0d0-Pold(ja,ia))) ! Pmat scaling with benzene as reference
              Api(ia,ja) = Api(ja,ia)
            end if
          end do

          apisave = Api
          call gfnffqmsolve(printlevel,Api,S,.false.,4000.0d0,npi,0,nelpi,dum,occ,eps,io,myunit)  !diagonalize, 4000 better than 300

          do i = 1,npi  ! save IP/EA
            if (occ(i) .gt. 0.5) then
              pisip(pis) = eps(i)   ! IP
              if (i+1 .lt. npi) pisea(pis) = eps(i+1) ! EA
            end if
          end do
          if (abs(dum-eold) .lt. 1.d-4) exit  ! end of iterations
          Pold = Api
          eold = dum
        end do
! end of iterative loop
        if (printlevel >= 2) then
          write (myunit,'(''Hueckel system :'',i3,'' charge : '',i3,'' ndim/Nel :'',2i5, &
      &         3x, ''eps(HOMO/LUMO)'',2f12.6)') pis,ipis(pis),npi,nelpi,pisip(pis),pisea(pis)
        end if
        if (pisip(pis) .gt. 0.40) then
          if (printlevel >= 1) then
            write (myunit,'(a,i0,a)') 'WARNING: probably wrong pi occupation for system ',pis,'. Second attempt with Nel=Nel-1!'
            do i = 1,nat
              if (piadr4(i) .ne. 0) write (myunit,*) 'at,nb,topo%hyb,Npiel:',i,pse(at(i)),sum(neigh%nb(neigh%numnb,i,:)),topo%hyb(i),piel(i)
            end do
          end if
          nelpi = nelpi-1
          Api = Apisave
          call gfnffqmsolve(printlevel,Api,S,.false.,4000.0d0,npi,0,nelpi,dum,occ,eps,io,myunit)  !diagonalize
          !call PREIG(6,occ,1.0d0,eps,1,npi)
          do i = 1,npi  ! save IP/EA
            if (occ(i) .gt. 0.5) then
              pisip(pis) = eps(i)   ! IP
              if (i+1 .lt. npi) pisea(pis) = eps(i+1) ! EA
            end if
          end do
          if (printlevel >= 2) then
            write (myunit,'(''Hueckel system :'',i3,'' charge : '',i3,'' ndim/Nel :'',2i5, &
        &         3x, ''eps(HOMO/LUMO)'',2f12.6)') pis,ipis(pis),npi,nelpi,pisip(pis),pisea(pis)
          end if
        end if
! save BO
        do i = 1,neigh%nbond
          jj = neigh%blist(1,i)
          ii = neigh%blist(2,i)
          ja = piadr4(jj)
          ia = piadr4(ii)
          if (ia .gt. 0.and.ja .gt. 0) then
            pibo(i) = Api(ja,ia)
            pbo(lin(ii,jj)) = Api(ja,ia)
            itmp(ii) = 1
            itmp(jj) = 1
          end if
        end do
        deallocate (Api,apisave,Pold,S,occ,eps)
      end do
! end of pi system loop
      piadr = itmp  ! array used for identifying pi atoms in following codes
      deallocate (pispop,pisip,pisea,ipis,pimvec,piel)
    end if
!----------- end Hueckel

  end subroutine hueckel_solve

!========================================================================================!
end module gfnff_topo_iniphases
