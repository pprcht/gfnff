! This file is part of xtb.
!
! Copyright (C) 2019-2020 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

module gfnff_engrad_module

  use iso_fortran_env,only:wp => real64,sp => real32,stdout => output_unit
  use gfnff_ini2
  use gfnff_data_types,only:TGFFData,TGFFNeighbourList,new,TGFFTopology, &
    &                       TCell,TDispersionData
  use gfnff_neighbor,only:TNeigh
  use gfnff_gbsa,only:TBorn
  use gfnff_param,only:sqrtZr4r2,gffVersion,gfnff_thresholds
  use gfnff_helpers
  use gfnff_cn
  use gfnff_rab
  use gfnff_gdisp0
  use gfnff_timer_mod,only:gfnff_timer
  use gfnff_math_wrapper
  implicit none
  private
  public :: gfnff_eg,gfnff_results

!&<
  !> results log
  type :: gfnff_results
    real(wp) :: e_total   = 0.0_wp
    real(wp) :: e_rep     = 0.0_wp
    real(wp) :: e_es      = 0.0_wp
    real(wp) :: e_disp    = 0.0_wp
    real(wp) :: e_xb      = 0.0_wp
    real(wp) :: g_born    = 0.0_wp
    real(wp) :: g_sasa    = 0.0_wp
    real(wp) :: g_hb      = 0.0_wp
    real(wp) :: g_shift   = 0.0_wp
    real(wp) :: dipole(3) = (/0.0_wp,0.0_wp,0.0_wp/)
    real(wp) :: g_solv    = 0.0_wp
    real(wp) :: gnorm     = 0.0_wp
    real(wp) :: e_bond    = 0.0_wp
    real(wp) :: e_angl    = 0.0_wp
    real(wp) :: e_tors    = 0.0_wp
    real(wp) :: e_hb      = 0.0_wp
    real(wp) :: e_batm    = 0.0_wp
    real(wp) :: e_ext     = 0.0_wp
  end type gfnff_results

  real(wp),private,parameter :: pi = 3.1415926535897932385_wp
  real(wp),private,parameter :: sqrtpi = 1.77245385091_wp
!&>

! ══════════════════════════════════════════════════════════════════════════════
contains  !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════

!---------------------------------------------------
! GFN-FF
! energy and analytical gradient for given xyz and
! charge ichrg
! requires D3 ini (rcov,r2r4,copyc6) as well as
! gfnff_ini call
!
! the total energy is
! ees + edisp + erep + ebond + eangl + etors + ehb + exb + ebatm + eext
!
! uses EEQ charge and D3 routines
! basic trigonometry for bending and torsion angles
! taken slightly modified from QMDFF code
! repulsion and rabguess from xtb GFN0 part
!
! requires setup of
!     integer,allocatable :: blist(:,:)
!     integer,allocatable :: alist(:,:)
!     integer,allocatable :: tlist(:,:)
!     integer,allocatable ::b3list(:,:)
!     real(wp),allocatable:: vbond(:,:)
!     real(wp),allocatable:: vangl(:,:)
!     real(wp),allocatable:: vtors(:,:)
!     chi,gam,alp,cnf
!     repa,repz,alphanb
!
!---------------------------------------------------

  subroutine gfnff_eg(printlevel,n,at,xyz,cell,sigma,ichrg,g,etot,res_gff, &
        & param,topo,neigh,nlist,efield,solvation,update,version,accuracy,printunit)
    !***********************************************************************
    !* Compute GFN-FF energy and analytical gradient.
    !* Input:
    !*   printlevel  - verbosity (0=silent, 1=minimal timing, 2=standard, 3=verbose)
    !*   n/at/xyz    - system definition
    !*   cell        - periodic cell (npbc=0 for molecular)
    !*   ichrg       - total charge
    !*   param/topo/neigh/nlist - GFN-FF data structures
    !*   efield      - external electric field
    !*   solvation   - GBSA/ALPB model (optional, allocated if active)
    !*   update      - rebuild HB/XB lists if .true.
    !*   version/accuracy - parameter version and threshold accuracy
    !* Output:
    !*   g           - gradient (Eh/Bohr)
    !*   etot        - total energy (Eh)
    !*   res_gff     - energy decomposition
    !*   sigma       - stress tensor (Eh, non-zero only for PBC)
    !*   printunit   - output unit (optional, default: stdout)
    !***********************************************************************
    implicit none

    character(len=*),parameter :: source = 'gfnff_eg'
    real(wp),allocatable :: rec_tVec(:,:)
    type(TNeigh),intent(inout) :: neigh ! main type for introducing PBC
    integer,intent(in)  :: n,ichrg,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TCell),intent(in) :: cell
    type(TDispersionData) :: disp_par,mcdisp_par
    type(gfnff_results),intent(out) :: res_gff
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TGFFNeighbourList),intent(inout) :: nlist
    type(TBorn),allocatable,intent(inout) :: solvation
    real(wp),intent(in) :: efield(3)
    logical,intent(in) :: update
    integer,intent(in) :: version
    real(wp),intent(in) :: accuracy
    integer,intent(in) :: printlevel  !< verbosity (0=silent,1=timing,2=info,3=verbose)
    integer,intent(in),optional :: printunit  !< output unit (default: stdout)

    real(wp),intent(out) :: sigma(3,3) ! stress tensor
    real(wp),intent(out) :: g(3,n)
    real(wp),intent(out) :: etot
    logical :: pr
    integer :: myunit
    logical:: exitRun

    real(wp) :: edisp,ees,ebond,eangl,etors,erep,ehb,exb,ebatm,eext
    real(wp) :: gsolv,gborn,ghb,gsasa,gshift

    integer :: i,j,k,l,m,ij,nd3,iTr,iTri,iTrj,iTrk,iTrl,iTrDum,wscAt,atnb,inb,nbb

    integer :: ati,atj,iat,jat
    integer :: hbA,hbB,nbk,nbnbk
    logical :: ex,require_update
    integer :: nhb1,nhb2,nxb
    real(wp) :: r2,rab,qq0,erff,dd,dum1,r3(3),t8,dum,t22,t39,vec(3)
    real(wp) :: dx,dy,dz,yy,t4,t5,t6,alpha,t20
    real(wp) :: repab,t16,t19,t26,t27,xa,ya,za,cosa,de,t28
    real(wp) :: gammij,eesinf,etmp,phi
    real(wp) :: rn,dr,g3tmp(3,3),g4tmp(3,4)
    real(wp) :: rij,drij(3,n),drijdcn(2)

    real(wp),allocatable :: rinf(:,:)
    real(wp),allocatable :: grab0(:,:,:),rab0(:),eeqtmp(:,:),rabdcn(:,:)
    real(wp),allocatable :: cn(:),dcn(:,:,:),dcndr(:,:,:),dcndL(:,:,:),qtmp(:),dEdcn(:)
    real(wp),allocatable :: hb_cn(:),hb_dcn(:,:,:),dhbcndL(:,:,:)
    real(wp),allocatable :: gTrans(:,:),rTrans(:,:) ! reciprocal, direct translation vector for ES
    real(wp),allocatable :: sqrab(:),srab(:)
    real(wp),allocatable :: g5tmp(:,:)
    integer,allocatable  :: d3list(:,:)
    real(wp),allocatable :: xtmp(:)
    type(gfnff_timer) :: timer
    real(wp) :: dispthr,cnthr,repthr,hbthr1,hbthr2
    real(wp) :: dist(n,n)
    real(wp) :: ds(3,3)
    real(wp) :: convF  !convergence factor alpha, aka ewald parameter
    logical,allocatable :: considered_ABH(:,:,:)
    real(wp) :: mcf_ees,mcf_ehb,mcf_nrep,mcf_s8
    pr = printlevel >= 2
    if (present(printunit)) then
      myunit = printunit
    else
      myunit = stdout
    end if

    if (version == gffVersion%mcgfnff2023) then
      mcf_nrep = 1.343608_wp
      mcf_ees = 0.800222_wp
      mcf_ehb = 0.727406_wp
      mcf_s8 = 2.858671_wp
      mcdisp_par = TDispersionData(s6=1.0_wp,s8=mcf_s8,a1=0.58_wp,a2=4.8_wp,s9=0.0_wp)
      disp_par = TDispersionData(s6=1.0_wp,s8=2.0_wp,a1=0.58_wp,a2=4.8_wp,s9=0.0_wp)
    else
      mcdisp_par = TDispersionData(s6=1.0_wp,s8=2.0_wp,a1=0.58_wp,a2=4.8_wp,s9=0.0_wp)
      disp_par = TDispersionData(s6=1.0_wp,s8=2.0_wp,a1=0.58_wp,a2=4.8_wp,s9=0.0_wp)
      mcf_nrep = 1.0_wp
      mcf_ees = 1.0_wp
      mcf_ehb = 1.0_wp
    end if

    call gfnff_thresholds(accuracy,dispthr,cnthr,repthr,hbthr1,hbthr2)

    vec = cell%lattice(:,1)+cell%lattice(:,2)
    ! get translation vectors within maximum cutoff, but at least central 27 cells (for 3D)
    neigh%oldCutOff = 0.0_wp
    call neigh%getTransVec(n,at,xyz,cell,60.0_wp)

    ! get Distances between atoms for repulsion
    call neigh%getTransVec(n,at,xyz,cell,sqrt(repthr))

!&<
    g(:,:)  = 0.0_wp
    exb     = 0.0_wp
    ehb     = 0.0_wp
    erep    = 0.0_wp
    ees     = 0.0_wp
    edisp   = 0.0_wp
    ebond   = 0.0_wp
    eangl   = 0.0_wp
    etors   = 0.0_wp
    ebatm   = 0.0_wp
    eext    = 0.0_wp

    gsolv   = 0.0d0
    gsasa   = 0.0d0
    gborn   = 0.0d0
    ghb     = 0.0d0
    gshift  = 0.0d0

    sigma(:,:) = 0.0_wp
    etot    = 0.0_wp
!&>

    allocate (sqrab(n*(n+1)/2),srab(n*(n+1)/2),qtmp(n),g5tmp(3,n), &
    &         eeqtmp(2,n*(n+1)/2),d3list(2,n*(n+1)/2),dcn(3,n,n),cn(n), &
    &         dcndr(3,n,n),dcndL(3,3,n),hb_dcn(3,n,n),hb_cn(n),dhbcndL(3,3,n))

    if (printlevel >= 2) then
      call timer%new(10+count([allocated(solvation)]))
    else if (printlevel == 1) then
      ! minimal timer for iteration time only
      call timer%new(1)
      call timer%measure(1,'iter. time')
    end if

    if (pr) call timer%measure(1,'distance/D3 list')
    nd3 = 0
    do i = 1,n
      ij = i*(i-1)/2
      do j = 1,i
        if (j .eq. i) cycle ! dont calc distance to self for non-periodic distances (below)
        k = ij+j
        sqrab(k) = (xyz(1,i)-xyz(1,j))**2+&
        &  (xyz(2,i)-xyz(2,j))**2+&
        &  (xyz(3,i)-xyz(3,j))**2
        if (sqrab(k) .lt. dispthr) then
          nd3 = nd3+1
          d3list(1,nd3) = i
          d3list(2,nd3) = j
        end if
        srab(k) = sqrt(sqrab(k))
      end do

      ! The loop above only runs over the off diagonal elements !
      ! This initializes the unitialized diagonal to zero but does not !
      ! add it to the dispersion list. !
      sqrab(ij+i) = 0.0d0
      srab(ij+i) = 0.0d0
    end do
    if (cell%npbc .ne. 0) then
      dist = 0.0_wp
      !$omp parallel do collapse(2) default(none) shared(dist,xyz) &
      !$omp private(i,j)
      do i = 1,n
        do j = 1,n
          dist(j,i) = NORM2(xyz(:,j)-xyz(:,i))
        end do
      end do
      !$omp end parallel do
    end if

    if (pr) call timer%measure(1)

    !----------!
    ! Setup HB !
    !----------!

    if (pr) call timer%measure(10,'HB/XB (incl list setup)')
    if (allocated(nlist%q)) then
      nlist%initialized = size(nlist%q) == n
    end if
    call gfnff_hbset0(n,at,xyz,topo,nhb1,nhb2,nxb,neigh,nlist,hbthr1,hbthr2)
    nlist%initialized = nlist%initialized.and.nhb1 <= nlist%nhb1 &
       & .and.nhb2 <= nlist%nhb2.and.nxb <= nlist%nxb
    require_update = .not.nlist%initialized
    if (.not.nlist%initialized) then
      if (printlevel >= 2) then
        write (myunit,'(10x,"Number of HB bonds (bound hydrogen)",5x,i0,x,i0,x,i0)') &
              & nhb1
        write (myunit,'(10x,"Number of HB bonds (unbound hydrogen)",3x,i0,x,i0,x,i0)') &
              & nhb2
        write (myunit,'(10x,"Number of XB bonds",22x,i0,x,i0,x,i0)') &
              & nxb
      end if
      call new(nlist,n,5*nhb1,5*nhb2,3*nxb)
      nlist%hbrefgeo(:,:) = xyz
    end if
    if (update.or.require_update) then
      call gfnff_hbset(n,at,xyz,topo,neigh,nlist,hbthr1,hbthr2)
    end if
    if (pr) call timer%measure(10)

    !------------!
    ! Setup GBSA !
    !------------!

    if (allocated(solvation)) then
      call timer%measure(11,"GBSA")
      call solvation%update(at,xyz)
      call timer%measure(11)
    end if

    !------------!
    ! REP part   !
    ! non-bonded !
    !------------!

    if (pr) call timer%measure(2,'non bonded repulsion')
    !$omp parallel do default(none) reduction(+:erep, g, sigma) &
    !$omp shared(n, at, xyz, srab, sqrab, repthr, &
    !$omp topo, param, neigh, mcf_nrep) &
    !$omp private(iat, jat, iTr, iTrDum, m, ij, ati, atj, rab, r2, r3, vec, t8, t16, t19, t26, t27)
    do iat = 1,n
      do jat = 1,iat
        do iTr = 1,neigh%nTrans

          !First calculate erep, g and sigma
          r2 = NORM2(xyz(:,iat)-xyz(:,jat)+neigh%transVec(:,iTr))**2

          ! cycle when above cut-off and when atom would interact with itself
          if (r2 .gt. repthr.OR.r2 .lt. 1.0e-8_wp) cycle

          ! bonded repulsion is calculated seperately and therefore cycled here
          if (iTr .le. neigh%numctr) then
            if (neigh%bpair(iat,jat,iTr) .eq. 1) cycle ! list avoided because of memory
          end if
          ati = at(iat)
          atj = at(jat)
          rab = sqrt(r2)
          t16 = r2**0.75
          t19 = t16*t16
          ! alphanb is the same for all iTr>numctr
          if (iTr .gt. neigh%numctr) then
            iTrDum = neigh%numctr+1
          else
            iTrDum = iTr
          end if
          t8 = t16*topo%alphanb(iat,jat,iTrDum)
          t26 = exp(-t8)*param%repz(ati)*param%repz(atj)*param%repscaln*mcf_nrep
          erep = erep+(t26/rab) !energy
          t27 = t26*(1.5d0*t8+1.0d0)/t19
          r3 = (xyz(:,iat)-xyz(:,jat)+neigh%transVec(:,iTr))*t27
          vec = xyz(:,iat)-xyz(:,jat)+neigh%transVec(:,iTr)
          sigma(:,1) = sigma(:,1)-r3(1)*vec
          sigma(:,2) = sigma(:,2)-r3(2)*vec
          sigma(:,3) = sigma(:,3)-r3(3)*vec
          g(:,iat) = g(:,iat)-r3
          g(:,jat) = g(:,jat)+r3
        end do
      end do
    end do
    !$omp end parallel do
    if (pr) call timer%measure(2)

    ! just a extremely crude mode for 2D-3D conversion !
    ! i.e. an harmonic potential with estimated Re !
    if (version == gffVersion%harmonic2020) then
      ebond = 0
      !$omp parallel do default(none) reduction(+:ebond, g) &
      !$omp shared(topo, param, xyz, at) private(i, iat, jat, rab, r2, r3, rn, dum)
      do i = 1,topo%nbond
        iat = topo%blist(1,i)
        jat = topo%blist(2,i)
        r3 = xyz(:,iat)-xyz(:,jat)
        rab = sqrt(sum(r3*r3))
        rn = 0.7*(param%rcov(at(iat))+param%rcov(at(jat)))
        r2 = rn-rab
        ebond = ebond+0.1d0*r2**2  ! fixfc = 0.1
        dum = 0.1d0*2.0d0*r2/rab
        g(:,jat) = g(:,jat)+dum*r3
        g(:,iat) = g(:,iat)-dum*r3
      end do
      !$omp end parallel do
      etot = ebond+erep
      return
    end if

    !------------------------------!
    ! erf CN and gradient for disp !
    !------------------------------!

    if (cell%npbc .eq. 0) then
      if (pr) call timer%measure(3,'dCN')
      call gfnff_dlogcoord(n,at,xyz,srab,cn,dcn,cnthr,param) ! new erf used in GFN0
      dhbcndL = 0.0_wp
      if (sum(neigh%nr_hb) .gt. 0) call dncoord_erf(n,at,xyz,param%rcov,hb_cn,hb_dcn,900.0d0,topo,neigh,dhbcndL) ! HB erf CN
      if (pr) call timer%measure(3)

    else
      if (pr) call timer%measure(3,'dCN')
      if (sum(neigh%nr_hb) .gt. 0) call dncoord_erf(n,at,xyz,param%rcov,hb_cn,hb_dcn,900.0d0,topo,neigh,dhbcndL) ! HB erf CN
      call getCoordinationNumber(n,at,xyz,neigh%nTrans,neigh%transVec,60.0_wp,5,cn,dcndr,dcndL,param)
      if (pr) call timer%measure(3)
    end if

    !-----!
    ! EEQ !
    !-----!

    if (cell%npbc .eq. 0) then
      if (pr) call timer%measure(4,'EEQ energy and q')
      call goed_gfnff(accuracy .gt. 1,n,at,sqrab,srab,&         ! modified version
      &                dfloat(ichrg),eeqtmp,cn,nlist%q,ees,solvation,param,topo)  ! without dq/dr
      if (pr) call timer%measure(4)
    else
      if (pr) call timer%measure(4,'EEQ energy and q')
      call goed_pbc_gfnff(accuracy .gt. 1,n,at,xyz,dist,cell, &
      & dfloat(ichrg),eeqtmp,cn,nlist%q,ees,solvation,param,topo,gTrans, &
      & rTrans,xtmp,convF)  ! without dq/dr
      ees = ees*mcf_ees
      if (pr) call timer%measure(4)
    end if

    !--------!
    ! D3(BJ) !
    !--------!

    if (cell%npbc .eq. 0) then
      if (pr) call timer%measure(5,'D3')
      if (nd3 .gt. 0) then
        call d3_gradient(topo%dispm,n,at,xyz,nd3,d3list,topo%zetac6, &
        & param%d3r0,sqrtZr4r2,4.0d0,param%dispscale,cn,dcn,edisp,g)
      end if
      deallocate (d3list)
      if (pr) call timer%measure(5)
    else
      ! use adjusted dispersion parameters for inter-molecular interactions
      ! calculate inter-molecular dispersion
      call d3_gradientPBC(topo%dispm,n,at,xyz,cell,topo%fraglist,neigh%nTrans,neigh%transVec,mcdisp_par,4.0_wp,topo%zetac6, &
           & param%d3r0,60.0_wp,.true.,cn,dcndr,dcndL,edisp,g,sigma)
      ! calculate intra-molecular dispersion
      call d3_gradientPBC(topo%dispm,n,at,xyz,cell,topo%fraglist,neigh%nTrans,neigh%transVec,disp_par,4.0_wp,topo%zetac6, &
           & param%d3r0,60.0_wp,.false.,cn,dcndr,dcndL,edisp,g,sigma)
    end if

    !---------!
    ! ES part !
    !---------!

    if (cell%npbc .eq. 0) then
      if (pr) call timer%measure(6,'EEQ gradient')
      !$omp parallel do default(none) reduction (+:g) &
      !$omp shared(nlist,n,sqrab,srab,eeqtmp,xyz,at) &
      !$omp private(i,j,k,ij,r3,r2,rab,gammij,erff,dd)
      do i = 1,n
        k = i*(i-1)/2
        do j = 1,i-1
          ij = k+j
          r2 = sqrab(ij)
          rab = srab(ij)
          gammij = eeqtmp(1,ij)
          erff = eeqtmp(2,ij)
          dd = (2.0d0*gammij*exp(-gammij**2*r2) &
                   & /(sqrtpi*r2)-erff/(rab*r2))*nlist%q(i)*nlist%q(j)
          r3 = (xyz(:,i)-xyz(:,j))*dd
          g(:,i) = g(:,i)+r3
          g(:,j) = g(:,j)-r3
        end do
      end do
      !$omp end parallel do
      if (.not.pr) deallocate (eeqtmp)

      if (allocated(solvation)) then
        call timer%measure(11,"GBSA")
        call solvation%addGradient(at,xyz,nlist%q,nlist%q,g)
        call solvation%getEnergyParts(nlist%q,nlist%q,gborn,ghb,gsasa, &
        & gshift)
        gsolv = gsasa+gborn+ghb+gshift
        call timer%measure(11)
      else
        gborn = 0.0d0
        ghb = 0.0d0
      end if

      do i = 1,n
        qtmp(i) = nlist%q(i)*param%cnf(at(i))/(2.0d0*sqrt(cn(i))+1.d-16)
      end do

      call gemv(dcn,qtmp,g,alpha=-1.0_wp,beta=1.0_wp)
      if (pr) call timer%measure(6)
    else ! periodic case
      if (pr) call timer%measure(6,'EEQ gradient')

      call es_grad_sigma(n,at,xyz,cell,topo,nlist,rTrans,gTrans,xtmp,convF, &
                 & sigma,g,mcf_ees)

      if (.not.pr) deallocate (eeqtmp)

      if (allocated(solvation)) then
        call timer%measure(11,"GBSA")
        call solvation%addGradient(at,xyz,nlist%q,nlist%q,g)
        call solvation%getEnergyParts(nlist%q,nlist%q,gborn,ghb,gsasa, &
        & gshift)
        gsolv = gsasa+gborn+ghb+gshift
        call timer%measure(11)
      else
        gborn = 0.0d0
        ghb = 0.0d0
      end if

      ! qtmp =  q * dXdcn | where X is the right-hand side
      do i = 1,n
        qtmp(i) = nlist%q(i)*param%cnf(at(i))/(2.0d0*sqrt(cn(i))+1.d-16)
      end do

      call gemv(dcndr,qtmp,g,alpha=-mcf_ees,beta=1.0_wp)
      call gemv(dcndL,qtmp,sigma,alpha=-mcf_ees,beta=1.0_wp)

      if (pr) call timer%measure(6)
    end if

    !-----------------!
    ! SRB bonded part !
    !-----------------!

    if (pr) call timer%measure(7,'bonds')
    if (neigh%nbond .gt. 0) then
      ! rab0 and vbond go over nbond not atoms, vbond setup with pbc in gfnff_ini.f90
      allocate (grab0(3,n,neigh%nbond),rab0(neigh%nbond),rabdcn(2,neigh%nbond))
      rab0(:) = neigh%vbond(1,:)
      if (cell%npbc .eq. 0) then
        call gfnffdrab(n,at,cn,dcn,neigh%nbond,neigh%blist,rab0,grab0,rabdcn)
      else
        call gfnffdrab(n,at,cn,dcndr,neigh%nbond,neigh%blist,rab0,grab0,rabdcn)
      end if
      deallocate (dcn,dcndr)
      allocate (dEdcn(n),source=0.0_wp)
      allocate (considered_ABH(topo%hb_mapNAB,topo%hb_mapNAB,topo%hb_mapNH),source=.false.)

      !$omp parallel do default(none) reduction(+:g, ebond, sigma, dEdcn) &
      !$omp shared(grab0, topo, neigh, param, rab0, rabdcn, xyz, at, hb_cn, hb_dcn, n, dhbcndL, considered_ABH) &
      !$omp private(i, k, iat, jat, ij, rab, rij, drij, drijdcn, t8, dr, dum, yy, &
      !$omp& dx, dy, dz, t4, t5, t6, ati, atj, iTr)
      do i = 1,neigh%nbond
        jat = neigh%blist(1,i)
        iat = neigh%blist(2,i)
        iTr = neigh%blist(3,i)
        if (iTr .gt. neigh%nTrans) cycle
        ati = at(iat)
        atj = at(jat)
        ij = iat*(iat-1)/2+jat
        rab = NORM2(xyz(:,jat)-xyz(:,iat)+neigh%transVec(:,iTr))
        rij = rab0(i)
        drij = grab0(:,:,i)
        drijdcn = rabdcn(:,i)
        if (neigh%nr_hb(i) .ge. 1) then
          call egbond_hb(i,iat,jat,iTr,rab,rij,drij,drijdcn,hb_cn,hb_dcn,n,at,xyz,&
               &ebond,g,sigma,param,topo,neigh,dEdcn,dhbcndL,considered_ABH)
        else
          call egbond(i,iat,jat,iTr,rab,rij,drij,drijdcn,n,at,xyz,ebond,g,sigma,neigh,rab0,dEdcn)
        end if
      end do
      !$omp end parallel do
      call gemv(dcndL,dEdcn,sigma,alpha=1.0_wp,beta=1.0_wp)

      deallocate (dcndL)
      deallocate (hb_dcn)

      ! bonded REP !

      !$omp parallel do default(none) reduction(+:erep, g, sigma) &
      !$omp shared(topo, param, at, xyz, neigh) &
      !$omp private(i, iTr, iat, jat, ij, xa, ya, za, dx, dy, dz, r2, rab, ati, atj, &
      !$omp& alpha, repab, t16, t19, t26, t27)
      do i = 1,neigh%nbond
        jat = neigh%blist(1,i)
        iat = neigh%blist(2,i)
        iTr = neigh%blist(3,i)
        xa = xyz(1,iat)
        ya = xyz(2,iat)
        za = xyz(3,iat)
        dx = xa-xyz(1,jat)-neigh%transVec(1,iTr)
        dy = ya-xyz(2,jat)-neigh%transVec(2,iTr)
        dz = za-xyz(3,jat)-neigh%transvec(3,iTr)
        rab = NORM2(xyz(:,iat)-(xyz(:,jat)+neigh%transVec(:,iTr)))
        r2 = rab**2
        ati = at(iat)
        atj = at(jat)
        alpha = sqrt(param%repa(ati)*param%repa(atj))
        repab = param%repz(ati)*param%repz(atj)*param%repscalb
        t16 = r2**0.75d0
        t19 = t16*t16
        t26 = exp(-alpha*t16)*repab
        erep = erep+t26/rab !energy
        t27 = t26*(1.5d0*alpha*t16+1.0d0)/t19
        g(1,iat) = g(1,iat)-dx*t27
        g(2,iat) = g(2,iat)-dy*t27
        g(3,iat) = g(3,iat)-dz*t27
        g(1,jat) = g(1,jat)+dx*t27
        g(2,jat) = g(2,jat)+dy*t27
        g(3,jat) = g(3,jat)+dz*t27
        sigma(1,1) = sigma(1,1)-1.0_wp*dx*t27*dx
        sigma(1,2) = sigma(1,2)-1.0_wp*dx*t27*dy
        sigma(1,3) = sigma(1,3)-1.0_wp*dx*t27*dz
        sigma(2,1) = sigma(2,1)-1.0_wp*dy*t27*dx
        sigma(2,2) = sigma(2,2)-1.0_wp*dy*t27*dy
        sigma(2,3) = sigma(2,3)-1.0_wp*dy*t27*dz
        sigma(3,1) = sigma(3,1)-1.0_wp*dz*t27*dx
        sigma(3,2) = sigma(3,2)-1.0_wp*dz*t27*dy
        sigma(3,3) = sigma(3,3)-1.0_wp*dz*t27*dz
      end do
      !$omp end parallel do
    end if ! if neigh%nbond.gt.0
    if (pr) call timer%measure(7)

    !------!
    ! bend !
    !------!

    if (pr) call timer%measure(8,'bend and torsion')
    if (topo%nangl .gt. 0) then
      !$omp parallel do default(none) reduction (+:eangl, g, sigma) &
      !$omp shared(n, at, xyz, topo, neigh, param) &
      !$omp private(m, j, i, k, etmp, g3tmp,ds)
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
    end if

    !---------!
    ! torsion !
    !---------!

    if (topo%ntors .gt. 0) then
      !$omp parallel do default(none) reduction(+:etors, g, sigma) &
      !$omp shared(param, topo, neigh, n, at, xyz) &
      !$omp private(m, i, j, k, l,iTrl,iTrj,iTrk, etmp, g4tmp,ds)
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
    end if

    !----------------------------------------!
    ! triple bonded carbon torsion potential !
    !----------------------------------------!

    if (allocated(topo%sTorsl)) then
      m = size(topo%sTorsl(1,:))
      if (m .ne. 0) then
        do i = 1,m
          call sTors_eg(m,n,xyz,topo,etmp,g5tmp)
          etors = etors+etmp
          g = g+g5tmp
        end do
      end if
    end if
    if (pr) call timer%measure(8)

    !------------!
    ! BONDED ATM !
    !------------!

    if (pr) call timer%measure(9,'bonded ATM')
    if (topo%nbatm .gt. 0) then
      !$omp parallel do default(none) reduction(+:ebatm, g, sigma) &
      !$omp shared(n, at, xyz, srab, sqrab, topo, neigh, param) &
      !$omp private(i, j, k, l, iTrk, iTrl, etmp, g3tmp, ds)
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
        sigma = sigma+ds
        ebatm = ebatm+etmp
      end do
      !$omp end parallel do
    end if
    if (pr) call timer%measure(9)

    !-----!
    ! EHB !
    !-----!

    ! get correct num of transVec for hb lists (e.g. hblist1/2)
    call neigh%getTransVec(n,at,xyz,cell,sqrt(hbthr2))
    if (pr) call timer%measure(10,'HB/XB (incl list setup)')
    if (update.or.require_update) then
      call gfnff_hbset(n,at,xyz,topo,neigh,nlist,hbthr1,hbthr2)
    end if

    if (nlist%nhb1 .gt. 0) then
      !$omp parallel do default(none) reduction(+:ehb, g, sigma) &
      !$omp shared(topo,nlist, neigh, param, n, at, xyz, sqrab, srab, mcf_ehb) &
      !$omp private(i, j, k, l, iTri, iTrj, etmp, g3tmp)
      do i = 1,nlist%nhb1
        j = nlist%hblist1(1,i)
        k = nlist%hblist1(2,i)
        l = nlist%hblist1(3,i)
        iTri = nlist%hblist1(4,i)
        iTrj = nlist%hblist1(5,i)
        if (iTri .gt. neigh%nTrans.or.iTrj .gt. neigh%nTrans) cycle
        call abhgfnff_eg1(n,j,k,l,iTri,iTrj,at,xyz,topo%qa,etmp,&
                & g3tmp,param,topo,neigh,sigma,mcf_ehb)
        g(1:3,j) = g(1:3,j)+g3tmp(1:3,1)*mcf_ehb
        g(1:3,k) = g(1:3,k)+g3tmp(1:3,2)*mcf_ehb
        g(1:3,l) = g(1:3,l)+g3tmp(1:3,3)*mcf_ehb
        ehb = ehb+etmp*mcf_ehb
      end do
      !$omp end parallel do
    end if
    exitRun = .false.
    if (nlist%nhb2 .gt. 0) then

      !$omp parallel do default(none) reduction(+:ehb, g, sigma) &
      !$omp shared(topo,nlist, neigh, param, n, at, xyz, sqrab, srab, exitRun, mcf_ehb) &
      !$omp private(i, j, k, l,iTrj,iTrk,iTrDum,iTr, nbb, nbk, nbnbk, atnb, etmp, g5tmp)
      do i = 1,nlist%nhb2
        j = nlist%hblist2(1,i)  !   A
        k = nlist%hblist2(2,i)  !   B
        l = nlist%hblist2(3,i)  !   H -> always in central cell
        iTrj = nlist%hblist2(4,i) ! iTrA
        iTrk = nlist%hblist2(5,i) ! iTrB

        if (iTrj .gt. neigh%nTrans.or.iTrk .gt. neigh%nTrans) cycle
        ! prepare variables for check in carbonyl/nitro case !
        ! Number neighbors of C/N should be > 1 for carbonyl/nitro !
        ! get needed values for Carbonyl or Nitro case
        nbnbk = 0
        if (at(k) .eq. 8.and.sum(neigh%nb(neigh%numnb,k,:)) .eq. 1) then
          nbk = 0
          iTr = 0 ! nbk is the first neighbor of k !
          atnb = 0
          call neigh%jth_nb(n,xyz,nbk,1,k,iTr)
          ! get iTrC and cycle if out of cutoff (the cutoff used in last getTransVec call)
          iTrDum = neigh%fTrSum(iTr,iTrk)
          if (iTrDum .eq. -1.or.iTrDum .gt. neigh%nTrans) cycle
          ! get number neighbors of neighbor of k !
          if (nbk .ne. 0) then
            nbnbk = sum(neigh%nb(neigh%numnb,nbk,:))
            atnb = at(nbk)
          end if
        end if

        ! Carbonyl case R-C=O...H_A !
        if (at(k) .eq. 8.and.sum(neigh%nb(neigh%numnb,k,:)) .eq. 1.and.atnb .eq. 6 &
              & .and.nbnbk .gt. 1) then
          call abhgfnff_eg3(n,j,k,l,iTrj,iTrk,nbk,iTrDum,at,xyz,topo%qa,sqrab,&
                  & srab,etmp,g5tmp,param,topo,neigh,sigma,exitRun,mcf_ehb)

          ! Nitro case R-N=O...H_A !
        else if (at(k) .eq. 8.and.sum(neigh%nb(neigh%numnb,k,:)) .eq. 1.and.atnb .eq. 7 &
              &  .and.nbnbk .gt. 1) then
          call abhgfnff_eg3(n,j,k,l,iTrj,iTrk,nbk,iTrDum,at,xyz,topo%qa,sqrab,&
                  & srab,etmp,g5tmp,param,topo,neigh,sigma,exitRun,mcf_ehb)

          ! N hetero aromat !
        else if (at(k) .eq. 7.and.sum(neigh%nb(neigh%numnb,k,:)) .eq. 2) then
          call abhgfnff_eg2_rnr(n,j,k,l,iTrj,iTrk,at,xyz,topo%qa,sqrab,&
                   & srab,etmp,g5tmp,param,topo,neigh,sigma,mcf_ehb)

          ! default !
        else
          nbb = sum(neigh%nb(neigh%numnb,k,:))
          call abhgfnff_eg2new(n,j,k,l,iTrj,iTrk,nbb,at,xyz,topo%qa,sqrab,srab, &
             & etmp,g5tmp,param,topo,neigh,sigma,mcf_ehb)
        end if
        g = g+g5tmp*mcf_ehb
        ehb = ehb+etmp*mcf_ehb
        nlist%hbe2(i) = etmp

      end do
      !$omp end parallel do
    end if

    !-----!
    ! EXB !
    !-----!

    if (nlist%nxb .gt. 0) then
      !$omp parallel do default(none) reduction(+:exb, g, sigma) &
      !$omp shared(topo, neigh, nlist, param, n, at, xyz) &
      !$omp private(i, j, k, l, iTrk, iTrl, etmp, g3tmp)
      do i = 1,nlist%nxb
        j = nlist%hblist3(1,i)   ! A in central cell
        k = nlist%hblist3(2,i)   ! B
        l = nlist%hblist3(3,i)   ! X
        iTrk = nlist%hblist3(4,i) !iTrB
        iTrl = nlist%hblist3(5,i) !iTrX
        if (iTrk .gt. neigh%nTrans.or.iTrl .gt. neigh%nTrans) cycle
        if (j .ne. 0.and.k .ne. 0) then
          call rbxgfnff_eg(n,j,k,l,iTrk,iTrl,at,xyz,topo%qa,etmp,g3tmp,param,neigh,sigma)
          g(1:3,j) = g(1:3,j)+g3tmp(1:3,1)
          g(1:3,k) = g(1:3,k)+g3tmp(1:3,2)
          g(1:3,l) = g(1:3,l)+g3tmp(1:3,3)
          exb = exb+etmp
          nlist%hbe3(i) = etmp
        end if
      end do
      !$omp end parallel do
    end if
    if (pr) call timer%measure(10)

    ! external stuff !
    if (sum(abs(efield)) .gt. 1d-6) then
      do i = 1,n
        r3(:) = -nlist%q(i)*efield(:)
        g(:,i) = g(:,i)+r3(:)
        eext = eext+r3(1)*(xyz(1,i)-topo%xyze0(1,i))+&
        &                    r3(2)*(xyz(2,i)-topo%xyze0(2,i))+&
        &                    r3(3)*(xyz(3,i)-topo%xyze0(3,i))
      end do
    end if

    !--------------!
    ! total energy !
    !--------------!

    etot = ees+edisp+erep+ebond &
    &           +eangl+etors+ehb+exb+ebatm+eext &
    &           +gsolv

    !----------!
    ! printout !
    !----------!
    if (printlevel >= 2) then

      call timer%write(myunit,'E+G')
      if (abs(sum(nlist%q)-ichrg) .gt. 1.d-1) then ! check EEQ only once
        if (printlevel >= 1) then
          write (myunit,*) nlist%q
          write (myunit,*) sum(nlist%q),ichrg
          write (myunit,'("**ERROR**",a,1x,a)') 'EEQ charge constrain error',source
        end if
        return
      end if
      r3 = 0
      do i = 1,n
        r3(:) = r3(:)+nlist%q(i)*xyz(:,i)
      end do

      ! just for fit De calc !
      sqrab = 1.d+12
      srab = 1.d+6
      cn = 0
      allocate (rinf(n,n),source=1.d+6)

      ! asymtotically for R=inf, Etot is the SIE contaminted EES !
      ! which is computed here to get the atomization energy De,n,at(n) !
      if (cell%npbc .eq. 0) then
        call goed_gfnff(.true.,n,at,sqrab,srab,dfloat(ichrg),eeqtmp,cn,qtmp,eesinf,solvation,param,topo)
      else
        if (allocated(gTrans)) deallocate (gTrans)
        if (allocated(rTrans)) deallocate (rTrans)
        call goed_pbc_gfnff(.true.,n,at,xyz,rinf,cell, &
        & dfloat(ichrg),eeqtmp,cn,qtmp,eesinf,solvation,param,topo,gTrans, &
        & rTrans,xtmp,convF)  ! without dq/dr
      end if
      de = -(etot-eesinf)

      ! if geometry optimization !
    else if (printlevel == 1) then
      call timer%measure(1)
      ! stop timer and add tag !
      !call timer%write_timing(myunit,1)
    end if

    ! write resusts to res type !
    res_gff%e_total = etot
    res_gff%gnorm = sqrt(sum(g**2))
    res_gff%e_bond = ebond
    res_gff%e_angl = eangl
    res_gff%e_tors = etors
    res_gff%e_es = ees
    res_gff%e_rep = erep
    res_gff%e_disp = edisp
    res_gff%e_hb = ehb
    res_gff%e_xb = exb
    res_gff%e_batm = ebatm
    res_gff%e_ext = eext
    res_gff%g_hb = ghb
    res_gff%g_born = gborn
    res_gff%g_solv = gsolv
    res_gff%g_shift = gshift
    res_gff%g_sasa = gsasa

    call gemv(xyz,nlist%q,res_gff%dipole)

  end subroutine gfnff_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egbond(i,iat,jat,iTr,rab,rij,drij,drijdcn,n,at,xyz,e,g,sigma,neigh,rab0,dEdcn)
    implicit none
    !Dummy
    !type(TGFFTopology), intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    integer,intent(in)   :: i
    integer,intent(in)   :: n
    integer,intent(in)   :: iat
    integer,intent(in)   :: jat
    integer,intent(in)   :: at(n)
    integer,intent(in)   :: iTr
    real(wp),intent(in)    :: rab
    real(wp),intent(in)    :: rij
    real(wp),intent(in)    :: drij(3,n)
    real(wp),intent(in)    :: drijdcn(2)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(inout) :: dEdcn(n)
    real(wp),intent(inout) :: e
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)
    !Stack
    integer j,k
    real(wp) dr,dum,ri(3),rj(3),rk(3),vrab(3)
    real(wp) dx,dy,dz,dg(3)
    real(wp) yy
    real(wp) t4,t5,t6,t8
    real(wp) rab0(neigh%nbond)

    t8 = neigh%vbond(2,i)
    dr = rab-rij
    dum = neigh%vbond(3,i)*exp(-t8*dr**2)
    e = e+dum                      ! bond energy
    yy = 2.0d0*t8*dr*dum
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
    sigma(:,1) = sigma(:,1)+dg(1)*vrab
    sigma(:,2) = sigma(:,2)+dg(2)*vrab
    sigma(:,3) = sigma(:,3)+dg(3)*vrab

    t4 = yy*(dx/rab)
    t5 = yy*(dy/rab)
    t6 = yy*(dz/rab)
    dg(1) = t4
    dg(2) = t5
    dg(3) = t6
    g(1,jat) = g(1,jat)+t4 ! to avoid if in loop below
    g(2,jat) = g(2,jat)+t5
    g(3,jat) = g(3,jat)+t6
    ! 3B sigma ! product rule, rij depends on rab through cn
    dEdcn(jat) = dEdcn(jat)+yy*drijdcn(2)

    do k = 1,n !3B gradient
      dg = drij(:,k)*yy
      g(:,k) = g(:,k)+dg
      sigma(:,1) = sigma(:,1)+dg(1)*vrab
      sigma(:,2) = sigma(:,2)+dg(2)*vrab
      sigma(:,3) = sigma(:,3)+dg(3)*vrab
    end do

  end subroutine egbond

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egbond_hb(i,iat,jat,iTr,rab,rij,drij,drijdcn,hb_cn,hb_dcn,n,at,xyz,e,&
                  &g,sigma,param,topo,neigh,dEdcn,dhbcndL,considered_ABH)
    implicit none
    !Dummy
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
    real(wp),intent(in)    :: drij(3,n)
    real(wp),intent(in)    :: drijdcn(2)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(in)    :: hb_cn(n)
    real(wp),intent(in)    :: hb_dcn(3,n,n)
    real(wp),intent(in) :: dhbcndL(3,3,n)
    logical,intent(inout)  :: considered_ABH(topo%hb_mapNAB,topo%hb_mapNAB,topo%hb_mapNH)! only consider ABH triplets once; indep of iTr
    real(wp),intent(inout) :: e
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)
    real(wp),intent(inout) :: dEdcn(n)
    real(wp) :: rbh(3),r2
    !Stack
    integer j,k
    integer jA,jH,iTrA,iTrH,iTrB
    integer hbH,hbB,hbA
    integer mapA,mapB,mapH
    real(wp) dr,dum
    real(wp) dx,dy,dz,vrab(3),dg(3)
    real(wp) yy,zz
    real(wp) t1,t4,t5,t6,t8

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
    t8 = (-t1*hb_cn(hbH)+1.0)*neigh%vbond(2,i)
    dr = rab-rij
    dum = neigh%vbond(3,i)*exp(-t8*dr**2)
    e = e+dum                      ! bond energy
    yy = 2.0d0*t8*dr*dum
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
    sigma(:,1) = sigma(:,1)+dg(1)*vrab
    sigma(:,2) = sigma(:,2)+dg(2)*vrab
    sigma(:,3) = sigma(:,3)+dg(3)*vrab

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

    do k = 1,n !3B gradient
      dg = drij(:,k)*yy
      g(:,k) = g(:,k)+dg
      sigma(:,1) = sigma(:,1)+dg(1)*vrab
      sigma(:,2) = sigma(:,2)+dg(2)*vrab
      sigma(:,3) = sigma(:,3)+dg(3)*vrab
    end do

    zz = dum*neigh%vbond(2,i)*dr**2*t1
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
          !
          sigma = sigma+zz*dhbcndL(:,:,iat)
        end if
        if (hbH .eq. jat) then ! only hbH cf. energy -> dum above
          !
          sigma = sigma+zz*dhbcndL(:,:,jat)
        end if
        do k = 1,topo%bond_hb_Bn(j)
          hbB = topo%bond_hb_B(1,k,j)
          iTrB = topo%bond_hb_B(2,k,j)
          ! only add gradient one time per ABH triple (independent of iTrB)
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

  subroutine dncoord_erf(nat,at,xyz,rcov,cn,dcn,thr,topo,neigh,dcndL)

    implicit none

    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    integer,intent(in)   :: nat
    integer,intent(in)   :: at(nat)
    real(wp),intent(in)  :: xyz(3,nat)
    real(wp),intent(in)  :: rcov(:)
    real(wp),intent(out) :: cn(nat)
    real(wp),intent(out) :: dcn(3,nat,nat)
    !> Derivative of the CN with respect to strain deformations.
    real(wp),intent(out) :: dcndL(:,:,:)
    real(wp),intent(in),optional :: thr
    real(wp) :: cn_thr

    integer  :: i,j,iTrB,iTrH
    integer  :: lin,linAH
    integer  :: iat,jat
    integer  :: jA,jH
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
        if (r2 .gt. thr) cycle
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egbend(m,j,i,k,n,at,xyz,e,g,ds,param,topo,neigh)
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
    real(wp) theta,deda(3),vp(3),et,dij,c1
    real(wp) term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
    real(wp) dampjl,damp2jl,rn
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
    ! loop over stress tensor dimensions
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

  end subroutine egbend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egbend_nci_mul(j,i,k,vTrB,vTrC,c0,fc,n,at,xyz,e,g)
    implicit none

    integer n,at(n)
    integer j,i,k ! B,C,H
    real(wp) vTrB(3),vTrC(3)
    real(wp) c0,fc
    real(wp) xyz(3,n),g(3,3),e

    real(wp) kijk,va(3),vb(3),vc(3),cosa
    real(wp) dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) theta,deda(3),vp(3),et,dij,c1
    real(wp) term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
    real(wp) dampjl,damp2jl,rn

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egbend_nci(j,i,k,c0,kijk,n,at,xyz,e,g,param)

    implicit none

    type(TGFFData),intent(in) :: param
    integer n,at(n)
    integer i,j,k
    real(wp) c0,kijk
    real(wp) xyz(3,n),g(3,3),e

    real(wp) va(3),vb(3),vc(3),cosa
    real(wp) dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) theta,deda(3),vp(3),et,dij,c1
    real(wp) term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
    real(wp) dampjl,damp2jl,rn

    va(1:3) = xyz(1:3,i)
    vb(1:3) = xyz(1:3,j)
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

    call gfnffdampa_nci(at(i),at(j),rab2,dampij,damp2ij,param)
    call gfnffdampa_nci(at(k),at(j),rcb2,dampjk,damp2jk,param)
    damp = dampij*dampjk

    ! linear !
    if (pi-c0 .lt. 1.d-6) then
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

  end subroutine egbend_nci

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egtors(m,i,j,k,l,iTrl,iTrj,iTrk,n,at,xyz,e,g,ds,param,topo,neigh)
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    integer :: m,n,at(n)
    integer :: i,j,k,l,iTrl,iTrj,iTrk,dim1,dim2
    real(wp) :: xyz(3,n),g(3,4),e,ds(3,3)
    real(wp) :: vTrl(3),vTrj(3),vTrk(3)

    real(wp) :: c0,kijk,va(3),vb(3),vc(3),vd(3),cosa
    real(wp) :: dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) :: term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) :: rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) :: theta,deda(3),vp(3),et,dij,c1
    real(wp) :: term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) :: ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) :: rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
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

      phi = valijklffPBC(1,n,xyz,i,j,k,l,vTrj,vTrk,vTrl)
      call dphidrPBC(2,n,xyz,i,j,k,l,vTrj,vTrk,vTrl,phi,dda,ddb,ddc,ddd)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> torsion without distance damping !!! damping is inherint in the HB term
  subroutine egtors_nci_mul(i,j,k,l,vTrR,vTrB,vTrC,rn,phi,phi0,tshift,n,at,xyz,e,g)
    implicit none
    !Dummy
    integer n,at(n)
    integer i,j,k,l
    integer rn
    real(wp) :: vTrR(3),vTrB(3),vTrC(3)
    real(wp) :: phi,phi0,tshift
    real(wp) :: xyz(3,n),g(3,4),e
    !Stack
    real(wp) :: c0,fc,kijk,va(3),vb(3),vc(3),cosa
    real(wp) :: dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) :: term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) :: rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) :: theta,deda(3),vp(3),et,dij,c1
    real(wp) :: term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) :: ddd(3),ddc(3),ddb(3),dda(3),rjl
    real(wp) :: rij,rijk,rkl,rjk,dampkl,damp2kl
    real(wp) :: dampjl,damp2jl

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egtors_nci(i,j,k,l,rn,phi0,fc,n,at,xyz,e,g,param)

    implicit none

    type(TGFFData),intent(in) :: param
    integer :: n,at(n)
    integer :: i,j,k,l
    integer :: rn
    real(wp) :: phi0,fc
    real(wp) :: xyz(3,n),g(3,4),e

    real(wp) :: c0,kijk,va(3),vb(3),vc(3),cosa
    real(wp) :: dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) :: term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) :: rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) :: theta,deda(3),vp(3),et,dij,c1
    real(wp) :: term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) :: ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) :: rij,rijk,rkl,rjk,dampkl,damp2kl
    real(wp) :: dampjl,damp2jl

    vab(1:3) = xyz(1:3,i)-xyz(1:3,j)
    vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
    vdc(1:3) = xyz(1:3,k)-xyz(1:3,l)

    rij = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
    rjk = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
    rkl = vdc(1)*vdc(1)+vdc(2)*vdc(2)+vdc(3)*vdc(3)

    call gfnffdampt_nci(at(i),at(j),rij,dampij,damp2ij,param)
    call gfnffdampt_nci(at(k),at(j),rjk,dampjk,damp2jk,param)
    call gfnffdampt_nci(at(k),at(l),rkl,dampkl,damp2kl,param)

    damp = dampjk*dampij*dampkl
    phi = valijklff(n,xyz,i,j,k,l)

    call dphidr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)

    dphi1 = phi-phi0
    c1 = rn*dphi1+pi
    x1cos = cos(c1)
    x1sin = sin(c1)
    et = (1.+x1cos)*fc
    dij = -rn*x1sin*fc*damp

    term1(1:3) = et*damp2ij*dampjk*dampkl*vab(1:3)
    term2(1:3) = et*damp2jk*dampij*dampkl*vcb(1:3)
    term3(1:3) = et*damp2kl*dampij*dampjk*vdc(1:3)

    g(1:3,1) = dij*dda(1:3)+term1
    g(1:3,2) = dij*ddb(1:3)-term1+term2
    g(1:3,3) = dij*ddc(1:3)+term3-term2
    g(1:3,4) = dij*ddd(1:3)-term3

    e = et*damp

  end subroutine egtors_nci

!cccccccccccccccccccccccccccccccccccccccccccccc
! damping of bend and torsion for long
! bond distances to allow proper dissociation
!cccccccccccccccccccccccccccccccccccccccccccccc

  subroutine gfnffdampa(ati,atj,r2,damp,ddamp,param)

    implicit none
    type(TGFFData),intent(in) :: param
    integer ati,atj
    real(wp) r2,damp,ddamp,rr,rcut

    rcut = param%atcuta*(param%rcov(ati)+param%rcov(atj))**2
    rr = (r2/rcut)**2
    damp = 1.0d0/(1.0d0+rr)
    ddamp = -2.d0*2*rr/(r2*(1.0d0+rr)**2)

  end subroutine gfnffdampa

  subroutine gfnffdampt(ati,atj,r2,damp,ddamp,param)

    implicit none
    type(TGFFData),intent(in) :: param
    integer ati,atj
    real(wp) r2,damp,ddamp,rr,rcut

    rcut = param%atcutt*(param%rcov(ati)+param%rcov(atj))**2
    rr = (r2/rcut)**2
    damp = 1.0d0/(1.0d0+rr)
    ddamp = -2.d0*2*rr/(r2*(1.0d0+rr)**2)

  end subroutine gfnffdampt

  subroutine gfnffdampa_nci(ati,atj,r2,damp,ddamp,param)

    implicit none
    type(TGFFData),intent(in) :: param
    integer ati,atj
    real(wp) r2,damp,ddamp,rr,rcut

    rcut = param%atcuta_nci*(param%rcov(ati)+param%rcov(atj))**2
    rr = (r2/rcut)**2
    damp = 1.0d0/(1.0d0+rr)
    ddamp = -2.d0*2*rr/(r2*(1.0d0+rr)**2)

  end subroutine gfnffdampa_nci

  subroutine gfnffdampt_nci(ati,atj,r2,damp,ddamp,param)

    implicit none
    type(TGFFData),intent(in) :: param
    integer ati,atj
    real(wp) r2,damp,ddamp,rr,rcut

    rcut = param%atcutt_nci*(param%rcov(ati)+param%rcov(atj))**2
    rr = (r2/rcut)**2
    damp = 1.0d0/(1.0d0+rr)
    ddamp = -2.d0*2*rr/(r2*(1.0d0+rr)**2)

  end subroutine gfnffdampt_nci

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ref.: S. Alireza Ghasemi, Albert Hofstetter, Santanu Saha, and Stefan Goedecker
!       PHYSICAL REVIEW B 92, 045131 (2015)
!       Interatomic potentials for ionic systems with density functional accuracy
!       based on charge densities obtained by a neural network
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine goed_gfnff(single,n,at,sqrab,r,chrg,eeqtmp,cn,q,es,gbsa,param,topo)

    implicit none

    character(len=*),parameter :: source = 'gfnff_eg_goed'
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo

    !> real*4 flag for solver
    logical,intent(in)  :: single

    !> number of atoms
    integer,intent(in)  :: n

    !> ordinal numbers
    integer,intent(in)  :: at(n)

    !> squared dist
    real(wp),intent(in)  :: sqrab(n*(n+1)/2)

    !> distance
    real(wp),intent(in)  :: r(n*(n+1)/2)

    !> total charge
    real(wp),intent(in)  :: chrg

    !> coordination number
    real(wp),intent(in)  :: cn(n)

    !> output charges
    real(wp),intent(out) :: q(n)

    !> ES term
    real(wp),intent(out) :: es

    !>intermediates
    real(wp),intent(out) :: eeqtmp(2,n*(n+1)/2)

    !> Solvation
    type(TBorn),allocatable,intent(in) :: gbsa

    integer  :: m,i,j,k,ii,ij,io1,io2
    integer,allocatable :: ipiv(:)
    real(wp) :: gammij,tsqrt2pi,r2,tmp
    real(wp),allocatable :: A(:,:),x(:)
    real(sp),allocatable :: A4(:,:),x4(:)
    parameter(tsqrt2pi=0.797884560802866_wp)
    logical :: exitRun

    ! # atoms + chrg constrain + frag constrain !
    m = n+topo%nfrag
    allocate (A(m,m),x(m))

    !  setup RHS !

    do i = 1,n
      x(i) = topo%chieeq(i)+param%cnf(at(i))*sqrt(cn(i))
    end do

    A = 0

    !  setup A matrix !

    !$omp parallel default(none) &
    !$omp shared(topo,n,sqrab,r,eeqtmp,A,at) &
    !$omp private(i,j,k,ij,gammij,tmp)
    !$omp do schedule(dynamic)
    do i = 1,n

      A(i,i) = tsqrt2pi/sqrt(topo%alpeeq(i))+topo%gameeq(i) ! J of i !
      k = i*(i-1)/2

      do j = 1,i-1

        ij = k+j
        gammij = 1./sqrt(topo%alpeeq(i)+topo%alpeeq(j)) ! squared above !
        tmp = erf(gammij*r(ij))
        eeqtmp(1,ij) = gammij
        eeqtmp(2,ij) = tmp
        A(j,i) = tmp/r(ij)
        A(i,j) = A(j,i)

      end do
    end do
    !$omp enddo
    !$omp end parallel

    !  fragment charge constrain !

    do i = 1,topo%nfrag
      x(n+i) = topo%qfrag(i)
      do j = 1,n
        if (topo%fraglist(j) .eq. i) then
          A(n+i,j) = 1
          A(j,n+i) = 1
        end if
      end do
    end do

    if (allocated(gbsa)) then
      A(:n,:n) = A(:n,:n)+gbsa%bornMat(:,:)
    end if

!     call prmat(6,A,m,m,'A eg')
    allocate (ipiv(m))

    if (single) then

      allocate (A4(m,m),x4(m))
      A4 = A
      x4 = x
      deallocate (A,x)
      call sytrf_wrap(a4,ipiv,io1)
      call sytrs_wrap(a4,x4,ipiv,io2)
      q(1:n) = x4(1:n)
      deallocate (A4,x4)

    else

      call sytrf_wrap(a,ipiv,io1)
      call sytrs_wrap(a,x,ipiv,io2)
      q(1:n) = x(1:n)
      deallocate (A,x)

    end if

    exitrun = (io1 /= 0).or.(io2 /= 0)
    if (exitRun) then
      write (stdout,'("**ERROR**",a,1x,a)') 'Solving linear equations failed',source
      return
    end if

    if (n .eq. 1) q(1) = chrg

    !  energy !

    es = 0.0_wp
    do i = 1,n

      ii = i*(i-1)/2

      do j = 1,i-1
        ij = ii+j
        tmp = eeqtmp(2,ij)
        es = es+q(i)*q(j)*tmp/r(ij)
      end do

      es = es-q(i)*(topo%chieeq(i)+param%cnf(at(i))*sqrt(cn(i))) &
      &        +q(i)*q(i)*0.5d0*(topo%gameeq(i)+tsqrt2pi/sqrt(topo%alpeeq(i)))

    end do

  end subroutine goed_gfnff

  subroutine goed_pbc_gfnff(single,n,at,xyz,r,cell,chrg,eeqtmp,cn,q,es,&
                        & gbsa,param,topo,gTrans,rTrans,x,cf)
    implicit none
    character(len=*),parameter :: source = 'gfnff_eg_goed'
    ! Molecular structure information
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    logical,intent(in)  :: single     ! real*4 flag for solver
    integer,intent(in)  :: n          ! number of atoms
    integer,intent(in)     :: at(n)   ! ordinal numbers
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(in)    :: r(n,n)  ! dist
    type(TCell),intent(in) :: cell
    real(wp),intent(in)  :: chrg       ! total charge on system
    real(wp),intent(in)  :: cn(n)      ! CN
    real(wp),intent(out) :: q(n)       ! output charges
    real(wp),intent(out) :: es         ! ES energyi
    real(wp),allocatable,intent(out) :: x(:)
    real(wp),intent(out) :: eeqtmp(2,n*(n+1)/2)    ! intermediates
    real(wp),intent(out) :: cf !convergence factor
    type(TBorn),allocatable,intent(in) :: gbsa

    ! local variables
    integer  :: m,i,j,k,ii,ij,io1,io2
    integer,allocatable :: ipiv(:)
    real(wp) :: gammij,tsqrt2pi,r2,tmp
    real(wp),allocatable :: A(:,:),x_right(:)
    real(sp),allocatable :: A4(:,:),x4(:)

    ! for calc of Coulomb matrix Amat
    real(wp),allocatable :: Amat(:,:),Amat_or(:,:)
    real(wp),allocatable,intent(out) :: rTrans(:,:)
    real(wp),allocatable,intent(out) :: gTrans(:,:)
    real(wp) :: vec(3)
    integer :: iRp,iG1,iG2,iG3,iT1,iT2,iT3
    integer,parameter :: ewaldCutD(3) = 2
    integer,parameter :: ewaldCutR(3) = 2
    real(wp) :: avgAlpeeq

    ! parameter
    parameter(tsqrt2pi=0.797884560802866_wp)
    logical :: exitRun

    m = n+topo%nfrag ! # atoms + chrg constrain + frag constrain

    allocate (Amat(m,m),Amat_or(m,m),x(m),x_right(m),source=0.0_wp) ! matrix contains constrains -> linear equations

    ! calculate rTrans, gTrans
    iRp = 0
    allocate (gTrans(3,product(2*ewaldCutR+1)-1),source=0.0_wp)
    do iG1 = -ewaldCutR(1),ewaldCutR(1)
      do iG2 = -ewaldCutR(2),ewaldCutR(2)
        do iG3 = -ewaldCutR(3),ewaldCutR(3)
          if (iG1 == 0.and.iG2 == 0.and.iG3 == 0) cycle
          iRp = iRp+1
          vec(:) = [iG1,iG2,iG3]
          gTrans(:,iRp) = matmul(cell%rec_lat,vec)
        end do
      end do
    end do

    iRp = 0
    allocate (rTrans(3,product(2*ewaldCutD+1)))
    do iT1 = -ewaldCutD(1),ewaldCutD(1)
      do iT2 = -ewaldCutD(2),ewaldCutD(2)
        do iT3 = -ewaldCutD(3),ewaldCutD(3)
          iRp = iRp+1
          vec(:) = [iT1,iT2,iT3]
          rTrans(:,iRp) = matmul(cell%lattice,vec)
        end do
      end do
    end do
!  calc eeqtmp
    !$omp parallel default(none) &
    !$omp shared(topo,n,r,eeqtmp) &
    !$omp private(i,j,k,ij,gammij,tmp)
    !$omp do schedule(dynamic)
    do i = 1,n
      k = i*(i-1)/2
      do j = 1,i-1
        ij = k+j
        gammij = 1./sqrt(topo%alpeeq(i)+topo%alpeeq(j)) ! squared above
        tmp = erf(gammij*r(i,j))
        eeqtmp(1,ij) = gammij
        eeqtmp(2,ij) = tmp
      end do
    end do
    !$omp enddo
    !$omp end parallel

    ! cf, aka ewald parameter
    avgAlpeeq = sum(topo%alpeeq)/n
    cf = get_cf(rTrans,gTrans,cell%volume,avgAlpeeq)

    ! build Ewald matrix
    call get_amat_3d(n,at,xyz,cell,topo,cf,rTrans,gTrans,Amat)

    !  setup RHS
    do i = 1,n
      x(i) = topo%chieeq(i)+param%cnf(at(i))*sqrt(cn(i))
    end do

!  fragment charge constrain
    do i = 1,topo%nfrag
      x(n+i) = topo%qfrag(i)
      do j = 1,n
        if (topo%fraglist(j) .eq. i) then
          Amat(n+i,j) = 1
          Amat(j,n+i) = 1
        end if
      end do
    end do

    if (allocated(gbsa)) then
      Amat(:n,:n) = Amat(:n,:n)+gbsa%bornMat(:,:)
    end if

    allocate (ipiv(m))

    Amat_or = Amat
    if (single) then
      allocate (A4(m,m),x4(m))
      A4 = Amat
      x4 = x
      x_right(1:n) = x(1:n)
      deallocate (Amat)
      call sytrf_wrap(a4,ipiv,io1)
      call sytrs_wrap(a4,x4,ipiv,io2)
      q(1:n) = x4(1:n)
      x = x4
      deallocate (A4,x4)
    else
      x_right(1:n) = x(1:n)
      call sytrf_wrap(Amat,ipiv,io1)
      call sytrs_wrap(Amat,x,ipiv,io2)
      q(1:n) = x(1:n)
      deallocate (Amat)
    end if

    exitRUn = (io1 /= 0).or.(io2 /= 0)
    if (exitRun) then
      write (stdout,'("**ERROR**",a,1x,a)') 'Solving linear equations failed',source
      return
    end if

    if (n .eq. 1) q(1) = chrg

    ! calculate E_es = q^T*(0.5*A*q-X) |  ^T is the transpose
    ! First calc bracket term with symv, result is saved in x_right
    ! now the x is used as q to save memory
    call symv(Amat_or,x,x_right,alpha=0.5_wp,beta=-1.0_wp)
    ! about symv similar to gemv
    ! y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y

    !Calc ES as dot product (=matrix product) since q and x are vectors
    !es = dot_product(x, x_right)
    es = dot(x,x_right)!*1.05_wp
  end subroutine goed_pbc_gfnff

  function get_cf(rTrans,gTrans,vol,avgAlp) result(cf)
    ! output the ewald splitting parameter
    real(wp) :: cf
    ! parameter from goed_pbc_gfnff ewaldCutD and ewaldCutR
    integer,parameter :: ewaldCutD(3) = 2
    integer,parameter :: ewaldCutR(3) = 2
    ! real space lattice vectors
    real(wp),intent(in) :: rTrans(:,:)
    ! reciprocal space lattice vectors
    real(wp),intent(in) :: gTrans(:,:)
    ! unit cell volume
    real(wp),intent(in) :: vol
    ! average alphaEEQ value
    real(wp),intent(in) :: avgAlp
    ! smallest reciprocal and real space Vectors
    real(wp) :: minG,minR
    ! approx Value of reciprocal and real part electrostatics
    real(wp) :: gPart,rPart
    !
    real(wp) :: gam
    ! current cf
    real(wp) :: cfCurr
    !
    real(wp) :: lenR,minTmp,gr_diff,diffMin
    real(wp) :: a(100)
    ! tolerance for golden-section search
    real(wp),parameter :: tol = 1.0e-4_wp
    ! golden ratio
    real(wp),parameter :: goldr = (1.0_wp+sqrt(2.0_wp))/2.0_wp
    ! evaluation points for gss
    real(wp) :: x1,x2,x3,x4
    !
    integer :: i,j,iter

    cf = 0.14999_wp ! default
    gam = 1.0_wp/sqrt(2*avgAlp)

    ! get smallest real and reciprocal vector
    minG = sqrt(minval(sum(gTrans(:,:)**2,dim=1)))
    minR = huge(1.0_wp)
    do i = 1,size(rTrans,dim=2)
      lenR = sqrt(sum(rTrans(:,i)**2))
      if (lenR .ne. 0.0_wp) then
        minR = min(minR,lenR)
      end if
    end do

    ! golden-section search algorithm for convergence factor
    iter = 0
    x1 = 1.0e-8_wp  ! left margin
    x4 = 2.0e+0_wp  ! right margin
    x2 = x4-(x4-x1)/goldr
    x3 = x1+(x4-x1)/goldr
    do while ((x4-x1) > tol)
      iter = iter+1
      if (grFct(x2,minG,minR,vol,gam) .lt. grFct(x3,minG,minR,vol,gam)) then
        x4 = x3
      else
        x1 = x2
      end if
      ! recalculate x2 and x3
      x2 = x4-(x4-x1)/goldr
      x3 = x1+(x4-x1)/goldr
    end do

    cf = (x1+x4)/2.0_wp
  end function get_cf

  function grFct(cfCurr,minG,minR,vol,gam) result(gr_diff)
    real(wp) :: gr_diff
    ! smallest reciprocal and real space Vectors
    real(wp),intent(in) :: cfCurr,minG,minR,vol,gam
    ! approx Value of reciprocal and real part electrostatics
    real(wp) :: gPart,rPart

    gPart = 4.0_wp*pi*exp(-minG**2/(4.0_wp*cfCurr**2))/(vol*minG**2)
    rPart = -erf(cfCurr*minR)/minR+erf(gam*minR)/minR
    gr_diff = abs(gPart-rPart)
  end function grFct

  subroutine es_grad_sigma(nat,at,xyz,cell,topo,nlist,rTrans,gTrans,xtmp,cf, &
           & sigma,gradient,mcf_ees)
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TCell),intent(in) :: cell
    type(TGFFTopology),intent(in) :: topo
    type(TGFFNeighbourList),intent(in) :: nlist
    real(wp),intent(in) :: cf
    !real(wp), intent(in) :: qvec(:)
    real(wp),intent(in) :: gTrans(:,:)
    real(wp),intent(in) :: rTrans(:,:)
    real(wp),intent(in) :: xtmp(nat+topo%nfrag)
    real(wp),intent(inout) :: sigma(3,3)
    real(wp),intent(inout) :: gradient(3,nat)
    real(wp),intent(in) :: mcf_ees
    real(wp),allocatable :: dXvecdr(:,:,:)
    real(wp),allocatable :: amatdr(:,:,:)
    real(wp),allocatable :: amatdL(:,:,:)
    real(wp),allocatable :: atrace(:,:)
    allocate (amatdr(3,nat,nat),amatdL(3,3,nat),source=0.0_wp)
    allocate (atrace(3,nat),source=0.0_wp)

    ! new routine from D4 for calculating derivatives
    call get_damat_3d(nat,at,xyz,cell,topo,cf,xtmp,rTrans,gTrans,amatdr,amatdL,atrace)

    ! y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y
    call gemv(amatdr,nlist%q,gradient,alpha=mcf_ees,beta=1.0_wp)

    ! Ees=q^T*(0.5*A*q-X)  here is the 0.5qA'q part. The -qX' part is below the es_grad_sigma call
    call gemv(amatdL,nlist%q,sigma,alpha=0.5_wp*mcf_ees,beta=1.0_wp)

  end subroutine es_grad_sigma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               code below taken from dftd4                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_amat_3d(nat,at,xyz,cell,topo,alpha,rTrans,gTrans,amat)
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TCell),intent(in) :: cell
    type(TGFFTopology),intent(in) :: topo
    real(wp),intent(in) :: alpha
    real(wp),intent(in) :: rTrans(:,:)
    real(wp),intent(in) :: gTrans(:,:)
    real(wp),intent(out) :: amat(:,:)

    integer :: iat,jat,izp,jzp,img
    real(wp) :: vec(3),gam,wsw,dtmp,rtmp,vol
    real(wp),parameter :: zero(3) = 0.0_wp
    real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
    real(wp),parameter :: sqrtpi = 1.772453850905516_wp

    amat(:,:) = 0.0_wp

    vol = cell%volume ! abs(matdet_3x3(cell%lattice))

    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:amat) shared(nat,at,xyz,cell, topo, rTrans, gTrans, alpha, vol) &
    !$omp private(iat, jat, gam, wsw, vec, dtmp, rtmp)
    do iat = 1,nat
      do jat = 1,iat-1
        gam = 1.0_wp/sqrt(topo%alpeeq(iat)+topo%alpeeq(jat))
        wsw = cell%wsc%w(jat,iat)
        do img = 1,cell%wsc%itbl(jat,iat)
          vec = xyz(:,iat)-xyz(:,jat) &
             & -(cell%lattice(:,1)*cell%wsc%lattr(1,img,jat,iat) &
             &  +cell%lattice(:,2)*cell%wsc%lattr(2,img,jat,iat) &
             &  +cell%lattice(:,3)*cell%wsc%lattr(3,img,jat,iat))
          call get_amat_dir_3d(vec,gam,alpha,rTrans,dtmp)
          call get_amat_rec_3d(vec,vol,alpha,gTrans,rtmp)
          amat(jat,iat) = amat(jat,iat)+(dtmp+rtmp)*wsw
          amat(iat,jat) = amat(iat,jat)+(dtmp+rtmp)*wsw
        end do
      end do

      gam = 1.0_wp/sqrt(2.0_wp*topo%alpeeq(iat))
      wsw = cell%wsc%w(iat,iat)
      do img = 1,cell%wsc%itbl(iat,iat)
        vec = zero

        call get_amat_dir_3d(vec,gam,alpha,rTrans,dtmp)
        call get_amat_rec_3d(vec,vol,alpha,gTrans,rtmp)
        amat(iat,iat) = amat(iat,iat)+(dtmp+rtmp)*wsw
      end do

      dtmp = topo%gameeq(iat)+sqrt2pi/sqrt(topo%alpeeq(iat))-2*alpha/sqrtpi
      amat(iat,iat) = amat(iat,iat)+dtmp
    end do
    !$omp end parallel do

    amat(nat+1,1:nat+1) = 1.0_wp
    amat(1:nat+1,nat+1) = 1.0_wp
    amat(nat+1,nat+1) = 0.0_wp

  end subroutine get_amat_3d

  subroutine get_amat_dir_3d(rij,gam,alp,trans,amat)
    real(wp),intent(in) :: rij(3)
    real(wp),intent(in) :: gam
    real(wp),intent(in) :: alp
    real(wp),intent(in) :: trans(:,:)
    real(wp),intent(out) :: amat

    integer :: itr
    real(wp) :: vec(3),r1,tmp
    real(wp),parameter :: eps = 1.0e-9_wp

    amat = 0.0_wp

    do itr = 1,size(trans,2)
      vec(:) = rij+trans(:,itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      tmp = erf(gam*r1)/r1-erf(alp*r1)/r1
      amat = amat+tmp
    end do

  end subroutine get_amat_dir_3d

  subroutine get_amat_rec_3d(rij,vol,alp,trans,amat)
    real(wp),intent(in) :: rij(3)
    real(wp),intent(in) :: vol
    real(wp),intent(in) :: alp
    real(wp),intent(in) :: trans(:,:)
    real(wp),intent(out) :: amat

    integer :: itr
    real(wp) :: fac,vec(3),g2,tmp
    real(wp),parameter :: eps = 1.0e-9_wp

    amat = 0.0_wp
    fac = 4*pi/vol

    do itr = 1,size(trans,2)
      vec(:) = trans(:,itr)
      g2 = dot_product(vec,vec)
      if (g2 < eps) cycle
      tmp = cos(dot_product(rij,vec))*fac*exp(-0.25_wp*g2/(alp*alp))/g2
      amat = amat+tmp
    end do

  end subroutine get_amat_rec_3d

  subroutine get_damat_3d(nat,at,xyz,cell,topo,alpha,qvec,rTrans,gTrans,dadr,dadL,atrace)
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TCell),intent(in) :: cell
    type(TGFFTopology),intent(in) :: topo
    real(wp),intent(in) :: alpha
    real(wp),intent(in) :: qvec(:)
    real(wp),intent(in) :: rTrans(:,:)
    real(wp),intent(in) :: gTrans(:,:)
    real(wp),intent(out) :: dadr(:,:,:)
    real(wp),intent(out) :: dadL(:,:,:)
    real(wp),intent(out) :: atrace(:,:)

    integer :: iat,jat,izp,jzp,img
    real(wp) :: vol,gam,wsw,vec(3),dG(3),dS(3,3)
    real(wp) :: dGd(3),dSd(3,3),dGr(3),dSr(3,3)
    real(wp),parameter :: zero(3) = 0.0_wp

    atrace(:,:) = 0.0_wp
    dadr(:,:,:) = 0.0_wp
    dadL(:,:,:) = 0.0_wp

    vol = cell%volume ! abs(matdet_3x3(cell%lattice))

    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:atrace, dadr, dadL) &
    !$omp shared(nat,at,xyz,cell, topo, alpha, vol, rTrans, gTrans, qvec) &
    !$omp private(iat, jat, img, gam, wsw, vec, dG, dS, &
    !$omp& dGr, dSr, dGd, dSd)
    do iat = 1,nat
      do jat = 1,iat-1
        dG(:) = 0.0_wp
        dS(:,:) = 0.0_wp
        gam = 1.0_wp/sqrt(topo%alpeeq(iat)+topo%alpeeq(jat))
        wsw = cell%wsc%w(jat,iat)
        do img = 1,cell%wsc%itbl(jat,iat)
          vec = xyz(:,iat)-xyz(:,jat) &
             & -(cell%lattice(:,1)*cell%wsc%lattr(1,img,jat,iat) &
             &  +cell%lattice(:,2)*cell%wsc%lattr(2,img,jat,iat) &
             &  +cell%lattice(:,3)*cell%wsc%lattr(3,img,jat,iat))
          call get_damat_dir_3d(vec,gam,alpha,rTrans,dGd,dSd)
          call get_damat_rec_3d(vec,vol,alpha,gTrans,dGr,dSr)
          dG = dG+(dGd+dGr)*wsw
          dS = dS+(dSd+dSr)*wsw
        end do
        atrace(:,iat) = +dG*qvec(jat)+atrace(:,iat)
        atrace(:,jat) = -dG*qvec(iat)+atrace(:,jat)
        dadr(:,iat,jat) = +dG*qvec(iat)+dadr(:,iat,jat)
        dadr(:,jat,iat) = -dG*qvec(jat)+dadr(:,jat,iat)
        dadL(:,:,jat) = +dS*qvec(iat)+dadL(:,:,jat)
        dadL(:,:,iat) = +dS*qvec(jat)+dadL(:,:,iat)
      end do

      dS(:,:) = 0.0_wp
      gam = 1.0_wp/sqrt(2.0_wp*topo%alpeeq(iat))
      wsw = cell%wsc%w(iat,iat)
      do img = 1,cell%wsc%itbl(iat,iat)
        vec = zero
        call get_damat_dir_3d(vec,gam,alpha,rTrans,dGd,dSd)
        call get_damat_rec_3d(vec,vol,alpha,gTrans,dGr,dSr)
        dS = dS+(dSd+dSr)*wsw
      end do
      dadL(:,:,iat) = +dS*qvec(iat)+dadL(:,:,iat)
    end do
    !$omp end parallel do

  end subroutine get_damat_3d

  subroutine get_damat_dir_3d(rij,gam,alp,trans,dg,ds)
    real(wp),intent(in) :: rij(3)
    real(wp),intent(in) :: gam
    real(wp),intent(in) :: alp
    real(wp),intent(in) :: trans(:,:)
    real(wp),intent(out) :: dg(3)
    real(wp),intent(out) :: ds(3,3)

    integer :: itr
    real(wp) :: vec(3),r1,r2,gtmp,atmp,gam2,alp2
    real(wp),parameter :: eps = 1.0e-9_wp
    real(wp),parameter :: sqrtpi = 1.772453850905516_wp

    dg(:) = 0.0_wp
    ds(:,:) = 0.0_wp

    gam2 = gam*gam
    alp2 = alp*alp

    do itr = 1,size(trans,2)
      vec(:) = rij+trans(:,itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1*r1
      gtmp = +2*gam*exp(-r2*gam2)/(sqrtpi*r2)-erf(r1*gam)/(r2*r1)
      atmp = -2*alp*exp(-r2*alp2)/(sqrtpi*r2)+erf(r1*alp)/(r2*r1)
      dg(:) = dg+(gtmp+atmp)*vec
      ds(:,1) = ds(:,1)+(gtmp+atmp)*vec(1)*vec
      ds(:,2) = ds(:,2)+(gtmp+atmp)*vec(2)*vec
      ds(:,3) = ds(:,3)+(gtmp+atmp)*vec(3)*vec
    end do

  end subroutine get_damat_dir_3d

  subroutine get_damat_rec_3d(rij,vol,alp,trans,dg,ds)
    real(wp),intent(in) :: rij(3)
    real(wp),intent(in) :: vol
    real(wp),intent(in) :: alp
    real(wp),intent(in) :: trans(:,:)
    real(wp),intent(out) :: dg(3)
    real(wp),intent(out) :: ds(3,3)

    integer :: itr
    real(wp) :: fac,vec(3),g2,gv,etmp,dtmp,alp2
    real(wp),parameter :: unity(3,3) = reshape(&
       & [1,0,0,0,1,0,0,0,1],shape(unity))
    real(wp),parameter :: eps = 1.0e-9_wp

    dg(:) = 0.0_wp
    ds(:,:) = 0.0_wp
    fac = 4*pi/vol
    alp2 = alp*alp

    do itr = 1,size(trans,2)
      vec(:) = trans(:,itr)
      g2 = dot_product(vec,vec)
      if (g2 < eps) cycle
      gv = dot_product(rij,vec)
      etmp = fac*exp(-0.25_wp*g2/alp2)/g2
      dtmp = -sin(gv)*etmp
      dg(:) = dg+dtmp*vec
      ds(:,1) = ds(:,1)+etmp*cos(gv) &
               & *((2.0_wp/g2+0.5_wp/alp2)*vec(1)*vec-unity(:,1))
      ds(:,2) = ds(:,2)+etmp*cos(gv) &
               & *((2.0_wp/g2+0.5_wp/alp2)*vec(2)*vec-unity(:,2))
      ds(:,3) = ds(:,3)+etmp*cos(gv) &
               & *((2.0_wp/g2+0.5_wp/alp2)*vec(3)*vec-unity(:,3))
    end do

  end subroutine get_damat_rec_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               code above taken from dftd4                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!cccccccccccccccccccccccccccccccccccc
! HB energy and analytical gradient
!cccccccccccccccccccccccccccccccccccc

!> Case 1: A...H...B
  subroutine abhgfnff_eg1(n,A,B,H,iTrA,iTrB,at,xyz,q,energy,gdr,param,topo,neigh,sigma,mcf_ehb)
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
    real(wp) xm,ym,zm
    real(wp) rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
    real(wp) drah(3),drbh(3),drab(3),drm(3)
    real(wp) dg(3),dga(3),dgb(3),dgh(3)
    real(wp) ga(3),gb(3),gh(3)
    real(wp) gi,denom,ratio,tmp,qhoutl,radab,rahprbh
    real(wp) ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo
    real(wp) bas,aci
    real(wp) eabh
    real(wp) aterm,rterm,dterm,sterm
    real(wp) qa,qb,qh
    real(wp) ca(2),cb(2)
    real(wp) gqa,gqb,gqh
    real(wp) caa,cbb
    real(wp) shortcut

    integer i,j!,iTrDum ! ij,lina

    gdr = 0
    energy = 0
    call hbonds(A,B,ca,cb,param,topo) ! get HB strength
!     A-B distance
    rab = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,B)+neigh%transVec(:,iTrB)))
    rab2 = rab**2
!     A-H distance
    rah = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,H)))
    rah2 = rah**2
!     B-H distance
    rbh = NORM2((xyz(:,B)+neigh%transVec(:,iTrB))-xyz(:,H))
    rbh2 = rbh**2

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))

!     out-of-line damp
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

!     long damping
    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

!     short damping
    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    rdamp = damp/rab2/rab

!     hydrogen charge scaled term
    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

!     hydrogen charge scaled term
    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

!     hydrogen charge scaled term
    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

!     donor-acceptor term
    rah4 = rah2*rah2
    rbh4 = rbh2*rbh2
    denom = 1.d0/(rah4+rbh4)

    caa = qa*ca(1)
    cbb = qb*cb(1)
    qhoutl = qh*outl

    bas = (caa*rah4+cbb*rbh4)*denom
    aci = (cb(2)*rah4+ca(2)*rbh4)*denom

!     energy
    rterm = -aci*rdamp*qhoutl
    energy = bas*rterm

!     gradient
    drah(1:3) = xyz(1:3,A)-xyz(1:3,H)+neigh%transVec(1:3,iTrA)
    drbh(1:3) = xyz(1:3,B)-xyz(1:3,H)+neigh%transVec(1:3,iTrB)
    drab(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))

    aterm = -aci*bas*rdamp*qh
    sterm = -rdamp*bas*qhoutl
    dterm = -aci*bas*qhoutl

    tmp = denom*denom*4.0d0
    dd24a = rah2*rbh4*tmp
    dd24b = rbh2*rah4*tmp

!     donor-acceptor part: bas
    gi = (caa-cbb)*dd24a*rterm
    ga(1:3) = gi*drah(1:3)
    gi = (cbb-caa)*dd24b*rterm
    gb(1:3) = gi*drbh(1:3)
    gh(1:3) = -ga(1:3)-gb(1:3)

!     donor-acceptor part: aci
    gi = (cb(2)-ca(2))*dd24a
    dga(1:3) = gi*drah(1:3)*sterm
    ga(1:3) = ga(1:3)+dga(1:3)

    gi = (ca(2)-cb(2))*dd24b
    dgb(1:3) = gi*drbh(1:3)*sterm
    gb(1:3) = gb(1:3)+dgb(1:3)

    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

!     damping part rab
    gi = rdamp*(-(2.d0*param%hbalp*ratio1/(1+ratio1))+(2.d0*param%hbalp*ratio3/(1+ratio3))-3.d0)/rab2
    dg(1:3) = gi*drab(1:3)*dterm
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!     out of line term: rab
    gi = aterm*2.d0*ratio2*expo*rahprbh/(1+ratio2)**2/(rahprbh-rab)/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!     out of line term: rah,rbh
    tmp = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    dga(1:3) = drah(1:3)*tmp/rah
    ga(1:3) = ga(1:3)+dga(1:3)
    dgb(1:3) = drbh(1:3)*tmp/rbh
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)
    ! sigma
    sigma(:,1) = sigma(:,1)+mcf_ehb*ga(1)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,2) = sigma(:,2)+mcf_ehb*ga(2)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,3) = sigma(:,3)+mcf_ehb*ga(3)*(xyz(:,A)+neigh%transVec(1:3,iTrA))
    sigma(:,1) = sigma(:,1)+mcf_ehb*gb(1)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,2) = sigma(:,2)+mcf_ehb*gb(2)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,3) = sigma(:,3)+mcf_ehb*gb(3)*(xyz(:,B)+neigh%transVec(1:3,iTrB))
    sigma(:,1) = sigma(:,1)+mcf_ehb*gh(1)*xyz(:,H)
    sigma(:,2) = sigma(:,2)+mcf_ehb*gh(2)*xyz(:,H)
    sigma(:,3) = sigma(:,3)+mcf_ehb*gh(3)*xyz(:,H)
!     move gradients into place
    gdr(1:3,1) = ga(1:3)
    gdr(1:3,2) = gb(1:3)
    gdr(1:3,3) = gh(1:3)

  end subroutine abhgfnff_eg1

!> Case 2: A-H...B including orientation of neighbors at B
  subroutine abhgfnff_eg2new(n,A,B,H,iTrA,iTrB,nbb,at,xyz,q,sqrab, &
                  & srab,energy,gdr,param,topo,neigh,sigma,mcf_ehb)
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
    real(wp) xm,ym,zm
    real(wp) rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
    real(wp) ranb(nbb),ranb2(nbb),rbnb(nbb),rbnb2(nbb)
    real(wp) drah(3),drbh(3),drab(3),drm(3)
    real(wp) dranb(3,nbb),drbnb(3,nbb)
    real(wp) dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
    real(wp) ga(3),gb(3),gh(3),gnb(3,nbb)
    real(wp) denom,ratio,qhoutl,radab
    real(wp) gi,gi_nb(nbb)
    real(wp) tmp1,tmp2(nbb)
    real(wp) rahprbh,ranbprbnb(nbb)
    real(wp) ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb(nbb)
    real(wp) eabh
    real(wp) aterm,dterm,nbterm
    real(wp) qa,qb,qh
    real(wp) ca(2),cb(2)
    real(wp) gqa,gqb,gqh
    real(wp) shortcut
    real(wp) const
    real(wp) outl_nb(nbb),outl_nb_tot
    real(wp) hbnbcut_save
    real(wp) vecDum(3)
    logical mask_nb(nbb)

!     proportion between Rbh und Rab distance dependencies
    real(wp) :: p_bh
    real(wp) :: p_ab

    integer i,j,inb,iTr!,iTrDum!,ij,lina
    p_bh = 1.d0+param%hbabmix
    p_ab = -param%hbabmix

    gdr = 0
    energy = 0

    call hbonds(A,B,ca,cb,param,topo)
!     Neighbours of B
    do i = 1,nbb
      inb = 0; iTr = 0 ! jth_nb output
      call neigh%jth_nb(n,xyz,inb,i,B,iTr) ! inb is the i-th nb of B when inb is shifted to iTr
!        compute distances
      vecDum = neigh%transVec(:,iTr)+neigh%transVec(:,iTrB)
      dranb(1:3,i) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,inb)+vecDum)
      drbnb(1:3,i) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-(xyz(1:3,inb)+vecDum)
!        A-nb(B) distance
      ranb2(i) = sum(dranb(1:3,i)**2)
      ranb(i) = sqrt(ranb2(i))
!        B-nb(B) distance
      rbnb2(i) = sum(drbnb(1:3,i)**2)
      rbnb(i) = sqrt(rbnb2(i))
    end do

!     A-B distance
    rab = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,B)+neigh%transVec(:,iTrB)))
    rab2 = rab**2
!     A-H distance
    rah = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,H)))
    rah2 = rah**2
!     B-H distance
    rbh = NORM2((xyz(:,B)+neigh%transVec(:,iTrB))-xyz(:,H))
    rbh2 = rbh**2

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))

!     out-of-line damp: A-H...B
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

!     out-of-line damp: A...nb(B)-B
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

!     long damping
    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

!     short damping
    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
    rbhdamp = damp*((p_bh/rbh2/rbh))
    rabdamp = damp*((p_ab/rab2/rab))
    rdamp = rbhdamp+rabdamp

!     hydrogen charge scaled term
    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

!     hydrogen charge scaled term
    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

!     hydrogen charge scaled term
    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

    qhoutl = qh*outl*outl_nb_tot

!     constant values, no gradient
    const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh
!     energy
    energy = -rdamp*qhoutl*const
!     gradient
    drah(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-xyz(1:3,H)
    drbh(1:3) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-xyz(1:3,H)
    drab(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))

    aterm = -rdamp*qh*outl_nb_tot*const
    nbterm = -rdamp*qh*outl*const
    dterm = -qhoutl*const

!------------------------------------------------------------------------------
!     damping part: rab
    gi = ((rabdamp+rbhdamp)*ddamp-3.d0*rabdamp)/rab2
    gi = gi*dterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = dg(1:3)
    gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     damping part: rbh
    gi = -3.d0*rbhdamp/rbh2
    gi = gi*dterm
    dg(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dg(1:3)
    gh(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     angular A-H...B term
!------------------------------------------------------------------------------
!     out of line term: rab
    tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    gi = -tmp1*rahprbh/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!     out of line term: rah,rbh
    gi = tmp1/rah
    dga(1:3) = gi*drah(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp1/rbh
    dgb(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

!------------------------------------------------------------------------------
!     angular A...nb(B)-B term
!------------------------------------------------------------------------------
!     out of line term: rab
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

!     out of line term: ranb,rbnb
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

!------------------------------------------------------------------------------
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

!------------------------------------------------------------------------------
!     move gradients into place
    gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
    gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
    gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
    do i = 1,nbb
      inb = 0; iTr = 0 ! jth_nb output
      call neigh%jth_nb(n,xyz,inb,i,B,iTr) ! inb is the i-th nb of B when shifted to iTr
      gdr(1:3,inb) = gdr(1:3,inb)+gnb(1:3,i)
    end do

    ! sigma according to gdr above
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

!> Case 2: A-H...B including LP position
  subroutine abhgfnff_eg2_rnr(n,A,B,H,iTrA,iTrB,at,xyz,q,sqrab,srab,energy,gdr,param,topo,neigh,sigma,mcf_ehb)
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
    real(wp) xm,ym,zm
    real(wp) rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
    real(wp) ranb(2),ranb2(2),rbnb(2),rbnb2(2)
    real(wp) drah(3),drbh(3),drab(3),drm(3),dralp(3),drblp(3)
    real(wp) dranb(3,2),drbnb(3,2)
    real(wp) dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
    real(wp) ga(3),gb(3),gh(3),gnb(3,2),gnb_lp(3),glp(3)
    real(wp) denom,ratio,qhoutl,radab
    real(wp) gi,gi_nb(2)
    real(wp) tmp1,tmp2(2),tmp3
    real(wp) rahprbh,ranbprbnb(2)
    real(wp) ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_lp,expo_nb(2)
    real(wp) eabh
    real(wp) aterm,dterm,nbterm,lpterm
    real(wp) qa,qb,qh
    real(wp) ca(2),cb(2)
    real(wp) gqa,gqb,gqh
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

!     proportion between Rbh und Rab distance dependencies
    real(wp) :: p_bh
    real(wp) :: p_ab
!     lone-pair out-of-line damping
    real(wp) hblpcut

    real(wp) :: vTrinb(3)
    integer i,j,nbb,inb,iTr!,iTrDum

    p_bh = 1.d0+param%hbabmix
    p_ab = -param%hbabmix

    gdr = 0
    energy = 0
    vector = 0
    lp_dist = 0.50-0.018*param%repz(at(B))
    hblpcut = 56

    call hbonds(A,B,ca,cb,param,topo)

    nbb = 2 ! given through if condition before call
!     Neighbours of B
    do i = 1,nbb
      inb = 0; iTr = 0 ! jth_nb output
      call neigh%jth_nb(n,xyz,inb,i,B,iTr) ! inb is the i-th nb of B when shifted to iTr
!        compute distances
      vTrinb = neigh%transVec(:,iTr)+neigh%transVec(:,iTrB)
      dranb(1:3,i) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,inb)+vTrinb)
      drbnb(1:3,i) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-(xyz(1:3,inb)+vTrinb)
!        A-nb(B) distance
      ranb2(i) = sum(dranb(1:3,i)**2)
      ranb(i) = sqrt(ranb2(i))
!        B-nb(B) distance
      rbnb2(i) = sum(drbnb(1:3,i)**2)
      rbnb(i) = sqrt(rbnb2(i))

      drnb(1:3,i) = (xyz(1:3,inb)+vTrinb)-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))
      vector = vector+drnb(1:3,i)
    end do

    vnorm = norm2(vector)
!     lonepair coordinates
    if (vnorm .gt. 1.d-10) then
      lp = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-lp_dist*(vector/vnorm)
    else
      lp = xyz(1:3,B)+neigh%transVec(1:3,iTrB)
      nbb = 0
    end if

!     A-B distance
    rab = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,B)+neigh%transVec(:,iTrB)))
    rab2 = rab**2
!     A-H distance
    rah = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,H)))
    rah2 = rah**2
!     B-H distance
    rbh = NORM2((xyz(:,B)+neigh%transVec(:,iTrB))-xyz(:,H))
    rbh2 = rbh**2

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))
!     out-of-line damp: A-H...B
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

!     out-of-line damp: A...LP-B
    rblp2 = sum(((xyz(1:3,B)+neigh%transVec(1:3,iTrB))-lp(1:3))**2)
    rblp = sqrt(rblp2)
    ralp2 = sum(((xyz(1:3,A)+neigh%transVec(1:3,iTrA))-lp(1:3))**2)
    ralp = sqrt(ralp2)
    ralpprblp = ralp+rblp+1.d-12
    expo_lp = (hblpcut/radab)*(ralpprblp/rab-1.d0)
    ratio2_lp = exp(expo_lp)
    outl_lp = 2.d0/(1.d0+ratio2_lp)

!     out-of-line damp: A...nb(B)-B
    do i = 1,nbb
      ranbprbnb(i) = ranb(i)+rbnb(i)+1.d-12
      expo_nb(i) = (param%hbnbcut/radab)*(ranbprbnb(i)/rab-1.d0)
      ratio2_nb(i) = exp(-expo_nb(i))**(1.0)
      outl_nb(i) = (2.d0/(1.d0+ratio2_nb(i)))-1.0d0
    end do
    outl_nb_tot = product(outl_nb)

!     long damping
    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

!     short damping
    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
    rbhdamp = damp*((p_bh/rbh2/rbh))
    rabdamp = damp*((p_ab/rab2/rab))
    rdamp = rbhdamp+rabdamp

!     hydrogen charge scaled term
    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

!     hydrogen charge scaled term
    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

!     hydrogen charge scaled term
    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

    qhoutl = qh*outl*outl_nb_tot*outl_lp

!     constant values, no gradient
    const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

!     energy
    energy = -rdamp*qhoutl*const
!     gradient
    drah(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-xyz(1:3,H)
    drbh(1:3) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-xyz(1:3,H)
    drab(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))
    dralp(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-lp(1:3)
    drblp(1:3) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-lp(1:3)

    aterm = -rdamp*qh*outl_nb_tot*outl_lp*const
    nbterm = -rdamp*qh*outl*outl_lp*const
    lpterm = -rdamp*qh*outl*outl_nb_tot*const
    dterm = -qhoutl*const

!------------------------------------------------------------------------------
!     damping part: rab
    gi = ((rabdamp+rbhdamp)*ddamp-3.d0*rabdamp)/rab2
    gi = gi*dterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = dg(1:3)
    gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     damping part: rbh
    gi = -3.d0*rbhdamp/rbh2
    gi = gi*dterm
    dg(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dg(1:3)
    gh(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     angular A-H...B term
!------------------------------------------------------------------------------
!     out of line term: rab
    tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    gi = -tmp1*rahprbh/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!     out of line term: rah,rbh
    gi = tmp1/rah
    dga(1:3) = gi*drah(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp1/rbh
    dgb(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

!!------------------------------------------------------------------------------
!!     angular A...LP-B term
!!------------------------------------------------------------------------------
!     out of line term: rab
    tmp3 = -2.d0*lpterm*ratio2_lp*expo_lp/(1+ratio2_lp)**2/(ralpprblp-rab)
    gi = -tmp3*ralpprblp/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!     out of line term: ralp,rblp
    gi = tmp3/ralp
    dga(1:3) = gi*dralp(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp3/(rblp+1.0d-12)
    dgb(1:3) = gi*drblp(1:3)
    gb(1:3) = gb(1:3)-dga(1:3)
    glp(1:3) = -dga(1:3)!-dgb(1:3)

!     neighbor part: LP
    unit_vec = 0
    do i = 1,3
      unit_vec(i) = -1
      gii(1:3,i) = -lp_dist*dble(nbb)*(unit_vec/vnorm+(vector*vector(i)/sum(vector**2)**(1.5d0)))
      unit_vec = 0
    end do
    gnb_lp = matmul(gii,glp)

!------------------------------------------------------------------------------
!     angular A...nb(B)-B term
!------------------------------------------------------------------------------
!     out of line term: rab
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

!     out of line term: ranb,rbnb
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

!------------------------------------------------------------------------------
!     move gradients into place
    gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
    gdr(1:3,B) = gdr(1:3,B)+gb(1:3)+gnb_lp(1:3)
    gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
    do i = 1,nbb
      inb = 0; iTr = 0 ! jth_nb output
      call neigh%jth_nb(n,xyz,inb,i,B,iTr) ! inb is the i-th nb of B when shifted to iTr
      gdr(1:3,inb) = gdr(1:3,inb)+gnb(1:3,i)-gnb_lp(1:3)/dble(nbb)
    end do

    ! sigma according to gdr above
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

!> case 3: A-H...B, B is 0=C including two in plane LPs at B
!> this is the multiplicative version of incorporationg etors and ebend
!> equal to abhgfnff_eg2_new multiplied by etors and eangl
  subroutine abhgfnff_eg3(n,A,B,H,iTrA,iTrB,C,iTrC,at,xyz,q,sqrab,srab,energy,&
                  & gdr,param,topo,neigh,sigma,exitRun,mcf_ehb)
    implicit none
    type(TGFFData),intent(in)     :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout)    :: neigh
    real(wp),intent(inout)        :: sigma(3,3)
    logical,intent(inout)         :: exitRun
    real(wp),intent(in)           :: mcf_ehb
    integer  :: A,B,H,iTrA,iTrB,n,at(n),C,iTrC
    real(wp) :: xyz(3,n),energy,gdr(3,n)
    real(wp) :: q(n)
    real(wp) :: sqrab(n*(n+1)/2)   ! squared dist
    real(wp) :: srab(n*(n+1)/2)    ! dist

    real(wp) :: outl,dampl,damps,rdamp,damp
    real(wp) :: ddamp,rabdamp,rbhdamp
    real(wp) :: ratio1,ratio2,ratio2_nb,ratio3
    real(wp) :: xm,ym,zm
    real(wp) :: rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
    real(wp) :: vTrR(3),vTrB(3),vTrC(3)
    real(wp) :: ranb,ranb2,rbnb,rbnb2
    real(wp) :: drah(3),drbh(3),drab(3),drm(3)
    real(wp) :: dranb(3),drbnb(3)
    real(wp) :: dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
    real(wp) :: ga(3),gb(3),gh(3),gnb(3)
    real(wp) :: phi,phi0,r0,t0,fc,tshift,bshift
    real(wp) :: eangl,etors,gangl(3,n),gtors(3,n)
    real(wp) :: etmp(20),g3tmp(3,3),g4tmp(3,4,20)
    real(wp) :: ratio,qhoutl,radab
    real(wp) :: gi,gi_nb
    real(wp) :: tmp1,tmp2
    real(wp) :: rahprbh,ranbprbnb
    real(wp) :: ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb
    real(wp) :: eabh
    real(wp) :: aterm,dterm,nbterm,bterm,tterm
    real(wp) :: qa,qb,qh
    real(wp) :: ca(2),cb(2)
    real(wp) :: gqa,gqb,gqh
    real(wp) :: shortcut
    integer :: tlist(6,sum(neigh%nb(neigh%numnb,C,:)))
    real(wp) :: vtors(2,sum(neigh%nb(neigh%numnb,C,:)))
    real(wp) :: const
    real(wp) :: outl_nb,outl_nb_tot
    logical :: mask_nb,t_mask(20)

!     proportion between Rbh und Rab distance dependencies
    real(wp) :: p_bh
    real(wp) :: p_ab

    integer :: D
    integer :: i,j,ii,jj,kk,ll,ij,lina,iTr,iTrDum,iTrR
    integer :: nbb,nbc
    integer :: ntors,rn
    real(wp) :: delr

    p_bh = 1.d0+param%hbabmix
    p_ab = -param%hbabmix

    gdr = 0
    energy = 0
    etors = 0
    gtors = 0
    eangl = 0
    gangl = 0
    call hbonds(A,B,ca,cb,param,topo)

    !Determine all neighbors for torsion term
    !  A
    !   \         tors:
    !    H        ll
    !     :        \
    !      O        jj
    !      ||       |
    !      C        kk
    !     / \       \
    !    R1  R2      ii
    !------------------------------------------

    nbb = 1 ! routine only called for nbb.eq.1
    nbc = sum(neigh%nb(neigh%numnb,C,:))
    ntors = nbc-nbb

!     Neighbours of B
!     compute distances
    dranb(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,C)+neigh%transVec(1:3,iTrC))
    drbnb(1:3) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-(xyz(1:3,C)+neigh%transVec(1:3,iTrC))
!     A-nb(B) distance
    ranb2 = sum(dranb(1:3)**2)
    ranb = sqrt(ranb2)
!     B-nb(B) distance
    rbnb2 = sum(drbnb(1:3)**2)
    rbnb = sqrt(rbnb2)
    rab = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,B)+neigh%transVec(:,iTrB)))
    rab2 = rab**2
!     A-H distance
    rah = NORM2((xyz(:,A)+neigh%transVec(:,iTrA))-(xyz(:,H)))
    rah2 = rah**2
!     B-H distance
    rbh = NORM2((xyz(:,B)+neigh%transVec(:,iTrB))-xyz(:,H))
    rbh2 = rbh**2

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))
!
!     out-of-line damp: A-H...B
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

!     out-of-line damp: A...nb(B)-B
    ! do i=1,nbb
    ranbprbnb = ranb+rbnb+1.d-12
    expo_nb = (param%hbnbcut/radab)*(ranbprbnb/rab-1.d0)
    ratio2_nb = exp(-expo_nb)
    outl_nb_tot = (2.d0/(1.d0+ratio2_nb))-1.0d0
    ! end do

!     long damping
    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

!     short damping
    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**6
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
    rbhdamp = damp*((p_bh/rbh2/rbh))
    rabdamp = damp*((p_ab/rab2/rab))
    rdamp = rbhdamp+rabdamp
    !Set up torsion paramter
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
        ! cycle iTr that are out of sqrt(hbtr2) cutoff
        ! vtors is set up for loop below, therefor the same cycle is included there aswell
        if (tlist(6,j) .le. 0.or.tlist(6,j) .gt. neigh%nTrans) cycle
        vtors(1,j) = pi/2.0
        vtors(2,j) = param%tors_hb
      end do
    end do
    !Calculate etors
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
    !Calculate gtors
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

    !Calculate eangl + gangl
    r0 = 120
    phi0 = r0*pi/180.
    bshift = param%bend_hb
    fc = 1.0d0-bshift
    call egbend_nci_mul(jj,kk,ll,vTrB,vTrC,phi0,fc,n,at,xyz,eangl,g3tmp)
    gangl(1:3,jj) = gangl(1:3,jj)+g3tmp(1:3,1)
    gangl(1:3,kk) = gangl(1:3,kk)+g3tmp(1:3,2)
    gangl(1:3,ll) = gangl(1:3,ll)+g3tmp(1:3,3)

!     hydrogen charge scaled term
    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

!     hydrogen charge scaled term
    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

!     hydrogen charge scaled term
    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

    qhoutl = qh*outl*outl_nb_tot

!     constant values, no gradient
    const = ca(2)*qa*cb(1)*qb*param%xhaci_coh
!     energy
    energy = -rdamp*qhoutl*eangl*etors*const
!     gradient
    drah(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-xyz(1:3,H)
    drbh(1:3) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-xyz(1:3,H)
    drab(1:3) = (xyz(1:3,A)+neigh%transVec(1:3,iTrA))-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))

    aterm = -rdamp*qh*outl_nb_tot*eangl*etors*const
    nbterm = -rdamp*qh*outl*eangl*etors*const
    dterm = -qhoutl*eangl*etors*const
    tterm = -rdamp*qhoutl*eangl*const
    bterm = -rdamp*qhoutl*etors*const

!------------------------------------------------------------------------------
!     damping part: rab
    gi = ((rabdamp+rbhdamp)*ddamp-3.d0*rabdamp)/rab2
    gi = gi*dterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = dg(1:3)
    gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     damping part: rbh
    gi = -3.d0*rbhdamp/rbh2
    gi = gi*dterm
    dg(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dg(1:3)
    gh(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     angular A-H...B term
!------------------------------------------------------------------------------
!     out of line term: rab
    tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    gi = -tmp1*rahprbh/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!     out of line term: rah,rbh
    gi = tmp1/rah
    dga(1:3) = gi*drah(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp1/rbh
    dgb(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

!------------------------------------------------------------------------------
!     angular A...nb(B)-B term
!------------------------------------------------------------------------------
!     out of line term: rab
    tmp2 = 2.d0*nbterm*ratio2_nb*expo_nb/(1+ratio2_nb)**2/(ranbprbnb-rab)
    gi_nb = -tmp2*ranbprbnb/rab2
    dg(1:3) = gi_nb*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!     out of line term: ranb,rbnb
    gi_nb = tmp2/ranb
    dga(1:3) = gi_nb*dranb(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi_nb = tmp2/rbnb
    dgb(1:3) = gi_nb*drbnb(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgnb(1:3) = -dga(1:3)-dgb(1:3)
    gnb(1:3) = dgnb(1:3)

!------------------------------------------------------------------------------
!     torsion term H...B=C<R1,R2
!------------------------------------------------------------------------------
    do i = 1,ntors
      ii = tlist(1,i)
      gdr(1:3,ii) = gdr(1:3,ii)+gtors(1:3,ii)*tterm
    end do
    gdr(1:3,jj) = gdr(1:3,jj)+gtors(1:3,jj)*tterm
    gdr(1:3,kk) = gdr(1:3,kk)+gtors(1:3,kk)*tterm
    gdr(1:3,ll) = gdr(1:3,ll)+gtors(1:3,ll)*tterm

!------------------------------------------------------------------------------
!     angle term H...B=C
!------------------------------------------------------------------------------
    gdr(1:3,jj) = gdr(1:3,jj)+gangl(1:3,jj)*bterm
    gdr(1:3,kk) = gdr(1:3,kk)+gangl(1:3,kk)*bterm
    gdr(1:3,ll) = gdr(1:3,ll)+gangl(1:3,ll)*bterm
    ! sigma
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
    ! torsion part
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
    ! jj, kk and ll same for every i in loop above (only ii or R changes)
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
    ! angle part
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

!------------------------------------------------------------------------------
!     move gradients into place
    gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
    gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
    gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
    ! B has one neighbor which is C in iTrC
    gdr(1:3,C) = gdr(1:3,C)+gnb(1:3)

  end subroutine abhgfnff_eg3

!> XB energy and analytical gradient
  subroutine rbxgfnff_eg(n,A,B,X,iTrB,iTrX,at,xyz,q,energy,gdr,param,neigh,sigma)

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
    real(wp) ::  rab,rax,rbx,rab2,rax2,rbx2,rax4,rbx4
    real(wp) ::  drax(3),drbx(3),drab(3),drm(3)
    real(wp) ::  dg(3),dga(3),dgb(3),dgx(3)
    real(wp) ::  gi,ga(3),gb(3),gx(3)
    real(wp) ::  ex1_a,ex2_a,ex1_b,ex2_b,ex1_x,ex2_x,expo
    real(wp) ::  aterm,dterm
    real(wp) ::  qa,qb,qx
    real(wp) ::  cx,cb
    real(wp) ::  gqa,gqb,gqx
    real(wp) ::  shortcut,const

    integer :: i,j

    gdr = 0
    energy = 0

    cb = 1. ! param%xhbas(at(B)) !
    cx = param%xbaci(at(X))

    ! compute distances;  A is in central cell !
    drax(1:3) = xyz(1:3,A)-(xyz(1:3,X)+neigh%transVec(1:3,iTrX))
    drbx(1:3) = (xyz(1:3,B)+neigh%transVec(1:3,iTrB))-(xyz(1:3,X)+neigh%transVec(1:3,iTrX))
    drab(1:3) = xyz(1:3,A)-(xyz(1:3,B)+neigh%transVec(1:3,iTrB))

    ! A-B distance !
    rab2 = sum(drab**2)
    rab = sqrt(rab2)

    ! A-X distance !
    rax2 = sum(drax**2)
    rax = sqrt(rax2)+1.d-12

    ! B-X distance !
    rbx2 = sum(drbx**2)
    rbx = sqrt(rbx2)+1.d-12

    ! out-of-line damp !
    expo = param%xbacut*((rax+rbx)/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow !
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

    ! long damping !
    ratio1 = (rbx2/param%hblongcut_xb)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

    ! short damping !
    shortcut = param%xbscut*(param%rad(at(A))+param%rad(at(B)))
    ratio3 = (shortcut/rbx2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    rdamp = damp/rbx2/rbx ! **2

    ! halogen charge scaled term !
    ex1_x = exp(param%xbst*q(X))
    ex2_x = ex1_x+param%xbsf
    qx = ex1_x/ex2_x

    ! donor charge scaled term !
    ex1_b = exp(-param%xbst*q(B))
    ex2_b = ex1_b+param%xbsf
    qb = ex1_b/ex2_b

    ! constant values, no gradient !
    const = cb*qb*cx*qx

    ! r^3 only slightly better than r^4 !
    aterm = -rdamp*const
    dterm = -outl*const
    energy = -rdamp*outl*const

    ! damping part rab !
    gi = rdamp*(-(2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3&
    &     /(1.d0+ratio3))-3.d0)/rbx2   ! 4,5,6 instead of 3 !
    gi = gi*dterm
    dg(1:3) = gi*drbx(1:3)
    gb(1:3) = dg(1:3)
    gx(1:3) = -dg(1:3)

    ! out of line term: rab !
    gi = 2.d0*ratio2*expo*(rax+rbx)/(1.d0+ratio2)**2/(rax+rbx-rab)/rab2
    gi = gi*aterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = +dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

    ! out of line term: rax,rbx !
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

    ! move gradients into place !
    gdr(1:3,1) = ga(1:3)
    gdr(1:3,2) = gb(1:3)
    gdr(1:3,3) = gx(1:3)
    ! sigma
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

!> taken from D3 ATM code
  subroutine batmgfnff_eg(n,iat,jat,kat,iTrj,iTrk,at,xyz,q,energy,g,ds,param,neigh)
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
    energy = c9*angr9 ! energy

    ! derivatives of each part w.r.t. r_ij,jk,ik
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

!> torsion term for rotation around triple bonded carbon
  subroutine sTors_eg(m,n,xyz,topo,energy,dg)
    integer,intent(in) :: m
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFTopology),intent(in) :: topo
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: dg(3,n)
    integer :: c1,c2,c3,c4
    integer :: i

    !> torsion angle between C1-C4
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

      ! dihedral angle in radians!
      phi = valijklff(n,xyz,c1,c2,c3,c4)
      call dphidr(n,xyz,c1,c2,c3,c4,phi,dp1,dp2,dp3,dp4)

      ! reference energy for torsion of 90° !
      ! calculated with DLPNO-CCSD(T) CBS on diphenylacetylene !
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

! ══════════════════════════════════════════════════════════════════════════════
! ══════════════════════════════════════════════════════════════════════════════
end module gfnff_engrad_module
