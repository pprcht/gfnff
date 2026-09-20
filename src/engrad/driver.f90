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

module gfnff_eg_driver

  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_topo_hbset,only:gfnff_hbset,gfnff_hbset0
  use gfnff_data_types,only:TGFFData,TGFFNeighbourList,new,TGFFTopology, &
    &                       TCell,TDispersionData
  use gfnff_neighbor,only:TNeigh
  use gfnff_alpb,only:TBorn
  use gfnff_param,only:sqrtZr4r2,gffVersion,gfnff_thresholds
  use gfnff_cn,only:getCoordinationNumber,gfnff_dlogcoord
  use gfnff_gdisp0,only:d3_gradient,d3_gradientPBC
  use gfnff_timing,only:gfnff_timer
  use gfnff_math_wrapper,only:gemv
  use gfnff_eg_terms,only:eg_efield
  use gfnff_eg_es,only:goed_gfnff,goed_pbc_gfnff,es_grad_sigma,es_grad_mol
  use gfnff_eg_hb,only:dncoord_erf, &
    &                  eg_hbonds_bound,eg_hbonds_unbound,eg_xbonds
  use gfnff_eg_bonded,only:eg_bonds,eg_bonds_harmonic,eg_angles,eg_torsions, &
    &                      eg_storsions,eg_batm
  use gfnff_eg_rep,only:eg_repulsion_nb,eg_repulsion_bonded
  implicit none
  private
  public :: gfnff_eg,gfnff_results
  public :: gff_term_rep,gff_term_es,gff_term_disp,gff_term_bond
  public :: gff_term_angl,gff_term_tors,gff_term_batm,gff_term_hb
  public :: gff_term_xb,gff_term_ext,gff_term_all

!&<
  !> Bitmask for the optional `terms` argument of gfnff_eg. It gates only the
  !> accumulation into etot and g; all prerequisites (lists, CNs, EEQ solve) still
  !> run, so a single term is available in isolation with unchanged charges and CNs.
  integer,parameter :: gff_term_rep  = 1     !> repulsion, non-bonded and bonded
  integer,parameter :: gff_term_es   = 2     !> EEQ electrostatics (+ solvation)
  integer,parameter :: gff_term_disp = 4     !> D3 dispersion
  integer,parameter :: gff_term_bond = 8     !> bond stretch
  integer,parameter :: gff_term_angl = 16    !> angle bending
  integer,parameter :: gff_term_tors = 32    !> torsion, incl. special torsion
  integer,parameter :: gff_term_batm = 64    !> bonded ATM three-body
  integer,parameter :: gff_term_hb   = 128   !> hydrogen bonds
  integer,parameter :: gff_term_xb   = 256   !> halogen bonds
  integer,parameter :: gff_term_ext  = 512   !> external electric field
  integer,parameter :: gff_term_all  = 1023  !> every term, the default
!&>

!&<
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

contains  !> MODULE PROCEDURES START HERE

  subroutine gfnff_eg(printlevel,n,at,xyz,cell,sigma,ichrg,g,etot,res_gff, &
        & param,topo,neigh,nlist,efield,solvation,update,version,accuracy,printunit, &
        & terms)
    !***********************************************************************
    !* GFN-FF energy and analytical gradient. Requires the D3 setup and a
    !* prior gfnff_ini call. Bend/torsion trigonometry adapted from QMDFF,
    !* repulsion and rabguess from the xtb GFN0 part.
    !* Input:
    !*   solvation - GBSA/ALPB model, allocated only if active
    !*   update    - rebuild the HB/XB lists
    !* Output:
    !*   g, etot   - gradient (Eh/Bohr) and total energy (Eh)
    !*   sigma     - stress tensor (Eh), non-zero only for PBC
    !*   res_gff   - energy decomposition
    !***********************************************************************
    implicit none

    character(len=*),parameter :: source = 'gfnff_eg'
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
    integer,intent(in),optional :: terms      !< gff_term_* mask (default: all)

    real(wp),intent(out) :: sigma(3,3) ! stress tensor
    real(wp),intent(out) :: g(3,n)
    real(wp),intent(out) :: etot
    logical :: pr
    integer :: myunit
    integer :: tmask

    real(wp) :: edisp,ees,ebond,eangl,etors,erep,ehb,exb,ebatm,eext
    real(wp) :: gsolv,gborn,ghb,gsasa,gshift

    integer :: i,nd3
    logical :: require_update

    real(wp),allocatable :: sqrab(:),srab(:)
    real(wp) :: dist(n,n)
    integer,allocatable  :: d3list(:,:)
    !> dcn is dCN/dr from gfnff_dlogcoord (molecular) or getCoordinationNumber (PBC)
    real(wp),allocatable :: cn(:),dcn(:,:,:),dcndL(:,:,:)
    real(wp),allocatable :: hb_cn(:),hb_dcn(:,:,:),dhbcndL(:,:,:)
    real(wp),allocatable :: eeqtmp(:,:),qtmp(:)
    real(wp),allocatable :: gTrans(:,:),rTrans(:,:) ! reciprocal, direct translation vector for ES
    real(wp),allocatable :: xtmp(:)
    real(wp) :: convF  !convergence factor alpha, aka ewald parameter

    type(gfnff_timer) :: timer
    real(wp) :: dispthr,cnthr,repthr,hbthr1,hbthr2
    !> mcGFN-FF term scaling factors; nrep/ees/ehb are 1.0 for standard GFN-FF
    real(wp) :: mcf_ees,mcf_ehb,mcf_nrep,mcf_s8

    pr = printlevel >= 2
    if (present(printunit)) then
      myunit = printunit
    else
      myunit = stdout
    end if
    if (present(terms)) then
      tmask = terms
    else
      tmask = gff_term_all
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

    !>-- translation vectors within the maximum cutoff, at least the central 27 cells (3D)
    neigh%oldCutOff = 0.0_wp
    call neigh%getTransVec(n,at,xyz,cell,60.0_wp)

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

    allocate (sqrab(n*(n+1)/2),srab(n*(n+1)/2),qtmp(n), &
    &         eeqtmp(2,n*(n+1)/2),d3list(2,n*(n+1)/2),dcn(3,n,n),cn(n), &
    &         dcndL(3,3,n),hb_dcn(3,n,n),hb_cn(n),dhbcndL(3,3,n))

    if (printlevel >= 2) then
      call timer%new(10+count([allocated(solvation)]))
    else if (printlevel == 1) then
      call timer%new(1)
      call timer%measure(1,'iter. time')
    end if

    if (pr) call timer%measure(1,'distance/D3 list')
    call build_distance_lists(n,xyz,dispthr,cell%npbc,sqrab,srab,nd3,d3list,dist)
    if (pr) call timer%measure(1)

    if (pr) call timer%measure(10,'HB/XB (incl list setup)')
    call setup_hb_lists(n,at,xyz,hbthr1,hbthr2,printlevel,myunit,topo,neigh, &
         & nlist,require_update)
    if (pr) call timer%measure(10)

    if (allocated(solvation)) then
      call timer%measure(11,"GBSA")
      call solvation%update(at,xyz)
      call timer%measure(11)
    end if

    if (pr) call timer%measure(2,'non bonded repulsion')
    if (iand(tmask,gff_term_rep) .ne. 0) then
      call eg_repulsion_nb(n,at,xyz,sqrab,repthr,mcf_nrep,param,topo,neigh,erep,g,sigma)
    end if
    if (pr) call timer%measure(2)

    !>-- crude mode for 2D-3D conversion: harmonic bonds with estimated Re.
    !>   It returns before the decomposition is filled at the end, so its two terms
    !>   are recorded here, else print_gfnff_results shows an all-zero breakdown.
    if (version == gffVersion%harmonic2020) then
      call eg_bonds_harmonic(n,at,xyz,param,neigh,ebond,g)
      etot = ebond+erep
      res_gff%e_bond = ebond
      res_gff%e_rep = erep
      res_gff%e_total = etot
      res_gff%gnorm = sqrt(sum(g**2))
      return
    end if

    if (cell%npbc .eq. 0) then
      if (pr) call timer%measure(3,'dCN')
      call gfnff_dlogcoord(n,at,xyz,srab,cn,dcn,cnthr,param) ! new erf used in GFN0
      dcndL = 0.0_wp
      dhbcndL = 0.0_wp
      if (sum(neigh%nr_hb) .gt. 0) call dncoord_erf(n,at,xyz,param%rcov,hb_cn,hb_dcn,900.0d0,topo,neigh,dhbcndL) ! HB erf CN
      if (pr) call timer%measure(3)

    else
      if (pr) call timer%measure(3,'dCN')
      if (sum(neigh%nr_hb) .gt. 0) call dncoord_erf(n,at,xyz,param%rcov,hb_cn,hb_dcn,900.0d0,topo,neigh,dhbcndL) ! HB erf CN
      call getCoordinationNumber(n,at,xyz,neigh%nTrans,neigh%transVec,60.0_wp,5,cn,dcn,dcndL,param)
      if (pr) call timer%measure(3)
    end if

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

    if (cell%npbc .eq. 0) then
      if (pr) call timer%measure(5,'D3')
      if (nd3 .gt. 0.and.iand(tmask,gff_term_disp) .ne. 0) then
        call d3_gradient(topo%dispm,n,at,xyz,nd3,d3list,topo%zetac6, &
        & param%d3r0,sqrtZr4r2,4.0d0,param%dispscale,cn,dcn,edisp,g)
      end if
      deallocate (d3list)
      if (pr) call timer%measure(5)
    else if (iand(tmask,gff_term_disp) .ne. 0) then
      !>-- inter-molecular dispersion with adjusted parameters, then intra-molecular.
      !>   NOTE: these two labels are the wrong way round, .true. selects
      !>   same-fragment pairs. Kept as-is, see the note in disp_gradient_latp.
      call d3_gradientPBC(topo%dispm,n,at,xyz,cell,topo%fraglist,neigh%nTrans,neigh%transVec,mcdisp_par,4.0_wp,topo%zetac6, &
           & param%d3r0,60.0_wp,.true.,cn,dcn,dcndL,edisp,g,sigma)
      call d3_gradientPBC(topo%dispm,n,at,xyz,cell,topo%fraglist,neigh%nTrans,neigh%transVec,disp_par,4.0_wp,topo%zetac6, &
           & param%d3r0,60.0_wp,.false.,cn,dcn,dcndL,edisp,g,sigma)
    end if

    !>-- the EEQ solve above always runs, its charges feed the HB/XB and field
    !>   terms; only the energy and the gradient are gated
    if (iand(tmask,gff_term_es) .eq. 0) ees = 0.0_wp

    if (iand(tmask,gff_term_es) .eq. 0) then
      deallocate (eeqtmp)
    else if (cell%npbc .eq. 0) then
      if (pr) call timer%measure(6,'EEQ gradient')
      call es_grad_mol(n,xyz,sqrab,srab,eeqtmp,nlist%q,g)
      deallocate (eeqtmp)

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

      deallocate (eeqtmp)

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

      !>-- qtmp = q * dX/dCN, where X is the right-hand side
      do i = 1,n
        qtmp(i) = nlist%q(i)*param%cnf(at(i))/(2.0d0*sqrt(cn(i))+1.d-16)
      end do

      call gemv(dcn,qtmp,g,alpha=-mcf_ees,beta=1.0_wp)
      call gemv(dcndL,qtmp,sigma,alpha=-mcf_ees,beta=1.0_wp)

      if (pr) call timer%measure(6)
    end if

    if (pr) call timer%measure(7,'bonds')
    if (neigh%nbond .gt. 0) then
      if (iand(tmask,gff_term_bond) .ne. 0) then
        call eg_bonds(n,at,xyz,cn,dcn,dcndL,hb_cn,hb_dcn,dhbcndL, &
             & param,topo,neigh,version,ebond,g,sigma)
      end if
      deallocate (dcn,dcndL,hb_dcn)

      if (iand(tmask,gff_term_rep) .ne. 0) then
        call eg_repulsion_bonded(n,at,xyz,param,neigh,erep,g,sigma)
      end if
    end if ! if neigh%nbond.gt.0
    if (pr) call timer%measure(7)

    if (pr) call timer%measure(8,'bend and torsion')
    if (iand(tmask,gff_term_angl) .ne. 0) then
      call eg_angles(n,at,xyz,param,topo,neigh,eangl,g,sigma)
    end if
    if (iand(tmask,gff_term_tors) .ne. 0) then
      call eg_torsions(n,at,xyz,param,topo,neigh,etors,g,sigma)
      call eg_storsions(n,xyz,topo,etors,g)
    end if
    if (pr) call timer%measure(8)

    if (pr) call timer%measure(9,'bonded ATM')
    if (iand(tmask,gff_term_batm) .ne. 0) then
      call eg_batm(n,at,xyz,param,topo,neigh,ebatm,g,sigma)
    end if
    if (pr) call timer%measure(9)

    !>-- correct number of translation vectors for the HB lists (hblist1/2)
    call neigh%getTransVec(n,at,xyz,cell,sqrt(hbthr2))
    if (pr) call timer%measure(10,'HB/XB (incl list setup)')
    if (update.or.require_update) then
      call gfnff_hbset(n,at,xyz,topo,neigh,nlist,hbthr1,hbthr2)
    end if

    if (iand(tmask,gff_term_hb) .ne. 0) then
      call eg_hbonds_bound(n,at,xyz,mcf_ehb,param,topo,neigh,nlist,ehb,g,sigma)
      call eg_hbonds_unbound(n,at,xyz,sqrab,srab,mcf_ehb, &
           & param,topo,neigh,nlist,ehb,g,sigma)
    end if

    if (iand(tmask,gff_term_xb) .ne. 0) then
      call eg_xbonds(n,at,xyz,param,topo,neigh,nlist,exb,g,sigma)
    end if
    if (pr) call timer%measure(10)

    if (iand(tmask,gff_term_ext) .ne. 0) then
      call eg_efield(n,xyz,efield,nlist%q,topo,eext,g)
    end if

    etot = ees+edisp+erep+ebond &
    &           +eangl+etors+ehb+exb+ebatm+eext &
    &           +gsolv

    if (printlevel >= 2) then

      call timer%write(myunit,'E+G')
      !>-- diagnostic only, the results below are filled either way
      if (abs(sum(nlist%q)-ichrg) .gt. 1.d-1) then
        write (myunit,*) nlist%q
        write (myunit,*) sum(nlist%q),ichrg
        write (myunit,'("**ERROR**",a,1x,a)') 'EEQ charge constrain error',source
      end if

    else if (printlevel == 1) then
      call timer%measure(1)
    end if

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

  subroutine build_distance_lists(n,xyz,dispthr,npbc,sqrab,srab,nd3,d3list,dist)
    !***********************************************************************
    !* Packed lower-triangle distances and the pairs within the dispersion
    !* cutoff; for PBC also the full matrix that the Ewald EEQ solver expects.
    !*   dispthr    - squared dispersion cutoff
    !*   sqrab/srab - squared and plain distances, diagonal zero
    !*   nd3/d3list - pairs inside the cutoff; only the first nd3 are defined
    !*   dist       - filled only for npbc /= 0
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,npbc
    real(wp),intent(in) :: xyz(3,n),dispthr
    real(wp),intent(out) :: sqrab(n*(n+1)/2),srab(n*(n+1)/2)
    integer,intent(out) :: nd3,d3list(2,n*(n+1)/2)
    real(wp),intent(out) :: dist(n,n)

    integer :: i,j,k,ij

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

      !>-- the loop above skips the diagonal: zero it without adding it to d3list
      sqrab(ij+i) = 0.0d0
      srab(ij+i) = 0.0d0
    end do

    if (npbc .ne. 0) then
      dist = 0.0_wp
      !$omp parallel do collapse(2) default(none) shared(dist,xyz,n) &
      !$omp private(i,j)
      do i = 1,n
        do j = 1,n
          dist(j,i) = NORM2(xyz(:,j)-xyz(:,i))
        end do
      end do
      !$omp end parallel do
    end if

  end subroutine build_distance_lists

  subroutine setup_hb_lists(n,at,xyz,hbthr1,hbthr2,printlevel,myunit,topo,neigh, &
        & nlist,require_update)
    !***********************************************************************
    !* Count the HB/XB candidates of the current geometry and (re)allocate
    !* nlist if it is missing or too small. The lists themselves are filled
    !* by the caller.
    !*   require_update - .true. if nlist was reallocated or the candidate
    !*                    count grew, so that the caller rebuilds the lists
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: hbthr1,hbthr2
    integer,intent(in) :: printlevel,myunit
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TGFFNeighbourList),intent(inout) :: nlist
    logical,intent(out) :: require_update

    integer :: nhb1,nhb2,nxb

    if (allocated(nlist%q)) then
      nlist%initialized = size(nlist%q) == n
    end if
    call gfnff_hbset0(n,at,xyz,topo,nhb1,nhb2,nxb,neigh,nlist,hbthr1,hbthr2)
    !>-- compare against the allocated bound, not nlist%nhb1: gfnff_hbset stores
    !>   the number of bonds found there, which would discard the headroom of
    !>   new(nlist,...) and force a reallocation as soon as one candidate appears
    if (nlist%initialized) then
      nlist%initialized = nhb1 <= size(nlist%hblist1,dim=2) &
         & .and.nhb2 <= size(nlist%hblist2,dim=2) &
         & .and.nxb <= size(nlist%hblist3,dim=2)
    end if
    require_update = .not.nlist%initialized
    !>-- more candidates than the last fill found: a pair came inside the cutoff.
    !>   Rebuild in place; resetting the reference geometry is what lets
    !>   gfnff_hbset pass its displacement check
    if (nlist%initialized) then
      if (nhb1 > nlist%nhb1.or.nhb2 > nlist%nhb2.or.nxb > nlist%nxb) then
        require_update = .true.
        nlist%hbrefgeo(:,:) = xyz
      end if
    end if
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
    !>-- the lists are not filled here: getTransVec has not run with the HB
    !>   cutoff yet, so they would be built against the repulsion translation set

  end subroutine setup_hb_lists

end module gfnff_eg_driver
