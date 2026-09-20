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
!> Cartesian nuclear Hessian of GFN-FF. Terms in `hess_terms_analytic` are
!> closed form; the rest are central finite differences of the analytic
!> gradient, restricted to exactly those terms via gfnff_eg's `terms` mask.
module gfnff_hess_driver

  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_data_types,only:TGFFData,TGFFTopology,TGFFNeighbourList,TCell
  use gfnff_neighbor,only:TNeigh
  use gfnff_alpb,only:TBorn
  use gfnff_eg_driver,only:gfnff_eg,gfnff_results, &
    &                          gff_term_rep,gff_term_es,gff_term_disp, &
    &                          gff_term_bond,gff_term_angl,gff_term_tors, &
    &                          gff_term_batm,gff_term_hb,gff_term_xb, &
    &                          gff_term_ext,gff_term_all
  use gfnff_hess_analysis,only:hess_symmetrize
  use gfnff_hess_rep,only:hess_repulsion_nb,hess_repulsion_bonded
  use gfnff_hess_bonded,only:hess_angles,hess_batm,hess_bonds, &
    &                          hess_torsions,hess_storsions, &
    &                          hess_torsions_available
  use gfnff_hess_es,only:hess_electrostatics
  use gfnff_hess_disp,only:hess_dispersion
  use gfnff_hess_hb,only:hess_xbonds,hess_hbonds_bound,hess_hbonds_unbound
  use gfnff_cn,only:gfnff_dlogcoord
  use gfnff_eg_hb,only:dncoord_erf
  use gfnff_param,only:gfnff_thresholds
  implicit none
  private

  public :: gfnff_hessian_core,hess_term_reference
  public :: hess_terms_analytic,hess_available_terms
  public :: hess_analytic_terms
  public :: hess_default_step

  !> Terms with a closed-form kernel. Whether a kernel applies to a given
  !> system is decided at runtime by hess_available_terms.
  integer,parameter :: hess_terms_analytic = gff_term_rep+gff_term_angl &
     & +gff_term_batm+gff_term_bond+gff_term_tors+gff_term_es &
     & +gff_term_disp+gff_term_xb+gff_term_hb+gff_term_ext

  !> Default central-difference step in bohr: O(h^2) truncation and O(eps/h)
  !> round-off of a ~1e-12 Eh/bohr gradient balance near 1e-8 Eh/bohr^2.
  real(wp),parameter :: hess_default_step = 5.0e-3_wp

contains  !> MODULE PROCEDURES START HERE

  subroutine gfnff_hessian_core(printlevel,n,at,xyz,cell,ichrg,param,topo,neigh, &
        & nlist,efield,solvation,version,accuracy,hess,energy,gradient,step, &
        & iostat,printunit)
    !***********************************************************************
    !* GFN-FF Cartesian Hessian d^2E/dR dR in Eh/bohr^2: closed-form kernels
    !* for hess_available_terms, central differences of the gradient for the
    !* rest. The lists (HB/XB, distance, D3) stay frozen while displacing.
    !* In:  cell (npbc /= 0 sends every term to the finite difference),
    !*      solvation (allocated if active), step (bohr, optional)
    !* Out: hess (3n,3n) with (c,A) -> 3*(A-1)+c; energy and gradient of all
    !*      terms at the reference geometry; iostat, always 0 at present
    !***********************************************************************
    implicit none
    integer,intent(in) :: printlevel
    integer,intent(in) :: n,at(n),ichrg
    real(wp),intent(in) :: xyz(3,n)
    type(TCell),intent(in) :: cell
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TGFFNeighbourList),intent(inout) :: nlist
    real(wp),intent(in) :: efield(3)
    integer,intent(in) :: version
    real(wp),intent(in) :: accuracy
    type(TBorn),allocatable,intent(inout) :: solvation
    real(wp),intent(out) :: hess(3*n,3*n)
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,n)
    real(wp),intent(in),optional :: step
    integer,intent(out),optional :: iostat
    integer,intent(in),optional :: printunit

    integer :: myunit,tfd,tana,io
    real(wp) :: h

    myunit = stdout
    if (present(printunit)) myunit = printunit
    h = hess_default_step
    if (present(step)) h = step
    io = 0
    hess = 0.0_wp

    !>-- reference energy and gradient, also primes the HB/XB lists
    call eg_at(printlevel,n,at,xyz,cell,ichrg,param,topo,neigh,nlist,efield, &
       & solvation,version,accuracy,gff_term_all,.true.,energy,gradient,myunit)

    tana = hess_available_terms(cell,neigh,topo,solvation,efield)
    tfd = gff_term_all-iand(gff_term_all,tana)
    if (tana .ne. 0) then
      call hess_analytic_terms(n,at,xyz,accuracy,version,tana,param,topo,neigh, &
         & nlist,hess,solvation)
    end if

    if (tfd .ne. 0) then
      call hess_finite_difference(printlevel,n,at,xyz,cell,ichrg,param,topo, &
         & neigh,nlist,efield,solvation,version,accuracy,tfd,h,hess,myunit)
    end if

    call hess_symmetrize(hess)

    if (present(iostat)) iostat = io

  end subroutine gfnff_hessian_core

  pure function hess_available_terms(cell,neigh,topo,solvation,efield) result(mask)
    !***********************************************************************
    !* Subset of hess_terms_analytic that is closed form for this system;
    !* the rest goes to the finite-difference pass. A periodic cell disables
    !* every kernel; special torsions or an unexpected torsion phase disable
    !* the torsion kernel; implicit solvation (shape correction included)
    !* disables nothing. A nonzero field disables gff_term_ext: eg_efield
    !* omits dq/dR from its gradient, whose exact derivative is then a
    !* non-symmetric Jacobian, so the finite difference mirrors the energy
    !* code instead.
    !***********************************************************************
    type(TCell),intent(in) :: cell
    type(TNeigh),intent(in) :: neigh
    type(TGFFTopology),intent(in) :: topo
    type(TBorn),allocatable,intent(in) :: solvation
    real(wp),intent(in) :: efield(3)
    integer :: mask

    mask = 0
    if (cell%npbc .ne. 0) return
    mask = hess_terms_analytic
    if (.not.hess_torsions_available(topo)) then
      mask = mask-iand(mask,gff_term_tors)
    end if
    if (sum(abs(efield)) .gt. 1.0e-6_wp) then
      mask = mask-iand(mask,gff_term_ext)
    end if

  end function hess_available_terms

  subroutine hess_analytic_terms(n,at,xyz,accuracy,version,terms,param,topo, &
        & neigh,nlist,hess,solvation)
    !***********************************************************************
    !* Add the closed-form Hessian of the terms in `terms` to hess (3n,3n),
    !* molecular systems only. Separate from the driver so the per-term
    !* tests run the same code path with a single bit set.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n),version,terms
    real(wp),intent(in) :: xyz(3,n),accuracy
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TGFFNeighbourList),intent(in) :: nlist
    type(TBorn),allocatable,intent(in),optional :: solvation
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: i,j,k,ncov
    logical :: havesolv
    logical :: needcn
    real(wp) :: dispthr,cnthr,repthr,hbthr1,hbthr2,mcf_nrep,mcf_ehb
    real(wp),allocatable :: sqrab(:),srab(:),cn(:),dcn(:,:,:)
    real(wp),allocatable :: hb_cn(:),hb_dcn(:,:,:),dhbcndL(:,:,:)

    call gfnff_thresholds(accuracy,dispthr,cnthr,repthr,hbthr1,hbthr2)
    mcf_nrep = mcgfnff_repulsion_scale(version)
    mcf_ehb = mcgfnff_hb_scale(version)

    !>-- packed squared and plain distances, the layout gfnff_eg builds
    allocate (sqrab(n*(n+1)/2),srab(n*(n+1)/2))
    do i = 1,n
      k = i*(i-1)/2
      do j = 1,i
        sqrab(k+j) = sum((xyz(:,i)-xyz(:,j))**2)
        srab(k+j) = sqrt(sqrab(k+j))
      end do
    end do

    needcn = iand(terms,gff_term_es+gff_term_bond+gff_term_disp) .ne. 0
    if (needcn) then
      allocate (cn(n),dcn(3,n,n))
      call gfnff_dlogcoord(n,at,xyz,srab,cn,dcn,cnthr,param)
    end if

    if (iand(terms,gff_term_rep) .ne. 0) then
      call hess_repulsion_nb(n,at,xyz,sqrab,repthr,mcf_nrep,param,topo,neigh,hess)
      if (neigh%nbond .gt. 0) then
        call hess_repulsion_bonded(n,at,xyz,param,neigh,hess)
      end if
    end if

    if (iand(terms,gff_term_angl) .ne. 0) then
      call hess_angles(n,at,xyz,param,topo,neigh,hess)
    end if

    if (iand(terms,gff_term_batm) .ne. 0) then
      call hess_batm(n,at,xyz,param,topo,neigh,hess)
    end if

    if (iand(terms,gff_term_tors) .ne. 0) then
      call hess_torsions(n,at,xyz,param,topo,neigh,hess)
      call hess_storsions(n,xyz,topo,hess)
    end if

    if (iand(terms,gff_term_es) .ne. 0) then
      havesolv = .false.
      if (present(solvation)) havesolv = allocated(solvation)
      if (havesolv) then
        call hess_electrostatics(n,at,xyz,srab,cnthr,cn,dcn,nlist%q,param,topo, &
           & hess,gbsa=solvation)
      else
        call hess_electrostatics(n,at,xyz,srab,cnthr,cn,dcn,nlist%q,param,topo, &
           & hess)
      end if
    end if

    if (iand(terms,gff_term_disp) .ne. 0) then
      call hess_dispersion(n,at,xyz,sqrab,srab,dispthr,cnthr,cn,dcn,param,topo,hess)
    end if

    if (iand(terms,gff_term_xb) .ne. 0) then
      call hess_xbonds(n,at,xyz,param,topo,neigh,nlist,hess)
    end if

    if (iand(terms,gff_term_hb) .ne. 0) then
      call hess_hbonds_bound(n,at,xyz,mcf_ehb,param,topo,neigh,nlist,hess)
      call hess_hbonds_unbound(n,at,xyz,mcf_ehb,param,topo,neigh,nlist,hess,ncov)
    end if

    !>-- gff_term_ext arrives only with a vanishing field: zero contribution

    if (iand(terms,gff_term_bond) .ne. 0) then
      allocate (hb_cn(n),hb_dcn(3,n,n),dhbcndL(3,3,n))
      hb_cn = 0.0_wp
      hb_dcn = 0.0_wp
      if (sum(neigh%nr_hb) .gt. 0) then
        call dncoord_erf(n,at,xyz,param%rcov,hb_cn,hb_dcn,900.0_wp,topo,neigh, &
           & dhbcndL)
      end if
      call hess_bonds(n,at,xyz,srab,cnthr,cn,dcn,hb_cn,hb_dcn, &
         & param,topo,neigh,version,hess)
      deallocate (hb_cn,hb_dcn,dhbcndL)
    end if

    if (needcn) deallocate (cn,dcn)

  end subroutine hess_analytic_terms

  subroutine hess_finite_difference(printlevel,n,at,xyz,cell,ichrg,param,topo, &
        & neigh,nlist,efield,solvation,version,accuracy,terms,h,hess,myunit)
    !***********************************************************************
    !* Add the central-difference Hessian of the terms in `terms` to hess:
    !*   H(:,b) += ( g(R + h e_b) - g(R - h e_b) ) / 2h,   b = (c,B).
    !* HB/XB lists stay fixed (update = .false.) to keep the function smooth.
    !***********************************************************************
    implicit none
    integer,intent(in) :: printlevel
    integer,intent(in) :: n,at(n),ichrg,version,terms,myunit
    real(wp),intent(in) :: xyz(3,n),efield(3),accuracy,h
    type(TCell),intent(in) :: cell
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TGFFNeighbourList),intent(inout) :: nlist
    type(TBorn),allocatable,intent(inout) :: solvation
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: iat,ic,ib,ndof
    real(wp) :: ep,em
    real(wp),allocatable :: xtmp(:,:),gp(:,:),gm(:,:)

    ndof = 3*n
    allocate (xtmp(3,n),gp(3,n),gm(3,n))

    do iat = 1,n
      do ic = 1,3
        ib = 3*(iat-1)+ic

        xtmp = xyz
        xtmp(ic,iat) = xyz(ic,iat)+h
        call eg_at(printlevel,n,at,xtmp,cell,ichrg,param,topo,neigh,nlist, &
           & efield,solvation,version,accuracy,terms,.false.,ep,gp,myunit)

        xtmp = xyz
        xtmp(ic,iat) = xyz(ic,iat)-h
        call eg_at(printlevel,n,at,xtmp,cell,ichrg,param,topo,neigh,nlist, &
           & efield,solvation,version,accuracy,terms,.false.,em,gm,myunit)

        hess(:,ib) = hess(:,ib)+reshape((gp-gm)/(2.0_wp*h), [ndof])
      end do
    end do

  end subroutine hess_finite_difference

  subroutine hess_term_reference(printlevel,n,at,xyz,cell,ichrg,param,topo, &
        & neigh,nlist,efield,solvation,version,accuracy,terms,hess,step,printunit)
    !***********************************************************************
    !* Finite-difference reference Hessian of the gff_term_* mask `terms`,
    !* used to validate the closed-form kernels; a single bit isolates one
    !* term. hess (3n,3n) is overwritten and symmetrised.
    !***********************************************************************
    implicit none
    integer,intent(in) :: printlevel
    integer,intent(in) :: n,at(n),ichrg,version,terms
    real(wp),intent(in) :: xyz(3,n),efield(3),accuracy
    type(TCell),intent(in) :: cell
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TGFFNeighbourList),intent(inout) :: nlist
    type(TBorn),allocatable,intent(inout) :: solvation
    real(wp),intent(out) :: hess(3*n,3*n)
    real(wp),intent(in),optional :: step
    integer,intent(in),optional :: printunit

    integer :: myunit
    real(wp) :: h,e0
    real(wp),allocatable :: g0(:,:)

    myunit = stdout
    if (present(printunit)) myunit = printunit
    h = hess_default_step
    if (present(step)) h = step
    hess = 0.0_wp

    !>-- one full call first so the HB/XB lists exist and stay frozen after
    allocate (g0(3,n))
    call eg_at(printlevel,n,at,xyz,cell,ichrg,param,topo,neigh,nlist,efield, &
       & solvation,version,accuracy,gff_term_all,.true.,e0,g0,myunit)

    call hess_finite_difference(printlevel,n,at,xyz,cell,ichrg,param,topo, &
       & neigh,nlist,efield,solvation,version,accuracy,terms,h,hess,myunit)

    call hess_symmetrize(hess)

  end subroutine hess_term_reference

  subroutine eg_at(printlevel,n,at,xyz,cell,ichrg,param,topo,neigh,nlist,efield, &
        & solvation,version,accuracy,terms,update,energy,gradient,myunit)
    !***********************************************************************
    !* gfnff_eg with the term mask, minus the stress tensor and results log.
    !***********************************************************************
    implicit none
    integer,intent(in) :: printlevel,n,at(n),ichrg,version,terms,myunit
    real(wp),intent(in) :: xyz(3,n),efield(3),accuracy
    type(TCell),intent(in) :: cell
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TGFFNeighbourList),intent(inout) :: nlist
    type(TBorn),allocatable,intent(inout) :: solvation
    logical,intent(in) :: update
    real(wp),intent(out) :: energy,gradient(3,n)

    real(wp) :: sigma(3,3)
    type(gfnff_results) :: res

    gradient = 0.0_wp
    call gfnff_eg(printlevel,n,at,xyz,cell,sigma,ichrg,gradient,energy,res, &
       & param,topo,neigh,nlist,efield,solvation,update,version,accuracy, &
       & printunit=myunit,terms=terms)

  end subroutine eg_at

  pure function mcgfnff_hb_scale(version) result(mcf_ehb)
    !***********************************************************************
    !* mcGFN-FF hydrogen bond scale, mirrors the top of gfnff_eg; else 1.0.
    !***********************************************************************
    use gfnff_param,only:gffVersion
    integer,intent(in) :: version
    real(wp) :: mcf_ehb
    if (version == gffVersion%mcgfnff2023) then
      mcf_ehb = 0.727406_wp
    else
      mcf_ehb = 1.0_wp
    end if
  end function mcgfnff_hb_scale

  pure function mcgfnff_repulsion_scale(version) result(mcf_nrep)
    !***********************************************************************
    !* mcGFN-FF non-bonded repulsion scale, mirrors the top of gfnff_eg;
    !* else 1.0.
    !***********************************************************************
    use gfnff_param,only:gffVersion
    integer,intent(in) :: version
    real(wp) :: mcf_nrep
    if (version == gffVersion%mcgfnff2023) then
      mcf_nrep = 1.343608_wp
    else
      mcf_nrep = 1.0_wp
    end if
  end function mcgfnff_repulsion_scale

end module gfnff_hess_driver
