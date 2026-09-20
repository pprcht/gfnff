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
!> Closed-form Hessian of the GFN-FF repulsion, counterpart of gfnff_eg_rep.
!> Non-bonded and bonded terms share one radial kernel with different parameters.
module gfnff_hess_rep

  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData,TGFFTopology
  use gfnff_neighbor,only:TNeigh
  use gfnff_hess_pair,only:pair_hess_block,scatter_pair_hessian,scatter_pair_row
  implicit none
  private

  public :: hess_repulsion_nb,hess_repulsion_bonded

contains  !> MODULE PROCEDURES START HERE

  pure subroutine rep_radial(zzij,alpha,r1,fp,fpp)
    !***********************************************************************
    !* fp = f'(r), fpp = f''(r) of f(r) = zzij*exp(-alpha*r^1.5)/r, built from
    !* the log derivative g = f'/f as f' = g f, f'' = (g' + g^2) f.
    !***********************************************************************
    real(wp),intent(in) :: zzij,alpha,r1
    real(wp),intent(out) :: fp,fpp

    real(wp) :: t,u,fval,gfac,gp

    t = r1*sqrt(r1)          !> r^1.5 without a libm pow call
    u = alpha*t
    fval = zzij*exp(-u)/r1
    gfac = -(1.5_wp*u+1.0_wp)/r1
    gp = (1.0_wp-0.75_wp*u)/(r1*r1)
    fp = gfac*fval
    fpp = (gp+gfac*gfac)*fval

  end subroutine rep_radial

  subroutine hess_repulsion_nb(n,at,xyz,sqrab,repthr,mcf_nrep,param,topo,neigh,hess)
    !***********************************************************************
    !* Non-bonded repulsion Hessian, counterpart of eg_repulsion_nb: all pairs
    !* within the cutoff that are not bonded in the central cell. Molecular
    !* systems only; hess (3n,3n) is incremented.
    !*   sqrab    - packed squared distances
    !*   repthr   - squared repulsion cutoff
    !*   mcf_nrep - mcGFN-FF scaling factor (1.0 for standard GFN-FF)
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: sqrab(n*(n+1)/2)
    real(wp),intent(in) :: repthr,mcf_nrep
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: iat,jat,ati,atj,ij,ihi,ilo
    real(wp) :: r2,r1,vec(3),ehat(3),xi(3),zzi,zzij,fp,fpp,alp
    real(wp) :: blk(3,3)

    !$omp parallel do default(none) shared(n, at, xyz, sqrab, repthr, &
    !$omp topo, param, neigh, mcf_nrep, hess) &
    !$omp private(iat, jat, ati, atj, ij, ihi, ilo, r1, r2, vec, ehat, xi, &
    !$omp zzi, zzij, fp, fpp, blk, alp) schedule(dynamic)
    do iat = 1,n
      ati = at(iat)
      xi = xyz(:,iat)
      zzi = param%repz(ati)*param%repscaln*mcf_nrep
      do jat = 1,n
        if (jat .eq. iat) cycle
        !>-- sqrab and alphanb hold the lower triangle only
        ihi = max(iat,jat)
        ilo = min(iat,jat)
        ij = ihi*(ihi-1)/2+ilo
        r2 = sqrab(ij)
        if (r2 .gt. repthr.or.r2 .lt. 1.0e-8_wp) cycle
        !>-- bonded pairs are handled in hess_repulsion_bonded
        if (neigh%bpair(iat,jat,1) .eq. 1) cycle
        atj = at(jat)
        zzij = zzi*param%repz(atj)
        r1 = sqrt(r2)
        !>-- H...H exponent scaled by bond separation; alphanb excludes that factor (see gfnff_ini)
        alp = topo%alphanb(ij)
        if (ati .eq. 1.and.atj .eq. 1) alp = alp*topo%hhrep(neigh%bpair(ihi,ilo,1))
        call rep_radial(zzij,alp,r1,fp,fpp)
        vec = xi-xyz(:,jat)
        ehat = vec/r1
        call pair_hess_block(fp,fpp,r1,ehat,blk)
        call scatter_pair_row(hess,iat,jat,blk)
      end do
    end do
    !$omp end parallel do

  end subroutine hess_repulsion_nb

  subroutine hess_repulsion_bonded(n,at,xyz,param,neigh,hess)
    !***********************************************************************
    !* Bonded repulsion Hessian, counterpart of eg_repulsion_bonded, over the
    !* bond list with exponent sqrt(repa_i repa_j) and scaling repscalb.
    !* Molecular only, serial (bond list is O(n)); hess (3n,3n) is incremented.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(inout) :: hess(3*n,3*n)

    integer :: i,iat,jat,ati,atj
    real(wp) :: vec(3),ehat(3),r1,r2,alpha,zzij,fp,fpp
    real(wp) :: blk(3,3)

    do i = 1,neigh%nbond
      jat = neigh%blist(1,i)
      iat = neigh%blist(2,i)
      vec = xyz(:,iat)-xyz(:,jat)
      r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
      r1 = sqrt(r2)
      ati = at(iat)
      atj = at(jat)
      alpha = sqrt(param%repa(ati)*param%repa(atj))
      zzij = param%repz(ati)*param%repz(atj)*param%repscalb
      call rep_radial(zzij,alpha,r1,fp,fpp)
      ehat = vec/r1
      call pair_hess_block(fp,fpp,r1,ehat,blk)
      call scatter_pair_hessian(hess,iat,jat,blk)
    end do

  end subroutine hess_repulsion_bonded

end module gfnff_hess_rep
