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
!> GFN-FF repulsion, split into a non-bonded part over all pairs within the
!> cutoff and a bonded part over the bond list. Both use exp(-alpha*r**1.5)/r.
module gfnff_eg_rep

  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData,TGFFTopology
  use gfnff_neighbor,only:TNeigh
  use gfnff_geometry,only:lin
  implicit none
  private

  public :: eg_repulsion_nb,eg_repulsion_bonded

contains  !> MODULE PROCEDURES START HERE

  subroutine eg_repulsion_nb(n,at,xyz,sqrab,repthr,mcf_nrep,param,topo,neigh, &
        & erep,g,sigma)
    !***********************************************************************
    !* Non-bonded repulsion over all atom pairs, periodic images included,
    !* within the squared cutoff repthr. Bonded pairs are left to
    !* eg_repulsion_bonded. sqrab holds the packed squared distances of the
    !* central cell, mcf_nrep is the mcGFN-FF scaling (1.0 for GFN-FF).
    !* erep, g and sigma are incremented.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: sqrab(n*(n+1)/2)
    real(wp),intent(in) :: repthr,mcf_nrep
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(inout) :: erep
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)

    integer :: iat,jat,iTr,iTrDum,ati,atj,ij
    real(wp) :: r2,rab,r3(3),vec(3),xi(3),t8,t16,t19,t26,t27,zzi,zzij,alp

    !$omp parallel do default(none) reduction(+:erep, g, sigma) &
    !$omp shared(n, at, xyz, sqrab, repthr, &
    !$omp topo, param, neigh, mcf_nrep) &
    !$omp private(iat, jat, iTr, iTrDum, ati, atj, ij, rab, r2, r3, vec, xi, t8, t16, t19, t26, t27, zzi, zzij, alp)
    do iat = 1,n
      ati = at(iat)
      ij = iat*(iat-1)/2
      xi = xyz(:,iat)
      zzi = param%repz(ati)*param%repscaln*mcf_nrep
      do jat = 1,iat
        atj = at(jat)
        zzij = zzi*param%repz(atj)
        do iTr = 1,neigh%nTrans

          !>-- skip pairs beyond the cutoff and an atom paired with itself; without
          !>   images r2 is tabulated and rejects before the vector is built
          if (neigh%nTrans .eq. 1) then
            r2 = sqrab(ij+jat)
            if (r2 .gt. repthr.OR.r2 .lt. 1.0e-8_wp) cycle
            vec = xi-xyz(:,jat)
          else
            vec = xi-xyz(:,jat)+neigh%transVec(:,iTr)
            r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
            if (r2 .gt. repthr.OR.r2 .lt. 1.0e-8_wp) cycle
          end if

          !>-- bonded pairs are handled in eg_repulsion_bonded
          if (iTr .le. neigh%numctr) then
            if (neigh%bpair(iat,jat,iTr) .eq. 1) cycle ! list avoided because of memory
          end if
          rab = sqrt(r2)
          t16 = rab*sqrt(rab)   ! r2**0.75, but without a libm pow call
          t19 = t16*t16
          !>-- H...H pairs scale the stored exponent by their bond count, which
          !>   is not available beyond the central cells
          alp = topo%alphanb(lin(jat,iat))
          if (ati .eq. 1.and.atj .eq. 1) then
            if (iTr .le. neigh%numctr) then
              alp = alp*topo%hhrep(neigh%bpair(iat,jat,iTr))
            else
              alp = alp*topo%hhrep(0)
            end if
          end if
          t8 = t16*alp
          t26 = exp(-t8)*zzij
          erep = erep+(t26/rab)
          t27 = t26*(1.5d0*t8+1.0d0)/t19
          r3 = vec*t27
          if (neigh%nTrans .ne. 1) then   !> stress only wanted for PBC
            sigma(:,1) = sigma(:,1)-r3(1)*vec
            sigma(:,2) = sigma(:,2)-r3(2)*vec
            sigma(:,3) = sigma(:,3)-r3(3)*vec
          end if
          g(:,iat) = g(:,iat)-r3
          g(:,jat) = g(:,jat)+r3
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine eg_repulsion_nb

  subroutine eg_repulsion_bonded(n,at,xyz,param,neigh,erep,g,sigma)
    !***********************************************************************
    !* Repulsion of the bonded pairs in neigh%blist, with repa and repscalb
    !* in place of alphanb/repscaln. erep, g and sigma are incremented.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFData),intent(in) :: param
    type(TNeigh),intent(in) :: neigh
    real(wp),intent(inout) :: erep
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: sigma(3,3)

    integer :: i,iat,jat,iTr,ati,atj
    real(wp) :: xa,ya,za,dx,dy,dz,r2,rab,alpha,repab,t16,t19,t26,t27

    !$omp parallel do default(none) reduction(+:erep, g, sigma) &
    !$omp shared(param, at, xyz, neigh) &
    !$omp private(i, iTr, iat, jat, xa, ya, za, dx, dy, dz, r2, rab, ati, atj, &
    !$omp& alpha, repab, t16, t19, t26, t27)
    !> Both atom indices come out of the bond list, so two iterations can
    !> hit the same entry of g: an atom takes part in several bonds. That is
    !> a scatter a vectoriser may only touch with conflict detection, and ifx
    !> gets it wrong at -O2 and above, silently dropping contributions. The
    !> loop is short and memory bound, so forbidding it costs nothing.
    !> The directive is a comment to compilers that do not know it.
!DIR$ NOVECTOR
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
      r2 = dx*dx+dy*dy+dz*dz
      rab = sqrt(r2)
      ati = at(iat)
      atj = at(jat)
      alpha = sqrt(param%repa(ati)*param%repa(atj))
      repab = param%repz(ati)*param%repz(atj)*param%repscalb
      t16 = r2**0.75d0
      t19 = t16*t16
      t26 = exp(-alpha*t16)*repab
      erep = erep+t26/rab
      t27 = t26*(1.5d0*alpha*t16+1.0d0)/t19
      g(1,iat) = g(1,iat)-dx*t27
      g(2,iat) = g(2,iat)-dy*t27
      g(3,iat) = g(3,iat)-dz*t27
      g(1,jat) = g(1,jat)+dx*t27
      g(2,jat) = g(2,jat)+dy*t27
      g(3,jat) = g(3,jat)+dz*t27
      if (neigh%nTrans .eq. 1) cycle   !> stress only wanted for PBC
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

  end subroutine eg_repulsion_bonded

end module gfnff_eg_rep
