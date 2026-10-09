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
!> Shared low-level helpers for the GFN-FF energy terms: the external field
!> coupling and the bend/torsion damping functions.
module gfnff_eg_terms

  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData,TGFFTopology
  implicit none
  private

  public :: gfnffdampa,gfnffdampt
  public :: gfnffdampa_nci,gfnffdampt_nci
  public :: eg_efield

contains  !> MODULE PROCEDURES START HERE

  subroutine eg_efield(n,xyz,efield,q,topo,eext,g)
    !************************************************************
    !* Interaction of the charges q with the external field efield,
    !* measured relative to the reference geometry topo%xyze0.
    !* Skipped for a numerically zero field. eext and g are
    !* incremented.
    !************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: efield(3)
    real(wp),intent(in) :: q(n)
    type(TGFFTopology),intent(in) :: topo
    real(wp),intent(inout) :: eext
    real(wp),intent(inout) :: g(3,n)

    integer :: i
    real(wp) :: r3(3)

    if (sum(abs(efield)) .le. 1d-6) return

    do i = 1,n
      r3(:) = -q(i)*efield(:)
      g(:,i) = g(:,i)+r3(:)
      eext = eext+r3(1)*(xyz(1,i)-topo%xyze0(1,i))+&
      &                    r3(2)*(xyz(2,i)-topo%xyze0(2,i))+&
      &                    r3(3)*(xyz(3,i)-topo%xyze0(3,i))
    end do

  end subroutine eg_efield

  subroutine gfnffdampa(ati,atj,r2,damp,ddamp,param)
    !************************************************************
    !* Bend damping, 1 at short and 0 at long distance, so that
    !* bonds can dissociate. Input are the atomic numbers ati/atj
    !* and the squared distance r2, cutoff from atcuta and rcov.
    !* Output: damp and ddamp = (1/r)*d(damp)/dr
    !************************************************************
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
    !************************************************************
    !* Torsion damping; gfnffdampa with the torsion cutoff atcutt.
    !************************************************************
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
    !************************************************************
    !* Bend damping function for the NCI variant (atcuta_nci).
    !************************************************************
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
    !************************************************************
    !* Torsion damping function for the NCI variant (atcutt_nci).
    !************************************************************
    implicit none
    type(TGFFData),intent(in) :: param
    integer ati,atj
    real(wp) r2,damp,ddamp,rr,rcut

    rcut = param%atcutt_nci*(param%rcov(ati)+param%rcov(atj))**2
    rr = (r2/rcut)**2
    damp = 1.0d0/(1.0d0+rr)
    ddamp = -2.d0*2*rr/(r2*(1.0d0+rr)**2)

  end subroutine gfnffdampt_nci
end module gfnff_eg_terms
