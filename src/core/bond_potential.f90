! ------------------------------------------------------------------------------
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
! along with gfnff.  If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------
!> GFN-FF bond stretch potential and its derivatives, shared by egbond,
!> egbond_hb and hess_bonds so the three call sites cannot drift apart.
module gfnff_bond_potential

  use iso_fortran_env,only:wp => real64
  use gfnff_param_tables,only:gffVersion
  implicit none
  private

  public :: bond_potential,bond_potential_hb

  !> exp(-1/2): the Gaussian's value at its inflection point, where the
  !> conformer variant leaves it. Written out rather than computed because
  !> exp() in a constant expression is not portable across compilers.
  real(wp),parameter :: einf = 0.60653065971263342426_wp

contains  !> MODULE PROCEDURES START HERE

  pure subroutine bond_potential(version,alpha,amp,d,e,ed,edd)
    !***********************************************************************
    !* Bond stretch energy e, ed = dE/dd and optional edd = d2E/dd2 for a
    !* bond outside a hydrogen bridge. amp = -D is negative when bonding,
    !* d = r - r0(CN), version is a gffVersion member. Shape in well().
    !***********************************************************************
    implicit none
    integer,intent(in) :: version
    real(wp),intent(in) :: alpha,amp,d
    real(wp),intent(out) :: e,ed
    real(wp),intent(out),optional :: edd

    real(wp) :: dum(4)

    call well(version,alpha,amp,d,e,ed,dum(1),dum(2),dum(3),dum(4))
    if (present(edd)) edd = dum(1)

  end subroutine bond_potential

  pure subroutine bond_potential_hb(version,al0,t1,hbcn,amp,d,e,ed,en,edd,enn,edn)
    !***********************************************************************
    !* Bond stretch energy for a bond in a hydrogen bridge, where the
    !* steepness alpha = (1 - t1 n) al0 depends on the hydrogen bond
    !* coordination number n = hbcn of the bridging H; t1 = 1 - vbond_scale.
    !* Returns e, ed = dE/dd, en = dE/dn and optionally the second
    !* derivatives edd, enn, edn. version, amp, d as in bond_potential.
    !***********************************************************************
    implicit none
    integer,intent(in) :: version
    real(wp),intent(in) :: al0,t1,hbcn,amp,d
    real(wp),intent(out) :: e,ed,en
    real(wp),intent(out),optional :: edd,enn,edn

    real(wp) :: alpha,dedd,ea,eaa,eda,dadn

    alpha = (-t1*hbcn+1.0_wp)*al0
    call well(version,alpha,amp,d,e,ed,dedd,ea,eaa,eda)
    if (present(edd)) edd = dedd

    !>-- The Gaussian branch keeps the closed forms of the published
    !>   implementation verbatim. They are the chain rule below written out,
    !>   but in a different order of multiplication, and floating point
    !>   multiplication is not associative: going through dE/dalpha instead
    !>   moves the hydrogen bonded gradients in the last bit.
    if (version /= gffVersion%conformer2020) then
      en = e*al0*d**2*t1
      if (present(enn)) enn = (al0*t1*d*d)**2*e
      if (present(edn)) edn = 2.0_wp*al0*t1*d*e*(1.0_wp-alpha*d*d)
      return
    end if

    dadn = -t1*al0
    en = ea*dadn
    if (present(enn)) enn = eaa*dadn*dadn
    if (present(edn)) edn = eda*dadn

  end subroutine bond_potential_hb

  pure subroutine well(version,alpha,amp,d,e,ed,edd,ea,eaa,eda)
    !***********************************************************************
    !* Bond well and its first and second partials in d and alpha.
    !* conformer2020 replaces the dissociative Gaussian tail past the
    !* inflection point, alpha d^2 = 1/2, by the tangent there. The join is
    !* C2 in d and in alpha, which bond_potential_hb needs.
    !***********************************************************************
    implicit none
    integer,intent(in) :: version
    real(wp),intent(in) :: alpha,amp,d
    real(wp),intent(out) :: e,ed,edd,ea,eaa,eda

    real(wp) :: sa,ec,ad,s

    if (version /= gffVersion%conformer2020.or.alpha*d**2 < 0.5_wp) then
      e = amp*exp(-alpha*d**2)
      ed = -2.0_wp*alpha*d*e
      edd = (-2.0_wp*alpha+4.0_wp*alpha*alpha*d*d)*e
      ea = -d*d*e
      eaa = d*d*d*d*e
      eda = (-2.0_wp*d+2.0_wp*alpha*d*d*d)*e
      return
    end if

    !>-- tangent past the inflection point, taken on |d| to stay even.
    !>   sign() would take the magnitude of its first argument, and amp is
    !>   negative, so the direction is carried explicitly
    sa = sqrt(alpha+alpha)
    ec = amp*einf
    ad = abs(d)
    s = sign(1.0_wp,d)
    e = ec*(2.0_wp-sa*ad)
    ed = -ec*sa*s
    edd = 0.0_wp
    ea = -ec*ad/sa
    eaa = ec*ad/(sa*(alpha+alpha))
    eda = -ec*s/sa

  end subroutine well

end module gfnff_bond_potential
