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

!> Vector, angle and torsion geometry with the analytic derivatives the
!> gradient needs, in both the molecular and the periodic form. Each PBC
!> variant takes the translation vectors of the atoms involved, so an angle
!> or torsion can span cell boundaries.
!>
!> lin() lives here too: it indexes the packed lower triangle that the
!> distance and pair arrays are stored in, which is the data structure all
!> of these routines are read alongside.
module gfnff_geometry
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  public :: lin
  public :: vlen,vsub,valijklff,valijklffPBC,torsPBC
  public :: omega,domegadr,dphidr,bangl,banglPBC,impsc
  public :: omegaPBC,domegadrPBC,dphidrPBC
  public :: crossproduct
  public :: crprod
  interface crprod
    module procedure crossprod
  end interface crprod

  real(wp),parameter :: eps = 1.0d-14
  real(wp),parameter :: pi = 3.1415926535897932384626433832795029d0

!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!

  subroutine vsub(a,b,c,n)
    implicit none
    integer n
    real(wp) :: a(n),b(n),c(n)
    integer :: i
    do i = 1,n
      c(i) = a(i)-b(i)
    end do
    return
  end subroutine vsub

! ══════════════════════════════════════════════════════════════════════════════

  real(wp) function vlen(a)
    implicit none !double precision (a-h,o-z)
    real(wp) :: a(3)
    real(wp) :: tot

    tot = a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
    vlen = 0.0d0
    if (tot .gt. 0.0d0) vlen = dsqrt(tot)

    return
  end function vlen

! ══════════════════════════════════════════════════════════════════════════════

  real(wp) function valijklff(natoms,xyz,i,j,k,l)
    implicit none
    integer :: ic,i,j,k,l,natoms
    real(wp) :: xyz(3,natoms)
    real(wp) :: ra(3),rb(3),rc(3),na(3),nb(3)
    real(wp) :: thab,thbc
    real(wp) :: nan,nbn,snanb,deter

    !>-- get torsion coordinate
    do ic = 1,3
      ra(ic) = xyz(ic,j)-xyz(ic,i)
      rb(ic) = xyz(ic,k)-xyz(ic,j)
      rc(ic) = xyz(ic,l)-xyz(ic,k)
    end do

    !>-- determinante of rb,ra,rc
    deter = ra(1)*(rb(2)*rc(3)-rb(3)*rc(2))  &
   &      -ra(2)*(rb(1)*rc(3)-rb(3)*rc(1)) &
   &      +ra(3)*(rb(1)*rc(2)-rb(2)*rc(1))

    thab = valijk(natoms,xyz,i,k,j)
    thbc = valijk(natoms,xyz,j,l,k)
    call crossprod(ra,rb,na)
    call crossprod(rb,rc,nb)
    nan = vecnorm(na,3,1)
    nbn = vecnorm(nb,3,1)

    snanb = 0.0d0
    do ic = 1,3
      snanb = snanb+na(ic)*nb(ic)
    end do
    if (abs(abs(snanb)-1.d0) .lt. eps) then
      snanb = sign(1.d0,snanb)
    end if

    valijklff = acos(snanb)
  end function valijklff
! ──────────────────────────────────────────────────────────────────────────────
  real(wp) function valijklffPBC(mo,natoms,xyz,i,j,k,l,vTrj,vTrk,vTrl)

    implicit none

    integer::  ic,i,j,k,l,natoms,mo

    real(wp) :: xyz(3,natoms),vTrj(3),vTrk(3),vTrl(3)
    real(wp) :: ra(3),rb(3),rc(3),na(3),nb(3)
    real(wp) :: nan,nbn,snanb

    !> get torsion coordinate
    if (mo .eq. 1) then ! egtors call -> j (=ii) in central cell
      do ic = 1,3
        ra(ic) = xyz(ic,j)-(xyz(ic,i)+vTrl(ic))
        rb(ic) = (xyz(ic,k)+vTrj(ic))-xyz(ic,j)
        rc(ic) = (xyz(ic,l)+vTrk(ic))-(xyz(ic,k)+vTrj(ic))
      end do
    else ! abhgfnff_eg3 call -> l (=H) in central cell
      do ic = 1,3
        ra(ic) = (xyz(ic,j)+vTrk(ic))-(xyz(ic,i)+vTrj(ic)) ! B - R
        rb(ic) = (xyz(ic,k)+vTrl(ic))-(xyz(ic,j)+vTrk(ic)) ! C - B
        rc(ic) = xyz(ic,l)-(xyz(ic,k)+vTrl(ic)) ! H - C
      end do
    end if

    ! NOTE: the triple product of ra,rb,rc and the two bond angles thab/thbc
    ! used to be evaluated here and then discarded; they never entered the
    ! result. Each thab/thbc cost two vector normalisations and an acos.

    call crossprod(ra,rb,na)
    call crossprod(rb,rc,nb)
    nan = vecnorm(na,3,1)
    nbn = vecnorm(nb,3,1)

    snanb = 0.0d0
    do ic = 1,3 ! scalar product of the crossproducts
      snanb = snanb+na(ic)*nb(ic)
    end do
    if (abs(abs(snanb)-1.d0) .lt. eps) then
      snanb = sign(1.d0,snanb)
    end if

    valijklffPBC = acos(snanb)
  end function valijklffPBC

! ══════════════════════════════════════════════════════════════════════════════

  real(wp) Function valijk(nat,xyz,j,k,i)
    implicit none
    integer :: nat,j,k,i,ic
    real(wp) :: ra(3),rb(3),rab,eps
    real(wp) :: xyz(3,nat),ran,rbn
    parameter(eps=1.d-14)

    do ic = 1,3
      ra(ic) = xyz(ic,j)-xyz(ic,i)
      rb(ic) = xyz(ic,k)-xyz(ic,i)
    end do

    ran = vecnorm(ra,3,1)
    rbn = vecnorm(rb,3,1)
    rab = 0.d0
    do ic = 1,3
      rab = rab+ra(ic)*rb(ic)
    end do

    if (abs(abs(rab)-1.d0) .lt. eps) then
      rab = sign(1.d0,rab)
    end if
    valijk = acos(rab)

  End Function valijk
! ──────────────────────────────────────────────────────────────────────────────
  real(wp) Function valijkPBC(mode,nat,xyz,j,k,i,vTr1,vTr2,vTr3)
    implicit none
    integer mode,nat,j,k,i,ic
    real(wp) :: ra(3),rb(3),rab
    real(wp) :: xyz(3,nat),ran,rbn,vTr1(3),vTr2(3),vTr3(3)

    if (mode .eq. 1) then ! here j=l,k=j,i=i are inserted, vTr1=vTrl, vTr2=vTrj
      do ic = 1,3
        ra(ic) = (xyz(ic,j)+vTr1(ic))-xyz(ic,i)
        rb(ic) = (xyz(ic,k)+vTr2(ic))-xyz(ic,i)
      end do
    elseif (mode .eq. 2) then ! here j=i,k=k,i=j are inserted vTr1=vTrk, vTr2=vTrj
      do ic = 1,3
        ra(ic) = xyz(ic,j)-(xyz(ic,i)+vTr2(ic))
        rb(ic) = (xyz(ic,k)+vTr1(ic))-(xyz(ic,i)+vTr2(ic))
      end do
    elseif (mode .eq. 3) then ! here j=R k=C i=B vTr1=vTrR vTr2=vTrC vTr3=vTrB
      do ic = 1,3
        ra(ic) = (xyz(ic,j)+vTr1(ic))-(xyz(ic,i)+vTr3(ic)) ! R - B
        rb(ic) = (xyz(ic,k)+vTr2(ic))-(xyz(ic,i)+vTr3(ic)) ! C - B
      end do
    elseif (mode .eq. 4) then ! here j=B k=H i=C vTr1=vTrB vTr2=vTrC
      do ic = 1,3
        ra(ic) = (xyz(ic,j)+vTr1(ic))-(xyz(ic,i)+vTr2(ic))
        rb(ic) = (xyz(ic,k))-(xyz(ic,i)+vTr2(ic))
      end do
    end if

    ran = vecnorm(ra,3,1)
    rbn = vecnorm(rb,3,1)
    rab = 0.d0
    do ic = 1,3
      rab = rab+ra(ic)*rb(ic)
    end do

    if (abs(abs(rab)-1.d0) .lt. eps) then
      rab = sign(1.d0,rab)
    end if
    valijkPBC = acos(rab)

  end function valijkPBC

! ══════════════════════════════════════════════════════════════════════════════

  subroutine crossprod(ra,rb,rab)
    implicit none
    real(wp) :: ra(3),rb(3),rab(3)
    rab(1) = ra(2)*rb(3)-ra(3)*rb(2)
    rab(2) = ra(3)*rb(1)-ra(1)*rb(3)
    rab(3) = ra(1)*rb(2)-ra(2)*rb(1)
  end Subroutine crossprod

  function crossproduct(ra,rb) result(rab)
    implicit none
    real(wp) :: ra(3),rb(3),rab(3)
    call crossprod(ra,rb,rab)
  end function crossproduct

  real(wp) Function vecnorm(r,n,inorm)
    implicit none
    integer :: i,n,inorm
    real(wp) :: r(n),or,sp,rn
    sp = 0.0_wp
    do i = 1,n
      sp = sp+r(i)*r(i)
    end do
    rn = sqrt(sp)
    if (inorm .gt. 0) then
      if (abs(rn) .gt. 1.d-14) then
        or = 1.0_wp/rn
        do i = 1,n
          r(i) = or*r(i)
        end do
      end if
    end if
    vecnorm = rn
  end function vecnorm

!========================================================================================!

  pure elemental integer function lin(i1,i2)
    integer,intent(in) :: i1,i2
    integer :: idum1,idum2
    idum1 = max(i1,i2)
    idum2 = min(i1,i2)
    lin = idum2+idum1*(idum1-1)/2
    return
  end function lin

!========================================================================================!

  real(wp) Function omega(nat,xyz,i,j,k,l)
    !>   Calculates the inversion angle
    implicit none
    integer :: ic,nat,i,j,k,l

    real(wp) :: xyz(3,nat)
    real(wp) :: rd(3),re(3),rn(3),rv(3),rnv
    real(wp) :: rnn,rvn

    do ic = 1,3
      re(ic) = xyz(ic,i)-xyz(ic,j)
      rd(ic) = xyz(ic,k)-xyz(ic,j)
      rv(ic) = xyz(ic,l)-xyz(ic,i)
    end do
    call crossprod(re,rd,rn)
    rnn = vecnorm(rn,3,1)
    rvn = vecnorm(rv,3,1)

    rnv = rn(1)*rv(1)+rn(2)*rv(2)+rn(3)*rv(3)
    omega = asin(rnv)

  End Function omega

  Subroutine domegadr(nat,xyz,i,j,k,l,omega, &
  &          domegadri,domegadrj,domegadrk,domegadrl)
    !> inversion derivatives
    implicit none
    integer  :: ic,i,j,k,l,nat
    real(wp) :: omega,sinomega
    real(wp) :: xyz(3,nat),onenner,rnn,rvn
    real(wp) :: rn(3),rv(3),rd(3),re(3),rdme(3),rve(3)
    real(wp) :: rne(3),rdv(3),rdn(3)
    real(wp) :: rvdme(3),rndme(3),nenner
    real(wp) :: domegadri(3),domegadrj(3),domegadrk(3),domegadrl(3),eps
    parameter(eps=1.d-14)

    sinomega = sin(omega)

    do ic = 1,3
      rv(ic) = xyz(ic,l)-xyz(ic,i)
      rd(ic) = xyz(ic,k)-xyz(ic,j)
      re(ic) = xyz(ic,i)-xyz(ic,j)

      rdme(ic) = rd(ic)-re(ic)
    end do

    call crossprod(re,rd,rn)
    rvn = vecnorm(rv,3,0)
    rnn = vecnorm(rn,3,0)

    call crossprod(rv,re,rve)
    call crossprod(rn,re,rne)
    call crossprod(rd,rv,rdv)
    call crossprod(rd,rn,rdn)
    call crossprod(rv,rdme,rvdme)
    call crossprod(rn,rdme,rndme)

    nenner = rnn*rvn*cos(omega)
    if (abs(nenner) .gt. eps) then
      onenner = 1.d0/nenner
      do ic = 1,3
! ... domega/dri
        domegadri(ic) = onenner*(rdv(ic)-rn(ic)- &
                                 sinomega*(rvn/rnn*rdn(ic)-rnn/rvn*rv(ic)))

! ... domega/drj
        domegadrj(ic) = onenner*(rvdme(ic)-sinomega*rvn/rnn*rndme(ic))

! ... domega/drk
        domegadrk(ic) = onenner*(rve(ic)-sinomega*rvn/rnn*rne(ic))

! ... domega/drl
        domegadrl(ic) = onenner*(rn(ic)-sinomega*rnn/rvn*rv(ic))
      end do
    else
      do ic = 1,3
        domegadri(ic) = 0.d0
        domegadrj(ic) = 0.d0
        domegadrk(ic) = 0.d0
        domegadrl(ic) = 0.d0
      end do
    end if

  End Subroutine domegadr

! ──────────────────────────────────────────────────────────────────────────────
  real(wp) Function omegaPBC(nat,xyz,i,j,k,l,vTr1,vTr2,vTr3)
    !   Calculates the inversion angle (with PBC)
    !  .....................................................................
    implicit none
    integer :: ic,nat,i,j,k,l

    real(wp) :: xyz(3,nat),vTr1(3),vTr2(3),vTr3(3),&
       &        rd(3),re(3),rn(3),rv(3),rnv,&
       &        rnn,rvn
    ! out-of-plane case from ini; atoms and iTr's sorted by distance to atom i
    ! i=central, j=1st nb, k=2nd, l=3rd
    do ic = 1,3
      re(ic) = xyz(ic,i)-(xyz(ic,j)+vTr2(ic))                ! Vec central to 1st nb
      rd(ic) = (xyz(ic,k)+vTr3(ic))-(xyz(ic,j)+vTr2(ic))            ! Vec 1st to 2nd nb
      rv(ic) = (xyz(ic,l)+vTr1(ic))-xyz(ic,i) ! Vec central to 3rd nb
    end do
    call crossprod(re,rd,rn)
    rnn = vecnorm(rn,3,1)
    rvn = vecnorm(rv,3,1)

    rnv = rn(1)*rv(1)+rn(2)*rv(2)+rn(3)*rv(3)
    omegaPBC = asin(rnv)

  End function omegaPBC
! ──────────────────────────────────────────────────────────────────────────────

  Subroutine domegadrPBC(nat,xyz,i,j,k,l,vTr1,vTr2,vTr3,omega,&
        &            domegadri,domegadrj,domegadrk,domegadrl)
    !     inversion derivatives (with PBC)
    !  .....................................................................
    implicit none
    integer :: ic,i,j,k,l,nat

    real(wp) ::  omega,sinomega,&
       &         vTr1(3),vTr2(3),vTr3(3), &
       &         xyz(3,nat),onenner,rnn,rvn,&
       &         rn(3),rv(3),rd(3),re(3),rdme(3),rve(3),&
       &         rne(3),rdv(3),rdn(3),&
       &         rvdme(3),rndme(3),nenner,&
       &         domegadri(3),domegadrj(3),domegadrk(3),domegadrl(3)

    sinomega = sin(omega)

    do ic = 1,3
      re(ic) = xyz(ic,i)-(xyz(ic,j)+vTr2(ic))            ! Vec central to 1st nb
      rd(ic) = (xyz(ic,k)+vTr3(ic))-(xyz(ic,j)+vTr2(ic))  ! Vec 1st to 2nd nb
      rv(ic) = (xyz(ic,l)+vTr1(ic))-xyz(ic,i)             ! Vec central to 3rd nb

      rdme(ic) = rd(ic)-re(ic)
    end do

    call crossprod(re,rd,rn)
    rvn = vecnorm(rv,3,0)
    rnn = vecnorm(rn,3,0)

    call crossprod(rv,re,rve)
    call crossprod(rn,re,rne)
    call crossprod(rd,rv,rdv)
    call crossprod(rd,rn,rdn)
    call crossprod(rv,rdme,rvdme)
    call crossprod(rn,rdme,rndme)

    nenner = rnn*rvn*cos(omega)
    if (abs(nenner) .gt. eps) then
      onenner = 1.d0/nenner
      do ic = 1,3
        ! ... domega/dri
        domegadri(ic) = onenner*(rdv(ic)-rn(ic)-&
           &                       sinomega*(rvn/rnn*rdn(ic)-rnn/rvn*rv(ic)))

        ! ... domega/drj
        domegadrj(ic) = onenner*(rvdme(ic)-sinomega*rvn/rnn*rndme(ic))

        ! ... domega/drk
        domegadrk(ic) = onenner*(rve(ic)-sinomega*rvn/rnn*rne(ic))

        ! ... domega/drl
        domegadrl(ic) = onenner*(rn(ic)-sinomega*rnn/rvn*rv(ic))
      end do
    else
      do ic = 1,3
        domegadri(ic) = 0.d0
        domegadrj(ic) = 0.d0
        domegadrk(ic) = 0.d0
        domegadrl(ic) = 0.d0
      end do
    end if
  end subroutine domegadrPBC

! ══════════════════════════════════════════════════════════════════════════════
  Subroutine dphidr(nat,xyz,i,j,k,l,phi, &
  &                dphidri,dphidrj,dphidrk,dphidrl)
    !> the torsion derivatives
    implicit none
    integer :: ic,i,j,k,l,nat
    real(wp) :: sinphi,cosphi,onenner
    real(wp) :: ra(3),rb(3),rc(3),rab(3),rac(3),rbc(3),rbb(3)
    real(wp) :: raa(3),rba(3),rapba(3),rapbb(3),rbpca(3),rbpcb(3)
    real(wp) :: rapb(3),rbpc(3),na(3),nb(3),nan,nbn
    real(wp) :: dphidri(3),dphidrj(3),dphidrk(3),dphidrl(3)
    real(wp) :: xyz(3,nat),phi,nenner,eps

    parameter(eps=1.d-14)

    cosphi = cos(phi)
    sinphi = sin(phi)
    do ic = 1,3
      ra(ic) = xyz(ic,j)-xyz(ic,i)
      rb(ic) = xyz(ic,k)-xyz(ic,j)
      rc(ic) = xyz(ic,l)-xyz(ic,k)

      rapb(ic) = ra(ic)+rb(ic)
      rbpc(ic) = rb(ic)+rc(ic)
    end do

    call crossprod(ra,rb,na)
    call crossprod(rb,rc,nb)
    nan = vecnorm(na,3,0)
    nbn = vecnorm(nb,3,0)

    nenner = nan*nbn*sinphi
    if (abs(nenner) .lt. eps) then
      dphidri = 0
      dphidrj = 0
      dphidrk = 0
      dphidrl = 0
      onenner = 1.0d0/(nan*nbn)
    else
      onenner = 1.d0/nenner
    end if
    call crossprod(na,rb,rab)
    call crossprod(nb,ra,rba)
    call crossprod(na,rc,rac)
    call crossprod(nb,rb,rbb)
    call crossprod(nb,rc,rbc)
    call crossprod(na,ra,raa)

    call crossprod(rapb,na,rapba)
    call crossprod(rapb,nb,rapbb)
    call crossprod(rbpc,na,rbpca)
    call crossprod(rbpc,nb,rbpcb)

! ... dphidri
    do ic = 1,3
      dphidri(ic) = onenner*(cosphi*nbn/nan*rab(ic)-rbb(ic))

! ... dphidrj
      dphidrj(ic) = onenner*(cosphi*(nbn/nan*rapba(ic) &
                                     +nan/nbn*rbc(ic)) &
                             -(rac(ic)+rapbb(ic)))
! ... dphidrk
      dphidrk(ic) = onenner*(cosphi*(nbn/nan*raa(ic) &
                                     +nan/nbn*rbpcb(ic)) &
                             -(rba(ic)+rbpca(ic)))
! ... dphidrl
      dphidrl(ic) = onenner*(cosphi*nan/nbn*rbb(ic)-rab(ic))
    end do

  End Subroutine dphidr
! ──────────────────────────────────────────────────────────────────────────────
  Subroutine dphidrPBC(mode,nat,xyz,i,j,k,l,vTrR,vTrB,vTrC,phi,&
        &                  dphidri,dphidrj,dphidrk,dphidrl)
    !     the torsion derivatives with PBC images
    implicit none

    integer :: mode,ic,i,j,k,l,nat

    real(wp) :: vTrR(3),vTrB(3),vTrC(3), &
    &           sinphi,cosphi,onenner,&
    &           ra(3),rb(3),rc(3),rab(3),rac(3),rbc(3),rbb(3),&
    &           raa(3),rba(3),rapba(3),rapbb(3),rbpca(3),rbpcb(3),&
    &           rapb(3),rbpc(3),na(3),nb(3),nan,nbn,&
    &           dphidri(3),dphidrj(3),dphidrk(3),dphidrl(3),&
    &           xyz(3,nat),phi,nenner

    cosphi = cos(phi)
    sinphi = sin(phi)
    if (mode .eq. 1) then
      do ic = 1,3
        ra(ic) = xyz(ic,j)+vTrB(ic)-xyz(ic,i)-vTrR(ic)
        rb(ic) = xyz(ic,k)+vTrC(ic)-xyz(ic,j)-vTrB(ic)
        rc(ic) = xyz(ic,l)-xyz(ic,k)-vTrC(ic)

        rapb(ic) = ra(ic)+rb(ic)
        rbpc(ic) = rb(ic)+rc(ic)
      end do
    elseif (mode .eq. 2) then
      do ic = 1,3
        ra(ic) = -(xyz(ic,i)+vTrC(ic))+xyz(ic,j)
        rb(ic) = -xyz(ic,j)+(xyz(ic,k)+vTrR(ic))
        rc(ic) = -(xyz(ic,k)+vTrR(ic))+(xyz(ic,l)+vTrB(ic))

        rapb(ic) = ra(ic)+rb(ic)
        rbpc(ic) = rb(ic)+rc(ic)
      end do
    end if
    call crossprod(ra,rb,na)
    call crossprod(rb,rc,nb)
    nan = vecnorm(na,3,0)
    nbn = vecnorm(nb,3,0)

    nenner = nan*nbn*sinphi
    if (abs(nenner) .lt. eps) then
      dphidri = 0
      dphidrj = 0
      dphidrk = 0
      dphidrl = 0
      if (abs(nan*nbn) .gt. eps) then
        onenner = 1.0d0/(nan*nbn)
      else
        onenner = 0.0d0
      end if
    else
      onenner = 1.d0/nenner
    end if

    call crossprod(na,rb,rab)
    call crossprod(nb,ra,rba)
    call crossprod(na,rc,rac)
    call crossprod(nb,rb,rbb)
    call crossprod(nb,rc,rbc)
    call crossprod(na,ra,raa)

    call crossprod(rapb,na,rapba)
    call crossprod(rapb,nb,rapbb)
    call crossprod(rbpc,na,rbpca)
    call crossprod(rbpc,nb,rbpcb)

    if (abs(onenner) .gt. eps) then
      do ic = 1,3
        ! ... dphidri
        dphidri(ic) = onenner*(cosphi*nbn/nan*rab(ic)-rbb(ic))
        ! ... dphidrj
        dphidrj(ic) = onenner*(cosphi*(nbn/nan*rapba(ic)&
           &                                +nan/nbn*rbc(ic))&
           &                        -(rac(ic)+rapbb(ic)))
        ! ... dphidrk
        dphidrk(ic) = onenner*(cosphi*(nbn/nan*raa(ic)&
           &                             +nan/nbn*rbpcb(ic))&
           &                        -(rba(ic)+rbpca(ic)))
        ! ... dphidrl
        dphidrl(ic) = onenner*(cosphi*nan/nbn*rbb(ic)-rab(ic))
      end do
    else
      dphidri = 0.0d0
      dphidrj = 0.0d0
      dphidrk = 0.0d0
      dphidrl = 0.0d0
    end if

  End subroutine dphidrPBC

! ──────────────────────────────────────────────────────────────────────────────
  subroutine torsPBC(mo,nat,xyz,i,j,k,l,vTr1,vTr2,vTr3,phi,&
        &            dphidri,dphidrj,dphidrk,dphidrl)
    !***********************************************************************
    !* Torsion angle and its Cartesian derivatives in one pass.
    !* Replaces a valijklffPBC + dphidrPBC pair, which each rebuilt the same
    !* three bond vectors, the same two normal vectors and the same two
    !* vector norms, and additionally round tripped the angle through
    !* acos -> cos/sin.
    !* Input:
    !*   mo        - 1: bonded torsion (vTr1/2/3 = vTrj,vTrk,vTrl)
    !*               2: HB torsion     (vTr1/2/3 = vTrR,vTrB,vTrC)
    !*   nat/xyz   - system definition
    !*   i,j,k,l   - the four atoms of the torsion
    !* Output:
    !*   phi       - torsion angle in radians, in [0,pi]
    !*   dphidr*   - derivative of phi w.r.t. each of the four positions
    !***********************************************************************
    implicit none
    integer,intent(in) :: mo,nat,i,j,k,l
    real(wp),intent(in) :: xyz(3,nat),vTr1(3),vTr2(3),vTr3(3)
    real(wp),intent(out) :: phi
    real(wp),intent(out) :: dphidri(3),dphidrj(3),dphidrk(3),dphidrl(3)

    integer :: ic
    real(wp) :: ra(3),rb(3),rc(3),rapb(3),rbpc(3),na(3),nb(3)
    real(wp) :: rab(3),rac(3),rbc(3),rbb(3),raa(3),rba(3)
    real(wp) :: rapba(3),rapbb(3),rbpca(3),rbpcb(3)
    real(wp) :: nan,nbn,snanb,cosphi,sinphi,nenner,onenner

    if (mo .eq. 1) then   !> bonded torsion, j in the central cell
      do ic = 1,3
        ra(ic) = xyz(ic,j)-(xyz(ic,i)+vTr3(ic))
        rb(ic) = (xyz(ic,k)+vTr1(ic))-xyz(ic,j)
        rc(ic) = (xyz(ic,l)+vTr2(ic))-(xyz(ic,k)+vTr1(ic))
      end do
    else                  !> HB torsion, l (=H) in the central cell
      do ic = 1,3
        ra(ic) = (xyz(ic,j)+vTr2(ic))-(xyz(ic,i)+vTr1(ic))
        rb(ic) = (xyz(ic,k)+vTr3(ic))-(xyz(ic,j)+vTr2(ic))
        rc(ic) = xyz(ic,l)-(xyz(ic,k)+vTr3(ic))
      end do
    end if
    do ic = 1,3
      rapb(ic) = ra(ic)+rb(ic)
      rbpc(ic) = rb(ic)+rc(ic)
    end do

    call crossprod(ra,rb,na)
    call crossprod(rb,rc,nb)
    nan = vecnorm(na,3,0)
    nbn = vecnorm(nb,3,0)

    !> cos(phi) straight from the normals; no normalise-then-dot round trip
    if (nan*nbn .gt. eps) then
      snanb = (na(1)*nb(1)+na(2)*nb(2)+na(3)*nb(3))/(nan*nbn)
    else
      snanb = 0.0d0
    end if
    snanb = min(1.0d0,max(-1.0d0,snanb))
    if (abs(abs(snanb)-1.d0) .lt. eps) snanb = sign(1.d0,snanb)
    phi = acos(snanb)

    cosphi = snanb
    !> phi is in [0,pi] so sin(phi) is non negative
    sinphi = sqrt(max(0.0d0,1.0d0-snanb*snanb))

    nenner = nan*nbn*sinphi
    if (abs(nenner) .lt. eps) then
      dphidri = 0
      dphidrj = 0
      dphidrk = 0
      dphidrl = 0
      if (abs(nan*nbn) .gt. eps) then
        onenner = 1.0d0/(nan*nbn)
      else
        onenner = 0.0d0
      end if
    else
      onenner = 1.d0/nenner
    end if

    if (abs(onenner) .le. eps) then
      dphidri = 0.0d0
      dphidrj = 0.0d0
      dphidrk = 0.0d0
      dphidrl = 0.0d0
      return
    end if

    call crossprod(na,rb,rab)
    call crossprod(nb,ra,rba)
    call crossprod(na,rc,rac)
    call crossprod(nb,rb,rbb)
    call crossprod(nb,rc,rbc)
    call crossprod(na,ra,raa)

    call crossprod(rapb,na,rapba)
    call crossprod(rapb,nb,rapbb)
    call crossprod(rbpc,na,rbpca)
    call crossprod(rbpc,nb,rbpcb)

    do ic = 1,3
      dphidri(ic) = onenner*(cosphi*nbn/nan*rab(ic)-rbb(ic))
      dphidrj(ic) = onenner*(cosphi*(nbn/nan*rapba(ic)&
         &                                +nan/nbn*rbc(ic))&
         &                        -(rac(ic)+rapbb(ic)))
      dphidrk(ic) = onenner*(cosphi*(nbn/nan*raa(ic)&
         &                             +nan/nbn*rbpcb(ic))&
         &                        -(rba(ic)+rbpca(ic)))
      dphidrl(ic) = onenner*(cosphi*nan/nbn*rbb(ic)-rab(ic))
    end do

  end subroutine torsPBC

! ══════════════════════════════════════════════════════════════════════════════

  pure subroutine bangl(xyz,i,j,k,angle)
    implicit none
    real(wp),intent(in)  :: xyz(3,*)
    real(wp),intent(out) :: angle
    integer,intent(in)  :: i,j,k

    real(wp) d2ij,d2jk,d2ik,xy,temp

    d2ij = sum((xyz(:,i)-xyz(:,j))**2)
    d2jk = sum((xyz(:,j)-xyz(:,k))**2)
    d2ik = sum((xyz(:,i)-xyz(:,k))**2)
    xy = sqrt(d2ij*d2jk)
    temp = 0.5d0*(d2ij+d2jk-d2ik)/xy
    if (temp .gt. 1.0d0) temp = 1.0d0
    if (temp .lt. -1.0d0) temp = -1.0d0
    angle = acos(temp)

  end subroutine bangl

  pure subroutine banglPBC(mode,xyz,i,j,k,iTr,iTr2,transVec,angle)
    implicit none
    real(wp),intent(in)  :: xyz(3,*)
    integer,intent(in)  :: mode,i,j,k,iTr,iTr2  ! j is in the middle
    real(wp),intent(in) :: transVec(:,:)
!    type(TNeigh),intent(in) :: neigh
    real(wp),intent(out) :: angle

    real(wp) :: d2ij,d2jk,d2ik,xy,temp,trV(3),trV2(3)
    !trV = neigh%transVec(:,iTr)
    !trV2 = neigh%transVec(:,iTr2)
    trV = transVec(:,iTr)
    trV2 = transVec(:,iTr2)
    if (mode .eq. 1) then
      d2ij = sum(((xyz(:,i)+trV)-xyz(:,j))**2)
      d2jk = sum((xyz(:,j)-(xyz(:,k)+trV2))**2)
      d2ik = sum(((xyz(:,i)+trV)-(xyz(:,k)+trV2))**2)
    end if
    if (mode .eq. 2) then
      d2ij = sum((xyz(:,i)-(xyz(:,j)+trV))**2)
      d2jk = sum(((xyz(:,j)+trV)-(xyz(:,k)+trV2))**2)
      d2ik = sum((xyz(:,i)-(xyz(:,k)+trV2))**2)
    end if
    xy = sqrt(d2ij*d2jk)
    temp = 0.5d0*(d2ij+d2jk-d2ik)/xy  ! the angle is between side dij and djk
    if (temp .gt. 1.0d0) temp = 1.0d0
    if (temp .lt. -1.0d0) temp = -1.0d0
    angle = acos(temp)

  end subroutine banglPBC
!========================================================================================!

  pure subroutine impsc(a,b,c)
    implicit none
    real(wp),intent(in)  :: a(3),b(3)
    real(wp),intent(out) :: c
    integer  :: i
    real(wp) :: rimp,al,bl

    rimp = 0.0d0

    do i = 1,3
      rimp = rimp+a(i)*b(i)
    end do

    al = norm2(a)
    bl = norm2(b)

    if (al .gt. 0.0d0.and.bl .gt. 0.0d0) then
      c = rimp/(al*bl)
    else
      c = 0.0d0
    end if

    return
  end subroutine impsc
!========================================================================================!
end module gfnff_geometry
