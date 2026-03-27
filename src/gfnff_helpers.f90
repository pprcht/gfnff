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

!> A collection of helper routines needed by gfnff
module gfnff_helpers
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  !> public from this module
  public :: lin
  public :: mrecgff,mrecgffPBC
  !public :: getring36,ssort
  public :: vlen,vsub,valijklff,valijklffPBC
  public :: omega,domegadr,dphidr,bangl,banglPBC,impsc
  public :: omegaPBC,domegadrPBC,dphidrPBC

  public :: crprod
  interface crprod
    module procedure crossprod
  end interface crprod
  public :: crossproduct

  public :: readl

  public :: bisectSearch
  interface bisectSearch
    module procedure :: bisectSearchReal
    module procedure :: bisectSearchInteger
  end interface bisectSearch

  public :: indexHeapSort

  real(wp),private,parameter :: eps = 1.0d-14
  real(wp),private,parameter :: pi = 3.1415926535897932384626433832795029d0

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine mrecgff(nat,nb,molcount,molvec)
    ! molcount: number of total fragments (increased during search)
    ! nat: overall number of atoms
    ! molvec: assignment vector of atom to fragment
    implicit none
    integer :: nat,molvec(nat),i,j,molcount,nb(20,nat)
    real(wp),allocatable :: bond(:,:)
    logical,allocatable :: taken(:)

    allocate (taken(nat),bond(nat,nat))
    bond = 0
    do i = 1,nat
      do j = 1,nb(20,i)
        bond(i,nb(j,i)) = 1
        bond(nb(j,i),i) = 1
      end do
    end do
    molvec = 0
    molcount = 1
    taken = .false.
    do i = 1,nat
      if (.not.taken(i)) then
        molvec(i) = molcount
        taken(i) = .true.
        call mrecgff2(nb,i,taken,nat,bond,molvec,molcount)
        molcount = molcount+1
      end if
    end do
    molcount = molcount-1
  end subroutine mrecgff
  recursive subroutine mrecgff2(nb,i,taken,nat,bond,molvec,molcnt)
    implicit none
    integer :: i,nat,molcnt,molvec(nat),j,icn,k,nb(20,nat)
    real(wp) :: bond(nat,nat)
    logical :: taken(nat)

    icn = nb(20,i)
    do k = 1,icn
      j = maxloc(bond(:,i),1)
      bond(j,i) = 0
      if (i .eq. j) cycle
      if (.not.taken(j)) then
        molvec(j) = molcnt
        taken(j) = .true.
        call mrecgff2(nb,j,taken,nat,bond,molvec,molcnt)
      end if
    end do
  end subroutine mrecgff2

!========================================================================================!

  subroutine mrecgffPBC(nat,numctr,numnb,nb,molcount,molvec)
    !**********************************************
    !* PBC-aware fragment search over the 3D
    !* neighbour array nb(numnb,nat,numctr).
    !* Atoms bonded across cell boundaries are
    !* placed in the same fragment.
    !**********************************************
    implicit none
    integer,intent(in)    :: nat,numctr,numnb,nb(numnb,nat,numctr)
    integer,intent(inout) :: molvec(nat),molcount
    integer :: i,j,iTr
    real(wp),allocatable :: bond(:,:,:)
    logical,allocatable  :: taken(:)

    allocate (taken(nat),bond(nat,nat,numctr))
    bond = 0.0_wp
    do i = 1,nat
      do iTr = 1,numctr
        do j = 1,nb(numnb,i,iTr)
          bond(nb(j,i,iTr),i,iTr) = 1.0_wp
        end do
      end do
    end do

    if (int(sum(bond)) .ne. sum(nb(numnb,:,:))) then
      write (*,*)
      write (*,'(a,2i10)') ' Warning (mrecgffPBC): bond sum mismatch', &
        & int(sum(bond)),sum(nb(numnb,:,:))
      write (*,*)
    end if

    molvec = 0
    molcount = 1
    taken = .false.
    do i = 1,nat
      if (.not.taken(i)) then
        molvec(i) = molcount
        taken(i) = .true.
        call mrecgff2PBC(numctr,numnb,nat,nb,i,taken,bond,molvec,molcount)
        molcount = molcount+1
      end if
    end do
    molcount = molcount-1
  end subroutine mrecgffPBC

  recursive subroutine mrecgff2PBC(numctr,numnb,nat,nb,i,taken,bond,molvec,molcnt)
    !**********************************************
    !* Recursive helper for mrecgffPBC.
    !**********************************************
    implicit none
    integer,intent(in)    :: numctr,numnb,nat,nb(numnb,nat,numctr),i,molcnt
    integer :: j,icn,k,iTr,j_iTr(2)
    real(wp),intent(inout) :: bond(nat,nat,numctr)
    integer,intent(inout)  :: molvec(nat)
    logical,intent(inout)  :: taken(nat)

    icn = sum(nb(numnb,i,:))
    do k = 1,icn
      j_iTr = maxloc(bond(:,i,:))
      j = j_iTr(1)
      iTr = j_iTr(2)
      bond(j,i,iTr) = 0.0_wp
      if (i .eq. j.and.iTr .eq. 1) cycle
      if (.not.taken(j)) then
        molvec(j) = molcnt
        taken(j) = .true.
        call mrecgff2PBC(numctr,numnb,nat,nb,j,taken,bond,molvec,molcnt)
      end if
    end do
  end subroutine mrecgff2PBC

! ══════════════════════════════════════════════════════════════════════════════

!  subroutine ssort(n,edum,ind)
!    implicit none
!    integer :: n,ii,k,j,i,sc1
!    real(wp) :: edum(n),pp
!    integer :: ind(n)
!
!    do ii = 2,n
!      i = ii-1
!      k = i
!      pp = edum(i)
!      do j = ii,n
!        if (edum(j) .gt. pp) exit
!        k = j
!        pp = edum(j)
!      end do
!      if (k .eq. i) exit
!      edum(k) = edum(i)
!      edum(i) = pp
!      sc1 = ind(i)
!      ind(i) = ind(k)
!      ind(k) = sc1
!    end do
!  end subroutine ssort

! ══════════════════════════════════════════════════════════════════════════════

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

    real(wp) :: xyz(3,natoms),vTrj(3),vTrk(3),vTrl(3),vDum(3)
    real(wp) :: ra(3),rb(3),rc(3),na(3),nb(3)
    real(wp) :: rab,rbc,thab,thbc
    real(wp) :: nan,nbn,rcn,snanb,deter

    vDum = 0.0 ! Dummy vector for valijkPBC function

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

    !> determinante of rb,ra,rc  (triple product)
    deter = ra(1)*(rb(2)*rc(3)-rb(3)*rc(2)) &
    &      -ra(2)*(rb(1)*rc(3)-rb(3)*rc(1)) &
    &      +ra(3)*(rb(1)*rc(2)-rb(2)*rc(1))

    if (mo .eq. 1) then ! not used
      thab = valijkPBC(1,natoms,xyz,i,k,j,vTrl,vTrj,vDum)
      thbc = valijkPBC(2,natoms,xyz,j,l,k,vTrk,vTrj,vDum)
    else
      thab = valijkPBC(3,natoms,xyz,i,k,j,vTrj,vTrl,vTrk) !    i=R       k=C       j=B
      ! vTrj=vTrR vTrl=vTrC vTrk=vTrB
      thbc = valijkPBC(4,natoms,xyz,j,l,k,vTrk,vTrl,vDum) !    j=B    l=H    k=C
      ! vTrk=vTrB     vTrl=vTrC
    end if
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
       &        rkjn,rljn,rnn,rvn
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
    &           sinphi,cosphi,onenner,thab,thbc,&
    &           ra(3),rb(3),rc(3),rab(3),rac(3),rbc(3),rbb(3),&
    &           raa(3),rba(3),rapba(3),rapbb(3),rbpca(3),rbpcb(3),&
    &           rapb(3),rbpc(3),na(3),nb(3),nan,nbn,&
    &           dphidri(3),dphidrj(3),dphidrk(3),dphidrl(3),&
    &           xyz(3,nat),phi,nenner,vz

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
    trV  = transVec(:,iTr)   
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
  subroutine readl(a1,x,n)
    character(len=*) ::  a1
    real(wp) :: x(:)
    integer,intent(out) :: n
    n = size(x,1) !will be overwrotten in getfloats
    call getfloats(a1,x,n)
  end subroutine readl
!========================================================================================!
!cuts the at blanks and tabstops and returns all floats and strings in order of occurence
  subroutine getfloats(line,floats,cf)
    implicit none
    real(wp),intent(inout) :: floats(*)
    character(len=*),intent(in) :: line
    integer,intent(inout) :: cf
    real(wp) :: num
    character(len=128) :: str,stmp
    character(len=80) strings(3)
    character(len=1) digit
    integer :: i,ty,cs,cfmax

    cfmax = cf !on input cf is the dimension of floats()
    stmp = ''
    cs = 0
    cf = 0
    strings(:) = ''
    do i = 1,len(trim(line))
      digit = line(i:i)
      if (digit .ne. ' '.and.digit .ne. char(9)) then  !should exclude tabstops and blanks, 9 is ascii code for tab
        stmp = trim(stmp)//trim(digit)
      elseif (stmp .ne. '') then
        call checktype(stmp,num,str,ty)      !get type of string, 0=number, 1=character
        if (ty .eq. 0) then
          cf = cf+1
          if (cf .le. cfmax) floats(cf) = num
        elseif (ty .eq. 1) then
          cs = cs+1
          !   strings(cs)=trim(str)
        else
          write (*,*) 'Problem in checktype, must abort'
          exit
        end if
        stmp = ''
      end if
      if (i .eq. len(trim(line))) then  !special case: end of line
        call checktype(stmp,num,str,ty)
        if (ty .eq. 0) then
          cf = cf+1
          if (cf .le. cfmax) floats(cf) = num
        elseif (ty .eq. 1) then
          cs = cs+1
          !   strings(cs)=trim(str)
        else
          write (*,*) 'Problem in checktype, must abort'
          exit
        end if
        stmp = ''
      end if
    end do
  contains
    subroutine checktype(field,num,str,ty)
      implicit none
      character(len=*) :: field,str
      real(wp) :: num
      integer :: e,ty
      logical :: is_num
      ty = 99
      str = ''
      is_num = .false.
      read (field,'(F10.5)',IOSTAT=e) num !cast string on real and get error code; 0 means success.
      if (e .eq. 0) is_num = .true.
      if (is_num) then
        if (index(field,'.') .ne. 0) then  !check for integer/real
          read (field,'(F30.16)') num
          ty = 0
        else                       !if integer, add .0 to string; otherwise cast to real does not work
          str = trim(field)//'.0'
          read (str,'(F30.16)') num
          str = ''
          ty = 0
        end if
      else
        str = trim(field)
        ty = 1
      end if
    end subroutine checktype
  end subroutine getfloats

!================================================================================!

!> Integer case for bisection search
  pure subroutine bisectSearchInteger(j,xx,x)

    !> Located element such that xx(j) <= x < xx(j+1)
    integer,intent(out) :: j

    !> Array of values in monotonic order to search through
    integer,intent(in) :: xx(:)

    !> Value to locate j for
    integer,intent(in) :: x

    integer :: n
    integer :: jlower,jupper,jcurr

    n = size(xx)
    if (n == 0) then
      j = 0
      return
    end if

    if (x < xx(1)) then
      j = 0
    else if (x == xx(1)) then
      j = 1
    else if (x == xx(n)) then
      j = n-1
    else if (x > xx(n)) then
      j = n
    else
      jlower = 0
      jcurr = n+1
      do while ((jcurr-jlower) > 1)
        jupper = (jcurr+jlower)/2
        if ((xx(n) >= xx(1)).eqv.(x >= xx(jupper))) then
          jlower = jupper
        else
          jcurr = jupper
        end if
      end do
      j = jlower
    end if

  end subroutine bisectSearchInteger

!> Real case for bisection search
  pure subroutine bisectSearchReal(j,xx,x,tol)

    !> Located element such that xx(j) <= x < xx(j+1)
    integer,intent(out) :: j

    !> Array of values in monotonic order to search through
    real(wp),intent(in) :: xx(:)

    !> Value to locate j for
    real(wp),intent(in) :: x

    !> Tolerance for equality comparision
    real(wp),intent(in),optional :: tol

    integer :: n
    integer :: jlower,jupper,jcurr
    real(wp) :: rTol
    logical :: ascending

    n = size(xx)
    if (n == 0) then
      j = 0
      return
    end if

    if (present(tol)) then
      rTol = tol
    else
      rTol = epsilon(0.0_wp)
    end if

    if (x < xx(1)-rTol) then
      j = 0
    else if (abs(x-xx(1)) <= rTol) then
      j = 1
    else if (abs(x-xx(n)) <= rTol) then
      j = n-1
    else if (x > xx(n)+rTol) then
      j = n
    else
      ascending = (xx(n) >= xx(1))
      jlower = 0
      jcurr = n+1
      do while ((jcurr-jlower) > 1)
        jupper = (jcurr+jlower)/2
        if (ascending.eqv.(x >= xx(jupper)+rTol)) then
          jlower = jupper
        else
          jcurr = jupper
        end if
      end do
      j = jlower
    end if

  end subroutine bisectSearchReal

! ──────────────────────────────────────────────────────────────────────────────

!> Real case heap sort returning an index.
!  Based on Numerical Recipes Software 1986-92
  pure subroutine indexHeapSort(indx,array,tolerance)
    !> Indexing array on return
    integer,intent(out) :: indx(:)
    !> Array of values to be sorted
    real(wp),intent(in) :: array(:)
    !> Tolerance for equality of two elements
    real(wp),intent(in),optional :: tolerance

    integer :: n,ir,ij,il,ii
    integer :: indxTmp
    real(wp) :: arrayTmp,tol

    !:ASSERT(size(array)==size(indx))

    if (present(tolerance)) then
      tol = tolerance
    else
      tol = epsilon(0.0_wp)
    end if

    do ii = 1,size(indx)
      indx(ii) = ii
    end do
    n = size(array)
    if (n <= 1) return
    il = n/2+1
    ir = n
    do
      if (il > 1) then
        il = il-1
        indxTmp = indx(il)
        arrayTmp = array(indxTmp)
      else
        indxTmp = indx(ir)
        arrayTmp = array(indxTmp)
        indx(ir) = indx(1)
        ir = ir-1
        if (ir < 1) then
          indx(1) = indxTmp
          return
        end if
      end if
      ii = il
      ij = 2*il
      do while (ij <= ir)
        if (ij < ir) then
          if (array(indx(ij)) < array(indx(ij+1))-tol) then
            ij = ij+1
          end if
        end if
        if (arrayTmp < array(indx(ij))-tol) then
          indx(ii) = indx(ij)
          ii = ij
          ij = 2*ij
        else
          ij = ir+1
        end if
      end do
      indx(ii) = indxTmp
    end do

  end subroutine indexHeapSort

!========================================================================================!
end module gfnff_helpers
