!================================================================================!
! This file is part of gfnff.
!
! Copyright (C) 2023 Philipp Pracht
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
!--------------------------------------------------------------------------------!
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert, Sebastian Spicher, Stefan Grimme
!> at https://github.com/grimme-lab/xtb
!================================================================================!
!> A collection of helper routines needed by gfnff
module gfnff_helpers
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  !> public from this module
  public :: lin
  public :: mrecgff
  public :: getring36,ssort
  public :: vlen,vsub,valijklff
  public :: omega,domegadr,dphidr,bangl,impsc

  public :: crprod
  interface crprod
    module procedure crossprod
  end interface crprod

  public :: readl

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
  subroutine getring36(n,at,nbin,a0_in,cout,irout)
    implicit none
    integer :: cout(10,20),irout(20)  ! output: atomlist, ringsize, # of rings in irout(20)
    integer :: at(n)
    integer :: n,nbin(20,n),a0,i,nb(20,n),a0_in
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
    integer :: n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
    integer :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
    integer :: list(n),m,mm,nn,c(10),cdum(10,600),iring
    integer :: adum1(0:n),adum2(0:n),kk,j,idum(600),same(600),k
    real(wp) :: w(n),av,sd

    if (n .le. 2) return

    nn = nbin(20,a0_in)

    cdum = 0
    kk = 0
    do m = 1,nn
!     if(nb(m,a0_in).eq.1) cycle
      nb = nbin
      do i = 1,n
        if (nb(20,i) .eq. 1) nb(20,i) = 0
      end do

      do mm = 1,nn
        w(mm) = dble(mm)
        list(mm) = mm
      end do
      w(m) = 0.0d0
      call ssort(nn,w,list)
      do mm = 1,nn
        nb(mm,a0_in) = nbin(list(mm),a0_in)
      end do

      iring = 0
      c = 0

      a0 = a0_in
      n0 = nb(20,a0)

      do i1 = 1,n0
        a1 = nb(i1,a0)
        if (a1 .eq. a0) cycle
        n1 = nb(20,a1)
        do i2 = 1,n1
          a2 = nb(i2,a1)
          if (a2 .eq. a1) cycle
          n2 = nb(20,a2)
          do i3 = 1,n2
            a3 = nb(i3,a2)
            n3 = nb(20,a3)
            if (a3 .eq. a2) cycle
            c(1) = a1
            c(2) = a2
            c(3) = a3
            if (a3 .eq. a0.and.chkrng(n,3,c)) then
              iring = 3
              kk = kk+1
              cdum(1:iring,kk) = c(1:iring)
              idum(kk) = iring
            end if
            do i4 = 1,n3
              a4 = nb(i4,a3)
              n4 = nb(20,a4)
              if (a4 .eq. a3) cycle
              c(4) = a4
              if (a4 .eq. a0.and.chkrng(n,4,c)) then
                iring = 4
                kk = kk+1
                cdum(1:iring,kk) = c(1:iring)
                idum(kk) = iring
              end if
              do i5 = 1,n4
                a5 = nb(i5,a4)
                n5 = nb(20,a5)
                if (a5 .eq. a4) cycle
                c(5) = a5
                if (a5 .eq. a0.and.chkrng(n,5,c)) then
                  iring = 5
                  kk = kk+1
                  cdum(1:iring,kk) = c(1:iring)
                  idum(kk) = iring
                end if
                do i6 = 1,n5
                  a6 = nb(i6,a5)
                  n6 = nb(20,a6)
                  if (a6 .eq. a5) cycle
                  c(6) = a6
                  if (a6 .eq. a0.and.chkrng(n,6,c)) then
                    iring = 6
                    kk = kk+1
                    cdum(1:iring,kk) = c(1:iring)
                    idum(kk) = iring
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
    end do

! compare
    same = 0
    do i = 1,kk
      do j = i+1,kk
        if (idum(i) .ne. idum(j)) cycle ! different ring size
        if (same(j) .eq. 1) cycle ! already double
        adum1 = 0
        adum2 = 0
        do m = 1,10
          i1 = cdum(m,i)
          i2 = cdum(m,j)
          adum1(i1) = 1
          adum2(i2) = 1
        end do
        if (sum(abs(adum1-adum2)) .ne. 0) then
          same(j) = 0
        else
          same(j) = 1
        end if
      end do
    end do

    m = 0
    do i = 1,kk
      if (same(i) .eq. 0) then
        m = m+1
        if (m .gt. 20) exit !stop 'too many rings'
        irout(m) = idum(i)     ! number of atoms in ring m
        nn = idum(i)
        cout(1:nn,m) = cdum(1:nn,i)
        i2 = 0
        do k = 1,nn            ! determine if its a hetereo
          i1 = at(cdum(k,i))
          i2 = i2+i1
        end do
        av = dble(i2)/dble(nn)
        sd = 0
        cout(m,19) = 0
        do k = 1,nn
          i1 = at(cdum(k,i))
          sd = sd+(av-dble(i1))**2
        end do
        if (sd .gt. 1.d-6) cout(m,19) = idint(1000.*sqrt(sd)/dble(nn))
      end if
    end do
    irout(20) = m  ! number of rings for this atom

    return
  end subroutine getring36

  logical function chkrng(nn,n,c)
    implicit none
    integer n,idum(nn),nn,c(10),i,j
    chkrng = .true.
    idum = 0
    do i = 1,n
      idum(c(i)) = idum(c(i))+1
    end do
    j = 0
    do i = 1,nn
      if (idum(i) .eq. 1) j = j+1
    end do
    if (j .ne. n) chkrng = .false.
  end function chkrng

  subroutine ssort(n,edum,ind)
    implicit none
    integer :: n,ii,k,j,m,i,sc1
    real(wp) :: edum(n),pp
    integer :: ind(n)

    do ii = 2,n
      i = ii-1
      k = i
      pp = edum(i)
      do j = ii,n
        if (edum(j) .gt. pp) exit
        k = j
        pp = edum(j)
      end do
      if (k .eq. i) exit
      edum(k) = edum(i)
      edum(i) = pp
      sc1 = ind(i)
      ind(i) = ind(k)
      ind(k) = sc1
    end do
  end subroutine ssort
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

!========================================================================================!

  real(wp) function vlen(a)
    implicit none !double precision (a-h,o-z)
    real(wp) :: a(3)
    real(wp) :: tot

    tot = a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
    vlen = 0.0d0
    if (tot .gt. 0.0d0) vlen = dsqrt(tot)

    return
  end function vlen

!========================================================================================!

  real(wp) function valijklff(natoms,xyz,i,j,k,l)
    implicit none
    integer :: ic,i,j,k,l,natoms
    real(wp) :: xyz(3,natoms)
    real(wp) :: ra(3),rb(3),rc(3),na(3),nb(3)
    real(wp) :: rab,rbc,thab,thbc
    real(wp) :: nan,nbn,rcn,snanb,deter

    real(wp),parameter :: eps = 1.0d-14
    real(wp),parameter :: pi = 3.1415926535897932384626433832795029d0

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

  subroutine crossprod(ra,rb,rab)
    implicit none
    real(wp) :: ra(3),rb(3),rab(3)
    rab(1) = ra(2)*rb(3)-ra(3)*rb(2)
    rab(2) = ra(3)*rb(1)-ra(1)*rb(3)
    rab(3) = ra(1)*rb(2)-ra(2)*rb(1)
  end Subroutine crossprod

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
    real(wp) :: rkjn,rljn,rnn,rvn

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

!========================================================================================!

  Subroutine dphidr(nat,xyz,i,j,k,l,phi, &
  &                dphidri,dphidrj,dphidrk,dphidrl)
    !> the torsion derivatives

    implicit none
    integer :: ic,i,j,k,l,nat
    real(wp) :: sinphi,cosphi,onenner,thab,thbc
    real(wp) :: ra(3),rb(3),rc(3),rab(3),rac(3),rbc(3),rbb(3)
    real(wp) :: raa(3),rba(3),rapba(3),rapbb(3),rbpca(3),rbpcb(3)
    real(wp) :: rapb(3),rbpc(3),na(3),nb(3),nan,nbn
    real(wp) :: dphidri(3),dphidrj(3),dphidrk(3),dphidrl(3)
    real(wp) :: xyz(3,nat),phi,nenner,eps,vz

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

!========================================================================================!

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
!========================================================================================!
end module gfnff_helpers
