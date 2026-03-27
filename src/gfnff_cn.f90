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
module gfnff_cn
  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData
  use gfnff_param,only: covalentRadD3,paulingEN
  implicit none
  private
  public :: gfnff_dlogcoord

  public :: cnType,getCoordinationNumber

  interface getCoordinationNumber
    module procedure :: getCoordinationNumberLP
  end interface getCoordinationNumber

  !> Possible counting functions for calculating coordination numbers
  type :: TCNTypeEnum

    !> Counting function not specified
    integer :: invalid = 0

    !> Original DFT-D3 coordination number
    integer :: exp = 1

    !> Faster decaying error function CN, better for dense systems
    integer :: erf = 2

    !> Error function CN with covalency correction
    integer :: cov = 3

    !> Particular long-ranged version of the DFT-D3 coordination number
    integer :: gfn = 4

    !> log CN, erf fct capped at max value cnmax
    integer :: log = 5

  end type TCNTypeEnum

  !> Enumerator for different coordination number types
  type(TCNTypeEnum),parameter :: cnType = TCNTypeEnum()

  abstract interface
    !> Abstract interface for the counting function (and its derivative)
    pure function countingFunction(k,r,r0)
      import :: wp

      !> Constant for counting function
      real(wp),intent(in) :: k

      !> Actual distance
      real(wp),intent(in) :: r

      !> Critical distance
      real(wp),intent(in) :: r0

      !> Value of the counting function in the range of [0,1]
      real(wp) :: countingFunction

    end function countingFunction
  end interface

  !> Parameter for electronegativity scaling
  real(wp),parameter :: k4 = 4.10451_wp

  !> Parameter for electronegativity scaling
  real(wp),parameter :: k5 = 19.08857_wp

  !> Parameter for electronegativity scaling
  real(wp),parameter :: k6 = 2*11.28174_wp**2

  !> pi
  real(wp),parameter :: pi = 3.1415926535897932385_wp 


!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CN routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> logCN derivative saved in dlogCN array
  subroutine gfnff_dlogcoord(n,at,xyz,rab,logCN,dlogCN,thr2,param)
    implicit none
    type(TGFFData),intent(in) :: param
    integer,intent(in)  :: n
    integer,intent(in)  :: at(n)
    real(wp),intent(in)  :: xyz(3,n)
    real(wp),intent(in)  :: rab(n*(n+1)/2)
    real(wp),intent(out) :: logCN(n)
    real(wp),intent(out) :: dlogCN(3,n,n)
    real(wp),intent(in)  :: thr2
    real(wp)             :: cn(n)
    !> counter
    integer  :: i,j,ij,ii
    !> distances
    real(wp) :: r
    real(wp) :: r0
    real(wp) :: dr
    real(wp),dimension(3) :: rij
    !> treshold
    real(wp) :: thr
    !> cn parameter
    real(wp) :: erfCN
    real(wp) :: dlogdcni
    real(wp) :: dlogdcnj
    real(wp) :: derivative
    !> local parameter
    real(wp),parameter :: kn = -7.5_wp
    !> initialize vaiiables
    cn = 0.0_wp
    logCN = 0.0_wp
    dlogCN = 0.0_wp
    dlogdcni = 0.0_wp
    dlogdcnj = 0.0_wp
    derivative = 0.0_wp
    r = 0.0_wp
    r0 = 0.0_wp
    rij = 0.0_wp
    thr = sqrt(thr2)
    !> create error function CN
    do i = 2,n
      ii = i*(i-1)/2
      do j = 1,i-1
        ij = ii+j
        r = rab(ij)
        if (r .gt. thr) cycle
        r0 = (param%rcov(at(i))+param%rcov(at(j)))
        dr = (r-r0)/r0
        !> hier kommt die CN funktion hin
!           erfCN = create_expCN(16.0d0,r,r0)
        erfCN = 0.5_wp*(1.0_wp+erf(kn*dr))
        cn(i) = cn(i)+erfCN
        cn(j) = cn(j)+erfCN
      end do
    end do
    !> create cutted logarithm CN + derivatives
    do i = 1,n
      ii = i*(i-1)/2
      logCN(i) = create_logCN(cn(i),param)
      !> get dlogCN/dCNi
      dlogdcni = create_dlogCN(cn(i),param)
      do j = 1,i-1
        ij = ii+j
        !> get dlogCN/dCNj
        dlogdcnj = create_dlogCN(cn(j),param)
        r = rab(ij)
        if (r .gt. thr) cycle
        r0 = (param%rcov(at(i))+param%rcov(at(j)))
        !> get derfCN/dRij
        derivative = create_derfCN(kn,r,r0)
!           derivative = create_dexpCN(16.0d0,r,r0)
        rij = derivative*(xyz(:,j)-xyz(:,i))/r
        !> project rij gradient onto cartesians
        dlogCN(:,j,j) = dlogdcnj*rij+dlogCN(:,j,j)
        dlogCN(:,i,j) = -dlogdcnj*rij
        dlogCN(:,j,i) = dlogdcni*rij
        dlogCN(:,i,i) = -dlogdcni*rij+dlogCN(:,i,i)
      end do
    end do

  contains

    pure elemental function create_logCN(cn,param) result(count)
      type(TGFFData),intent(in) :: param
      real(wp),intent(in) :: cn
      real(wp) :: count
      count = log(1+exp(param%cnmax))-log(1+exp(param%cnmax-cn))
    end function create_logCN

    pure elemental function create_dlogCN(cn,param) result(count)
      type(TGFFData),intent(in) :: param
      real(wp),intent(in) :: cn
      real(wp) :: count
      count = exp(param%cnmax)/(exp(param%cnmax)+exp(cn))
    end function create_dlogCN

!    pure elemental function create_erfCN(k,r,r0) result(count)
!      real(wp),intent(in) :: k
!      real(wp),intent(in) :: r
!      real(wp),intent(in) :: r0
!      real(wp) :: count
!      real(wp) :: dr
!      dr = (r-r0)/r0
!      count = 0.5_wp*(1.0_wp+erf(kn*dr))
!    end function create_erfCN

    pure elemental function create_derfCN(k,r,r0) result(count)
      real(wp),intent(in) :: k
      real(wp),intent(in) :: r
      real(wp),intent(in) :: r0
      real(wp),parameter :: sqrtpi = 1.77245385091_wp
      real(wp) :: count
      real(wp) :: dr
      dr = (r-r0)/r0
      count = k/sqrtpi*exp(-k**2*dr*dr)/r0
    end function create_derfCN

!    pure elemental function create_expCN(k,r,r0) result(count)
!      real(wp),intent(in) :: k
!      real(wp),intent(in) :: r
!      real(wp),intent(in) :: r0
!      real(wp) :: count
!      count = 1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))
!    end function create_expCN

!    pure elemental function create_dexpCN(k,r,r0) result(count)
!      real(wp),intent(in) :: k
!      real(wp),intent(in) :: r
!      real(wp),intent(in) :: r0
!      real(wp) :: count
!      real(wp) :: expterm
!      expterm = exp(-k*(r0/r-1._wp))
!      count = (-k*r0*expterm)/(r**2*((expterm+1._wp)**2))
!    end function create_dexpCN

  end subroutine gfnff_dlogcoord

! ──────────────────────────────────────────────────────────────────────────────

!> Geometric fractional coordination number, supports both error function
!  and exponential counting functions.
  subroutine getCoordinationNumberLP(nat,at,xyz,ntrans,trans,cutoff,cf,cn,dcndr,dcndL,param)

    !> Molecular structure information
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)

    !> Number of lattice points
    integer,intent(in) :: ntrans

    !> Lattice points
    real(wp),intent(in) :: trans(:,:)

    !> Real space cutoff
    real(wp),intent(in) :: cutoff

    !> Coordination number type (by counting function)
    integer,intent(in) :: cf

    ! parameter type, holds cnmax
    type(TGFFData),intent(in) :: param

    !> Error function coordination number
    real(wp),intent(out) :: cn(:)

    !> Derivative of the CN with respect to the Cartesian coordinates
    real(wp),intent(out) :: dcndr(:,:,:)

    !> Derivative of the CN with respect to strain deformations
    real(wp),intent(out) :: dcndL(:,:,:)

    real(wp),parameter :: kcn_exp = 16.0_wp
    real(wp),parameter :: kcn_erf = 7.5_wp
    real(wp),parameter :: kcn_gfn = 10.0_wp

    select case (cf)
    case (cnType%exp)
      call ncoordLatP(nat,at,xyz,ntrans,trans,cutoff,kcn_exp,expCount,dexpCount, &
         & .false.,covalentRadD3,paulingEN,cn,dcndr,dcndL)
    case (cnType%erf)
      call ncoordLatP(nat,at,xyz,ntrans,trans,cutoff,kcn_erf,erfCount,derfCount, &
         & .false.,covalentRadD3,paulingEN,cn,dcndr,dcndL)
    case (cnType%cov)
      call ncoordLatP(nat,at,xyz,ntrans,trans,cutoff,kcn_erf,erfCount,derfCount, &
         & .true.,covalentRadD3,paulingEN,cn,dcndr,dcndL)
    case (cnType%gfn)
      call ncoordLatP(nat,at,xyz,ntrans,trans,cutoff,kcn_gfn,gfnCount,dgfnCount, &
         & .false.,covalentRadD3,paulingEN,cn,dcndr,dcndL)
    case (cnType%log)
      call ncoordLatP(nat,at,xyz,ntrans,trans,cutoff,kcn_erf,erfCount,derfCount, &
         & .false.,covalentRadD3,paulingEN,cn,dcndr,dcndL)
      ! calculate modified CN and derivatives
      call cutCoordinationNumber(nat,cn,dcndr,dcndL,param%cnmax)

    end select

  end subroutine getCoordinationNumberLP

!> Actual implementation of the coordination number, takes a generic counting
!  function to return the respective CN.
  subroutine ncoordLatP(nat,at,xyz,ntrans,trans,cutoff,kcn,cfunc,dfunc,enscale, &
        & rcov,en,cn,dcndr,dcndL)

    !> Molecular structure information
    integer,intent(in) :: nat,at(nat) 
    real(wp),intent(in) :: xyz(3,nat) 

    !> Number of lattice points
    integer,intent(in) :: ntrans

    !> Lattice points
    real(wp),intent(in) :: trans(:,:)

    !> Real space cutoff
    real(wp),intent(in) :: cutoff

    !> Function implementing the counting function
    procedure(countingFunction) :: cfunc

    !> Function implementing the derivative of counting function w.r.t. distance
    procedure(countingFunction) :: dfunc

    !> Use a covalency criterium by Pauling EN's
    logical,intent(in) :: enscale

    !> Steepness of counting function
    real(wp),intent(in) :: kcn

    !> Covalent radius
    real(wp),intent(in) :: rcov(:)

    !> Electronegativity
    real(wp),intent(in) :: en(:)

    !> Error function coordination number.
    real(wp),intent(out) :: cn(:)

    !> Derivative of the CN with respect to the Cartesian coordinates.
    real(wp),intent(out) :: dcndr(:,:,:)

    !> Derivative of the CN with respect to strain deformations.
    real(wp),intent(out) :: dcndL(:,:,:)

    integer :: iat,jat,ati,atj,itr
    real(wp) :: r2,r1,rc,rij(3),countf,countd(3),stress(3,3),den,cutoff2

    cn = 0.0_wp
    dcndr = 0.0_wp
    dcndL = 0.0_wp
    cutoff2 = cutoff**2
    !$omp parallel do default(none) private(den) shared(enscale, rcov, en)&
    !$omp reduction(+:cn, dcndr, dcndL) shared(nat,at,xyz, kcn, trans, cutoff2, ntrans) &
    !$omp private(jat, itr, ati, atj, r2, rij, r1, rc, countf, countd, stress)
    do iat = 1,nat
      ati = at(iat)
      do jat = 1,iat
        atj = at(jat)

        if (enscale) then
          den = k4*exp(-(abs(en(ati)-en(atj))+k5)**2/k6)
        else
          den = 1.0_wp
        end if

        do itr = 1,ntrans
          rij = xyz(:,iat)-(xyz(:,jat)+trans(:,itr))
          r2 = sum(rij**2)
          if (r2 > cutoff2.or.r2 < 1.0e-12_wp) cycle
          r1 = sqrt(r2)

          rc = rcov(ati)+rcov(atj)

          countf = den*cfunc(kcn,r1,rc)
          countd = den*dfunc(kcn,r1,rc)*rij/r1

          cn(iat) = cn(iat)+countf
          if (iat .ne. jat.or.itr .ne. 1) then
            cn(jat) = cn(jat)+countf
          end if

          dcndr(:,iat,iat) = dcndr(:,iat,iat)+countd
          dcndr(:,jat,jat) = dcndr(:,jat,jat)-countd
          dcndr(:,iat,jat) = dcndr(:,iat,jat)+countd
          dcndr(:,jat,iat) = dcndr(:,jat,iat)-countd

          stress(:,1) = countd(1)*rij
          stress(:,2) = countd(2)*rij
          stress(:,3) = countd(3)*rij

          dcndL(:,:,iat) = dcndL(:,:,iat)+stress
          if (iat .ne. jat.or.itr .ne. 1) then
            dcndL(:,:,jat) = dcndL(:,:,jat)+stress
          end if
        end do
      end do
    end do

    !$omp end parallel do

  end subroutine ncoordLatP

!> Error function counting function for coordination number contributions.
  pure function erfCount(k,r,r0) result(count)

    !> Steepness of the counting function.
    real(wp),intent(in) :: k

    !> Current distance.
    real(wp),intent(in) :: r

    !> Cutoff radius.
    real(wp),intent(in) :: r0

    real(wp) :: count

    count = 0.5_wp*(1.0_wp+erf(-k*(r-r0)/r0))

  end function erfCount

!> Derivative of the counting function w.r.t. the distance.
  pure function derfCount(k,r,r0) result(count)

    !> Steepness of the counting function.
    real(wp),intent(in) :: k

    !> Current distance.
    real(wp),intent(in) :: r

    !> Cutoff radius.
    real(wp),intent(in) :: r0

    real(wp),parameter :: sqrtpi = sqrt(pi)

    real(wp) :: count

    count = -k/sqrtpi/r0*exp(-k**2*(r-r0)**2/r0**2)

  end function derfCount

!> Exponential counting function for coordination number contributions.
  pure function expCount(k,r,r0) result(count)

    !> Steepness of the counting function.
    real(wp),intent(in) :: k

    !> Current distance.
    real(wp),intent(in) :: r

    !> Cutoff radius.
    real(wp),intent(in) :: r0

    real(wp) :: count

    count = 1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))

  end function expCount

!> Derivative of the counting function w.r.t. the distance.
  pure function dexpCount(k,r,r0) result(count)

    !> Steepness of the counting function.
    real(wp),intent(in) :: k

    !> Current distance.
    real(wp),intent(in) :: r

    !> Cutoff radius.
    real(wp),intent(in) :: r0

    real(wp) :: count
    real(wp) :: expterm

    expterm = exp(-k*(r0/r-1._wp))

    count = (-k*r0*expterm)/(r**2*((expterm+1._wp)**2))

  end function dexpCount

!> Exponential counting function for coordination number contributions.
  pure function gfnCount(k,r,r0) result(count)

    !> Steepness of the counting function.
    real(wp),intent(in) :: k

    !> Current distance.
    real(wp),intent(in) :: r

    !> Cutoff radius.
    real(wp),intent(in) :: r0

    real(wp) :: count

    count = expCount(k,r,r0)*expCount(2*k,r,r0+2)

  end function gfnCount

!> Derivative of the counting function w.r.t. the distance.
  pure function dgfnCount(k,r,r0) result(count)

    !> Steepness of the counting function.
    real(wp),intent(in) :: k

    !> Current distance.
    real(wp),intent(in) :: r

    !> Cutoff radius.
    real(wp),intent(in) :: r0

    real(wp) :: count

    count = dexpCount(k,r,r0)*expCount(2*k,r,r0+2) &
       &  +expCount(k,r,r0)*dexpCount(2*k,r,r0+2)

  end function dgfnCount

!> Cutoff function for large coordination numbers
  pure subroutine cutCoordinationNumber(nAtom,cn,dcndr,dcndL,maxCN)

    !> number of atoms
    integer,intent(in) :: nAtom

    !> on input coordination number, on output modified CN
    real(wp),intent(inout) :: cn(:)

    !> on input derivative of CN w.r.t. cartesian coordinates,
    !> on output derivative of modified CN
    real(wp),intent(inout),optional :: dcndr(:,:,:)

    !> on input derivative of CN w.r.t. strain deformation,
    !> on output derivative of modified CN
    real(wp),intent(inout),optional :: dcndL(:,:,:)

    !> maximum CN (not strictly obeyed)
    real(wp),intent(in),optional :: maxCN

    real(wp) :: cnmax
    integer :: iAt

    if (present(maxCN)) then
      cnmax = maxCN
    else
      cnmax = 4.5_wp
    end if

    if (cnmax <= 0.0_wp) return

    if (present(dcndL)) then
      do iAt = 1,nAtom
        dcndL(:,:,iAt) = dcndL(:,:,iAt)*dCutCN(cn(iAt),cnmax)
      end do
    end if

    if (present(dcndr)) then
      do iAt = 1,nAtom
        dcndr(:,:,iAt) = dcndr(:,:,iAt)*dCutCN(cn(iAt),cnmax)
      end do
    end if

    do iAt = 1,nAtom
      cn(iAt) = cutCN(cn(iAt),cnmax)
    end do

  end subroutine cutCoordinationNumber

!> Cutting function for the coordination number.
  elemental function cutCN(cn,cut) result(cnp)

    !> Current coordination number.
    real(wp),intent(in) :: cn

    !> Cutoff for the CN, this is not the maximum value.
    real(wp),intent(in) :: cut

    !> Cuting function vlaue
    real(wp) :: cnp

    cnp = log(1.0_wp+exp(cut))-log(1.0_wp+exp(cut-cn))

  end function cutCN

!> Derivative of the cutting function w.r.t. coordination number
  elemental function dCutCN(cn,cut) result(dcnpdcn)

    !> Current coordination number.
    real(wp),intent(in) :: cn

    !> Cutoff for the CN, this is not the maximum value.
    real(wp),intent(in) :: cut

    !> Derivative of the cutting function
    real(wp) :: dcnpdcn

    dcnpdcn = exp(cut)/(exp(cut)+exp(cn))

  end function dCutCn

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

end module gfnff_cn

