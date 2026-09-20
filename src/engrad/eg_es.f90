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
!> GFN-FF electrostatics: EEQ charges and energy from a direct (molecular) or
!> Ewald-summed (periodic) solve, the ES gradient and the periodic stress tensor.
module gfnff_eg_es

  use iso_fortran_env,only:wp => real64,sp => real32,stdout => output_unit
  use gfnff_data_types,only:TGFFData,TGFFNeighbourList,TGFFTopology,TCell
  use gfnff_wsc,only:wsc_images,wsc_maxcells
  use gfnff_alpb,only:TBorn
  use gfnff_math_wrapper,only:dot,symv,sytrf_wrap,sytrs_wrap
  implicit none
  private

  public :: goed_gfnff,goed_pbc_gfnff
  public :: es_grad_mol,es_grad_sigma

  real(wp),private,parameter :: pi = 3.1415926535897932385_wp
  real(wp),private,parameter :: sqrtpi = 1.77245385091_wp

contains  !> MODULE PROCEDURES START HERE

! Ref.: S. Alireza Ghasemi, Albert Hofstetter, Santanu Saha, and Stefan Goedecker
!       PHYSICAL REVIEW B 92, 045131 (2015)
!       Interatomic potentials for ionic systems with density functional accuracy
!       based on charge densities obtained by a neural network

  subroutine es_grad_mol(n,xyz,sqrab,srab,eeqtmp,q,g)
    !***********************************************************************
    !* Gradient of the molecular EEQ electrostatic energy at fixed charges,
    !* added to g. The charge response enters later through the CN chain rule.
    !* eeqtmp holds the packed gamma_ij and erf(gamma_ij*r_ij) from goed_gfnff.
    !***********************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: sqrab(n*(n+1)/2),srab(n*(n+1)/2)
    real(wp),intent(in) :: eeqtmp(2,n*(n+1)/2)
    real(wp),intent(in) :: q(n)
    real(wp),intent(inout) :: g(3,n)

    integer :: i,j,k,ij
    real(wp) :: r2,rab,gammij,erff,dd,r3(3)

    !$omp parallel do default(none) reduction (+:g) &
    !$omp shared(q,n,sqrab,srab,eeqtmp,xyz) &
    !$omp private(i,j,k,ij,r3,r2,rab,gammij,erff,dd)
    do i = 1,n
      k = i*(i-1)/2
      do j = 1,i-1
        ij = k+j
        r2 = sqrab(ij)
        rab = srab(ij)
        gammij = eeqtmp(1,ij)
        erff = eeqtmp(2,ij)
        dd = (2.0d0*gammij*exp(-gammij**2*r2) &
                 & /(sqrtpi*r2)-erff/(rab*r2))*q(i)*q(j)
        r3 = (xyz(:,i)-xyz(:,j))*dd
        g(:,i) = g(:,i)+r3
        g(:,j) = g(:,j)-r3
      end do
    end do
    !$omp end parallel do

  end subroutine es_grad_mol

  subroutine goed_gfnff(single,n,at,sqrab,r,chrg,eeqtmp,cn,q,es,gbsa,param,topo)
    !***********************************************************************
    !* EEQ charges q and electrostatic energy es from a direct solve (real32
    !* if single) of the linear system with one charge constraint per fragment.
    !* eeqtmp returns the packed gamma_ij and erf(gamma_ij*r_ij) for es_grad_mol.
    !***********************************************************************

    implicit none

    character(len=*),parameter :: source = 'gfnff_eg_goed'
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo

    logical,intent(in)  :: single

    integer,intent(in)  :: n

    integer,intent(in)  :: at(n)

    real(wp),intent(in)  :: sqrab(n*(n+1)/2)

    real(wp),intent(in)  :: r(n*(n+1)/2)

    real(wp),intent(in)  :: chrg

    real(wp),intent(in)  :: cn(n)

    real(wp),intent(out) :: q(n)

    real(wp),intent(out) :: es

    real(wp),intent(out) :: eeqtmp(2,n*(n+1)/2)

    type(TBorn),allocatable,intent(in) :: gbsa

    integer  :: m,i,j,k,ii,ij,io1,io2
    integer,allocatable :: ipiv(:)
    real(wp) :: gammij,tsqrt2pi,tmp
    real(wp),allocatable :: A(:,:),x(:)
    real(sp),allocatable :: A4(:,:),x4(:)
    parameter(tsqrt2pi=0.797884560802866_wp)
    logical :: exitRun

    m = n+topo%nfrag
    allocate (A(m,m),x(m))

    do i = 1,n
      x(i) = topo%chieeq(i)+param%cnf(at(i))*sqrt(cn(i))
    end do

    !>-- clear only the constraint border, the loop below writes the full
    !>   n x n block
    if (m > n) then
      A(:,n+1:m) = 0.0_wp
      A(n+1:m,:) = 0.0_wp
    end if

    !$omp parallel default(none) &
    !$omp shared(topo,n,sqrab,r,eeqtmp,A,at) &
    !$omp private(i,j,k,ij,gammij,tmp)
    !$omp do schedule(dynamic)
    do i = 1,n

      A(i,i) = tsqrt2pi/sqrt(topo%alpeeq(i))+topo%gameeq(i) ! J of i
      k = i*(i-1)/2

      do j = 1,i-1

        ij = k+j
        gammij = 1./sqrt(topo%alpeeq(i)+topo%alpeeq(j)) ! squared above
        tmp = erf(gammij*r(ij))
        eeqtmp(1,ij) = gammij
        eeqtmp(2,ij) = tmp
        A(j,i) = tmp/r(ij)
        A(i,j) = A(j,i)

      end do
    end do
    !$omp enddo
    !$omp end parallel

    do i = 1,topo%nfrag
      x(n+i) = topo%qfrag(i)
      do j = 1,n
        if (topo%fraglist(j) .eq. i) then
          A(n+i,j) = 1
          A(j,n+i) = 1
        end if
      end do
    end do

    if (allocated(gbsa)) then
      A(:n,:n) = A(:n,:n)+gbsa%bornMat(:,:)
    end if

    allocate (ipiv(m))

    if (single) then

      allocate (A4(m,m),x4(m))
      A4 = A
      x4 = x
      deallocate (A,x)
      call sytrf_wrap(a4,ipiv,io1)
      call sytrs_wrap(a4,x4,ipiv,io2)
      q(1:n) = x4(1:n)
      deallocate (A4,x4)

    else

      call sytrf_wrap(a,ipiv,io1)
      call sytrs_wrap(a,x,ipiv,io2)
      q(1:n) = x(1:n)
      deallocate (A,x)

    end if

    exitrun = (io1 /= 0).or.(io2 /= 0)
    if (exitRun) then
      write (stdout,'("**ERROR**",a,1x,a)') 'Solving linear equations failed',source
      return
    end if

    if (n .eq. 1) q(1) = chrg

    !>-- E_es = q^T*(0.5*A*q - X), summed over pairs

    es = 0.0_wp
    do i = 1,n

      ii = i*(i-1)/2

      do j = 1,i-1
        ij = ii+j
        tmp = eeqtmp(2,ij)
        es = es+q(i)*q(j)*tmp/r(ij)
      end do

      es = es-q(i)*(topo%chieeq(i)+param%cnf(at(i))*sqrt(cn(i))) &
      &        +q(i)*q(i)*0.5d0*(topo%gameeq(i)+tsqrt2pi/sqrt(topo%alpeeq(i)))

    end do

  end subroutine goed_gfnff

  subroutine goed_pbc_gfnff(single,n,at,xyz,r,cell,chrg,eeqtmp,cn,q,es,&
                        & gbsa,param,topo,gTrans,rTrans,x,cf)
    !***********************************************************************
    !* Periodic EEQ charges q and electrostatic energy es from the Ewald-summed
    !* Coulomb matrix. Also returns the reciprocal (without G=0) and real space
    !* translations gTrans/rTrans, two cells in each direction, the solution
    !* vector x including the constraint entries and the Ewald parameter cf.
    !***********************************************************************
    implicit none
    character(len=*),parameter :: source = 'gfnff_eg_goed'
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    logical,intent(in)  :: single     ! real*4 flag for solver
    integer,intent(in)  :: n          ! number of atoms
    integer,intent(in)     :: at(n)   ! ordinal numbers
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(in)    :: r(n,n)  ! dist
    type(TCell),intent(in) :: cell
    real(wp),intent(in)  :: chrg       ! total charge on system
    real(wp),intent(in)  :: cn(n)      ! CN
    real(wp),intent(out) :: q(n)       ! output charges
    real(wp),intent(out) :: es         ! ES energy
    real(wp),allocatable,intent(out) :: x(:)
    real(wp),intent(out) :: eeqtmp(2,n*(n+1)/2)    ! intermediates
    real(wp),intent(out) :: cf !convergence factor
    type(TBorn),allocatable,intent(in) :: gbsa

    integer  :: m,i,j,k,ij,io1,io2
    integer,allocatable :: ipiv(:)
    real(wp) :: gammij,tsqrt2pi,tmp
    real(wp),allocatable :: x_right(:)
    real(sp),allocatable :: A4(:,:),x4(:)

    real(wp),allocatable :: Amat(:,:),Amat_or(:,:)
    real(wp),allocatable,intent(out) :: rTrans(:,:)
    real(wp),allocatable,intent(out) :: gTrans(:,:)
    real(wp) :: vec(3)
    integer :: iRp,iG1,iG2,iG3,iT1,iT2,iT3
    integer,parameter :: ewaldCutD(3) = 2
    integer,parameter :: ewaldCutR(3) = 2
    real(wp) :: avgAlpeeq

    parameter(tsqrt2pi=0.797884560802866_wp)
    logical :: exitRun

    m = n+topo%nfrag ! # atoms + chrg constrain + frag constrain

    allocate (Amat(m,m),Amat_or(m,m),x(m),x_right(m),source=0.0_wp) ! matrix contains constrains -> linear equations

    iRp = 0
    allocate (gTrans(3,product(2*ewaldCutR+1)-1),source=0.0_wp)
    do iG1 = -ewaldCutR(1),ewaldCutR(1)
      do iG2 = -ewaldCutR(2),ewaldCutR(2)
        do iG3 = -ewaldCutR(3),ewaldCutR(3)
          if (iG1 == 0.and.iG2 == 0.and.iG3 == 0) cycle
          iRp = iRp+1
          vec(:) = [iG1,iG2,iG3]
          gTrans(:,iRp) = matmul(cell%rec_lat,vec)
        end do
      end do
    end do

    iRp = 0
    allocate (rTrans(3,product(2*ewaldCutD+1)))
    do iT1 = -ewaldCutD(1),ewaldCutD(1)
      do iT2 = -ewaldCutD(2),ewaldCutD(2)
        do iT3 = -ewaldCutD(3),ewaldCutD(3)
          iRp = iRp+1
          vec(:) = [iT1,iT2,iT3]
          rTrans(:,iRp) = matmul(cell%lattice,vec)
        end do
      end do
    end do
    !$omp parallel default(none) &
    !$omp shared(topo,n,r,eeqtmp) &
    !$omp private(i,j,k,ij,gammij,tmp)
    !$omp do schedule(dynamic)
    do i = 1,n
      k = i*(i-1)/2
      do j = 1,i-1
        ij = k+j
        gammij = 1./sqrt(topo%alpeeq(i)+topo%alpeeq(j)) ! squared above
        tmp = erf(gammij*r(i,j))
        eeqtmp(1,ij) = gammij
        eeqtmp(2,ij) = tmp
      end do
    end do
    !$omp enddo
    !$omp end parallel

    avgAlpeeq = sum(topo%alpeeq)/n
    cf = get_cf(rTrans,gTrans,cell%volume,avgAlpeeq)

    call get_amat_3d(n,at,xyz,cell,topo,cf,rTrans,gTrans,Amat)

    do i = 1,n
      x(i) = topo%chieeq(i)+param%cnf(at(i))*sqrt(cn(i))
    end do

    !>-- get_amat_3d fills row/column n+1 for a single total-charge constraint;
    !>   clear it so that every fragment constrains its own atoms only
    Amat(n+1:m,:) = 0.0_wp
    Amat(:,n+1:m) = 0.0_wp
    do i = 1,topo%nfrag
      x(n+i) = topo%qfrag(i)
      do j = 1,n
        if (topo%fraglist(j) .eq. i) then
          Amat(n+i,j) = 1
          Amat(j,n+i) = 1
        end if
      end do
    end do

    if (allocated(gbsa)) then
      Amat(:n,:n) = Amat(:n,:n)+gbsa%bornMat(:,:)
    end if

    allocate (ipiv(m))

    Amat_or = Amat
    if (single) then
      allocate (A4(m,m),x4(m))
      A4 = Amat
      x4 = x
      x_right(1:n) = x(1:n)
      deallocate (Amat)
      call sytrf_wrap(a4,ipiv,io1)
      call sytrs_wrap(a4,x4,ipiv,io2)
      q(1:n) = x4(1:n)
      x = x4
      deallocate (A4,x4)
    else
      x_right(1:n) = x(1:n)
      call sytrf_wrap(Amat,ipiv,io1)
      call sytrs_wrap(Amat,x,ipiv,io2)
      q(1:n) = x(1:n)
      deallocate (Amat)
    end if

    exitRUn = (io1 /= 0).or.(io2 /= 0)
    if (exitRun) then
      write (stdout,'("**ERROR**",a,1x,a)') 'Solving linear equations failed',source
      return
    end if

    if (n .eq. 1) q(1) = chrg

    !>-- E_es = q^T*(0.5*A*q - X); symv (y := alpha*A*x + beta*y) leaves the
    !>   bracket in x_right, x serves as q
    call symv(Amat_or,x,x_right,alpha=0.5_wp,beta=-1.0_wp)

    es = dot(x,x_right)
  end subroutine goed_pbc_gfnff

  function get_cf(rTrans,gTrans,vol,avgAlp) result(cf)
    !***********************************************************************
    !* Ewald splitting parameter: the cf in [1e-8,2] that minimises grFct,
    !* found by a section search. avgAlp is the average EEQ alpha.
    !***********************************************************************
    real(wp) :: cf
    integer,parameter :: ewaldCutD(3) = 2
    integer,parameter :: ewaldCutR(3) = 2
    real(wp),intent(in) :: rTrans(:,:)
    real(wp),intent(in) :: gTrans(:,:)
    real(wp),intent(in) :: vol
    real(wp),intent(in) :: avgAlp
    real(wp) :: minG,minR
    real(wp) :: gam
    real(wp) :: lenR
    real(wp),parameter :: tol = 1.0e-4_wp
    real(wp),parameter :: goldr = (1.0_wp+sqrt(2.0_wp))/2.0_wp
    real(wp) :: x1,x2,x3,x4
    integer :: i,iter

    cf = 0.14999_wp ! default
    gam = 1.0_wp/sqrt(2*avgAlp)

    minG = sqrt(minval(sum(gTrans(:,:)**2,dim=1)))
    minR = huge(1.0_wp)
    do i = 1,size(rTrans,dim=2)
      lenR = sqrt(sum(rTrans(:,i)**2))
      if (lenR .ne. 0.0_wp) then
        minR = min(minR,lenR)
      end if
    end do

    !>-- golden-section type search; goldr is (1+sqrt(2))/2, not the golden ratio
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
      x2 = x4-(x4-x1)/goldr
      x3 = x1+(x4-x1)/goldr
    end do

    cf = (x1+x4)/2.0_wp
  end function get_cf

  function grFct(cfCurr,minG,minR,vol,gam) result(gr_diff)
    !***********************************************************************
    !* Mismatch between the leading reciprocal and real space Ewald terms at
    !* splitting parameter cfCurr, estimated from the shortest translations.
    !***********************************************************************
    real(wp) :: gr_diff
    real(wp),intent(in) :: cfCurr,minG,minR,vol,gam
    real(wp) :: gPart,rPart

    gPart = 4.0_wp*pi*exp(-minG**2/(4.0_wp*cfCurr**2))/(vol*minG**2)
    rPart = -erf(cfCurr*minR)/minR+erf(gam*minR)/minR
    gr_diff = abs(gPart-rPart)
  end function grFct

  subroutine es_grad_sigma(nat,at,xyz,cell,topo,nlist,rTrans,gTrans,xtmp,cf, &
           & sigma,gradient,mcf_ees)
    !***********************************************************************
    !* Periodic ES contribution to gradient and sigma from the contracted
    !* Ewald A-matrix derivative, scaled by mcf_ees. xtmp is the vector A is
    !* differentiated against; the charges nlist%q are contracted with it.
    !***********************************************************************
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TCell),intent(in) :: cell
    type(TGFFTopology),intent(in) :: topo
    type(TGFFNeighbourList),intent(in) :: nlist
    real(wp),intent(in) :: cf
    real(wp),intent(in) :: gTrans(:,:)
    real(wp),intent(in) :: rTrans(:,:)
    real(wp),intent(in) :: xtmp(nat+topo%nfrag)
    real(wp),intent(inout) :: sigma(3,3)
    real(wp),intent(inout) :: gradient(3,nat)
    real(wp),intent(in) :: mcf_ees
    real(wp) :: dEdL(3,3)
    real(wp),allocatable :: dEdr(:,:)
    allocate (dEdr(3,nat),source=0.0_wp)
    dEdL = 0.0_wp

    call get_damat_3d(nat,at,xyz,cell,topo,cf,xtmp,nlist%q,rTrans,gTrans,dEdr,dEdL)

    gradient(:,:) = gradient(:,:)+mcf_ees*dEdr(:,:)

    !>-- only the 0.5*q*A'*q part of Ees = q^T*(0.5*A*q - X); the caller adds -q*X'
    sigma(:,:) = sigma(:,:)+0.5_wp*mcf_ees*dEdL(:,:)

  end subroutine es_grad_sigma

! code from here to the end of the module taken from dftd4
  subroutine get_wsc_reference(cell,nat,xyz,wref)
    !***********************************************************************
    !* Geometry the Wigner-Seitz images are searched against: the one the
    !* cell was set up with (cell%wsc_xyz), else xyz. The image assignment is
    !* a discrete choice the gradient does not differentiate, so letting it
    !* move with the atoms would put small steps into the energy surface.
    !***********************************************************************
    implicit none
    type(TCell),intent(in) :: cell
    integer,intent(in) :: nat
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),allocatable,intent(out) :: wref(:,:)

    if (allocated(cell%wsc_xyz)) then
      if (size(cell%wsc_xyz,dim=2) == nat) then
        wref = cell%wsc_xyz
        return
      end if
    end if
    wref = xyz

  end subroutine get_wsc_reference

  subroutine get_gfactors(gTrans,vol,alp,gfac)
    !***********************************************************************
    !* Ewald prefactor 4*pi/V * exp(-G^2/4a^2) / G^2 for every reciprocal
    !* lattice vector. It is pair independent, hence built once per matrix.
    !* Vanishing G gets a zero factor, which the consumers skip (G=0 term).
    !***********************************************************************
    implicit none
    real(wp),intent(in) :: gTrans(:,:)
    real(wp),intent(in) :: vol,alp
    real(wp),allocatable,intent(out) :: gfac(:)

    integer :: itr
    real(wp) :: g2,fac
    real(wp),parameter :: eps = 1.0e-9_wp

    allocate (gfac(size(gTrans,2)))
    fac = 4*pi/vol
    do itr = 1,size(gTrans,2)
      g2 = gTrans(1,itr)**2+gTrans(2,itr)**2+gTrans(3,itr)**2
      if (g2 < eps) then
        gfac(itr) = 0.0_wp
      else
        gfac(itr) = fac*exp(-0.25_wp*g2/(alp*alp))/g2
      end if
    end do

  end subroutine get_gfactors

  subroutine get_amat_3d(nat,at,xyz,cell,topo,alpha,rTrans,gTrans,amat)
    !***********************************************************************
    !* Ewald-summed EEQ Coulomb matrix with self terms and the total charge
    !* constraint in row/column nat+1. alpha is the Ewald splitting parameter.
    !***********************************************************************
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TCell),intent(in) :: cell
    type(TGFFTopology),intent(in) :: topo
    real(wp),intent(in) :: alpha
    real(wp),intent(in) :: rTrans(:,:)
    real(wp),intent(in) :: gTrans(:,:)
    real(wp),intent(out) :: amat(:,:)

    integer :: iat,jat,img,wc,lattr(3,wsc_maxcells)
    real(wp) :: vec(3),gam,wsw,dtmp,rtmp,vol
    real(wp),allocatable :: gfac(:),wref(:,:)
    real(wp),parameter :: zero(3) = 0.0_wp
    real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
    real(wp),parameter :: sqrtpi = 1.772453850905516_wp

    amat(:,:) = 0.0_wp

    vol = cell%volume ! abs(matdet_3x3(cell%lattice))

    call get_gfactors(gTrans,vol,alpha,gfac)
    call get_wsc_reference(cell,nat,xyz,wref)

    !>-- no reduction on amat: each element is written by exactly one iteration
    !$omp parallel do default(none) schedule(runtime) &
    !$omp shared(amat) shared(nat,at,xyz,cell, topo, rTrans, gTrans, alpha, vol, gfac, wref) &
    !$omp private(iat, jat, img, gam, wsw, vec, dtmp, rtmp, wc, lattr)
    do iat = 1,nat
      do jat = 1,iat-1
        gam = 1.0_wp/sqrt(topo%alpeeq(iat)+topo%alpeeq(jat))
        !>-- the reciprocal sum is lattice periodic and the image weights add
        !>   to one, so evaluate it once for the plain pair vector
        vec = xyz(:,iat)-xyz(:,jat)
        call get_amat_rec_3d(vec,gfac,gTrans,rtmp)
        amat(jat,iat) = amat(jat,iat)+rtmp
        amat(iat,jat) = amat(iat,jat)+rtmp
        !>-- the direct sum reaches only two cells and must be centred on the
        !>   nearest image; a 27 cell scan replaces a table for every pair
        call wsc_images(nat,wref,iat,jat,cell%lattice,cell%pbc,lattr,wc)
        wsw = 1.0_wp/real(wc,wp)
        do img = 1,wc
          vec = xyz(:,iat)-xyz(:,jat) &
             & -(cell%lattice(:,1)*lattr(1,img) &
             &  +cell%lattice(:,2)*lattr(2,img) &
             &  +cell%lattice(:,3)*lattr(3,img))
          call get_amat_dir_3d(vec,gam,alpha,rTrans,dtmp)
          amat(jat,iat) = amat(jat,iat)+dtmp*wsw
          amat(iat,jat) = amat(iat,jat)+dtmp*wsw
        end do
      end do

      !>-- self term, image independent since vec is zero
      gam = 1.0_wp/sqrt(2.0_wp*topo%alpeeq(iat))
      vec = zero
      call get_amat_dir_3d(vec,gam,alpha,rTrans,dtmp)
      call get_amat_rec_3d(vec,gfac,gTrans,rtmp)
      amat(iat,iat) = amat(iat,iat)+dtmp+rtmp

      dtmp = topo%gameeq(iat)+sqrt2pi/sqrt(topo%alpeeq(iat))-2*alpha/sqrtpi
      amat(iat,iat) = amat(iat,iat)+dtmp
    end do
    !$omp end parallel do

    amat(nat+1,1:nat+1) = 1.0_wp
    amat(1:nat+1,nat+1) = 1.0_wp
    amat(nat+1,nat+1) = 0.0_wp

  end subroutine get_amat_3d

  subroutine get_amat_dir_3d(rij,gam,alp,trans,amat)
    !***********************************************************************
    !* Real space Ewald part of one A-matrix element: sum over translations
    !* of (erf(gam*r)-erf(alp*r))/r with r = |rij+trans|.
    !***********************************************************************
    real(wp),intent(in) :: rij(3)
    real(wp),intent(in) :: gam
    real(wp),intent(in) :: alp
    real(wp),intent(in) :: trans(:,:)
    real(wp),intent(out) :: amat

    integer :: itr
    real(wp) :: vec(3),r1,tmp
    real(wp),parameter :: eps = 1.0e-9_wp

    amat = 0.0_wp

    do itr = 1,size(trans,2)
      vec(:) = rij+trans(:,itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      tmp = erf(gam*r1)/r1-erf(alp*r1)/r1
      amat = amat+tmp
    end do

  end subroutine get_amat_dir_3d

  subroutine get_amat_rec_3d(rij,gfac,trans,amat)
    !***********************************************************************
    !* Reciprocal space Ewald part of one A-matrix element: sum over G of
    !* gfac(G)*cos(G*rij), gfac from get_gfactors.
    !***********************************************************************
    real(wp),intent(in) :: rij(3)
    real(wp),intent(in) :: gfac(:)
    real(wp),intent(in) :: trans(:,:)
    real(wp),intent(out) :: amat

    integer :: itr

    amat = 0.0_wp

    do itr = 1,size(trans,2)
      if (gfac(itr) == 0.0_wp) cycle
      amat = amat+cos(rij(1)*trans(1,itr)+rij(2)*trans(2,itr) &
         &            +rij(3)*trans(3,itr))*gfac(itr)
    end do

  end subroutine get_amat_rec_3d

  subroutine get_damat_3d(nat,at,xyz,cell,topo,alpha,qvec,qcon,rTrans,gTrans,dEdr,dEdL)
    !***********************************************************************
    !* Derivative of the Ewald A-matrix, contracted on the fly so that no
    !* (3,nat,nat) array has to be stored or reduced over.
    !* Input:
    !*   qvec  - vector A is differentiated against (the EEQ right-hand side)
    !*   qcon  - vector the derivative is contracted with (the charges)
    !* Output, both accumulated over pairs, unscaled:
    !*   dEdr  - qvec(i)*sum_j dA(:,i,j)*qcon(j), the Cartesian part
    !*   dEdL  - sum_ij qvec(i)*dA(:,:,i,j)*qcon(j), the strain part
    !***********************************************************************
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TCell),intent(in) :: cell
    type(TGFFTopology),intent(in) :: topo
    real(wp),intent(in) :: alpha
    real(wp),intent(in) :: qvec(:)
    real(wp),intent(in) :: qcon(:)
    real(wp),intent(in) :: rTrans(:,:)
    real(wp),intent(in) :: gTrans(:,:)
    real(wp),intent(out) :: dEdr(:,:)
    real(wp),intent(out) :: dEdL(:,:)

    integer :: iat,jat,img,wc,lattr(3,wsc_maxcells)
    real(wp) :: vol,gam,wsw,vec(3),dG(3),dS(3,3),wij
    real(wp) :: dGd(3),dSd(3,3),dGr(3),dSr(3,3)
    real(wp),allocatable :: gfac(:),wref(:,:)
    real(wp),parameter :: zero(3) = 0.0_wp

    dEdr(:,:) = 0.0_wp
    dEdL(:,:) = 0.0_wp

    vol = cell%volume ! abs(matdet_3x3(cell%lattice))

    call get_gfactors(gTrans,vol,alpha,gfac)
    call get_wsc_reference(cell,nat,xyz,wref)

    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:dEdr, dEdL) &
    !$omp shared(nat,at,xyz,cell, topo, alpha, vol, rTrans, gTrans, qvec, qcon, gfac, wref) &
    !$omp private(iat, jat, img, gam, wsw, wij, vec, dG, dS, &
    !$omp& dGr, dSr, dGd, dSd, wc, lattr)
    do iat = 1,nat
      do jat = 1,iat-1
        gam = 1.0_wp/sqrt(topo%alpeeq(iat)+topo%alpeeq(jat))
        !>-- reciprocal part once per pair, direct part per image, see get_amat_3d
        vec = xyz(:,iat)-xyz(:,jat)
        call get_damat_rec_3d(vec,gfac,alpha,gTrans,dGr,dSr)
        dG(:) = dGr
        dS(:,:) = dSr
        call wsc_images(nat,wref,iat,jat,cell%lattice,cell%pbc,lattr,wc)
        wsw = 1.0_wp/real(wc,wp)
        do img = 1,wc
          vec = xyz(:,iat)-xyz(:,jat) &
             & -(cell%lattice(:,1)*lattr(1,img) &
             &  +cell%lattice(:,2)*lattr(2,img) &
             &  +cell%lattice(:,3)*lattr(3,img))
          call get_damat_dir_3d(vec,gam,alpha,rTrans,dGd,dSd)
          dG = dG+dGd*wsw
          dS = dS+dSd*wsw
        end do
        dEdr(:,iat) = dEdr(:,iat)+dG*(qvec(iat)*qcon(jat))
        dEdr(:,jat) = dEdr(:,jat)-dG*(qvec(jat)*qcon(iat))
        wij = qvec(iat)*qcon(jat)+qvec(jat)*qcon(iat)
        dEdL(:,:) = dEdL(:,:)+dS*wij
      end do

      !>-- self term, image independent since vec is zero
      gam = 1.0_wp/sqrt(2.0_wp*topo%alpeeq(iat))
      vec = zero
      call get_damat_dir_3d(vec,gam,alpha,rTrans,dGd,dSd)
      call get_damat_rec_3d(vec,gfac,alpha,gTrans,dGr,dSr)
      dS = dSd+dSr
      dEdL(:,:) = dEdL(:,:)+dS*(qvec(iat)*qcon(iat))
    end do
    !$omp end parallel do

  end subroutine get_damat_3d

  subroutine get_damat_dir_3d(rij,gam,alp,trans,dg,ds)
    !***********************************************************************
    !* Real space part of the A-matrix derivative for one pair: dg w.r.t.
    !* rij, ds the strain part (outer product with the translated vector).
    !***********************************************************************
    real(wp),intent(in) :: rij(3)
    real(wp),intent(in) :: gam
    real(wp),intent(in) :: alp
    real(wp),intent(in) :: trans(:,:)
    real(wp),intent(out) :: dg(3)
    real(wp),intent(out) :: ds(3,3)

    integer :: itr
    real(wp) :: vec(3),r1,r2,gtmp,atmp,gam2,alp2
    real(wp),parameter :: eps = 1.0e-9_wp
    real(wp),parameter :: sqrtpi = 1.772453850905516_wp

    dg(:) = 0.0_wp
    ds(:,:) = 0.0_wp

    gam2 = gam*gam
    alp2 = alp*alp

    do itr = 1,size(trans,2)
      vec(:) = rij+trans(:,itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1*r1
      gtmp = +2*gam*exp(-r2*gam2)/(sqrtpi*r2)-erf(r1*gam)/(r2*r1)
      atmp = -2*alp*exp(-r2*alp2)/(sqrtpi*r2)+erf(r1*alp)/(r2*r1)
      dg(:) = dg+(gtmp+atmp)*vec
      ds(:,1) = ds(:,1)+(gtmp+atmp)*vec(1)*vec
      ds(:,2) = ds(:,2)+(gtmp+atmp)*vec(2)*vec
      ds(:,3) = ds(:,3)+(gtmp+atmp)*vec(3)*vec
    end do

  end subroutine get_damat_dir_3d

  subroutine get_damat_rec_3d(rij,gfac,alp,trans,dg,ds)
    !***********************************************************************
    !* Reciprocal space counterpart of get_damat_dir_3d.
    !***********************************************************************
    real(wp),intent(in) :: rij(3)
    real(wp),intent(in) :: gfac(:)
    real(wp),intent(in) :: alp
    real(wp),intent(in) :: trans(:,:)
    real(wp),intent(out) :: dg(3)
    real(wp),intent(out) :: ds(3,3)

    integer :: itr
    real(wp) :: vec(3),g2,gv,etmp,dtmp,ctmp,pref,alp2
    real(wp),parameter :: unity(3,3) = reshape(&
       & [1,0,0,0,1,0,0,0,1],shape(unity))

    dg(:) = 0.0_wp
    ds(:,:) = 0.0_wp
    alp2 = alp*alp

    do itr = 1,size(trans,2)
      etmp = gfac(itr)
      if (etmp == 0.0_wp) cycle
      vec(:) = trans(:,itr)
      g2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
      gv = rij(1)*vec(1)+rij(2)*vec(2)+rij(3)*vec(3)
      dtmp = -sin(gv)*etmp
      ctmp = etmp*cos(gv)
      pref = 2.0_wp/g2+0.5_wp/alp2
      dg(:) = dg+dtmp*vec
      ds(:,1) = ds(:,1)+ctmp*(pref*vec(1)*vec-unity(:,1))
      ds(:,2) = ds(:,2)+ctmp*(pref*vec(2)*vec-unity(:,2))
      ds(:,3) = ds(:,3)+ctmp*(pref*vec(3)*vec-unity(:,3))
    end do

  end subroutine get_damat_rec_3d

end module gfnff_eg_es
