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
module gfnff_qm_setup
  use iso_fortran_env,only:wp => real64,sp => real64, stdout=>output_unit
  use gfnff_data_types,only:TGFFData,init,TGFFGenerator,TGFFTopology
  implicit none
  private
  public :: gfnffqmsolve

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

! solve QM Hamiltonian A in overlap basis S (if ovlp=.true., ZDO otherwise)
! and return density matrix in A (and energy weighted density in S if ovlp=.true.)
! ndim is the dimension of the problem for nel electrons with nopen more alpha than beta
  subroutine gfnffqmsolve(pr,A,S,ovlp,et,ndim,nopen,nel,eel,focc,e,io)
    implicit none
    character(len=*),parameter :: source = 'gfnffqmsolve()'
    integer :: ndim      ! # basis
    integer :: nopen     ! # of open shells
    integer :: nel       ! # of electrons
    logical :: ovlp      ! in overlap basis?
    logical :: pr        !
    real(wp) :: et       ! electronic temp Fermi smear
    real(wp) :: eel      ! electronic energy = sum_occ nocc*eps
    real(wp) :: focc(ndim)! occupations
    real(wp) :: e(ndim)  ! eigenvalues
    real(wp) :: A(ndim,ndim)
    real(wp) :: S(ndim,ndim)
    integer,intent(out) :: io ! output status
    integer  :: ihomoa,ihomob,i,liwork,info,lwork
    real(wp) :: ga,gb,efa,efb,nfoda,nfodb
    real(wp),allocatable :: X(:,:)
    real(wp),allocatable :: focca(:),foccb(:),aux(:)
    integer,allocatable  :: iwork(:),ifail(:)
    !> LAPACK
    external :: dsyev,dsygvd

    io = 0

    allocate (focca(ndim),foccb(ndim))

! ZDO case
    if (.not.ovlp) then

      lwork = 1+6*ndim+2*ndim**2
!     allocate(aux4(lwork),e4(ndim),A4(ndim,ndim))
!     A4 = A
!     call ssyev ('V','U',ndim,A4,ndim,e4,aux4,lwork,info)
!     A = A4
!     e = e4
!     deallocate(aux4,A4,e4)
      allocate (aux(lwork))
      call dsyev('V','U',ndim,A,ndim,e,aux,lwork,info)

    else
! with overlap case

      allocate (aux(1),iwork(1),ifail(ndim))
      call DSYGVD(1,'V','U',ndim,A,ndim,S,ndim,e,aux, &  !workspace query
     &               -1,IWORK,LIWORK,INFO)
      lwork = int(aux(1))
      liwork = iwork(1)
      deallocate (aux,iwork)

      allocate (aux(lwork),iwork(liwork))              !do it
      call DSYGVD(1,'V','U',ndim,A,ndim,S,ndim,e,aux, &
     &            LWORK,IWORK,LIWORK,INFO)

    end if
! DONE

! scale energy  so that gaps a roughly ok (in eV)

    e = e*0.1_wp*27.2113957_wp

    if (info .ne. 0) then
      write (stdout,*) 'diag error in ',source
      io = info 
    end if

    if (et .gt. 1.d-3) then
! Fermi smearing, convert restricted occ first to alpha/beta
      call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      if (ihomoa .le. ndim.and.ihomoa .gt. 0) then
        call fermismear(.false.,ndim,ihomoa,et,e,focca,nfoda,efa,ga)
      else
        focca = 0
      end if
      if (ihomob .le. ndim.and.ihomob .gt. 0) then
        call fermismear(.false.,ndim,ihomob,et,e,foccb,nfodb,efb,gb)
      else
        foccb = 0
      end if
      focc = focca+foccb
      if (ihomoa+1 .le. ndim) then
        if (abs(focc(ihomoa)-focc(ihomoa+1)) .lt. 1.d-4) then ! a perfect birad is anit-aromatic
          focc = 0                                        ! and hence we break the sym
          do i = 1,nel/2
            focc(i) = 2.0d0
          end do
          if(pr)then
          write (stdout,*) 'perfect biradical detected at FT-HMO level. Breaking the symmetry'
          write (stdout,*) 'because its assumed to be an anit-aromatic system like COT or CB.'
          endif
        end if
      end if
    else
      focc = 0
      do i = 1,nel/2
        focc(i) = 2.0d0
      end do
      if (2*(nel/2) .ne. nel) focc(nel/2+1) = 1.0
    end if

    focca = focc*e
    eel = sum(focca)
! density matrix
    S = A
    call dmat(ndim,focc,S,A)

    if (ovlp) then
      allocate (X(ndim,ndim))
! energy weighted density on S
      call dmat(ndim,focca,S,X)
      S = X
    end if

    deallocate (focca,foccb)

  end subroutine gfnffqmsolve
!========================================================================================!

  subroutine fermismear(prt,norbs,nel,t,eig,occ,fod,e_fermi,s)
    integer,intent(in)  :: norbs
    integer,intent(in)  :: nel
    real(wp),intent(in)  :: eig(norbs)
    real(wp),intent(out) :: occ(norbs)
    real(wp),intent(in)  :: t
    real(wp),intent(out) :: fod
    real(wp),intent(out) :: e_fermi
    logical,intent(in)  :: prt

    real(wp) :: bkt,occt,total_number
    real(wp) :: total_dfermi,dfermifunct,fermifunct,s,change_fermi

    real(wp),parameter :: autoev = 27.21138505_wp
    real(wp),parameter :: kB = 3.166808578545117e-06_wp
    real(wp),parameter :: boltz = kB*autoev
    real(wp),parameter :: thr = 1.d-9
    integer :: ncycle,i,j,m,k,i1,i2

    bkt = boltz*t

    e_fermi = 0.5*(eig(nel)+eig(nel+1))
    occt = nel

    do ncycle = 1,200
      total_number = 0.0
      total_dfermi = 0.0
      do i = 1,norbs
        fermifunct = 0.0
        if ((eig(i)-e_fermi)/bkt .lt. 50) then
          fermifunct = 1.0/(exp((eig(i)-e_fermi)/bkt)+1.0)
          dfermifunct = exp((eig(i)-e_fermi)/bkt)/ &
          &       (bkt*(exp((eig(i)-e_fermi)/bkt)+1.0)**2)
        else
          dfermifunct = 0.0
        end if
        occ(i) = fermifunct
        total_number = total_number+fermifunct
        total_dfermi = total_dfermi+dfermifunct
      end do
      change_fermi = (occt-total_number)/total_dfermi
      e_fermi = e_fermi+change_fermi
      if (abs(occt-total_number) .le. thr) exit
    end do

    fod = 0
    s = 0
    do i = 1,norbs
      if (occ(i) .gt. thr.and.1.0d00-occ(i) .gt. thr) &
      &   s = s+occ(i)*log(occ(i))+(1.0d0-occ(i))*log(1.0d00-occ(i))
      if (eig(i) .lt. e_fermi) then
        fod = fod+1.0d0-occ(i)
      else
        fod = fod+occ(i)
      end if
    end do
    s = s*kB*t

    if (prt) then
      write (*,'('' t,e(fermi),nfod : '',2f10.3,f10.6)') t,e_fermi,fod
    end if

  end subroutine fermismear
!========================================================================================!
  subroutine occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
    integer  :: nel
    integer  :: nopen
    integer  :: ndim
    integer  :: ihomoa
    integer  :: ihomob
    real(wp) :: focca(ndim)
    real(wp) :: foccb(ndim)
    integer  :: focc(ndim)
    integer  :: i,na,nb,ihomo

    focc = 0
    focca = 0
    foccb = 0
!>--- even nel
    if (mod(nel,2) .eq. 0) then
      ihomo = nel/2
      do i = 1,ihomo
        focc(i) = 2
      end do
      if (2*ihomo .ne. nel) then
        ihomo = ihomo+1
        focc(ihomo) = 1
        if (nopen .eq. 0) nopen = 1
      end if
      if (nopen .gt. 1) then
        do i = 1,nopen/2
          focc(ihomo-i+1) = focc(ihomo-i+1)-1
          focc(ihomo+i) = focc(ihomo+i)+1
        end do
      end if
!>--- odd nel
    else
      na = nel/2+(nopen-1)/2+1
      nb = nel/2-(nopen-1)/2
      do i = 1,na
        focc(i) = focc(i)+1
      end do
      do i = 1,nb
        focc(i) = focc(i)+1
      end do
    end if

    do i = 1,ndim
      if (focc(i) .eq. 2) then
        focca(i) = 1.0d0
        foccb(i) = 1.0d0
      end if
      if (focc(i) .eq. 1) focca(i) = 1.0d0
    end do

    ihomoa = 0
    ihomob = 0
    do i = 1,ndim
      if (focca(i) .gt. 0.99) ihomoa = i
      if (foccb(i) .gt. 0.99) ihomob = i
    end do

  end subroutine occu
!========================================================================================!
  subroutine dmat(ndim,focc,C,P)
    !> density matrix
    !> C: MO coefficient
    !> P  dmat
    use gfnff_math_wrapper,only:gemm
    integer,intent(in)  :: ndim
    real(wp),intent(in)  :: focc(:)
    real(wp),intent(in)  :: C(:,:)
    real(wp),intent(out) :: P(:,:)
    integer :: i,m
    real(wp),allocatable :: Ptmp(:,:)

    allocate (Ptmp(ndim,ndim))
    Ptmp = 0.0_wp

    do m = 1,ndim
      do i = 1,ndim
        Ptmp(i,m) = C(i,m)*focc(m)
      end do
    end do
    call gemm(C,Ptmp,P,transb='t')

    deallocate (Ptmp)
  end subroutine dmat

!========================================================================================!
end module gfnff_qm_setup
