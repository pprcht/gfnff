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
! along with gfnff. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert, Sebastian Spicher, Stefan Grimme
!> at https://github.com/grimme-lab/xtb
! ------------------------------------------------------------------------------

!> Bond-path distance matrix, packed lower triangle: pair(i,j) is 1, 2 or 3
!> bonds apart (5 if farther, 0 only on the diagonal). Drives non-bonded
!> screening; one byte per pair keeps the matrix small at large atom counts.
module gfnff_topo_bondmat
  use iso_fortran_env,only:int8
  use gfnff_neighbor,only:TNeigh
  use gfnff_geometry,only:lin
  implicit none
  private

  public :: nbondmat,nbondmat_pbc
  public :: pairsbond,pairsbond_pbc
  public :: countf

contains  !> MODULE PROCEDURES START HERE

  subroutine nbondmat(n,numnb,numctr,nb,pair)
    !***********************************
    !* Builds the packed bond-path distance matrix (see module header).
    !***********************************
    implicit none
    integer,intent(in)  :: n,numnb,numctr
    integer,intent(in)  :: nb(numnb,n,numctr)
    integer,intent(out) ::  pair(n*(n+1)/2)
    integer :: i,ni,newi,j,newatom,tag,d,i1,ni1,iii,ii,k
    integer,allocatable :: list(:,:),nlist(:,:),nnn(:),nn(:)
    logical :: da

    allocate (nnn(n),nn(n),list(5*n,n),nlist(5*n,n))

    nn(1:n) = nb(numnb,1:n,1)

    pair = 0
    list = 0
    do i = 1,n
      ni = nn(i)
      list(1:ni,i) = nb(1:ni,i,1)
    end do

    nlist = list

    pair = 0
    do i = 1,n
      do j = 1,nb(numnb,i,1)
        k = nb(j,i,1)
        pair(lin(k,i)) = 1
      end do
    end do

    tag = 1

    do d = 1,2 ! expand the neighbor shell twice, tagging 2- and 3-bond paths

      do i = 1,n
        ni = nn(i)
        newi = ni
        do ii = 1,ni
          i1 = list(ii,i)
          ni1 = nb(numnb,i1,1)
          do iii = 1,ni1
            newatom = nb(iii,i1,1)
            da = .false.
            do j = 1,newi
              if (newatom .eq. list(j,i)) da = .true.
            end do
            if (.not.da) then
              newi = newi+1
              nlist(newi,i) = newatom
            end if
          end do
        end do
        nnn(i) = newi
      end do

      list = nlist
      nn = nnn

      tag = tag+1
      call pairsbond(n,nn,list,pair,tag)

    end do
    do i = 1,n
      do j = 1,i
        if (i .ne. j.and.pair(lin(j,i)) .eq. 0) pair(lin(j,i)) = 5
      end do
    end do

  end subroutine nbondmat

  subroutine pairsbond(n,nn,list,pair,tag)
    !***********************************
    !* Sets unset pair(i,j) = tag where i,j are mutual list(1:nn,:) members.
    !***********************************
    implicit none
    integer n,nn(n),list(5*n,n),tag
    integer i,j,k,ni,nj,ii,jj,ij
    integer pair(n*(n+1)/2)
    logical dai,daj

    do i = 1,n
      ni = nn(i)
      ij = i*(i-1)/2
      do j = 1,i-1
        k = ij+j
        nj = nn(j)
        dai = .false.
        daj = .false.
        do ii = 1,ni
          if (list(ii,i) .eq. j) daj = .true.
        end do
        do jj = 1,nj
          if (list(jj,j) .eq. i) dai = .true.
        end do
        if (dai.and.daj.and.pair(k) .eq. 0) then
          pair(k) = tag
        end if
      end do
    end do

  end subroutine pairsbond

  subroutine nbondmat_pbc(n,numnb,numctr,nb,neigh,pair)
    !***********************************
    !* Periodic counterpart of nbondmat: builds pair(j,i,iTr) using the
    !* same convention; detects only paths within two neighboring cells.
    !***********************************
    implicit none
    type(TNeigh),intent(in) :: neigh ! for locating neighbor
    integer,intent(in)  :: n,numnb,numctr
    integer,intent(in)  :: nb(numnb,n,numctr)
    integer(int8),allocatable,intent(out) ::  pair(:,:,:)
    integer :: nnbi,nbi(2,numnb),cval
    integer :: i,inew,j,inb,ixnb,iTr,iTrnew,sumiTr,k,iTr2,l
    integer  :: nbr(numnb,n,numctr) ! reduced neighbor list no unpaired bonds
    logical :: hasnb ! true if atom i and k have a paired bond (both have each other as nb)
    integer :: tmpp(3,10*n),nt
    !>-- tmpp mirrors nbondmat's setup: for mindless03, bonds in bpair
    !>   between A and B need not sum to bpair(A,B). Only paired bonds are
    !>   kept here (metals can have one-sided bonds the partner does not
    !>   list back); tmpp remembers them for restoration below.
    tmpp = 0 ! tmpp_usage
    nt = 0     ! tmpp_usage
    nbr = nb
    do i = 1,n
      do iTr = 1,neigh%numctr
        do j = 1,nb(neigh%numnb,i,iTr)
          k = nb(j,i,iTr)
          hasnb = .false.
          do iTr2 = 1,neigh%numctr
            do l = 1,nb(neigh%numnb,k,iTr2)
              if (nb(l,k,iTr2) .eq. i) then
                hasnb = .true.
              end if
            end do
          end do
          if (.not.hasnb) then
            do l = j,neigh%numnb-2 ! fails if an atom has numnb-1 neighbors !
              nbr(l,i,iTr) = nb(l+1,i,iTr)
            end do
            nbr(numnb,i,iTr) = nbr(numnb,i,iTr)-1
            nt = nt+1       ! tmpp_usage
            tmpp(1,nt) = k    ! tmpp_usage
            tmpp(2,nt) = i    ! tmpp_usage
            tmpp(3,nt) = iTr  ! tmpp_usage
          end if
        end do
      end do
    end do

    allocate (pair(n,n,numctr),source=0_int8)
    do i = 1,n
      do iTr = 1,numctr
        do inb = 1,nbr(numnb,i,iTr)
          j = nbr(inb,i,iTr)
          pair(j,i,iTr) = 1_int8
        end do
      end do
      do cval = 1,2 ! 2nd/3rd neighbors; paths beyond two cells go undetected
        call countf(n,numctr,numnb,pair,i,cval,nnbi,nbi)
        do ixnb = 1,nnbi
          inew = nbi(1,ixnb)
          iTrnew = nbi(2,ixnb)
          if (iTrnew .eq. 1) then
            do iTr = 1,numctr
              do inb = 1,nbr(numnb,inew,iTr)
                j = nbr(inb,inew,iTr)
                if (pair(j,i,iTr) .ne. 0) cycle ! take shortest path only
                if (j .eq. i.AND.iTr .eq. 1) cycle ! dont count to-self-bonds
                pair(j,i,iTr) = int(cval+1,int8)
              end do
            end do
          else
            !>-- same-cell neighbor indices match across cells: use nb(numnb,inew,iTr=1) first
            do inb = 1,nbr(numnb,inew,1)
              j = nbr(inb,inew,1)
              if (pair(j,i,iTrnew) .ne. 0) cycle ! take shortest path only
              pair(j,i,iTrnew) = int(cval+1,int8) ! to-self-bonds not possible since iTr.ne.1
            end do
            do iTr = 2,numctr
              do inb = 1,nbr(numnb,inew,iTr)
                j = nbr(inb,inew,iTr)
                sumiTr = neigh%fTrSum(iTr,iTrnew) ! combined translation, for pair's index
                if (sumiTr .eq. -1.or.(sumiTr .eq. 1.and.j .eq. i)) cycle ! -1: outside 27 cells; 1: to-self
                if (sumiTr .gt. 27) cycle
                if (pair(j,i,sumiTr) .ne. 0) cycle ! take shortest path only

                pair(j,i,sumiTr) = int(cval+1,int8)
              end do
            end do
          end if
        end do ! ixnb
      end do ! cval
    end do ! i
    do i = 1,n
      do j = 1,n
        do iTr = 1,numctr
          if (pair(j,i,iTr) .eq. 0.and.j .ne. i) pair(j,i,iTr) = 5_int8
        end do
      end do
    end do

    !>-- tmpp_usage: restore the deleted one-sided bonds as bonds again
    do l = 1,10*n
      if (tmpp(1,l) .eq. 0) exit
      k = tmpp(1,l)
      i = tmpp(2,l)
      iTr = tmpp(3,l)
      pair(k,i,iTr) = 1_int8
      pair(i,k,neigh%iTrNeg(iTr)) = 1_int8
    end do

  end subroutine nbondmat_pbc

  subroutine pairsbond_pbc(n,nn,list,pair,tag)
    !***********************************
    !* Periodic-path counterpart of pairsbond (see there).
    !***********************************
    implicit none
    integer n,nn(n),list(5*n,n),tag
    integer i,j,k,ni,nj,ii,jj,ij
    integer pair(n*(n+1)/2)
    logical dai,daj

    do i = 1,n
      ni = nn(i)
      ij = i*(i-1)/2
      do j = 1,i-1
        k = ij+j
        nj = nn(j)
        dai = .false.
        daj = .false.
        do ii = 1,ni
          if (list(ii,i) .eq. j) daj = .true.
        end do
        do jj = 1,nj
          if (list(jj,j) .eq. i) dai = .true.
        end do
        if (dai.and.daj.and.pair(k) .eq. 0) then
          pair(k) = tag
        end if
      end do
    end do

  end subroutine pairsbond_pbc

  subroutine countf(n,numctr,numnb,pair,i,cval,nnbi,nbi)
    !***********************************
    !* Collects atom i's neighbors exactly cval bonds away (nbi holds index/cell).
    !***********************************
    integer,intent(in) :: n,numctr,numnb,i,cval
    integer(int8),intent(in) :: pair(n,n,numctr)
    integer,intent(inout) :: nnbi,nbi(2,numnb)
    integer :: l,m

    nnbi = 0
    nbi = 0
    do l = 1,n
      do m = 1,numctr
        if (pair(l,i,m) .eq. cval) then
          nnbi = nnbi+1
          nbi(1,nnbi) = l  ! idx of nb,  j
          nbi(2,nnbi) = m  ! which cell, iTr
        end if
      end do
    end do
  end subroutine countf

end module gfnff_topo_bondmat
