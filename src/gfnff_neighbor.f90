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
!> PBC-aware neighbour list type TNeigh.
!> Ported from xtb's xtb_gfnff_neighbor module; all xtb-specific types replaced
!> by plain Fortran arrays and integer iostat flags.
module gfnff_neighbor
  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData
  implicit none
  private

  public :: TNeigh

  !> PBC-aware neighbour list
  type :: TNeigh

    !> nb(j,i,iTr): atom j is a neighbour of atom i shifted by translation
    !> vector iTr.  nb(numnb,i,iTr) = count of neighbours of i in cell iTr.
    integer,allocatable :: nb(:,:,:)

    !> Full neighbour list (including metals)
    integer,allocatable :: nbf(:,:,:)

    !> Non-metal / normally-coordinated neighbour list
    integer,allocatable :: nbm(:,:,:)

    !> Bond pair count: bpair(j,i,iTr) = number of bonds between i and j
    !> when j is translated by iTr
    integer,allocatable :: bpair(:,:,:)
    integer,allocatable :: molbpair(:)

    !> Maximum neighbours per atom (index numnb stores the count)
    integer :: numnb = 42

    !> Number of central cells (1 for molecular, 27 for 3D-PBC)
    integer :: numctr = 1

    !> iTrSum(lin(i,j)): local index of the translation vector that is the
    !> sum of the vectors at unique IDs i and j
    integer,allocatable :: iTrSum(:)

    !> iTrNeg(i): local index of the negated translation vector of i
    integer,allocatable :: iTrNeg(:)

    !> iTrUsed(iTr): .true. if bonds to this cell should be counted
    !> (avoids double-counting)
    logical,allocatable :: iTrUsed(:)

    integer :: nbond = 0
    integer :: nbond_blist = 0

    !> Bonded atom pairs
    integer,allocatable :: blist(:,:)

    !> Number of H-bonds per bond
    integer,allocatable :: nr_hb(:)

    !> Bond potential parameters (shift, steepness, prefactor)
    real(wp),allocatable :: vbond(:,:)

    !> Number of active translation vectors within cutoff
    integer :: nTrans = 0

    !> Dimension of iTrSum/iTrNeg (= size of trVecInt along dim 2)
    integer :: iTrDim = 0

    !> idTr(k): unique ID of the k-th active translation vector
    integer,allocatable :: idTr(:)

    !> Trid(id): local index k such that idTr(k) == id; 0 if outside cutoff
    integer,allocatable :: Trid(:)

    !> Integer coefficients of lattice-vector linear combinations
    integer,allocatable :: trVecInt(:,:)

    !> Cartesian translation vectors (3 x nTrans), in Bohr
    real(wp),allocatable :: transVec(:,:)

    !> Cached cutoff from last getTransVec call (no-op if unchanged)
    real(wp) :: oldCutOff = 0.0_wp

  contains

    procedure :: init_n
    procedure :: get_nb
    procedure :: getTransVec
    procedure :: nbLoc
    procedure :: jth_nb
    procedure :: id2v
    procedure :: fTrSum

  end type TNeigh

!========================================================================================!
!========================================================================================!
contains
!========================================================================================!
!========================================================================================!

  subroutine init_n(self,nat,npbc,io)
    !**********************************************
    !* Initialise TNeigh: set numctr, fill lookup
    !* tables iTrSum/iTrNeg/iTrUsed.
    !* Input: nat (not used here but kept for
    !*   future allocation), npbc (0 or 3;
    !*   other values emit a warning via io).
    !**********************************************
    class(TNeigh),intent(inout) :: self
    integer,intent(in)  :: nat
    integer,intent(in)  :: npbc
    integer,intent(out) :: io

    integer :: iTrDim

    io = 0

    select case (npbc)
    case (3)
      self%numctr = 27
    case (0)
      self%numctr = 1
    case default
      write(*,'(a,i0,a)') 'Warning (gfnff_neighbor/init_n): npbc=',npbc, &
        & ' — PBC only fully tested for 3D systems.'
      if (npbc > 0) then
        self%numctr = 27
      else
        self%numctr = 1
      end if
    end select

    ! Allocate trVecInt with a safe default size so that filliTrSum can run.
    ! getTransVec will resize it properly for the actual cutoff.
    if (.not.allocated(self%trVecInt)) then
      allocate(self%trVecInt(3,729),source=0)
    end if
    iTrDim = min(size(self%trVecInt,dim=2),729)
    self%iTrDim = iTrDim
    call filliTrSum(iTrDim,self%trVecInt,self%iTrSum,self%iTrNeg)
    call filliTrUsed(self)
    self%oldCutOff = 0.0_wp

  end subroutine init_n

!========================================================================================!

  subroutine get_nb(self,nat,at,xyz,lattice,rtmp,mchar,icase,f_in,f2_in,param)
    !**********************************************
    !* Compute all pairwise distances over the
    !* central numctr cells and fill the nb/nbf/nbm
    !* arrays accordingly via fillnb.
    !* icase selects which list to fill:
    !*   1 = nbf, 2 = nb, 3 = nbm
    !**********************************************
    class(TNeigh),intent(inout) :: self
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: lattice(3,3)
    real(wp),intent(in) :: rtmp(nat*(nat+1)/2)
    real(wp),intent(in) :: mchar(nat)
    integer,intent(in)  :: icase
    real(wp),intent(in) :: f_in, f2_in
    type(TGFFData),intent(in) :: param

    real(wp),allocatable :: dist(:,:,:)
    integer :: i,j,iTr

    allocate(dist(nat,nat,self%numctr),source=0.0_wp)
    !$omp parallel do collapse(3) default(none) &
    !$omp shared(dist,nat,xyz,self) private(iTr,i,j)
    do iTr = 1,self%numctr
      do i = 1,nat
        do j = 1,nat
          dist(j,i,iTr) = norm2(xyz(:,i) - (xyz(:,j) + self%transVec(:,iTr)))
        end do
      end do
    end do
    !$omp end parallel do
    call fillnb(self,nat,at,rtmp,dist,mchar,icase,f_in,f2_in,param)

  end subroutine get_nb

!========================================================================================!

  subroutine getTransVec(self,nat,npbc,pbc,lattice,cutoff,io)
    !**********************************************
    !* Build the list of translation vectors within
    !* the given cutoff.  Caches result; no-op if
    !* cutoff is unchanged from last call.
    !* For periodic systems the central 27 cells
    !* are always included regardless of cutoff.
    !**********************************************
    class(TNeigh),intent(inout) :: self
    integer,intent(in)  :: nat
    integer,intent(in)  :: npbc
    logical,intent(in)  :: pbc(3)
    real(wp),intent(in) :: lattice(3,3)
    real(wp),intent(in) :: cutoff
    integer,intent(out) :: io

    integer :: ix,iy,iz,i,nLat,nShell,d,indx,id
    integer,allocatable :: idTrTmp(:)
    real(wp),allocatable :: sVec(:),vec(:),trVecTmp(:,:)
    real(wp) :: sh

    io = 0

    if (self%oldCutOff .eq. cutoff) return
    self%oldCutOff = cutoff

    if (npbc .eq. 0) then
      if (.not.allocated(self%idTr)) allocate(self%idTr(1))
      if (.not.allocated(self%Trid)) allocate(self%Trid(1))
      self%nTrans = 1
      self%idTr = 1
      self%Trid = 1
      if (.not.allocated(self%trVecInt)) allocate(self%trVecInt(3,1),source=0)
      if (.not.allocated(self%transVec)) allocate(self%transVec(3,1),source=0.0_wp)
      return
    end if

    nLat = npbc

    ! Find shortest lattice vector
    allocate(sVec(3))
    sVec = 1.0e15_wp
    do i = 1,3
      if (.not.pbc(i)) cycle
      if (norm2(lattice(:,i)) < norm2(sVec)) sVec = lattice(:,i)
    end do
    if (norm2(sVec) .eq. 0.0_wp) then
      write(*,'(a)') 'Error (gfnff_neighbor/getTransVec): zero lattice vector on periodic axis.'
      io = 1
      return
    end if

    nShell = ceiling(cutoff / norm2(sVec))
    d = (2*nShell + 1)**nLat

    if (allocated(self%trVecInt)) then
      if (size(self%trVecInt,dim=2) < d) deallocate(self%trVecInt)
    end if

    if (.not.allocated(self%trVecInt)) then
      allocate(self%trVecInt(3,d),source=0)
      indx = 1
      if (nLat .eq. 3) then
        do i = 1,nShell
          do iz = -i,i
            do iy = -i,i
              do ix = -i,i
                if (abs(ix).ne.i .and. abs(iy).ne.i .and. abs(iz).ne.i) cycle
                indx = indx + 1
                self%trVecInt(:,indx) = [ix,iy,iz]
              end do
            end do
          end do
        end do
      elseif (nLat .eq. 2) then
        do i = 1,nShell
          do iz = -i,i
            do iy = -i,i
              do ix = -i,i
                if (.not.pbc(1)) then
                  if (abs(iy).ne.i .and. abs(iz).ne.i) cycle
                else if (.not.pbc(2)) then
                  if (abs(ix).ne.i .and. abs(iz).ne.i) cycle
                else if (.not.pbc(3)) then
                  if (abs(ix).ne.i .and. abs(iy).ne.i) cycle
                end if
                indx = indx + 1
                if (pbc(1)) self%trVecInt(1,indx) = ix
                if (pbc(2)) self%trVecInt(2,indx) = iy
                if (pbc(3)) self%trVecInt(3,indx) = iz
              end do
            end do
          end do
        end do
      elseif (nLat .eq. 1) then
        do i = 1,nShell
          do iz = -i,i
            do iy = -i,i
              do ix = -i,i
                if (pbc(1)) then
                  if (abs(ix).ne.i) cycle
                else if (pbc(2)) then
                  if (abs(iy).ne.i) cycle
                else if (pbc(3)) then
                  if (abs(iz).ne.i) cycle
                end if
                indx = indx + 1
                if (pbc(1)) self%trVecInt(1,indx) = ix
                if (pbc(2)) self%trVecInt(2,indx) = iy
                if (pbc(3)) self%trVecInt(3,indx) = iz
              end do
            end do
          end do
        end do
      else
        write(*,'(a)') 'Error (gfnff_neighbor/getTransVec): unsupported nLat.'
        io = 1
        return
      end if
    end if

    if (allocated(self%Trid)) deallocate(self%Trid)
    allocate(self%Trid(d),source=0)
    allocate(idTrTmp(d),source=0)
    allocate(vec(3),source=0.0_wp)
    indx = 0

    if (nLat .eq. 3) then
      allocate(trVecTmp(3,d))
      do id = 1,d
        vec = self%trVecInt(1,id)*lattice(:,1) &
            + self%trVecInt(2,id)*lattice(:,2) &
            + self%trVecInt(3,id)*lattice(:,3)
        if (id > 27) then
          if (id .eq. 28)   sh = 2.0_wp
          if (id .eq. 126)  sh = 3.0_wp
          if (id .eq. 344)  sh = 4.0_wp
          if (id .eq. 730)  sh = 5.0_wp
          if (id .eq. 1332) sh = 6.0_wp
          if (id .eq. 2198) sh = 7.0_wp
          if (id .eq. 3376) sh = 8.0_wp
          if (id .eq. 4914) sh = 9.0_wp
          if (id .eq. 6860) sh = 10.0_wp
          if (id .eq. 9262) sh = 11.0_wp
          if (id .eq. 12168) sh = 12.0_wp
          if (id .eq. 15626) sh = 13.0_wp
          if (id .eq. 19684) sh = 14.0_wp
          if (norm2(vec)*(1.0_wp - 1.0_wp/sh) > cutoff) cycle
        end if
        indx = indx + 1
        idTrTmp(indx) = id
        self%Trid(id) = indx
        trVecTmp(:,indx) = vec
      end do
      if (allocated(self%transVec)) deallocate(self%transVec)
      if (allocated(self%idTr))    deallocate(self%idTr)
      allocate(self%transVec(3,indx))
      allocate(self%idTr(indx),source=0)
      self%transVec = trVecTmp(:,1:indx)
      self%idTr     = idTrTmp(1:indx)
    elseif (nLat .eq. 2) then
      allocate(trVecTmp(3,d))
      do id = 1,d
        vec = self%trVecInt(1,id)*lattice(:,1) &
            + self%trVecInt(2,id)*lattice(:,2)
        if (id > 9) then
          if (id .eq. 10) sh = 2.0_wp
          if (norm2(vec) > cutoff) cycle
        end if
        indx = indx + 1
        idTrTmp(indx) = id
        trVecTmp(:,indx) = vec
      end do
      if (allocated(self%transVec)) deallocate(self%transVec)
      if (allocated(self%idTr))    deallocate(self%idTr)
      allocate(self%transVec(3,indx))
      allocate(self%idTr(indx),source=0)
      self%transVec = trVecTmp(:,1:indx)
      self%idTr     = idTrTmp(1:indx)
    elseif (nLat .eq. 1) then
      allocate(trVecTmp(3,d))
      do id = 1,d
        vec = self%trVecInt(1,id)*lattice(:,1)
        if (indx > 27) then
          if (norm2(vec) > cutoff) cycle
        end if
        indx = indx + 1
        idTrTmp(indx) = id
        trVecTmp(:,indx) = vec
      end do
      if (allocated(self%transVec)) deallocate(self%transVec)
      if (allocated(self%idTr))    deallocate(self%idTr)
      allocate(self%transVec(3,indx))
      allocate(self%idTr(indx),source=0)
      self%transVec = trVecTmp(:,1:indx)
      self%idTr     = idTrTmp(1:indx)
    end if

    self%nTrans = size(self%transVec,dim=2)

  end subroutine getTransVec

!========================================================================================!

  function id2v(self,lattice,id) result(vector)
    !**********************************************
    !* Convert a unique translation-vector ID to
    !* the corresponding Cartesian vector (Bohr).
    !**********************************************
    class(TNeigh),intent(inout) :: self
    real(wp),intent(in) :: lattice(3,3)
    integer,intent(in)  :: id
    real(wp) :: vector(3)

    vector = 0.0_wp
    if (id <= size(self%trVecInt,dim=2)) then
      vector = self%trVecInt(1,id)*lattice(:,1) &
             + self%trVecInt(2,id)*lattice(:,2) &
             + self%trVecInt(3,id)*lattice(:,3)
    end if

  end function id2v

!========================================================================================!

  subroutine nbLoc(self,n,nblist,indexi,locarr)
    !**********************************************
    !* Return a compact array of cells that contain
    !* neighbours of atom indexi.
    !* locarr(:,k): neighbour atom indices in cell k.
    !* locarr(numnb-1,k) = count, locarr(numnb,k) = iTr.
    !**********************************************
    class(TNeigh),intent(in) :: self
    integer,intent(in)  :: n
    integer,intent(in)  :: nblist(self%numnb,n,self%numctr)
    integer,intent(in)  :: indexi
    integer,allocatable,intent(out) :: locarr(:,:)

    integer :: j,iTr,counter,counter2

    counter = 0
    do iTr = 1,self%numctr
      if (nblist(self%numnb,indexi,iTr) .ne. 0) counter = counter + 1
    end do
    allocate(locarr(self%numnb,counter),source=0)
    counter2 = 0
    do iTr = 1,self%numctr
      if (nblist(self%numnb,indexi,iTr) .eq. 0) cycle
      counter = 0
      counter2 = counter2 + 1
      do j = 1,self%numnb - 2
        if (nblist(j,indexi,iTr) .ne. 0) then
          counter = counter + 1
          locarr(counter,counter2) = nblist(j,indexi,iTr)
          locarr(self%numnb,counter2)   = iTr
          locarr(self%numnb-1,counter2) = nblist(self%numnb,indexi,iTr)
        end if
      end do
    end do

  end subroutine nbLoc

!========================================================================================!

  subroutine jth_nb(self,n,xyz,j,jth,i,iTr)
    !**********************************************
    !* Return index j and cell iTr of the jth-
    !* nearest neighbour of atom i (sorted by
    !* distance).
    !**********************************************
    class(TNeigh),intent(inout) :: self
    integer,intent(in)    :: n
    real(wp),intent(in)   :: xyz(3,n)
    integer,intent(inout) :: j,iTr
    integer,intent(in)    :: jth,i

    integer :: k,l,iTrdum,sizenb,jdum
    integer,allocatable  :: indx(:)
    real(wp),allocatable :: distnb(:)

    sizenb = sum(self%nb(self%numnb,i,:))
    allocate(distnb(sizenb),source=0.0_wp)
    allocate(indx(sizenb),source=0)

    j   = 0
    iTr = 0
    l   = 0
    do iTrdum = 1,self%numctr
      do k = 1,self%nb(self%numnb,i,iTrdum)
        jdum = self%nb(k,i,iTrdum)
        l = l + 1
        distnb(l) = norm2(xyz(:,i) - (xyz(:,jdum) + self%transVec(:,iTrdum)))
      end do
    end do

    ! Simple insertion sort to get the sorted index array
    call indexSort(indx,distnb,sizenb)

    l = 0
    do iTrdum = 1,self%numctr
      do k = 1,self%nb(self%numnb,i,iTrdum)
        jdum = self%nb(k,i,iTrdum)
        l = l + 1
        if (l .eq. indx(jth)) then
          iTr = iTrdum
          j   = jdum
        end if
      end do
    end do

  end subroutine jth_nb

!========================================================================================!

  function fTrSum(self,indx1,indx2) result(indxSum)
    !**********************************************
    !* Return the local index of the translation
    !* vector that is the sum of the vectors at
    !* local indices indx1 and indx2.
    !* Returns -1 if the sum falls outside cutoff.
    !**********************************************
    class(TNeigh) :: self
    integer :: indx1,indx2
    integer :: indxSum

    integer :: id1,id2,k,idOfSum

    indxSum = -1
    if (abs(indx1) > self%nTrans .or. abs(indx2) > self%nTrans) return
    id1 = self%idTr(indx1)
    id2 = self%idTr(indx2)
    k = lin(id1,id2)
    if (k <= size(self%iTrSum)) then
      idOfSum = self%iTrSum(k)
      if (idOfSum <= 0) return
      if (idOfSum > size(self%Trid)) return
      if (self%Trid(idOfSum) > 0 .and. self%Trid(idOfSum) <= self%nTrans) then
        indxSum = self%Trid(idOfSum)
      end if
    end if

  end function fTrSum

!========================================================================================!
!> MODULE-PRIVATE HELPERS
!========================================================================================!

  subroutine fillnb(self,n,at,rad,dist,mchar,icase,f,f2,param)
    !**********************************************
    !* Inner loop: fill nb/nbf/nbm from the
    !* precomputed distance array dist(n,n,numctr).
    !* icase: 1=nbf, 2=nb, 3=nbm.
    !**********************************************
    implicit none
    type(TGFFData),intent(in)  :: param
    type(TNeigh),intent(inout) :: self
    integer,intent(in)  :: n,at(n),icase
    real(wp),intent(in) :: rad(n*(n+1)/2)
    real(wp),intent(in) :: dist(n,n,self%numctr)
    real(wp),intent(in) :: mchar(n),f,f2

    integer :: i,j,k,iTr,nn,hc_crit,nnfi,nnfj
    integer,allocatable :: tag(:,:,:)
    real(wp) :: rco,fm

    allocate(tag(n,n,self%numctr),source=0)

    if (icase .eq. 1) then
      if (.not.allocated(self%nbf)) then
        allocate(self%nbf(self%numnb,n,self%numctr),source=0)
      else
        self%nbf = 0
      end if
    end if
    if (icase .eq. 2) then
      if (.not.allocated(self%nb)) then
        allocate(self%nb(self%numnb,n,self%numctr),source=0)
      else
        self%nb = 0
      end if
    end if
    if (icase .eq. 3) then
      if (.not.allocated(self%nbm)) then
        allocate(self%nbm(self%numnb,n,self%numctr),source=0)
      else
        self%nbm = 0
      end if
    end if

    do iTr = 1,self%numctr
      do i = 1,n
        do j = 1,n
          k = lin(j,i)
          if (dist(j,i,iTr) <= 0.0_wp) cycle
          fm = 1.0_wp
          if (icase .eq. 1) then
            if (param%metal(at(i)) .eq. 2) fm = fm*f2
            if (param%metal(at(j)) .eq. 2) fm = fm*f2
            if (param%metal(at(i)) .eq. 1) fm = fm*(f2 + 0.025_wp)
            if (param%metal(at(j)) .eq. 1) fm = fm*(f2 + 0.025_wp)
          end if
          if (icase .eq. 2) then
            nnfi = sum(self%nbf(self%numnb,i,:))
            nnfj = sum(self%nbf(self%numnb,j,:))
            hc_crit = 6
            if (param%group(at(i)) <= 2) hc_crit = 4
            if (nnfi > hc_crit) cycle
            hc_crit = 6
            if (param%group(at(j)) <= 2) hc_crit = 4
            if (nnfj > hc_crit) cycle
          end if
          if (icase .eq. 3) then
            nnfi = sum(self%nbf(self%numnb,i,:))
            nnfj = sum(self%nbf(self%numnb,j,:))
            if (mchar(i) > 0.25_wp .or. param%metal(at(i)) > 0) cycle
            if (mchar(j) > 0.25_wp .or. param%metal(at(j)) > 0) cycle
            if (nnfi > param%normcn(at(i)) .and. at(i) > 10) cycle
            if (nnfj > param%normcn(at(j)) .and. at(j) > 10) cycle
          end if
          rco = rad(k)
          if (dist(j,i,iTr) < fm*f*rco) tag(j,i,iTr) = 1
        end do
      end do
    end do

    do iTr = 1,self%numctr
      do i = 1,n
        nn = 0
        do j = 1,n
          if (tag(j,i,iTr) .eq. 1) then
            nn = nn + 1
            if (icase .eq. 1) self%nbf(nn,i,iTr) = j
            if (icase .eq. 2) self%nb(nn,i,iTr)  = j
            if (icase .eq. 3) self%nbm(nn,i,iTr) = j
          end if
        end do
        if (icase .eq. 1) self%nbf(self%numnb,i,iTr) = nn
        if (icase .eq. 2) self%nb(self%numnb,i,iTr)  = nn
        if (icase .eq. 3) self%nbm(self%numnb,i,iTr) = nn
      end do
    end do

  end subroutine fillnb

!========================================================================================!

  subroutine filliTrSum(iTrDim,trVecInt,iTrSum,iTrNeg)
    !**********************************************
    !* Precompute iTrSum and iTrNeg lookup tables.
    !**********************************************
    integer,intent(in) :: iTrDim
    integer,intent(in) :: trVecInt(3,iTrDim)
    integer,allocatable,intent(out) :: iTrSum(:),iTrNeg(:)

    integer :: i,j,s,k,kk
    integer :: vec1(3),vec2(3),vsum(3)

    allocate(iTrSum(iTrDim*(iTrDim+1)/2),source=-1)
    allocate(iTrNeg(iTrDim),source=-1)

    if (iTrDim .eq. 1) then
      ! nothing to do for the molecular case
    else
      do i = 1,iTrDim
        kk = i*(i-1)/2
        do j = 1,i
          k    = kk + j
          vec1 = trVecInt(:,i)
          vec2 = trVecInt(:,j)
          vsum = vec1 + vec2
          do s = 1,iTrDim
            if (j .eq. 1) then
              if (vec1(1) .eq. -trVecInt(1,s) .and. &
                  vec1(2) .eq. -trVecInt(2,s) .and. &
                  vec1(3) .eq. -trVecInt(3,s)) then
                iTrNeg(i) = s
                iTrNeg(s) = i
              end if
            end if
            if (vsum(1) .eq. trVecInt(1,s) .and. &
                vsum(2) .eq. trVecInt(2,s) .and. &
                vsum(3) .eq. trVecInt(3,s)) then
              iTrSum(k) = s
            end if
          end do
        end do
      end do
    end if
    iTrNeg(1) = 1
    iTrSum(1) = 1

  end subroutine filliTrSum

!========================================================================================!

  subroutine filliTrUsed(neigh)
    !**********************************************
    !* Flag which cell directions to use to avoid
    !* double-counting bonded pairs.
    !**********************************************
    class(TNeigh),intent(inout) :: neigh

    integer :: iTr,iTrF

    if (allocated(neigh%iTrUsed)) deallocate(neigh%iTrUsed)
    allocate(neigh%iTrUsed(neigh%numctr))
    neigh%iTrUsed = .true.
    do iTr = neigh%numctr,2,-1
      if (neigh%iTrUsed(iTr)) then
        iTrF = neigh%iTrNeg(iTr)
        if (iTrF > 0 .and. iTrF <= neigh%numctr) neigh%iTrUsed(iTrF) = .false.
      end if
    end do

  end subroutine filliTrUsed

!========================================================================================!

  subroutine indexSort(indx,arr,n)
    !**********************************************
    !* Simple O(n^2) insertion sort returning the
    !* permutation index array.  Suitable for the
    !* small neighbour lists encountered here.
    !**********************************************
    integer,intent(in)  :: n
    real(wp),intent(in) :: arr(n)
    integer,intent(out) :: indx(n)

    integer :: i,j,tmp
    do i = 1,n
      indx(i) = i
    end do
    do i = 2,n
      j = i
      do while (j > 1)
        if (arr(indx(j)) < arr(indx(j-1))) then
          tmp = indx(j); indx(j) = indx(j-1); indx(j-1) = tmp
        else
          exit
        end if
        j = j - 1
      end do
    end do

  end subroutine indexSort

!========================================================================================!

  integer function lin(i1,i2)
    !**********************************************
    !* Linear index for upper triangular storage.
    !**********************************************
    integer,intent(in) :: i1,i2
    integer :: i,j
    if (i1 > i2) then
      i = i1; j = i2
    else
      i = i2; j = i1
    end if
    lin = i*(i-1)/2 + j
  end function lin

!========================================================================================!
end module gfnff_neighbor
