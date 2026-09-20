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

!> Hydrogen and halogen bond list construction.
!>
!> Two entry points with different jobs: gfnff_hbset0 counts candidates so the
!> lists can be allocated with headroom, and gfnff_hbset fills them. They are
!> kept in step deliberately; a screening change in one must be mirrored in
!> the other or the fill will overrun what the count reserved.
!>
!> The bond_hb* routines build the separate list of A-H...B units attached to
!> a bond, which the bonded hydrogen bond term needs.
module gfnff_topo_hbset
  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData,TGFFTopology,TGFFNeighbourList
  use gfnff_neighbor,only:TNeigh
!$ use omp_lib
  implicit none
  private

  public :: gfnff_hbset,gfnff_hbset0
  public :: bond_hbset,bond_hbset0
  public :: bond_hb_AHB_set,bond_hb_AHB_set1,bond_hb_AHB_set0
  public :: hbonds

contains  !> MODULE PROCEDURES START HERE

  subroutine gfnff_hbset(n,at,xyz,topo,neigh,nlist,hbthr1,hbthr2)
    !***********************************
    !* Fill nlist%hblist1/2/3 with HB/XB
    !* candidates (OpenMP parallel scan).
    !* Counting companion: gfnff_hbset0.
    !***********************************
    implicit none
    integer,intent(in) :: n
    integer,intent(in) :: at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    type(TGFFNeighbourList),intent(inout) :: nlist
    real(wp),intent(in) :: hbthr1,hbthr2

    integer :: i,j,k,nh,ix
    integer :: iTri,iTrj,iTrDum
    real(wp) :: rmsd,rab,rih,rjh,xi(3),xj(3)
    real(wp),allocatable :: rihl(:)
    logical,allocatable :: ibnd(:)
    logical :: ijnonbond
    integer :: nhb1,nhb2,nxb
!$  integer,parameter :: N_MAX_LIST = 800 !< keep approx. 32 kb of integer(int64)
!$  integer,allocatable :: hblist1(:,:),hblist2(:,:),hblist3(:,:)

    !>-- refresh only on first call or after a substantial move
    rmsd = sqrt(sum((xyz-nlist%hbrefgeo)**2))/dble(n)
    if (.not. (rmsd < 1.d-6.or.rmsd > 0.3d0)) return

    nlist%nhb1 = 0
    nlist%nhb2 = 0
    nlist%nxb = 0

    !$omp parallel default(none) &
    !$omp shared(topo, neigh, nlist, xyz, hbthr1, hbthr2) &
    !$omp private(iTri, iTrj, iTrDum, ix, i, j, k, nh, rab, rih, rjh) &
    !$omp private(xi, xj, rihl, ibnd) &
    !$omp private(ijnonbond, hblist1, hblist2, hblist3, nhb1, nhb2, nxb)

!$  allocate (hblist1(5,N_MAX_LIST),source=0)
!$  allocate (hblist2(5,N_MAX_LIST),source=0)
!$  allocate (hblist3(5,N_MAX_LIST),source=0)

    nhb1 = 0
    nhb2 = 0
    nxb = 0

    !>-- the A-B scan needs a hydrogen to find anything
    if (topo%nathbH > 0) then

      allocate (rihl(topo%nathbH),ibnd(topo%nathbH))

      !$omp do schedule(dynamic, 8)
      do ix = 1,topo%nathbAB
        i = topo%hbatABl(1,ix)
        j = topo%hbatABl(2,ix)
        do iTri = 1,neigh%nTrans ! go through i shifts
          xi(:) = xyz(1:3,i)+neigh%transVec(1:3,iTri)
          !>-- A...H distance and A-H bond flag are independent of the j
          !>   shift, so build them once per i shift instead of nTrans times
          do k = 1,topo%nathbH
            nh = topo%hbatHl(1,k)
            rihl(k) = sum((xyz(1:3,nh)-xi(:))**2)
            ibnd(k) = .false.
            if (iTri <= neigh%numctr) ibnd(k) = neigh%bpair(i,nh,iTri) == 1
          end do
          do iTrj = 1,neigh%nTrans ! go through j shifts
            xj(:) = xyz(1:3,j)+neigh%transVec(1:3,iTrj)
            rab = sum((xi(:)-xj(:))**2)
            if (rab > hbthr1) cycle
            ! combined shift for neigh% distances/bpair of the two shifted atoms
            iTrDum = neigh%fTrSum(neigh%iTrNeg(iTri),iTrj)
            if (iTrDum > neigh%nTrans.or.iTrDum < -1.or.iTrDum == 0) cycle ! invalid shift
            if (iTrDum <= neigh%numctr.and.iTrDum > 0) then
              ijnonbond = neigh%bpair(j,i,iTrDum) /= 1
            else
              ! i and j are not in neighboring cells
              ijnonbond = .true.
            end if
            do k = 1,topo%nathbH
              nh = topo%hbatHl(1,k) ! nh always in central cell
              ! i is the bonded A
              if (ibnd(k).and.ijnonbond) then
                nhb2 = nhb2+1
                hblist2(1,nhb2) = i
                hblist2(2,nhb2) = j
                hblist2(3,nhb2) = nh
                hblist2(4,nhb2) = iTri
                hblist2(5,nhb2) = iTrj
!$              if (nhb2 == N_MAX_LIST) call update_hblist2(nlist,nhb2,hblist2)
                cycle ! already counted, excluded from nhb1
              end if
              ! j is the bonded A
              if (iTrj <= neigh%numctr) then
                if (neigh%bpair(j,nh,iTrj) == 1.and.ijnonbond) then
                  nhb2 = nhb2+1
                  hblist2(1,nhb2) = j
                  hblist2(2,nhb2) = i
                  hblist2(3,nhb2) = nh
                  hblist2(4,nhb2) = iTrj
                  hblist2(5,nhb2) = iTri
!$                if (nhb2 == N_MAX_LIST) call update_hblist2(nlist,nhb2,hblist2)
                  cycle ! already counted, excluded from nhb1
                end if
              end if
              ! neither A nor B is covalently bonded to H
              rjh = sum((xyz(1:3,nh)-xj(:))**2)
              if (rab+rihl(k)+rjh < hbthr2) then
                nhb1 = nhb1+1
                hblist1(1,nhb1) = i
                hblist1(2,nhb1) = j
                hblist1(3,nhb1) = nh
                hblist1(4,nhb1) = iTri
                hblist1(5,nhb1) = iTrj
!$              if (nhb1 == N_MAX_LIST) call update_hblist1(nlist,nhb1,hblist1)
              end if
            end do ! k
          end do ! iTrj
        end do ! iTri
      end do ! ix
      !$omp end do nowait

      deallocate (rihl,ibnd)

    end if

    !>-- for the nxb list only i is not shifted
    !$omp do schedule(dynamic, 16)
    do ix = 1,topo%natxbAB
      i = topo%xbatABl(1,ix) ! A
      j = topo%xbatABl(2,ix) ! B
      iTrj = topo%xbatABl(4,ix) ! iTrB
      if (iTrj > neigh%nTrans.or.iTrj < -1.or.iTrj == 0) cycle ! invalid shift
      rab = sum((xyz(1:3,j)-xyz(1:3,i)+neigh%transVec(1:3,iTrj))**2)
      if (rab > hbthr2) cycle
      nxb = nxb+1
      hblist3(1,nxb) = i ! A
      hblist3(2,nxb) = j ! B
      hblist3(3,nxb) = topo%xbatABl(3,ix) ! X
      hblist3(4,nxb) = iTrj
      hblist3(5,nxb) = topo%xbatABl(5,ix) ! iTrX
!$    if (nxb == N_MAX_LIST) call update_hblist3(nlist,nxb,hblist3)
    end do
    !$omp end do nowait

    call update_hblist1(nlist,nhb1,hblist1)
    call update_hblist2(nlist,nhb2,hblist2)
    call update_hblist3(nlist,nxb,hblist3)

!$  deallocate (hblist1,hblist2,hblist3)

    !$omp end parallel

    nlist%hbrefgeo = xyz

  contains

    subroutine update_hblist1(neigh_list,nhb,hblist)
      type(TGFFNeighbourList),intent(inout) :: neigh_list
      integer,intent(inout) :: nhb
      integer,intent(inout) :: hblist(5,nhb)
      !$omp critical (list1)
!$    neigh_list%hblist1(:,nlist%nhb1+1:nlist%nhb1+nhb) = hblist(:,1:nhb)
      neigh_list%nhb1 = neigh_list%nhb1+nhb
      !$omp end critical (list1)
!$    nhb = 0
!$    hblist = 0
    end subroutine update_hblist1

    subroutine update_hblist2(neigh_list,nhb,hblist)
      type(TGFFNeighbourList),intent(inout) :: neigh_list
      integer,intent(inout) :: nhb
      integer,intent(inout) :: hblist(5,nhb)
      !$omp critical (list2)
!$    neigh_list%hblist2(:,nlist%nhb2+1:nlist%nhb2+nhb) = hblist(:,1:nhb)
      neigh_list%nhb2 = neigh_list%nhb2+nhb
      !$omp end critical (list2)
!$    nhb = 0
!$    hblist = 0
    end subroutine update_hblist2

    subroutine update_hblist3(neigh_list,nxb,hblist)
      type(TGFFNeighbourList),intent(inout) :: neigh_list
      integer,intent(inout) :: nxb
      integer,intent(inout) :: hblist(5,nxb)
      !$omp critical (list3)
!$    neigh_list%hblist3(:,nlist%nxb+1:nlist%nxb+nxb) = hblist(:,1:nxb)
      neigh_list%nxb = neigh_list%nxb+nxb
      !$omp end critical (list3)
!$    nxb = 0
!$    hblist = 0
    end subroutine update_hblist3

  end subroutine gfnff_hbset

  subroutine bond_hbset(n,at,xyz,npbc,bond_hbn,bond_hbl,topo,neigh,hbthr1,hbthr2)
    !***********************************
    !* Fill bond_hbl with the A-H...B
    !* triplets attached to each bond.
    !* Columns: A,B,H,iTrA,iTrB,iTrH.
    !* Counting companion: bond_hbset0.
    !***********************************
    implicit none
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    integer,intent(in) :: n
    integer,intent(in) :: at(n)
    integer,intent(in) :: bond_hbn
    integer,intent(out) :: bond_hbl(6,bond_hbn)
    real(wp),intent(in) :: xyz(3,n)
    integer,intent(in) :: npbc
    real(wp),intent(in) :: hbthr1,hbthr2

    integer i,j,k,nh,ix,iTri,iTrj,iTrDum
    integer bond_nr
    real(wp) rab
    logical ijnonbond
    real(wp) :: vab(3)

    bond_nr = 0
    bond_hbl = 0
    do ix = 1,topo%nathbAB
      !>-- get i and j from the AB list; which is A and which is B is not
      !>   decided yet
      i = topo%hbatABl(1,ix)
      j = topo%hbatABl(2,ix)
      do iTri = 1,neigh%nTrans
        do iTrj = 1,neigh%nTrans
          iTrDum = neigh%fTrSum(neigh%iTrNeg(iTri),iTrj)
          vab = (xyz(:,i)+neigh%transVec(:,iTri))-(xyz(:,j)+neigh%transVec(:,iTrj))
          rab = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
          if (rab .gt. hbthr1) cycle
          if (iTrDum .le. neigh%numctr.and.iTrDum .gt. 0) then
            ijnonbond = neigh%bpair(j,i,iTrDum) .ne. 1
          else
            ! i and j are not in neighboring cells
            ijnonbond = .true.
          end if
          do k = 1,topo%nathbH
            nh = topo%hbatHl(1,k)  ! ALWAYS in central cell
            !>-- i is the A; only possible if A is adjacent to the central
            !>   cell, since H always is
            if (iTri .le. neigh%numctr) then
              if (neigh%bpair(i,nh,iTri) .eq. 1.and.ijnonbond) then
                bond_nr = bond_nr+1
                bond_hbl(1,bond_nr) = i     ! 1=A
                bond_hbl(2,bond_nr) = j     ! 2=B
                bond_hbl(3,bond_nr) = nh    ! 3=H
                bond_hbl(4,bond_nr) = iTri  ! 4=iTrA
                bond_hbl(5,bond_nr) = iTrj  ! 5=iTrB
                bond_hbl(6,bond_nr) = 1     ! 6=iTrH
              end if
            end if
            ! j is the A
            if (iTrj .le. neigh%numctr) then
              if (neigh%bpair(j,nh,iTrj) .eq. 1.and.ijnonbond) then
                bond_nr = bond_nr+1
                bond_hbl(1,bond_nr) = j     ! 1=A
                bond_hbl(2,bond_nr) = i     ! 2=B
                bond_hbl(3,bond_nr) = nh    ! 3=H
                bond_hbl(4,bond_nr) = iTrj  ! 4=iTrA
                bond_hbl(5,bond_nr) = iTri  ! 5=iTrB
                bond_hbl(6,bond_nr) = 1     ! 6=iTrH
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine bond_hbset

  subroutine bond_hbset0(n,at,xyz,npbc,bond_hbn,topo,neigh,hbthr1,hbthr2)
    !***********************************
    !* Count the A-H...B triplets bond_hbset
    !* will fill, to size bond_hbl. Branch
    !* logic mirrors bond_hbset.
    !***********************************
    implicit none
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh
    integer,intent(in) :: n
    integer,intent(in) :: at(n)
    integer,intent(out) :: bond_hbn
    real(wp),intent(in) :: xyz(3,n)
    integer,intent(in) :: npbc
    real(wp),intent(in) :: hbthr1,hbthr2

    integer i,j,k,nh,ix,iTri,iTrj,iTrDum
    real(wp) rab
    logical ijnonbond
    real(wp) :: vab(3)

    bond_hbn = 0
    do ix = 1,topo%nathbAB
      i = topo%hbatABl(1,ix)
      j = topo%hbatABl(2,ix)
      do iTri = 1,neigh%nTrans
        do iTrj = 1,neigh%nTrans
          iTrDum = neigh%fTrSum(neigh%iTrNeg(iTri),iTrj)
          if (iTrDum .eq. -1.or.iTrDum .gt. neigh%numctr) then
            vab = (xyz(:,i)+neigh%transVec(:,iTri))-(xyz(:,j)+neigh%transVec(:,iTrj))
            rab = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
            if (rab .gt. hbthr1) cycle
            ijnonbond = .true. ! not in neighboring or same cell for iTrDum=-1
          else
            vab = xyz(:,i)-(xyz(:,j)+neigh%transVec(:,iTrDum))
            rab = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
            if (rab .gt. hbthr1) cycle
            ijnonbond = neigh%bpair(j,i,iTrDum) .ne. 1
          end if
          do k = 1,topo%nathbH
            nh = topo%hbatHl(1,k)  ! ALWAYS in central cell
            ! i is the A
            if (iTri .le. neigh%numctr) then
              if (neigh%bpair(i,nh,iTri) .eq. 1.and.ijnonbond) then
                bond_hbn = bond_hbn+1
              end if
            end if
            ! j is the A
            if (iTrj .le. neigh%numctr) then
              if (neigh%bpair(j,nh,iTrj) .eq. 1.and.ijnonbond) then
                bond_hbn = bond_hbn+1
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine bond_hbset0

  subroutine bond_hb_AHB_set(n,at,numbond,bond_hbn,bond_hbl,tot_AHB_nr,lin_AHB,topo,neigh)
    !***********************************
    !* Group the A-H...B triplets from
    !* bond_hbset by shared A-H unit into
    !* topo%bond_hb_AH/B/Bn. Sizes come
    !* from bond_hb_AHB_set0 and set1.
    !***********************************
    implicit none
    type(TNeigh),intent(inout) :: neigh
    type(TGFFTopology),intent(inout) :: topo
    integer,intent(in)  :: n
    integer,intent(in)  :: numbond
    integer,intent(in)  :: at(n)
    integer,intent(in)  :: bond_hbn
    integer,intent(in)  :: bond_hbl(6,bond_hbn)
    integer,intent(in)  :: tot_AHB_nr
    integer,intent(inout) :: lin_AHB(4,0:tot_AHB_nr)
    integer :: i,j
    integer :: ii,jj,iTr,iTrA,iTrH,iTrB
    integer :: ia,ja
    integer :: hbH,hbA
    integer :: Hat,Aat
    integer :: Bat,atB
    integer :: tot_count,t
    integer :: AH_count
    integer :: B_count
    integer :: lin_diff
    logical :: iTr_same

    tot_count = 0
    AH_count = 0
    B_count = 0
    lin_diff = 0

    do i = 1,numbond
      jj = neigh%blist(1,i)
      ii = neigh%blist(2,i)
      iTr = neigh%blist(3,i)
      ia = at(ii)
      ja = at(jj)
      if (ia .eq. 1) then
        hbH = ii
        hbA = jj
      else if (ja .eq. 1) then
        hbH = jj
        hbA = ii
      else
        cycle
      end if
      if (at(hbA) .eq. 7.or.at(hbA) .eq. 8) then
        do j = 1,bond_hbn
          Bat = bond_hbl(2,j)
          iTrB = bond_hbl(5,j)
          atB = at(Bat)
          Aat = bond_hbl(1,j)
          iTrA = bond_hbl(4,j)
          Hat = bond_hbl(3,j)
          iTrH = bond_hbl(6,j) ! always 1
          if ((hbA .eq. Aat.and.hbH .eq. Hat.and.iTr .eq. iTrA.and.iTrH .eq. 1).or.&
              &hbA .eq. Aat.and.hbH .eq. Hat.and.iTr .eq. iTrH.and.iTrA .eq. 1) then
            if (atB .eq. 7.or.atB .eq. 8) then
              tot_count = tot_count+1
              t = tot_count
              lin_AHB(1,tot_count) = hbA
              lin_AHB(2,tot_count) = hbH
              lin_AHB(3,tot_count) = iTrA
              lin_AHB(4,tot_count) = iTrH
              if (lin_AHB(1,tot_count)-lin_AHB(1,tot_count-1) .eq. 0.and.&
                 &lin_AHB(2,tot_count)-lin_AHB(2,tot_count-1) .eq. 0) then
                lin_diff = 0
              else
                lin_diff = 1
              end if
              iTr_same = iTrA .eq. lin_AHB(3,tot_count-1).and.iTrH .eq. lin_AHB(4,tot_count-1)
              if (lin_diff .eq. 0.and.iTr_same) B_count = B_count+1
              if (lin_diff .ne. 0.or..not.iTr_same) then
                AH_count = AH_count+1
                topo%bond_hb_AH(1,AH_count) = hbA
                topo%bond_hb_AH(2,AH_count) = hbH
                topo%bond_hb_AH(3,AH_count) = iTrA
                topo%bond_hb_AH(4,AH_count) = iTrH
                B_count = 1
              end if
              topo%bond_hb_Bn(AH_count) = B_count
              topo%bond_hb_B(1,B_count,AH_count) = Bat
              topo%bond_hb_B(2,B_count,AH_count) = iTrB
            end if
          else
            cycle
          end if
        end do
      end if
    end do

  end subroutine bond_hb_AHB_set

  subroutine bond_hb_AHB_set1(n,at,numbond,bond_hbn,bond_hbl,tot_AHB_nr,lin_AHB,AH_count,bmax,topo,neigh)
    !***********************************
    !* Sizing pass for bond_hb_AHB_set:
    !* same matching, returns AH_count and
    !* bmax to allocate topo%bond_hb_AH/B/Bn.
    !* Also sets topo%isABH and neigh%nr_hb.
    !***********************************
    implicit none
    type(TNeigh),intent(inout) :: neigh
    type(TGFFTopology),intent(inout) :: topo
    integer,intent(in)  :: n
    integer,intent(in)  :: numbond
    integer,intent(in)  :: at(n)
    integer,intent(in)  :: bond_hbn
    integer,intent(in)  :: bond_hbl(6,bond_hbn)
    integer,intent(in)  :: tot_AHB_nr
    integer,intent(inout) :: lin_AHB(4,0:tot_AHB_nr)
    integer,intent(out) :: AH_count
    integer,intent(out) :: bmax
    integer :: i,j
    integer :: ii,jj,iTr,iTrA,iTrH
    integer :: ia,ja
    integer :: hbH,hbA
    integer :: Hat,Aat
    integer :: Bat,atB
    integer :: tot_count
    integer :: B_count
    integer :: lin_diff
    logical :: iTr_same

    tot_count = 0
    AH_count = 0
    B_count = 1
    bmax = 1
    lin_diff = 0

    do i = 1,numbond ! order must follow blist, indexes neigh%nr_hb(i)
      jj = neigh%blist(1,i)
      ii = neigh%blist(2,i)
      iTr = neigh%blist(3,i)
      ia = at(ii)
      ja = at(jj)
      if (ia .eq. 1) then
        hbH = ii
        hbA = jj
      else if (ja .eq. 1) then
        hbH = jj
        hbA = ii
      else
        cycle
      end if
      if (at(hbA) .eq. 7.or.at(hbA) .eq. 8) then
        do j = 1,bond_hbn
          Bat = bond_hbl(2,j)
          atB = at(Bat)
          Aat = bond_hbl(1,j)
          iTrA = bond_hbl(4,j)
          Hat = bond_hbl(3,j)
          iTrH = bond_hbl(6,j) ! always 1
          if ((hbA .eq. Aat.and.hbH .eq. Hat.and.iTr .eq. iTrA.and.iTrH .eq. 1).or.&
              &(hbA .eq. Aat.and.hbH .eq. Hat.and.iTr .eq. iTrH.and.iTrA .eq. 1)) then
            if (atB .eq. 7.or.atB .eq. 8) then
              tot_count = tot_count+1
              lin_AHB(1,tot_count) = hbA
              lin_AHB(2,tot_count) = hbH
              lin_AHB(3,tot_count) = iTrA
              lin_AHB(4,tot_count) = iTrH
              topo%isABH(Bat) = .true.
              if (lin_AHB(1,tot_count)-lin_AHB(1,tot_count-1) .eq. 0.and.&
                 &lin_AHB(2,tot_count)-lin_AHB(2,tot_count-1) .eq. 0) then
                lin_diff = 0
              else
                lin_diff = 1
              end if
              iTr_same = iTrA .eq. lin_AHB(3,tot_count-1).and.iTrH .eq. lin_AHB(4,tot_count-1)
              if (lin_diff .eq. 0.and.iTr_same) B_count = B_count+1
              if (lin_diff .ne. 0.or..not.iTr_same) then
                AH_count = AH_count+1
                B_count = 1
              end if
              if (B_count .gt. bmax) bmax = B_count
            end if
          else
            cycle
          end if
        end do
        topo%isABH(hbA) = .true.
        topo%isABH(hbH) = .true.
        neigh%nr_hb(i) = B_count
      end if
    end do

  end subroutine bond_hb_AHB_set1

  subroutine bond_hb_AHB_set0(n,at,numbond,bond_hbn,bond_hbl,tot_AHB_nr,neigh)
    !***********************************
    !* Count the A-H...B matches that
    !* bond_hb_AHB_set1/set will walk,
    !* to size lin_AHB.
    !***********************************
    implicit none
    type(TNeigh),intent(inout) :: neigh
    integer,intent(in)  :: n
    integer,intent(in)  :: numbond
    integer,intent(in)  :: at(n)
    integer,intent(in)  :: bond_hbn
    integer,intent(in)  :: bond_hbl(6,bond_hbn)
    integer,intent(out) :: tot_AHB_nr
    integer :: i,j
    integer :: ii,jj,iTr,iTrA,iTrH
    integer :: ia,ja
    integer :: hbH,hbA
    integer :: Hat,Aat
    integer :: Bat,atB

    tot_AHB_nr = 0
    do i = 1,numbond
      jj = neigh%blist(1,i)
      ii = neigh%blist(2,i)
      iTr = neigh%blist(3,i)
      ia = at(ii)
      ja = at(jj)
      if (ia .eq. 1) then
        hbH = ii
        hbA = jj
      else if (ja .eq. 1) then
        hbH = jj
        hbA = ii
      else
        cycle
      end if
      if (at(hbA) .eq. 7.or.at(hbA) .eq. 8) then
        do j = 1,bond_hbn
          Bat = bond_hbl(2,j)
          atB = at(Bat)
          Aat = bond_hbl(1,j)
          iTrA = bond_hbl(4,j)
          Hat = bond_hbl(3,j)
          iTrH = bond_hbl(6,j) ! always 1
          if ((hbA .eq. Aat.and.hbH .eq. Hat.and.iTr .eq. iTrA.and.iTrH .eq. 1).or.&
              &hbA .eq. Aat.and.hbH .eq. Hat.and.iTr .eq. iTrH.and.iTrA .eq. 1) then
            if (atB .eq. 7.or.atB .eq. 8) then
              tot_AHB_nr = tot_AHB_nr+1
            end if
          else
            cycle
          end if
        end do
      end if
    end do

  end subroutine bond_hb_AHB_set0

  subroutine gfnff_hbset0(n,at,xyz,topo,nhb1,nhb2,nxb,neigh,nlist,hbthr1,hbthr2)
    !***********************************
    !* Count the HB/XB candidates
    !* gfnff_hbset would find, without
    !* filling any lists. Used to check
    !* whether hblist1/2/3 still fit.
    !***********************************
    implicit none
    integer,intent(in) :: n
    integer,intent(in) :: at(n)
    real(wp),intent(in) :: xyz(3,n)
    type(TGFFTopology),intent(in) :: topo
    integer,intent(out) :: nhb1
    integer,intent(out) :: nhb2
    integer,intent(out) :: nxb
    type(TNeigh),intent(in) :: neigh
    type(TGFFNeighbourList),intent(in) :: nlist
    real(wp),intent(in) :: hbthr1,hbthr2

    integer :: i,j,k,nh,ix
    integer :: iTri,iTrj,iTrDum
    logical :: ijnonbond
    real(wp) :: rab,rih,rjh,xi(3),xj(3)
    real(wp),allocatable :: rihl(:)
    logical,allocatable :: ibnd(:)

    nhb1 = 0
    nhb2 = 0
    nxb = 0

    !$omp parallel default(none) &
    !$omp reduction(+:nhb1, nhb2, nxb) &
    !$omp shared(topo, neigh, xyz, hbthr1, hbthr2) &
    !$omp private(iTri, iTrj, iTrDum, ix, i, j, k, nh, rab, rih, rjh, ijnonbond) &
    !$omp private(xi, xj, rihl, ibnd)

    !>-- the A-B scan needs a hydrogen to find anything
    if (topo%nathbH > 0) then

      allocate (rihl(topo%nathbH),ibnd(topo%nathbH))

      !$omp do schedule(dynamic, 8)
      do ix = 1,topo%nathbAB
        i = topo%hbatABl(1,ix)
        j = topo%hbatABl(2,ix)
        do iTri = 1,neigh%nTrans ! go through i shifts
          xi(:) = xyz(1:3,i)+neigh%transVec(1:3,iTri)
          !>-- A...H distance and A-H bond flag are independent of the j
          !>   shift, so build them once per i shift instead of nTrans times
          do k = 1,topo%nathbH
            nh = topo%hbatHl(1,k)
            rihl(k) = sum((xyz(1:3,nh)-xi(:))**2)
            ibnd(k) = .false.
            if (iTri <= neigh%numctr) ibnd(k) = neigh%bpair(i,nh,iTri) == 1
          end do
          do iTrj = 1,neigh%nTrans ! go through j shifts
            xj(:) = xyz(1:3,j)+neigh%transVec(1:3,iTrj)
            rab = sum((xi(:)-xj(:))**2)
            if (rab > hbthr1) cycle
            ! combined shift for neigh% distances/bpair of the two shifted atoms
            iTrDum = neigh%fTrSum(neigh%iTrNeg(iTri),iTrj)
            if (iTrDum > neigh%nTrans.or.iTrDum < -1.or.iTrDum == 0) cycle
            if (iTrDum <= neigh%numctr.and.iTrDum > 0) then
              ijnonbond = neigh%bpair(j,i,iTrDum) /= 1
            else
              ! i and j are not in neighboring cells
              ijnonbond = .true.
            end if
            do k = 1,topo%nathbH
              ! i is the bonded A
              if (ibnd(k).and.ijnonbond) then
                nhb2 = nhb2+1
                cycle
              end if
              nh = topo%hbatHl(1,k) ! nh always in central cell
              ! j is the bonded A
              if (iTrj <= neigh%numctr) then
                if (neigh%bpair(j,nh,iTrj) == 1.and.ijnonbond) then
                  nhb2 = nhb2+1
                  cycle
                end if
              end if
              ! neither A nor B is covalently bonded to H
              rjh = sum((xyz(1:3,nh)-xj(:))**2)
              if (rab+rihl(k)+rjh < hbthr2) then
                nhb1 = nhb1+1
              end if
            end do ! k
          end do ! iTrj
        end do ! iTri
      end do ! ix
      !$omp end do nowait

      deallocate (rihl,ibnd)

    end if

    !$omp do schedule(dynamic)
    do ix = 1,topo%natxbAB
      i = topo%xbatABl(1,ix)
      j = topo%xbatABl(2,ix)
      iTrj = topo%xbatABl(4,ix)
      if (iTrj > neigh%nTrans.or.iTrj < -1.or.iTrj == 0) cycle
      rab = sum((xyz(1:3,j)-(xyz(1:3,i)+neigh%transVec(1:3,iTrj)))**2)
      if (rab > hbthr2) cycle
      nxb = nxb+1
    end do
    !$omp end do

    !$omp end parallel

  end subroutine gfnff_hbset0

  subroutine hbonds(i,j,ci,cj,param,topo)
    !***********************************
    !* HB strength for atoms i,j:
    !* (1)=basicity, (2)=acidity.
    !***********************************
    implicit none
    type(TGFFTopology),intent(in) :: topo
    type(TGFFData),intent(in) :: param
    integer i,j
    real(wp) ci(2),cj(2)
    ci(1) = topo%hbbas(i)
    cj(1) = topo%hbbas(j)
    ci(2) = topo%hbaci(i)
    cj(2) = topo%hbaci(j)
  end subroutine hbonds

end module gfnff_topo_hbset
