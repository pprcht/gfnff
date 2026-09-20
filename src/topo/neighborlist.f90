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

!> Builds the neighbour list the whole topology is derived from.
!>
!> gfnff_neigh applies the covalent-radius criterion twice: once with plain
!> radii, once with CN-corrected radii, so a bond from an over-coordinated
!> first guess can still be dropped. getnb does the criterion test for one
!> case; nn_nearest_noM finds the nearest non-metal neighbour for the metal
!> special cases.
module gfnff_topo_neighborlist
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_data_types,only:TGFFData,TGFFTopology,TCell
  use gfnff_neighbor,only:TNeigh
  use gfnff_geometry,only:lin,banglPBC
  use gfnff_rab,only:gfnffrab
!$ use omp_lib
  implicit none
  private

  public :: gfnff_neigh,getnb,nn_nearest_noM

contains  !> MODULE PROCEDURES START HERE

  subroutine gfnff_neigh(makeneighbor,natoms,at,xyz,cell,rab,fq,f_in,f2_in,lintr, &
                        & mchar,hyb,itag,param,topo,neigh,nb_call,printlevel,printunit)
    !***********************************************************************
    !* Determine hybridisation states and neighbour lists.
    !* printlevel >= 1: print warnings/errors (bond, hybridisation, CN)
    !* printunit: output unit (default: stdout)
    !***********************************************************************
    implicit none
    character(len=*),parameter :: source = 'gfnff_topo_neighborlist'
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(inout) :: topo
    type(TNeigh),intent(inout) :: neigh ! contains nb, nbf and nbm
    logical,intent(in) :: makeneighbor,nb_call
    integer,intent(in) :: natoms
    integer,intent(in) :: at(natoms)
    integer :: hyb(natoms)
    integer :: itag(natoms)
    real(wp) :: rab(natoms*(natoms+1)/2)
    real(wp) :: xyz(3,natoms)
    type(TCell),intent(in) :: cell
    real(wp) :: mchar(natoms)
    real(wp) :: fq
    real(wp) :: f_in,f2_in               ! radius scaling for atoms/metal atoms recpectively
    real(wp) :: lintr                    ! threshold for linearity
    integer,intent(in),optional :: printlevel  !< verbosity (0=silent,1=errors,2=info,3=verbose)
    integer,intent(in),optional :: printunit   !< output unit (default: stdout)

    integer :: mylevel,myunit
    logical :: etacoord
    integer,allocatable :: nbdum(:,:,:),nbdum2(:,:),locarr(:,:)
    real(wp),allocatable :: cn(:),rtmp(:)
    integer :: i,j,k,jj,kk,ll,ati,nb20i,nbdiff,nbmdiff,nni,nh,nm
    integer :: ai,aj,nn,im,ncm,l,no,iTr,iTr2,numnbf,numnbm,numnb,idxdum,idxdum2,numctr
    integer :: nat
    real(wp) :: f1,phi,f2,rco,fat(103)
    real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
    data fat/103*1.0d0/

!>-- hand-tuned per-element radius scale factors
    fat(1) = 1.02
    fat(4) = 1.03
    fat(5) = 1.02
    fat(8) = 1.02
    fat(9) = 1.05
    fat(10) = 1.10
    fat(11) = 1.01
    fat(12) = 1.02
    fat(15) = 0.97
    fat(18) = 1.10
    fat(19) = 1.02
    fat(20) = 1.02
    fat(38) = 1.02
    fat(34) = 0.99
    fat(50) = 1.01
    fat(51) = 0.99
    fat(52) = 0.95
    fat(53) = 0.98
    fat(56) = 1.02
    fat(76) = 1.02
    fat(82) = 1.06
    fat(83) = 0.95

    mylevel = 0
    if (present(printlevel)) mylevel = printlevel
    if (present(printunit)) then
      myunit = printunit
    else
      myunit = stdout
    end if

    nat = natoms
    allocate (cn(natoms),rtmp(natoms*(natoms+1)/2),nbdum2(20,natoms))
    rtmp = 0.0

!>-- determine the neighbor list
    if (makeneighbor) then

      do i = 1,natoms
        cn(i) = dble(param%normcn(at(i)))
      end do
      call gfnffrab(natoms,at,cn,rtmp) ! guess RAB based on "normal" CN
      do i = 1,natoms
        ai = at(i)
        f1 = fq
        if (param%metal(ai) > 0) f1 = f1*2.0d0
        do j = 1,i
          f2 = fq
          aj = at(j)
          if (param%metal(aj) > 0) f2 = f2*2.0d0
          k = lin(j,i)
          rco = rtmp(k)
          rtmp(k) = rtmp(k)-topo%qa(i)*f1-topo%qa(j)*f2 ! change radius of atom i and j with charge
          rtmp(k) = rtmp(k)*fat(ai)*fat(aj)
        end do
      end do

      call neigh%get_nb(nat,at,xyz,cell,rab,rtmp,mchar,1,f_in,f2_in,param) ! nbf
      !>-- neigh%nb is only used for hyb states here; nbf overwrites it later
      call neigh%get_nb(nat,at,xyz,cell,rab,rtmp,mchar,2,f_in,f2_in,param) ! nb
      call neigh%get_nb(nat,at,xyz,cell,rab,rtmp,mchar,3,f_in,f2_in,param) ! nbm

      !>-- reuse the caller-supplied neighbor list instead of rebuilding it
    else

      neigh%nbf = neigh%nb
      neigh%nbm = neigh%nb

    end if

    itag = 0 ! save special hyb info
    numctr = neigh%numctr ! number of central cells considered (e.g. 1 for molec case)

!>-- tag atoms in nb(19,i) if they belong to a cluster (avoids the ring search)
    do i = 1,natoms
      if (sum(neigh%nbf(neigh%numnb,i,:)) .eq. 0.and.param%group(at(i)) .ne. 8) then
        if (mylevel >= 1) then
          write (myunit,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'')')
          write (myunit,'(''  warning: no bond partners for atom'',i4)') i
          write (myunit,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'')')
        end if
      end if
      if (at(i) .lt. 11.and.sum(neigh%nbf(neigh%numnb,i,:)) .gt. 2) then
        do iTr = 1,numctr
          do k = 1,neigh%nbf(neigh%numnb,i,iTr)
            kk = neigh%nbf(k,i,iTr)
            if (param%metal(at(kk)) .ne. 0.or.sum(neigh%nb(neigh%numnb,kk,:)) .gt. 4) then
              neigh%nb(neigh%numnb-1,i,1) = 1  ! ring search is limited to unit cell.
              neigh%nbf(neigh%numnb-1,i,1) = 1  ! Assumption: If the conditions are true
              neigh%nbm(neigh%numnb-1,i,1) = 1  ! in one cell they are true in all cells
            end if
          end do
        end do
      end if
    end do

!>-- hybridization states
    if (.not.allocated(nbdum)) &
      & allocate (nbdum(neigh%numnb,nat,numctr),source=0)
    do i = 1,natoms
      ati = at(i)
      numnbf = sum(neigh%nbf(neigh%numnb,i,:))
      numnbm = sum(neigh%nbm(neigh%numnb,i,:))
!>-- detect pi bonding to a metal, so hyb is taken from the reduced
!>   (metal-free) neighbor list
      etacoord = .false.
      if (ati .le. 10) then
        if (ati .eq. 6.and.numnbf .ge. 4.and.numnbm .eq. 3) etacoord = .true.  ! CP case
        if (ati .eq. 6.and.numnbf .eq. 3.and.numnbm .eq. 2) etacoord = .true.  ! alkyne case
        nm = 0
        do iTr = 1,numctr
          do k = 1,neigh%nbf(neigh%numnb,i,iTr)  ! how many metals ? and which
            kk = neigh%nbf(k,i,iTr)
            if (param%metal(at(kk)) .ne. 0) then
              nm = nm+1
              im = kk
            end if
          end do
        end do
        if (nm .eq. 0) then
          etacoord = .false.  ! etacoord makes no sense without metals!
        elseif (nm .eq. 1) then  ! distinguish M-CR2-R i.e. not an eta coord.
          ncm = 0
          do iTr = 1,numctr
            do k = 1,neigh%nbf(neigh%numnb,i,iTr)  !
              if (neigh%nbf(k,i,iTr) .ne. im) then ! all neighbors that are not the metal im
                kk = neigh%nbf(k,i,iTr)
                do l = 1,sum(neigh%nbf(neigh%numnb,kk,:))
                  if (neigh%nbf(l,kk,iTr) .eq. im) ncm = ncm+1 ! ncm=1 is alkyne, =2 is cp
                end do
              end if
            end do
          end do
          if (ncm .eq. 0) etacoord = .false.
        end if
      end if
      if (etacoord) then
        itag(i) = -1
        nbdum(:,i,:) = neigh%nbm(:,i,:)
      else
        nbdum(:,i,:) = neigh%nbf(:,i,:) ! take full set of neighbors by default
      end if
    end do

    do i = 1,natoms
      ati = at(i)
      hyb(i) = 0    ! don't know it
      numnbm = sum(neigh%nbm(neigh%numnb,i,:))
      numnbf = sum(neigh%nbf(neigh%numnb,i,:))
      numnb = sum(neigh%nb(neigh%numnb,i,:))
      nbdiff = numnbf-numnb
      nbmdiff = numnbf-numnbm
      nb20i = sum(nbdum(neigh%numnb,i,:))
      !>-- count H and O neighbors of i
      nh = 0
      no = 0
      do iTr = 1,numctr
        do j = 1,nb20i
          if (nbdum(j,i,iTr) .eq. 0) cycle
          if (at(nbdum(j,i,iTr)) .eq. 1) nh = nh+1
          if (at(nbdum(j,i,iTr)) .eq. 8) no = no+1
        end do
      end do
!>-- H
      if (param%group(ati) .eq. 1) then
        if (nb20i .eq. 2) hyb(i) = 1 ! bridging H
        if (nb20i .gt. 2) hyb(i) = 3 ! M+ tetra coord
        if (nb20i .gt. 4) hyb(i) = 0 ! M+ HC
      end if
!>-- Be
      if (param%group(ati) .eq. 2) then
        if (nb20i .eq. 2) hyb(i) = 1 ! bridging M
        if (nb20i .gt. 2) hyb(i) = 3 ! M+ tetra coord
        if (nb20i .gt. 4) hyb(i) = 0 !
      end if
!>-- B
      if (param%group(ati) .eq. 3) then
        if (nb20i .gt. 4) hyb(i) = 3
        if (nb20i .gt. 4.and.ati .gt. 10.and.nbdiff .eq. 0) hyb(i) = 5
        if (nb20i .eq. 4) hyb(i) = 3
        if (nb20i .eq. 3) hyb(i) = 2
        if (nb20i .eq. 2) hyb(i) = 1
      end if
!>-- C
      if (param%group(ati) .eq. 4) then
        if (nb20i .ge. 4) hyb(i) = 3
        if (nb20i .gt. 4.and.ati .gt. 10.and.nbdiff .eq. 0) hyb(i) = 5
        if (nb20i .eq. 3) hyb(i) = 2
        if (nb20i .eq. 2) then
          !>-- locate the two neighbors for the bangl call
          call neigh%nbLoc(natoms,nbdum,i,locarr)
          if (size(locarr,dim=2) .eq. 1) then
            idxdum = locarr(1,1)
            idxdum2 = locarr(2,1)
            iTr = locarr(neigh%numnb,1)
            iTr2 = locarr(neigh%numnb,1)
            deallocate (locarr)
          elseif (size(locarr,dim=2) .eq. 2) then
            idxdum = locarr(1,1)
            idxdum2 = locarr(1,2)
            iTr = locarr(neigh%numnb,1)
            iTr2 = locarr(neigh%numnb,2)
            deallocate (locarr)
          else
            if (mylevel >= 1) write (myunit,'("**ERROR**",a,1x,a)') ' Hybridization failed. Neighbors could not be located.',source
          end if
          call banglPBC(1,xyz,idxdum,i,idxdum2,iTr,iTr2,neigh%transVec,phi)
          if (phi*180./pi .lt. 150.0) then                         ! geometry dep. setup! GEODEP
            hyb(i) = 2  ! otherwise, carbenes will not be recognized
            itag(i) = 1  ! tag for Hueckel and HB routines
          else
            hyb(i) = 1  ! linear triple bond etc
          end if
          if (topo%qa(i) .lt. -0.4) then
            hyb(i) = 2
            itag(i) = 0  ! tag for Hueckel and HB routines
          end if
        end if
        if (nb20i .eq. 1) hyb(i) = 1  ! CO
      end if
!>-- N
      if (param%group(ati) .eq. 5) then
        if (nb20i .ge. 4) hyb(i) = 3
        if (nb20i .gt. 4.and.ati .gt. 10.and.nbdiff .eq. 0) hyb(i) = 5
        if (nb20i .eq. 3) hyb(i) = 3
        if (nb20i .eq. 3.and.ati .eq. 7) then
          kk = 0
          ll = 0
          nn = 0
          do iTr = 1,numctr
            do j = 1,3
              jj = nbdum(j,i,iTr)
              if (jj .eq. 0) exit ! if there is no 1st nb there is no nb at all
              if (at(jj) .eq. 8.and.sum(neigh%nb(neigh%numnb,jj,:)) .eq. 1) kk = kk+1 ! check for NO2 or R2-N=O
              if (at(jj) .eq. 5.and.sum(neigh%nb(neigh%numnb,jj,:)) .eq. 4) ll = ll+1 ! check for B-N, if the CN(B)=4 the N is loosely bound and sp2
              if (at(jj) .eq. 16.and.sum(neigh%nb(neigh%numnb,jj,:)) .eq. 4) nn = nn+1 ! check for N-SO2-
            end do
          end do
          if (nn .eq. 1.and.ll .eq. 0.and.kk .eq. 0) hyb(i) = 3
          if (ll .eq. 1.and.nn .eq. 0) hyb(i) = 2
          if (kk .ge. 1) then
            hyb(i) = 2
            itag(i) = 1  ! tag for Hueckel with no el. for the N in NO2
          end if
          if (nbmdiff .gt. 0.and.nn .eq. 0) hyb(i) = 2  ! pyridin N coord. to heavy atom
        end if
        if (nb20i .eq. 2) then
          hyb(i) = 2
          !>-- locate the two neighbors for the bangl call
          call neigh%nbLoc(natoms,nbdum,i,locarr)
          if (size(locarr,dim=2) .eq. 1) then
            idxdum = locarr(1,1)
            idxdum2 = locarr(2,1)
            iTr = locarr(neigh%numnb,1)
            iTr2 = locarr(neigh%numnb,1)
            deallocate (locarr)
          elseif (size(locarr,dim=2) .eq. 2) then
            idxdum = locarr(1,1)
            idxdum2 = locarr(1,2)
            iTr = locarr(neigh%numnb,1)
            iTr2 = locarr(neigh%numnb,2)
            deallocate (locarr)
          else
            if (mylevel >= 1) write (myunit,'("**ERROR**",a,1x,a)') ' Hybridization failed. Neighbors could not be located.',source
          end if
          call banglPBC(1,xyz,idxdum,i,idxdum2,iTr,iTr2,neigh%transVec,phi)
          jj = idxdum
          kk = idxdum2
          if (sum(nbdum(neigh%numnb,jj,:)) .eq. 1.and.at(jj) .eq. 6) hyb(i) = 1  ! R-N=C
          if (sum(nbdum(neigh%numnb,kk,:)) .eq. 1.and.at(kk) .eq. 6) hyb(i) = 1  ! R-N=C
          if (sum(nbdum(neigh%numnb,jj,:)) .eq. 1.and.at(jj) .eq. 7) hyb(i) = 1  ! R-N=N in e.g. diazomethane
          if (sum(nbdum(neigh%numnb,kk,:)) .eq. 1.and.at(kk) .eq. 7) hyb(i) = 1  ! R-N=N in e.g. diazomethane
          if (idxdum .gt. 0.and.param%metal(at(idxdum)) .gt. 0) hyb(i) = 1 ! M-NC-R in e.g. nitriles
          if (idxdum2 .gt. 0.and.param%metal(at(idxdum2)) .gt. 0) hyb(i) = 1 ! M-NC-R in e.g. nitriles
          if (at(jj) .eq. 7.and.at(kk) .eq. 7.and. &
&         sum(nbdum(neigh%numnb,jj,:)) .le. 2.and.sum(nbdum(neigh%numnb,kk,:)) .le. 2) hyb(i) = 1  ! N=N=N
          if (phi*180./pi .gt. lintr) hyb(i) = 1  ! geometry dep. setup! GEODEP
        end if
        if (nb20i .eq. 1) hyb(i) = 1
      end if
!>-- O
      if (param%group(ati) .eq. 6) then
        if (nb20i .ge. 3) hyb(i) = 3
        if (nb20i .gt. 3.and.ati .gt. 10.and.nbdiff .eq. 0) hyb(i) = 5
        if (nb20i .eq. 2) hyb(i) = 3
        if (nb20i .eq. 2.and.nbmdiff .gt. 0) then
          call nn_nearest_noM(i,natoms,at,xyz,neigh,rab,j,param) ! CN of closest non-M atom
          if (j .eq. 3) hyb(i) = 2 ! M-O-X konj
          if (j .eq. 4) hyb(i) = 3 ! M-O-X non
        end if
        if (nb20i .eq. 1) hyb(i) = 2
        if (nb20i .eq. 1.and.nbdiff .eq. 0) then
          call neigh%nbLoc(natoms,neigh%nb,i,locarr)
          iTr = locarr(neigh%numnb,1)
          deallocate (locarr)
          if (sum(neigh%nb(neigh%numnb,neigh%nb(1,i,iTr),:)) .eq. 1) hyb(i) = 1 ! CO
        end if
      end if
!>-- F
      if (param%group(ati) .eq. 7) then
        if (nb20i .eq. 2) hyb(i) = 1
        if (nb20i .gt. 2.and.ati .gt. 10) hyb(i) = 5
      end if
!>-- Ne
      if (param%group(ati) .eq. 8) then
        hyb(i) = 0
        if (nb20i .gt. 0.and.ati .gt. 2) hyb(i) = 5
      end if
!>-- done with main groups
      if (param%group(ati) .le. 0) then ! TMs
        nni = nb20i
        if (nh .ne. 0.and.nh .ne. nni) nni = nni-nh ! don't count Hs
        if (nni .le. 2) hyb(i) = 1
        if (nni .le. 2.and.param%group(ati) .le. -6) hyb(i) = 2
        if (nni .eq. 3) hyb(i) = 2
        if (nni .eq. 4.and.param%group(ati) .gt. -7) hyb(i) = 3  ! early TM, tetrahedral
        if (nni .eq. 4.and.param%group(ati) .le. -7) hyb(i) = 3  ! late TM, square planar
        if (nni .eq. 5.and.param%group(ati) .eq. -3) hyb(i) = 3  ! Sc-La are tetrahedral CN=5
      end if
    end do

    neigh%nb = nbdum ! list is complete but hyb determination is based only on reduced (without metals) list
    deallocate (nbdum)

    j = 0
    do i = 1,natoms
      numnb = sum(neigh%nb(neigh%numnb,i,:))
      if (numnb .gt. 12) j = j+1
      do iTr = 1,neigh%numctr
        do k = 1,neigh%nb(neigh%numnb,i,iTr)
          kk = neigh%nb(k,i,iTr)
          if (at(kk) .eq. 6.and.at(i) .eq. 6.and.itag(i) .eq. 1.and.itag(kk) .eq. 1) then ! check the very special situation of
            itag(i) = 0                                                           ! two carbene C bonded which is an arine
            itag(kk) = 0
          end if
        end do
      end do
    end do
    if (dble(j)/dble(natoms) .gt. 0.3.and.nb_call) then
      if (mylevel >= 1) write (myunit,'("**ERROR**",a,1x,a)') ' too many atoms with extreme high CN',source
    end if

  end subroutine gfnff_neigh

  subroutine getnb(n,at,rad,r,mchar,icase,f,f2,nbf,nb,param)
    !***********************************************************************
    !* Apply the bond-distance criterion to fill nb (CN in nb(20,:)).
    !* icase=1: full list. icase=2: excludes over-coordinated (non-)metal
    !* atoms. icase=3: excludes metals and unusually coordinated atoms.
    !***********************************************************************
    implicit none
    type(TGFFData),intent(in) :: param
    integer n,at(n),nbf(20,n),nb(20,n)
    real(wp) rad(n*(n+1)/2),r(n*(n+1)/2),mchar(n),f,f2

    integer :: i,j,k,nn,icase,hc_crit,nnfi,nnfj
    integer :: tag(n*(n+1)/2)
    real(wp) :: rco,fm

    nb = 0 ! resulting array (nbf is full from first call)
    tag = 0
    do i = 1,n
      nnfi = nbf(20,i)                  ! full CN of i, only valid for icase > 1
      do j = 1,i-1
        nnfj = nbf(20,j)               ! full CN of j
        fm = 1.0d0
        if (icase .eq. 1) then
          if (param%metal(at(i)) .eq. 2) fm = fm*f2 !change radius for metal atoms
          if (param%metal(at(j)) .eq. 2) fm = fm*f2
          if (param%metal(at(i)) .eq. 1) fm = fm*(f2+0.025)
          if (param%metal(at(j)) .eq. 1) fm = fm*(f2+0.025)
        end if
        if (icase .eq. 2) then
          hc_crit = 6
          if (param%group(at(i)) .le. 2) hc_crit = 4
          if (nnfi .gt. hc_crit) cycle
          hc_crit = 6
          if (param%group(at(j)) .le. 2) hc_crit = 4
          if (nnfj .gt. hc_crit) cycle
        end if
        if (icase .eq. 3) then
          if (mchar(i) .gt. 0.25.or.param%metal(at(i)) .gt. 0) cycle   ! metal case
          if (mchar(j) .gt. 0.25.or.param%metal(at(j)) .gt. 0) cycle   ! metal case
          if (nnfi .gt. param%normcn(at(i)).and.at(i) .gt. 10) cycle   ! HC case
          if (nnfj .gt. param%normcn(at(j)).and.at(j) .gt. 10) cycle   ! HC case
        end if
        k = lin(j,i)
        rco = rad(k) !(rad(i)+rad(j))/0.5291670d0
        if (r(k) .lt. fm*f*rco) tag(k) = 1  ! r: actual distance; fm*f*rco: threshold
      end do
    end do

    do i = 1,n
      nn = 0
      do j = 1,n
        if (tag(lin(j,i)) .eq. 1.and.i .ne. j) then
          nn = nn+1
          nb(nn,i) = j
        end if
      end do
      nb(20,i) = nn
    end do

  end subroutine getnb

  subroutine nn_nearest_noM(ii,n,at,xyz,neigh,r,nn,param)
    !***********************************************************************
    !* CN of the non-metal neighbor of atom ii that is closest in space.
    !***********************************************************************
    implicit none
    type(TGFFData),intent(in) :: param
    type(TNeigh),intent(in) :: neigh
    integer,intent(in) :: ii,n,at(n)
    real(wp),intent(in) :: xyz(3,n)
    integer,intent(inout) :: nn
    real(wp),intent(in) :: r(n*(n+1)/2)

    integer jmin,j,jj,numnb,iTr
    real(wp) :: dist
    real(wp) rmin

    numnb = neigh%numnb
    nn = 0
    rmin = 1.d+42
    jmin = 0
    dist = 0.0
    do iTr = 1,neigh%numctr
      do j = 1,neigh%nb(numnb,ii,iTr)
        jj = neigh%nb(j,ii,iTr)
        if (param%metal(at(jj)) .ne. 0) cycle
        dist = NORM2(xyz(:,ii)-(xyz(:,jj)+neigh%transVec(:,iTr)))
        if (dist .lt. rmin) then
          rmin = dist
          jmin = jj
        end if
      end do
    end do

    if (jmin .gt. 0) nn = sum(neigh%nb(numnb,jmin,:))

  end subroutine nn_nearest_noM

end module gfnff_topo_neighborlist
