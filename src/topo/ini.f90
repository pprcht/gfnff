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
module gfnff_topo_ini
  use iso_fortran_env,only:wp => real64,sp => real32,stdout => output_unit,int8

  use gfnff_param,only:gfnff_thresholds,pse,gffVersion
  use gfnff_data_types,only:TGFFData,TGFFTopology,TGFFGenerator,TCell
  use gfnff_neighbor,only:TNeigh
  use gfnff_latticepoint,only:TLatticePoint,init_l

  use gfnff_topo_neighborlist,only:gfnff_neigh
  use gfnff_topo_eeq,only:goedeckera,qheavy
  use gfnff_topo_bondmat,only:nbondmat_pbc
  use gfnff_topo_iniphases,only:set_hb_xb_lists,perceive_rings, &
    &                           set_bonded_triples,set_pair_exponents, &
    &                           set_bonded_parameters,hueckel_solve
  use gfnff_topo_predicates,only:pilist,nofs,amide,amideH
  use gfnff_cn,only:gfnff_dlogcoord,getCoordinationNumber
  use gfnff_geometry,only:lin,omegaPBC
  use gfnff_fragments,only:mrecgffPBC
  !> (type,external) and not the bare form: bare leaves a call to an undeclared
  !> subroutine legal, which surfaces only at link time. With external the
  !> compiler names the file and line instead.
  implicit none(type,external)
  private

  public :: gfnff_ini

contains   !> MODULE PROCEDURES START HERE

  subroutine gfnff_ini(printlevel,makeneighbor,nat,at,xyz,ichrg, &
    &                  gen,param,topo,neigh,cell,efield, &
    &                  accuracy,version,io,printunit)
    !***********************************************************************
    !* Main GFN-FF topology and parameter initialization.
    !* Input:
    !*   printlevel  - verbosity (0=silent, 1=errors, 2=info, 3=verbose)
    !*   makeneighbor - rebuild neighbour list if .true.
    !*   nat/at/xyz  - system definition
    !*   ichrg       - total charge
    !*   gen/param/topo/neigh/cell - GFN-FF data structures
    !*   efield      - external electric field
    !*   accuracy    - threshold accuracy parameter
    !*   version     - force field version; selects which rule system
    !*                 assigns the bonded term parameters
    !* Output:
    !*   io          - error status (0 = success)
    !*   printunit   - output unit (optional, default: stdout)
    !***********************************************************************
    implicit none
    character(len=*),parameter :: source = 'gfnff_ini'
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: ichrg
    type(TNeigh),intent(inout) :: neigh ! main type for introducing PBC
    type(TGFFTopology),intent(inout) :: topo
    type(TGFFGenerator),intent(in) :: gen
    type(TGFFData),intent(in) :: param
    type(TCell),intent(in) :: cell
    real(wp),intent(in) :: efield(3)
    real(wp),intent(in) :: accuracy
    integer,intent(in) :: version

    integer,intent(in) :: printlevel    !< verbosity (0=silent,1=errors,2=info,3=verbose)
    logical,intent(in) :: makeneighbor  !< rebuild neighbour list if .true.
    integer,intent(out) :: io
    integer,intent(in),optional :: printunit  !< output unit (default: stdout)

    integer :: ati,atj,i,j,k,l,nn,ii,jj,kk,ll,m,ij,idum,ip,ji
    integer :: pis,nh
    integer :: picount,npiall
    integer :: nm
    integer :: niel(103)
    integer :: qloop_count,ifrag
    integer :: iTr,iTrj,iTrk,iTrl,iTrtmp
    real(wp) :: vTrl(3),vTrj(3),vTrk(3),vec(3),MaxCutOff

    real(wp) :: ff,phi
    real(wp) :: dum,dum1,dum2
    real(wp) :: ees
    real(wp) :: bohr

    parameter(bohr=1.0_wp/0.52917726_wp)

    real(wp),parameter :: rabd_cutoff = 13.0_wp

    logical :: picon
    logical :: piat,ex

    integer,allocatable  :: btyp(:),imetal(:),nbm(:,:),nbf(:,:)
    integer,allocatable  :: itag(:)
    integer,allocatable  :: piadr(:),piadr2(:),piadr3(:),piadr4(:)
    integer,allocatable  :: itmp(:),sring(:,:),cring(:,:,:)
    integer,allocatable  :: ipis(:),pimvec(:),nbpi(:,:,:)
    integer,allocatable  :: bdum(:,:,:),cdum(:,:,:)
    real(wp),allocatable :: rab(:)
    real(wp),allocatable :: sqrab(:)
    real(wp),allocatable :: cn(:)
    real(wp),allocatable :: dcn(:,:,:),dcndL(:,:,:)
    real(wp),allocatable :: dgam(:),dxi(:)
    real(wp),allocatable :: mchar(:)
    real(wp),allocatable :: rtmp(:)
    real(wp),allocatable :: pbo(:)
    real(wp),allocatable :: qtmp(:),dqa(:),qah(:)
    real(wp),allocatable :: pibo(:)
    real(sp),allocatable :: rabd(:,:)
    real(wp),allocatable :: transVec(:,:)
    type(TLatticePoint)  :: latPoint
    real(wp) :: lattice(3,3)
    integer :: boundaryCondition = 0
    integer,allocatable:: locarr(:,:)

    integer  :: ich,err,myunit
    real(wp) :: dispthr,cnthr,repthr,hbthr1,hbthr2
    logical :: exitRun,nb_call,adjLnAn,pr,pr2

    real(wp),parameter :: pi = 3.1415926535897932385_wp

    io = 0
    exitRun = .false.
    if (present(printunit)) then
      myunit = printunit
    else
      myunit = stdout
    end if
    pr = printlevel >= 2
    pr2 = printlevel >= 3

    call gfnff_thresholds(accuracy,dispthr,cnthr,repthr,hbthr1,hbthr2)

    lattice(:,:) = cell%lattice
    boundaryCondition = cell%npbc

    if (printlevel >= 2) then
      write (myunit,*)
      write (myunit,'(10x,"entering GFN-FF setup routine... ",i0)') nat
    end if

    if (printlevel >= 2) then
      write (myunit,*)
      write (myunit,'(10x,"==================== Thresholds ====================")')
      write (myunit,'(10x,"CN  :",f12.5)') cnthr
      write (myunit,'(10x,"rep :",f12.5)') repthr
      write (myunit,'(10x,"disp:",f12.5)') dispthr
      write (myunit,'(10x,"HB1 :",f12.5)') hbthr1
      write (myunit,'(10x,"HB2 :",f12.5)') hbthr2
      write (myunit,*)
    end if

    allocate (rab(nat*(nat+1)/2),source=0.0d0)
    allocate (cn(nat),source=0.0d0)
    allocate (sqrab(nat*(nat+1)/2),source=0.0d0)
    allocate (topo%hyb(nat),source=0)
    allocate (rtmp(nat*(nat+1)/2),source=0.0d0)
    allocate (pbo(nat*(nat+1)/2),source=0.0d0)
    allocate (piadr(nat),source=0)
    allocate (piadr2(nat),source=0)
    allocate (itmp(nat),source=0)
    allocate (itag(nat),source=0)
    allocate (sring(20,nat),source=0)
    allocate (cring(10,20,nat),source=0)
    allocate (piadr3(nat),source=0)
    allocate (piadr4(nat),source=0)
    allocate (qtmp(nat),source=0.0d0)
    allocate (dxi(nat),source=0.0d0)
    allocate (dgam(nat),source=0.0d0)
    allocate (topo%chieeq(nat),source=0.0d0)
    allocate (topo%gameeq(nat),source=0.0d0)
    allocate (topo%alpeeq(nat),source=0.0d0)
    allocate (topo%qa(nat),source=0.0d0)
    allocate (dqa(nat),source=0.0d0)
    allocate (qah(nat),source=0.0d0)
    allocate (nbm(20,nat),source=0)
    allocate (mchar(nat),source=0.0d0)
    allocate (imetal(nat),source=0)
    allocate (topo%zetac6(nat*(nat+1)/2),source=0.0d0)
    allocate (topo%xyze0(3,nat),source=0.0d0)
    allocate (nbf(20,nat),source=0)

    niel = 0
    do i = 1,nat
      niel(at(i)) = niel(at(i))+1
    end do

    if (printlevel >= 2) then
      write (myunit,'(10x,"Pauling EN used:")')
      do i = 1,103
        if (niel(i) .gt. 0) write (myunit,'(10x,"Z :",i2,"  EN :",f6.2)') i,param%en(i)
      end do
      dum = sqrt(sum(efield**2))
      write (myunit,'(10x,"electric field strengths (au):",f6.3)') dum
      write (myunit,*)
      write (myunit,'(10x," ------------------------------------------------- ")')
      write (myunit,'(10x,"|           Force Field Initialization            |")')
      write (myunit,'(10x," ------------------------------------------------- ")')
      write (myunit,*)
    else
      dum = sqrt(sum(efield**2))
    end if

    !>-- translation vectors within the maximum cutoff (at least the central 27 cells)
    call neigh%getTransVec(nat,at,xyz,cell,60.0_wp)  ! needed for neigh%init_n -> filliTrSum
    call neigh%init_n(nat,at,xyz,cell)

    !>-- bond pair matrix and non-bonded pair exponents
    allocate (neigh%bpair(nat,nat,neigh%numctr),source=0_int8)
    allocate (topo%alphanb(nat*(nat+1)/2),source=0.0d0)

    !>-- distances and bonds
    topo%xyze0 = xyz ! initial geom

    if (printlevel >= 2) write (myunit,'(10x,"distances ...")')
    pbo = 0
    rab = 0
    sqrab = 0
    !>-- redo the translation vectors with the HB cutoff sqrt(hbthr2)
    call neigh%getTransVec(nat,at,xyz,cell,sqrt(hbthr2))

    do i = 1,nat
      ati = at(i)
      kk = i*(i-1)/2
      do j = 1,i-1
        atj = at(j)
        k = kk+j
        rab(k) = NORM2(xyz(:,i)-xyz(:,j))
        sqrab(k) = rab(k)**2
        if (rab(k) .lt. 1.d-3) then
          if (printlevel >= 1) then
            write (myunit,*) i,j,ati,atj,rab(k)
            write (myunit,'("**ERROR** ",a,1x,a)') "Particular close distance present",source
          end if
          exitRun = .true.
          exit
        end if
      end do
    end do

    if (exitRun) then
      io = -1
      return
    end if

    allocate (dcn(3,nat,nat),source=0.0d0)
    if (boundaryCondition .eq. 0) then
      call gfnff_dlogcoord(nat,at,xyz,rab,cn,dcn,cnthr,param) ! dcn needed

    else
      vec = lattice(:,1)+lattice(:,2)
      MaxCutOff = sqrt(norm2(vec)**2+norm2(lattice(:,3))**2 &
        & -2*dot_product(vec,lattice(:,3)))+1.0_wp !
      MaxCutOff = max(MaxCutoff,60.0_wp) ! at least 60

      call init_l(latPoint,nat,at,xyz,cell%lattice,cell%npbc,MaxCutOff)
      call latPoint%getLatticepoints(transVec,MaxCutOff)
      latPoint%ntrans = size(transVec,dim=2)
      allocate (dcndL(3,3,nat),source=0.0d0)
      call getCoordinationNumber(nat,at,xyz,latPoint%nTrans,transVec,40.0_wp,5,cn,dcn,dcndL,param)
      deallocate (dcndL)
    end if
    do i = 1,nat
      dum2 = 0
      do j = 1,nat
        dum2 = dum2+sqrt(dcn(1,j,i)**2+dcn(2,j,i)**2+dcn(3,j,i)**2)
      end do
      !>-- estimated metallic character: ratio of av. dCN and CN times an EN cut-off
      !>   function, used in the neighbor routine and for the BS estimate
      mchar(i) = exp(-0.005d0*param%en(at(i))**8)*dum2/(cn(i)+1.0d0)
    end do
    deallocate (dcn)

    !>-- neighbor list, hyb and ring info
    !>-- Ln or An present?
    adjLnAn = .false.
    do i = 1,nat
      if ((at(i) .ge. 57.and.at(i) .le. 71).or.(at(i) .ge. 89.and.at(i) .le. 103)) then
        adjLnAn = .true.
        exit
      end if
    end do

    topo%qa = 0
    qloop_count = 0
    nb_call = .false. ! gfnff_neigh was (already) called

    !>-- charge loop: repeated only if rqshrink is significant (or Ln/An present)
    do while ((qloop_count .lt. 2.and.gen%rqshrink .gt. 1.d-3).or.adjLnAn)

      if (printlevel >= 2) then
        write (myunit,'(10x,"----------------------------------------")')
        write (myunit,'(10x,"generating topology and atomic info file ...")')
      end if
      call gfnff_neigh(makeneighbor,nat,at,xyz,cell,rab,gen%rqshrink, &
         & gen%rthr,gen%rthr2,gen%linthr,mchar,topo%hyb,itag,param,topo,neigh,nb_call, &
         & printlevel,myunit)
      nb_call = .true.

      !>-- H bound to Ln or An: needs the charges of the first pass
      if (adjLnAn.and.allocated(topo%qa).and.qloop_count .ne. 0) then
        call adjust_NB_LnH_AnH(param,nat,at,xyz,topo,neigh)
        adjLnAn = .false.
      end if

      do i = 1,nat
        imetal(i) = param%metal(at(i))
        !>-- Sn,Pb,Bi with small CN are better described as non-metals.
        !>   The number of neighbors can only decrease from first to second qloop.
        if (sum(neigh%nb(neigh%numnb,i,:)) .le. 4.and.param%group(at(i)) .gt. 3) imetal(i) = 0
      end do

      !>-- number of bonds, each pair counted once. For bonds to other cells only
      !>   translation vectors with "positive" sign (along the axes) are considered.
      allocate (bdum(nat,nat,neigh%numctr),source=0)
      allocate (cdum(neigh%numnb,nat,neigh%numctr),source=0)
      k = 0
      do i = 1,nat
        do iTr = 1,neigh%numctr
          do j = 1,neigh%nb(neigh%numnb,i,iTr)
            l = neigh%nb(j,i,iTr)
            if (bdum(l,i,iTr) .eq. 0) then
              bdum(l,i,iTr) = 1
              bdum(i,l,neigh%iTrNeg(iTr)) = 1
              cdum(j,i,iTr) = 1
              k = k+1
            end if
          end do
        end do
      end do
      neigh%nbond = k
      neigh%nbond_blist = neigh%nbond
      topo%nbond_blist = neigh%nbond
      allocate (btyp(neigh%nbond),source=0)
      allocate (pibo(neigh%nbond),source=0.0d0)
      allocate (neigh%blist(3,neigh%nbond),source=0) !first dim now 3 for saving iTr
      k = 0
      do i = 1,nat
        do iTr = 1,neigh%numctr
          do j = 1,neigh%nb(neigh%numnb,i,iTr)
            if (cdum(j,i,iTr) .eq. 1) then
              k = k+1
              neigh%blist(1,k) = neigh%nb(j,i,iTr)
              neigh%blist(2,k) = i
              neigh%blist(3,k) = iTr
            end if
          end do
        end do
      end do
      if (allocated(bdum)) deallocate (bdum)
      if (allocated(cdum)) deallocate (cdum)
      if (k .ne. neigh%nbond) then
        if (printlevel >= 1) write (myunit,'("**WARNING** ",a,1x,a)') "Setup of blist not as expected, check your results.",source
      end if

      !>-- Hueckel setup for first-row sp2 and sp atoms: list of all possible pi atoms
      k = 0  ! counts number of possible pi atoms
      piadr = 0
      piadr2 = 0
      do i = 1,nat
        piat = (topo%hyb(i) .eq. 1.or.topo%hyb(i) .eq. 2).and.pilist(at(i)) ! sp or sp2 and CNOFS
        kk = 0
        do iTr = 1,neigh%numctr
          do j = 1,neigh%nb(neigh%numnb,i,iTr)
            jj = neigh%nb(j,i,iTr)
            if (at(i) .eq. 8.and.at(jj) .eq. 16.and.topo%hyb(jj) .eq. 5) then
              piat = .false.
              cycle        ! SO3   is not a pi
            end if
            if (topo%hyb(jj) .eq. 1.or.topo%hyb(jj) .eq. 2) kk = kk+1         ! attached to sp2 or sp
          end do
        end do
        picon = kk .gt. 0.and.nofs(at(i))                     ! an N,O,F (sp3) on sp2
        if (at(i) .eq. 7.and.sum(neigh%nb(neigh%numnb,i,:)) .gt. 3) cycle           ! NR3-X is not a pi
        if (at(i) .eq. 16.and.topo%hyb(i) .eq. 5) cycle           ! SO3   is not a pi
        if (picon.or.piat) then
          k = k+1
          piadr(k) = i
          piadr2(i) = k
        end if
      end do
      npiall = k
      !>-- pi neighbor list
      allocate (nbpi(neigh%numnb,npiall,neigh%numctr),pimvec(npiall),source=0)
      nbpi = 0
      do i = 1,nat
        if (piadr2(i) .eq. 0) cycle
        ii = piadr2(i)
        nbpi(neigh%numnb,ii,:) = 0
        do iTr = 1,neigh%numctr
          do j = 1,neigh%nb(neigh%numnb,i,iTr)
            k = neigh%nb(j,i,iTr)
            if (piadr2(k) .gt. 0) then
              nbpi(neigh%numnb,ii,iTr) = nbpi(neigh%numnb,ii,iTr)+1
              nbpi(nbpi(neigh%numnb,ii,iTr),ii,iTr) = piadr2(k)
            end if
          end do
        end do
      end do

      !>-- assign pi atoms to fragments
      call mrecgffPBC(npiall,neigh%numctr,neigh%numnb,nbpi,picount,pimvec)
      deallocate (nbpi)

      !>-- xi correction for EEQ

      dxi = 0 ! default none
      do i = 1,nat
        ati = at(i)
        nn = sum(neigh%nb(neigh%numnb,i,:))
        if (nn .eq. 0) cycle
        ip = piadr2(i)
        call neigh%jth_nb(nat,xyz,ji,1,i,iTrtmp)  ! ji is the first nb of i in cell iTr
        nh = 0
        nm = 0
        do iTr = 1,neigh%numctr
          do j = 1,neigh%nb(neigh%numnb,i,iTr)
            if (at(neigh%nb(j,i,iTr)) .eq. 1) nh = nh+1
            if (imetal(neigh%nb(j,i,iTr)) .ne. 0) nm = nm+1
          end do
        end do
        if (ati .eq. 5) dxi(i) = dxi(i)+nh*0.015
        if (ati .eq. 6.and.nn .eq. 2.and.itag(i) .eq. 1) dxi(i) = -0.15 ! make carbene more negative
     if (ati .eq. 6.and.nn .eq. 1.and.at(ji) .eq. 8.and.neigh%nb(neigh%numnb,ji,iTrtmp) .eq. 1) then
          dxi(ji) = 0.15! free CO
        end if
if (ati .eq. 8.and.nn .eq. 1.and.ip .ne. 0.and.at(ji) .eq. 7.and.piadr2(ji) .ne. 0) dxi(i) = 0.05    ! nitro oxygen, otherwise NO2 HBs are too strong
        if (ati .eq. 8.and.nn .eq. 2.and.nh .eq. 2) dxi(i) = -0.02    ! H2O
        if (param%group(ati) .eq. 6.and.nn .gt. 2) dxi(i) = dxi(i)+nn*0.005! good effect
        if (ati .eq. 8.or.ati .eq. 16) dxi(i) = dxi(i)-nh*0.005
        if (param%group(ati) .eq. 7.and.ati .gt. 9.and.nn .gt. 1) then ! polyvalent Cl,Br ...
          if (nm .eq. 0) then
            dxi(i) = dxi(i)-nn*0.021! good effect
          else
            dxi(i) = dxi(i)+nn*0.05 ! good effect for TMs
          end if
        end if
      end do

      !>-- atomic EEQ xi, here for the non-geom. dep. charges qa with CN = nb
      do i = 1,nat
        ati = at(i)
        dum = min(dble(sum(neigh%nb(neigh%numnb,i,:))),gen%cnmax)  ! limits it
        topo%chieeq(i) = -param%chi(ati)+dxi(i)+param%cnf(ati)*sqrt(dum)
        topo%gameeq(i) = param%gam(ati)
        !>-- true TM charges are small, so make the metals less electronegative for
        !>   the non-geom. dep. charges: more q+ reflects the true polarity better,
        !>   which is used for guessing various potential terms. Big positive effect.
        if (imetal(i) .eq. 2) then
          topo%chieeq(i) = topo%chieeq(i)-gen%mchishift
        end if
        topo%alpeeq(i) = param%alp(ati)**2
      end do

      !>-- topology based charges

      if (printlevel >= 2) write (myunit,'(10x,"pair mat ...")')
      !>-- number of covalent bonds between atoms (up to 4 bonds), with PBC
      if (qloop_count .eq. 1) then
        call nbondmat_pbc(nat,neigh%numnb,neigh%numctr,neigh%nb,neigh,neigh%bpair)
      end if

      if (printlevel >= 2) write (myunit,'(10x,"computing topology distances matrix with Floyd-Warshall algo ...")')
      allocate (rabd(nat,nat),source=0.0e0_sp)
      rabd = rabd_cutoff
      !>-- topological distances by Floyd-Warshall, used in the EEQ for the
      !>   approximate topology charges qa
      do i = 1,nat
        rabd(i,i) = 0.0
        do iTr = 1,neigh%numctr
          do k = 1,neigh%nb(neigh%numnb,i,iTr)
            j = neigh%nb(k,i,iTr)
            rabd(j,i) = param%rad(at(i))+param%rad(at(j))
            rabd(i,j) = rabd(j,i)
          end do
        end do
      end do
      !>-- 1,x distances from the direct-neighbor ones above
      do k = 1,nat
        do i = 1,nat
          if (rabd(i,k) > gen%tdist_thr) cycle
          do j = 1,nat
            if (rabd(k,j) > gen%tdist_thr) cycle !tdist_thr = 12.0
            if (rabd(i,j) > (rabd(i,k)+rabd(k,j))) then
              rabd(i,j) = rabd(i,k)+rabd(k,j) ! get at least 1,4 distances
            end if
          end do
        end do
      end do

      do i = 1,nat
        do j = 1,i-1
          ij = lin(j,i)
          if (rabd(j,i) .gt. gen%tdist_thr) rabd(j,i) = rabd_cutoff ! values not properly considered
          rtmp(ij) = gen%rfgoed1*rabd(j,i)/0.52917726d0
        end do
      end do
      deallocate (rabd)

      if (printlevel >= 2) write (myunit,'(10x,"making topology EEQ charges ...")')
      if (topo%nfrag .le. 1) then                           ! nothing is known
        !>-- determine fragments
        call mrecgffPBC(nat,neigh%numctr,neigh%numnb,neigh%nbf,topo%nfrag,topo%fraglist)
        if (printlevel >= 2) write (myunit,'(10x,"#fragments for EEQ constrain: ",i0)') topo%nfrag
        !>-- reference charges from file, if topo%refcharges is set
        if (allocated(topo%refcharges)) then
          inquire (file=trim(topo%refcharges),exist=ex)
          if (.not.ex) then
            if (printlevel >= 1) write (myunit,'("**ERROR** ",a,1x,a)') &
              & 'reference charge file '//trim(topo%refcharges)//' not found',source
            io = -1
            return
          end if
        else
          ex = .false.
        end if
        if (ex) then
          if (printlevel >= 2) write (myunit,'(10x,a)') &
            & trim(topo%refcharges)//" file detected, attempting to read ..."
          open (newunit=ich,file=trim(topo%refcharges),action='read')
          qtmp = 0
          err = 0
          i = 0
          do while (err == 0)
            read (ich,*,iostat=err) dum
            if (err /= 0) exit
            if (i < nat) then
              i = i+1
              qtmp(topo%fraglist(i)) = qtmp(topo%fraglist(i))+dum
            else
              if (printlevel >= 1) write (myunit,'("**WARNING** ",a,1x,a)') &
                & "More charges than atoms present, assuming missmatch",source
              err = 1
            end if
          end do
          if (is_iostat_end(err).and.i == nat) err = 0
          close (ich)
          if (err == 0) then
            if (i < nat.or.abs(sum(qtmp)-ichrg) > 1.0e-3_wp) then
              if (printlevel >= 1) write (myunit,'("**WARNING** ",a,1x,a)') &
                & "Rejecting external charges input due to missmatch",source
            else
              topo%qfrag = dnint(qtmp)
         if (printlevel >= 2) write (myunit,'(10x,"fragment charges from <",a,"> :",10(1x,F7.3))') &
                     & trim(topo%refcharges),topo%qfrag(1:topo%nfrag)
            end if
          else
            if (printlevel >= 1) write (myunit,'("**ERROR** ",a,1x,a)') &
              & "Could not initialize fragment charges from file",source
            io = -1
            return
          end if
        end if
        !>-- host-supplied atomic reference charges (e.g. CEH) passed in memory: summed
        !>   per fragment to the integer net-charge constraint qfrag of the EEQ model
        if (allocated(topo%refq)) then
          if (size(topo%refq) == nat) then
            qtmp = 0.0d0
            do i = 1,nat
              qtmp(topo%fraglist(i)) = qtmp(topo%fraglist(i))+topo%refq(i)
            end do
            if (abs(sum(qtmp(1:topo%nfrag))-ichrg) > 1.0e-1_wp) then
              if (printlevel >= 1) write (myunit,'("**WARNING** ",a,1x,a)') &
                & "Reference charges do not sum to total charge, ignoring them",source
            else
              topo%qfrag(1:topo%nfrag) = dnint(qtmp(1:topo%nfrag))
              if (printlevel >= 2) write (myunit,'(10x,"fragment charges from reference charges :",10(1x,F7.3))') &
                & topo%qfrag(1:topo%nfrag)
            end if
          end if
        end if
        if (nat .lt. 100.and.topo%nfrag .gt. 2.and.ichrg .ne. 0.and.sum(topo%qfrag(2:topo%nfrag)) .gt. 999) then
          itmp = 0
          do i = 1,nat
            itmp(topo%fraglist(i)) = itmp(topo%fraglist(i))+1
          end do
          if (printlevel >= 1) then
            do i = 1,topo%nfrag
              write (myunit,*) i,itmp(i)
            end do
            write (myunit,'("**ERROR** ",a,1x,a)') 'fragment charge input required',source
          end if
          io = -1
          return
        end if
        if (nat .ge. 100.and.topo%nfrag .gt. 2.and.ichrg .ne. 0.and.sum(topo%qfrag(2:topo%nfrag)) .gt. 999) then
          topo%qfrag(1) = ichrg
          topo%qfrag(2:topo%nfrag) = 0
        end if
        if (topo%nfrag .eq. 2.and.ichrg .ne. 0.and.sum(topo%qfrag(2:topo%nfrag)) .gt. 999) then
          if (printlevel >= 2) write (myunit,*) 'trying auto detection of charge on 2 fragments:'
          topo%qfrag(1) = 0
          topo%qfrag(2) = ichrg
          call goedeckera(nat,at,rtmp,topo%qa,dum1,topo,printlevel,myunit,exitRun)
          if (exitRun) then
     if (printlevel >= 1) write (myunit,'("**ERROR** ",a,1x,a)') "Failed to generate charges",source
            io = -1
            return
          end if
          topo%qfrag(2) = 0
          topo%qfrag(1) = ichrg
          call goedeckera(nat,at,rtmp,topo%qa,dum2,topo,printlevel,myunit,exitRun)
          if (exitRun) then
     if (printlevel >= 1) write (myunit,'("**ERROR** ",a,1x,a)') "Failed to generate charges",source
            io = -1
            return
          end if
          if (dum1 .lt. dum2) then
            topo%qfrag(1) = 0
            topo%qfrag(2) = ichrg
          end if
          if (printlevel >= 2) then
            write (myunit,*) 'dEes      :',dum1-dum2
            write (myunit,*) 'charge 1/2:',topo%qfrag(1:2)
          end if
        end if
      end if

      !>-- topology-only EEQ charges from the rabd values, with the "right" fragment charge
      call goedeckera(nat,at,rtmp,topo%qa,ees,topo,printlevel,myunit,exitRun)
      if (exitRun) then
     if (printlevel >= 1) write (myunit,'("**ERROR** ",a,1x,a)') "Failed to generate charges",source
        io = -1
        return
      end if

      !>-- estimate how much of the fragment charge sits on the pi subsystems
      if (picount .gt. 0.and.qloop_count .gt. 0) then
        allocate (ipis(picount),source=0)
        qtmp = topo%qa ! save the "right" ones
        qah = topo%qa
        !>-- heavy atoms only, i.e. H condensed onto its neighbor
        call qheavy(nat,at,neigh%numnb,neigh%numctr,neigh%nb,qah)
        do pis = 1,picount
          do k = 1,npiall
            if (pimvec(k) .eq. pis) then
              kk = piadr(k)
              ifrag = topo%fraglist(kk) !the pi atom of this pi fragment is in EEQ fragment ifrag
              exit
            end if
          end do
          dum2 = topo%qfrag(ifrag) ! save
          topo%qfrag(ifrag) = 0 ! make only this EEQ fragment neutral
          call goedeckera(nat,at,rtmp,topo%qa,ees,topo,printlevel,myunit,exitRun) ! for neutral
          if (exitRun) then
   if (printlevel >= 1) write (myunit,'("**ERROR** ",a,1x,a)') "Failed to generate charges",source
            io = -1
            return
          end if
          topo%qfrag(ifrag) = dum2 ! back
          call qheavy(nat,at,neigh%numnb,neigh%numctr,neigh%nb,topo%qa)
          dqa = qah-topo%qa ! difference charges upon ionization
          dum1 = 0
          dum = 0
          do k = 1,npiall
            if (pimvec(k) .eq. pis) dum = dum+dqa(piadr(k)) ! only pi atoms
          end do
          dum = dum*1.1 !charges tend to be slightly too small 1.1-1.2
          ipis(pis) = idnint(dum)
          dum1 = dum1+dum
        end do
        topo%qa = qtmp ! put "right" charges used in FF construction and for HB/XB in place
      end if

      if (qloop_count .eq. 0) then
        do i = 1,nat
          itmp(i) = sum(neigh%nb(neigh%numnb,i,:))
        end do
      end if
      qloop_count = qloop_count+1
      if ((qloop_count .lt. 2.and.gen%rqshrink .gt. 1.d-3).or.adjLnAn) then  ! do the loop only if factor is significant
        deallocate (btyp,pibo,pimvec,neigh%blist)
      end if
    end do

    !>-- change EEQ J with the estimated q, a kind of third-order term

    do i = 1,nat
      ff = 0                           ! do nothing
      if (at(i) .eq. 1) ff = -0.08 ! H
      if (at(i) .eq. 5) ff = -0.05 ! B
      if (at(i) .eq. 6) then
        ff = -0.27 ! C
        if (topo%hyb(i) .lt. 3) ff = -0.45 ! unsat
        if (topo%hyb(i) .lt. 2) ff = -0.34 ! unsat
      end if
      if (at(i) .eq. 7) then
        ff = -0.13 ! N
        if (piadr(i) .ne. 0) ff = -0.14
        if (amide(nat,at,topo%hyb,neigh%numnb,neigh%numctr,neigh%nb,piadr,i)) ff = -0.16
      end if
      if (at(i) .eq. 8) then
        ff = -0.15 ! O
        if (topo%hyb(i) .lt. 3) ff = -0.08 ! unsat
      end if
      if (at(i) .eq. 9) ff = 0.10 ! F
      if (at(i) .gt. 10) ff = -0.02 ! heavy
      if (at(i) .eq. 17) ff = -0.02 ! Cl
      if (at(i) .eq. 35) ff = -0.11 ! Br
      if (at(i) .eq. 53) ff = -0.07 ! I
      if (imetal(i) .eq. 1) ff = -0.08 ! M maing
      if (imetal(i) .eq. 2) ff = -0.9  ! M TM    ??? too large
      if (param%group(at(i)) .eq. 8) ff = 0.0  ! RG
      dgam(i) = topo%qa(i)*ff
    end do

    !>-- true EEQ parameters; they are atomic, not element specific
    do i = 1,nat
      topo%chieeq(i) = -param%chi(at(i))+dxi(i)
      topo%gameeq(i) = param%gam(at(i))+dgam(i)
      if (amideH(nat,at,topo%hyb,neigh%numnb,neigh%numctr,neigh%nb,piadr2,i,neigh)) topo%chieeq(i) = topo%chieeq(i)-0.02
      ff = 0
      if (at(i) .eq. 6) ff = 0.09
      if (at(i) .eq. 7) ff = -0.21
      if (param%group(at(i)) .eq. 6) ff = -0.03
      if (param%group(at(i)) .eq. 7) ff = 0.50
      if (imetal(i) .eq. 1) ff = 0.3
      if (imetal(i) .eq. 2) ff = -0.1
      topo%alpeeq(i) = (param%alp(at(i))+ff*topo%qa(i))**2
    end do
    deallocate (dgam,dxi)

    !>-- graph-given harmonic2020 path: it reads only the bond list, the topological
    !>   charges and the non-bonded pair exponents, all built from the graph and the
    !>   elements alone. Everything from here to set_bonded_parameters is perception
    !>   that reads the geometry, which is meaningless in this mode and could invent
    !>   rings and pi systems out of an unstructured cloud of atoms.
    if (version == gffVersion%harmonic2020.and.allocated(neigh%user_bondmat)) then
      call set_pair_exponents(nat,at,param,gen,neigh,topo)
      call set_hb_xb_empty(nat,topo)
      deallocate (rtmp)
      if (printlevel >= 2) then
        write (myunit,*)
        write (myunit,*) 'GFN-FF setup done (graph-defined, perception skipped).'
        write (myunit,*)
      end if
      return
    end if

    !>-- ring perception (smallest ring size)
    call perceive_rings(nat,at,xyz,cell,printlevel,myunit,neigh,sring,cring)

    !>-- bonded triples not covered by the bend and torsion terms
    call set_bonded_triples(nat,printlevel,myunit,neigh,topo,idum)
    if (idum /= 0) then
      io = idum
      return
    end if

    !>-- non-bonded pair exponents and the D4 zeta scaling
    call set_pair_exponents(nat,at,param,gen,neigh,topo)

    !>-- hydrogen and halogen bond participant lists
    call set_hb_xb_lists(nat,at,xyz,cell,param,gen,itag,piadr2,hbthr2, &
       & printlevel,myunit,neigh,topo)

    !>-- iterative Hueckel over the pi subsystems
    call hueckel_solve(nat,at,xyz,param,gen,rab,itag,npiall,picount,pimvec, &
       & ipis,piadr3,piadr4,itmp,piadr,pibo,pbo,printlevel,myunit,neigh,topo,idum)
    if (idum /= 0) then
      io = idum
      return
    end if
    !>-- modify hyb due to pi assignment, and output

    do i = 1,nat
      if (topo%hyb(i) .eq. 2.and.piadr(i) .eq. 0.and.sum(neigh%nb(neigh%numnb,i,:)) .eq. 3 &
              &.and.param%group(at(i)) .eq. 4) then ! C,Si,Ge... CN=3, no pi
        call neigh%nbLoc(nat,neigh%nb,i,locarr)
        jj = 0
        kk = 0
        ll = 0
        if (size(locarr,dim=2) .eq. 1) then     ! all 3 in same cell
          jj = locarr(1,1)
          kk = locarr(2,1)
          ll = locarr(3,1)
          iTrj = locarr(neigh%numnb,1)
          iTrk = locarr(neigh%numnb,1)
          iTrl = locarr(neigh%numnb,1)
        elseif (size(locarr,dim=2) .eq. 2) then ! in 2 cells
          jj = locarr(1,1)
          kk = locarr(1,2)
          ll = locarr(2,1)             !  guess its here
          iTrj = locarr(neigh%numnb,1)
          iTrk = locarr(neigh%numnb,2)
          iTrl = locarr(neigh%numnb,1)
          if (ll .eq. 0) then  ! if ll is still zero then it is in the other cell
            ll = locarr(2,2)
            iTrl = locarr(neigh%numnb,2)
          end if
        else                                  ! in 3 cells
          jj = locarr(1,1)
          kk = locarr(1,2)
          ll = locarr(1,3)
          iTrj = locarr(neigh%numnb,1)
          iTrk = locarr(neigh%numnb,2)
          iTrl = locarr(neigh%numnb,3)
        end if
        deallocate (locarr)
        vTrl = neigh%transVec(:,iTrl)
        vTrj = neigh%transVec(:,iTrj)
        vTrk = neigh%transVec(:,iTrk)
        phi = omegaPBC(nat,xyz,i,jj,kk,ll,vTrl,vTrj,vTrk)  ! the shitty second geom. dep. term GEODEP
        if (abs(phi)*180./pi .gt. 40.d0) topo%hyb(i) = 3  ! change to sp^3
      end if
    end do

    if (printlevel >= 2) then
      write (myunit,*)
     write (myunit,'(2x,"atom   neighbors  erfCN metchar sp-hybrid imet pi  qest     coordinates")')
      do i = 1,nat
        j = topo%hyb(i)
        if (amide(nat,at,topo%hyb,neigh%numnb,neigh%numctr,neigh%nb,piadr,i)) j = -topo%hyb(i)
        if (at(i) .eq. 6.and.itag(i) .eq. 1) j = -topo%hyb(i)
        write (myunit,'(i5,2x,a2,3x,i4,3x,f5.2,2x,f5.2,8x,i2,3x,i2,3x,i2,2x,f6.3,3f12.6)') &
    &             i,pse(at(i)),sum(neigh%nb(neigh%numnb,i,:)),cn(i),mchar(i),j,imetal(i),piadr(i),topo%qa(i),xyz(1:3,i)
      end do
    end if

    !>-- fragments and their charges for output (check for CT)
    if (printlevel >= 2) then
      write (myunit,'(/,''molecular fragment  # atoms  topo charge'')')
      do i = 1,topo%nfrag
        dum = 0
        m = 0
        do k = 1,nat
          if (topo%fraglist(k) .eq. i) then
            m = m+1
            dum = dum+topo%qa(k)
          end if
        end do
        write (myunit,'(5x,i3,10x,i4,10x,f8.3)') i,m,dum
      end do
      write (myunit,*)
    end if

    !>-- bonded-term parameters: which rule system assigns them. Everything above is
    !>   perception; downstream code consumes the arrays filled here (neigh%vbond,
    !>   topo%vangl, topo%vtors, topo%alphanb) over lists that are already built.
    !>   Another force field is a sibling routine writing those arrays plus a case here.
    !>   An unregistered version is an error, not a silent fallback to the GFN-FF rules.
    select case (version)
    case (gffVersion%angewChem2020,gffVersion%angewChem2020_1, &
       &  gffVersion%angewChem2020_2,gffVersion%harmonic2020, &
       &  gffVersion%mcgfnff2023,gffVersion%conformer2020)
      call set_bonded_parameters(nat,at,xyz,cell,param,gen,cn,rab,rtmp,mchar, &
         & pbo,pibo,piadr,imetal,itag,btyp,sring,cring,hbthr1,hbthr2, &
         & printlevel,myunit,neigh,topo,idum)
    case default
      if (printlevel >= 1) write (myunit,'("**ERROR** ",a,i0,1x,a)') &
         & 'no bonded-parameter rule system registered for version ',version,source
      io = 1
      return
    end select
    if (idum /= 0) then
      io = idum
      return
    end if
    deallocate (rtmp)

    if (printlevel >= 2) then
      write (myunit,*)
      write (myunit,*) 'GFN-FF setup done.'
      write (myunit,*)
    end if

  end subroutine gfnff_ini

  subroutine set_hb_xb_empty(nat,topo)
    !***********************************************************************
    !* Put the hydrogen- and halogen-bond lists into a defined empty state,
    !* for the paths that skip set_hb_xb_lists.
    !* The arrays are normally allocated inside set_hb_xb_lists. Here the lists
    !* get length zero, so that gfnff_hbset0 has nothing to walk by intent
    !* and not by accident.
    !***********************************************************************
    implicit none
    integer,intent(in) :: nat
    type(TGFFTopology),intent(inout) :: topo

    topo%nathbH = 0
    topo%nathbAB = 0
    topo%natxbAB = 0
    if (.not.allocated(topo%hbbas)) allocate (topo%hbbas(nat),source=1.0d0)
    if (.not.allocated(topo%hbaci)) allocate (topo%hbaci(nat),source=1.0d0)
    if (.not.allocated(topo%hbatHl)) allocate (topo%hbatHl(2,0),source=0)
    if (.not.allocated(topo%hbatABl)) allocate (topo%hbatABl(2,0),source=0)
    if (.not.allocated(topo%xbatABl)) allocate (topo%xbatABl(5,0),source=0)

  end subroutine set_hb_xb_empty

  subroutine adjust_NB_LnH_AnH(param,nat,at,xyz,topo,neigh)
    !***********************************************************************
    !* Remove H from the neighbor lists of Ln and An atoms, and the metal from
    !* the list of that H, if the topology charge topo%qa(H) exceeds
    !* qthr = -0.0281. Modifies neigh%nb in place.
    !***********************************************************************
    type(TGFFData),intent(in) :: param
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TGFFTopology),intent(in) :: topo
    type(TNeigh),intent(inout) :: neigh ! main type for introducing PBC
    integer nb_tmp(neigh%numnb)
    integer :: i,iTr,iTrH,idx,inb,l,k,m,nnb,count_idx
    real(wp),parameter :: qthr = -0.0281_wp ! charge threshold

    do i = 1,nat
      !>-- only Ln and An
      if ((at(i) .ge. 57.and.at(i) .le. 71).or.(at(i) .ge. 89.and.at(i) .le. 103)) then
        do iTr = 1,neigh%numctr
          nnb = neigh%nb(neigh%numnb,i,iTr)
          idx = 0
          do count_idx = 1,nnb
            idx = idx+1
            inb = neigh%nb(idx,i,iTr) ! atom index of neighbor of i
            if (inb .eq. 0) cycle
            if (at(inb) .eq. 1) then
              if (topo%qa(inb) .gt. qthr) then
                !>-- compact the neighbor list of this Ln/An without the H
                nb_tmp = neigh%nb(:,i,iTr)
                nb_tmp(neigh%numnb) = neigh%nb(neigh%numnb,i,iTr)-1
                nb_tmp(idx) = 0
                neigh%nb(1:neigh%numnb-2,i,iTr) = 0
                neigh%nb(neigh%numnb,i,iTr) = nb_tmp(neigh%numnb) ! number neighbors
                l = 0
                do k = 1,neigh%numnb-2
                  if (nb_tmp(k) .ne. 0) then
                    l = l+1
                    neigh%nb(l,i,iTr) = nb_tmp(k)
                  end if
                end do
                !>-- same for the neighbor list of this H
                do iTrH = 1,neigh%numctr
                  do k = 1,neigh%nb(neigh%numnb,inb,iTrH)
                    if (neigh%nb(k,inb,iTrH) .eq. i) then
                      nb_tmp = neigh%nb(:,inb,iTrH)
                      nb_tmp(neigh%numnb) = neigh%nb(neigh%numnb,inb,iTrH)-1
                      nb_tmp(k) = 0
                      neigh%nb(1:neigh%numnb-2,inb,iTrH) = 0
                      neigh%nb(neigh%numnb,inb,iTrH) = nb_tmp(neigh%numnb) ! number neighbors
                      l = 0
                      do m = 1,neigh%numnb-2
                        if (nb_tmp(m) .ne. 0) then
                          l = l+1
                          neigh%nb(l,inb,iTrH) = nb_tmp(m)
                        end if
                      end do
                    end if
                  end do
                end do
                !>-- an element was removed from neigh%nb, so step idx back
                idx = idx-1
              end if
            end if
          end do
        end do
      end if
    end do

  end subroutine adjust_NB_LnH_AnH

end module gfnff_topo_ini
