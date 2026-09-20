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

!> Element and environment predicates used while the topology is perceived
!> (pi system? amide? halogen-bond donor?), that the parametrisation
!> branches on. Cheap and side-effect free; called from several setup sites.
module gfnff_topo_predicates
  use iso_fortran_env,only:stdout => output_unit
  use gfnff_neighbor,only:TNeigh
  implicit none
  private

  public :: pilist,nofs,xatom
  public :: ctype,alphaCO,amide,amideH

contains  !> MODULE PROCEDURES START HERE

  logical function pilist(ati)
    !***********************************
    !* True if element ati (B,C,N,O,F,S,Cl) can carry a pi system.
    !***********************************
    integer ati
    pilist = .false.
    if (ati .eq. 5.or.ati .eq. 6.or.ati .eq. 7.or.ati .eq. 8.or.ati .eq. 9.or.ati .eq. 16.or.ati .eq. 17) pilist = .true.
  end function pilist

  logical function nofs(ati)
    !***********************************
    !* True if element ati is N, O, F, S or Cl.
    !***********************************
    integer ati
    nofs = .false.
    if (ati .eq. 7.or.ati .eq. 8.or.ati .eq. 9.or.ati .eq. 16.or.ati .eq. 17) nofs = .true.
  end function nofs

  logical function xatom(ati)
    !***********************************
    !* True if ati can be the X in a halogen bond A-X...B.
    !***********************************
    integer ati
    xatom = .false.
    if (ati .eq. 17.or.ati .eq. 35.or.ati .eq. 53.or.&
   &   ati .eq. 16.or.ati .eq. 34.or.ati .eq. 52.or.&
   &   ati .eq. 15.or.ati .eq. 33.or.ati .eq. 51) xatom = .true.
  end function xatom

  integer function ctype(n,at,numnb,numctr,nb,pi,a)
    !***********************************
    !* Returns 1 if atom a is a carbonyl (C=O) carbon, 0 otherwise.
    !***********************************
    integer n,a,at(n),numnb,numctr,nb(numnb,n,numctr),pi(n)
    integer i,no,j,iTr

    ctype = 0 ! don't know

    no = 0
    do iTr = 1,numctr
      do i = 1,nb(numnb,a,iTr)
        j = nb(i,a,iTr)
        if (at(j) .eq. 8.and.pi(j) .ne. 0) no = no+1
      end do
    end do

    if (no .eq. 1.and.pi(a) .ne. 0) ctype = 1 ! a C=O carbon

  end function ctype

  logical function alphaCO(n,at,hyb,numnb,numctr,nb,pi,a,b)
    !***********************************
    !* True if the sp3-C/pi-C bond a-b is alpha to a carbonyl.
    !***********************************
    integer n,a,b,at(n),hyb(n),numnb,numctr,nb(numnb,n,numctr),pi(n)
    integer i,j,no,iTr

    alphaCO = .false.
    if (pi(a) .ne. 0.and.hyb(b) .eq. 3.and.at(a) .eq. 6.and.at(b) .eq. 6) then
      no = 0
      do iTr = 1,numctr
        do i = 1,nb(numnb,a,iTr)
          j = nb(i,a,iTr)
          if (at(j) .eq. 8.and.pi(j) .ne. 0.and.sum(nb(numnb,j,:)) .eq. 1) no = no+1 ! a pi =O on the C?
        end do
      end do
      if (no .eq. 1) then
        alphaCO = .true.
        return
      end if
    end if
    if (pi(b) .ne. 0.and.hyb(a) .eq. 3.and.at(b) .eq. 6.and.at(a) .eq. 6) then
      no = 0
      do iTr = 1,numctr
        do i = 1,nb(numnb,b,iTr)
          j = nb(i,b,iTr)
          if (at(j) .eq. 8.and.pi(j) .ne. 0.and.sum(nb(numnb,j,:)) .eq. 1) no = no+1 ! a pi =O on the C?
        end do
      end do
      if (no .eq. 1) then
        alphaCO = .true.
        return
      end if
    end if

  end function alphaCO

  logical function amide(n,at,hyb,numnb,numctr,nb,pi,a)
    !***********************************
    !* True if N atom a is an amide nitrogen.
    !***********************************
    integer n,a,at(n),hyb(n),numnb,numctr,nb(numnb,n,numctr),pi(n)
    integer i,j,no,nc,ic,iTr

    amide = .false. ! don't know
    if (pi(a) .eq. 0.or.hyb(a) .ne. 3.or.at(a) .ne. 7) return

    nc = 0
    no = 0
    do iTr = 1,numctr
      do i = 1,nb(numnb,a,iTr)
        j = nb(i,a,iTr)
        if (at(j) .eq. 6.and.pi(j) .ne. 0) then  ! a pi C on N?
          nc = nc+1
          ic = j
        end if
      end do
    end do

    if (nc .eq. 1) then
      do iTr = 1,numctr
        do i = 1,nb(numnb,ic,iTr)
          j = nb(i,ic,iTr)
          if (at(j) .eq. 8.and.pi(j) .ne. 0.and.nb(numnb,j,iTr) .eq. 1) no = no+1 ! a pi =O on the C?
        end do
      end do
    end if

    if (no .eq. 1) amide = .true.

  end function amide

  logical function amideH(n,at,hyb,numnb,numctr,nb,pi,a,neigh)
    !***********************************
    !* True if H atom a sits on an amide nitrogen with an sp3-C substituent.
    !***********************************
    type(TNeigh),intent(in) :: neigh ! for locating neighbor
    integer n,a,at(n),hyb(n),numnb,numctr,nb(numnb,n,numctr),pi(n)
    integer,allocatable :: locarr(:,:)
    integer i,j,nc,nn,iTr

    amideH = .false. ! don't know
    if (sum(nb(numnb,a,:)) .ne. 1) return
    call neigh%nbLoc(n,nb,a,locarr) ! locarr gives iTr of cell with the neighbor
    if (size(locarr,dim=2) .gt. 1) write (stdout,*) 'WARNING: Neighbors in more cells than expected! source: topo/predicates, amideH'
    nn = nb(1,a,locarr(numnb,1))       ! the N
    deallocate (locarr)
    if (.not.amide(n,at,hyb,numnb,numctr,nb,pi,nn)) return

    nc = 0
    do iTr = 1,numctr
      do i = 1,nb(numnb,nn,iTr)
        j = nb(i,nn,iTr)
        if (at(j) .eq. 6.and.hyb(j) .eq. 3) then  ! a sp3 C on N?
          nc = nc+1
        end if
      end do
    end do

    if (nc .eq. 1) amideH = .true.

  end function amideH

end module gfnff_topo_predicates
