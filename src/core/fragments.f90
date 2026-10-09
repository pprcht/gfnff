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

!> Splits a system into disconnected fragments by a recursive walk over the
!> neighbour list. The result feeds the EEQ fragment charge constraint.
module gfnff_fragments
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  public :: mrecgff,mrecgffPBC

contains  !> MODULE PROCEDURES START HERE

  subroutine mrecgff(nat,nb,molcount,molvec)
    !**********************************************
    !* Fragment search: molvec(i) is the fragment
    !* of atom i, molcount the number of fragments.
    !* nb(20,i) is the neighbour count of atom i.
    !**********************************************
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
    !**********************************************
    !* Recursive helper for mrecgff: assigns every
    !* atom reachable from i to fragment molcnt.
    !**********************************************
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

  subroutine mrecgffPBC(nat,numctr,numnb,nb,molcount,molvec)
    !**********************************************
    !* Fragment search over the periodic neighbour
    !* array nb(numnb,nat,numctr). Atoms bonded
    !* across cell boundaries share a fragment.
    !**********************************************
    implicit none
    integer,intent(in)    :: nat,numctr,numnb,nb(numnb,nat,numctr)
    integer,intent(inout) :: molvec(nat),molcount
    integer :: i,j,iTr
    real(wp),allocatable :: bond(:,:,:)
    logical,allocatable  :: taken(:)

    allocate (taken(nat),bond(nat,nat,numctr))
    bond = 0.0_wp
    do i = 1,nat
      do iTr = 1,numctr
        do j = 1,nb(numnb,i,iTr)
          bond(nb(j,i,iTr),i,iTr) = 1.0_wp
        end do
      end do
    end do

    if (int(sum(bond)) .ne. sum(nb(numnb,:,:))) then
      write (*,*)
      write (*,'(a,2i10)') ' Warning (mrecgffPBC): bond sum mismatch', &
        & int(sum(bond)),sum(nb(numnb,:,:))
      write (*,*)
    end if

    molvec = 0
    molcount = 1
    taken = .false.
    do i = 1,nat
      if (.not.taken(i)) then
        molvec(i) = molcount
        taken(i) = .true.
        call mrecgff2PBC(numctr,numnb,nat,nb,i,taken,bond,molvec,molcount)
        molcount = molcount+1
      end if
    end do
    molcount = molcount-1
  end subroutine mrecgffPBC

  recursive subroutine mrecgff2PBC(numctr,numnb,nat,nb,i,taken,bond,molvec,molcnt)
    !**********************************************
    !* Recursive helper for mrecgffPBC.
    !**********************************************
    implicit none
    integer,intent(in)    :: numctr,numnb,nat,nb(numnb,nat,numctr),i,molcnt
    integer :: j,icn,k,iTr,j_iTr(2)
    real(wp),intent(inout) :: bond(nat,nat,numctr)
    integer,intent(inout)  :: molvec(nat)
    logical,intent(inout)  :: taken(nat)

    icn = sum(nb(numnb,i,:))
    do k = 1,icn
      j_iTr = maxloc(bond(:,i,:))
      j = j_iTr(1)
      iTr = j_iTr(2)
      bond(j,i,iTr) = 0.0_wp
      if (i .eq. j.and.iTr .eq. 1) cycle
      if (.not.taken(j)) then
        molvec(j) = molcnt
        taken(j) = .true.
        call mrecgff2PBC(numctr,numnb,nat,nb,j,taken,bond,molvec,molcnt)
      end if
    end do
  end subroutine mrecgff2PBC
end module gfnff_fragments
