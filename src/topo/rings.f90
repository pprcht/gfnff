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

!> Ring perception and ring-membership queries used by the parametrization.
!>
!> getring36 finds the smallest ring (3 to 6 members, plus one larger ring
!> path) through a given atom by walking the neighbour list. The rings*
!> routines then report, for an atom, bond, angle or torsion, the size of
!> the smallest ring shared by all of its members; that size selects the
!> ring-dependent force constants.
module gfnff_topo_rings
  use iso_fortran_env,only:wp => real64
  use gfnff_neighbor,only:TNeigh
  use gfnff_geometry,only:banglPBC
  implicit none
  private

  public :: ringsatom,ringsbond,ringsbend,ringstors,ringstorl
  public :: chktors,chkrng
  public :: getring36,ssort

contains  !> MODULE PROCEDURES START HERE

  subroutine ringsatom(n,i,c,s,rings)
    !***********************************
    !* Smallest ring size containing atom i, from ring lists c/s.
    !***********************************
    implicit none
    integer n,i,k,c(10,20,n),s(20,n),rings

    rings = 99
    do k = 1,s(20,i)    ! all rings of atom i
      if (s(k,i) .lt. rings) rings = s(k,i)
    end do

  end subroutine ringsatom

  subroutine ringsbond(n,i,j,c,s,rings)
    !***********************************
    !* Smallest ring containing bond i-j; 0 if none.
    !***********************************
    implicit none
    integer n,i,j,k,l,c(10,20,n),s(20,n),rings,rings1,rings2

    rings1 = 99
    rings2 = 99
    do k = 1,s(20,i)    ! all rings of atom i
      do l = 1,s(k,i)  ! all atoms of ring k
        if (c(l,k,i) .eq. j.and.s(k,i) .lt. rings1) then
          rings1 = s(k,i)
        end if
      end do
    end do
    do k = 1,s(20,j)    ! all rings of atom j
      do l = 1,s(k,j)  ! all atoms of ring k
        if (c(l,k,j) .eq. i.and.s(k,j) .lt. rings2) then
          rings2 = s(k,j)
        end if
      end do
    end do
    continue
    rings = min(rings1,rings2)
    if (rings .eq. 99) rings = 0

  end subroutine ringsbond

  subroutine ringsbend(n,i,j,k,c,s,rings)
    !***********************************
    !* Smallest ring containing angle i-j-k; 0 if none.
    !***********************************
    implicit none
    integer n,i,j,k,rings
    integer c(10,20,n),s(20,n)
    integer itest,rings1,rings2,rings3,m,l

    if (s(20,i) .eq. 0.or.s(20,j) .eq. 0.or.s(20,k) .eq. 0) then
      rings = 0
      return
    end if

    rings1 = 99
    rings2 = 99
    rings3 = 99

    do m = 1,s(20,i)    ! all rings of atom i
      itest = 0
      do l = 1,s(m,i)  ! all atoms of ring m
        if (c(l,m,i) .eq. j.or.c(l,m,i) .eq. k) itest = itest+1
      end do
      if (itest .eq. 2.and.s(m,i) .lt. rings1) rings1 = s(m,i)
    end do
    do m = 1,s(20,j)    ! all rings of atom j
      itest = 0
      do l = 1,s(m,j)  ! all atoms of ring m
        if (c(l,m,j) .eq. i.or.c(l,m,j) .eq. k) itest = itest+1
      end do
      if (itest .eq. 2.and.s(m,j) .lt. rings2) rings2 = s(m,j)
    end do
    do m = 1,s(20,k)    ! all rings of atom k
      itest = 0
      do l = 1,s(m,k)  ! all atoms of ring m
        if (c(l,m,k) .eq. i.or.c(l,m,k) .eq. j) itest = itest+1
      end do
      if (itest .eq. 2.and.s(m,j) .lt. rings3) rings3 = s(m,k)
    end do

    rings = min(rings1,rings2,rings3)
    if (rings .eq. 99) rings = 0

  end subroutine ringsbend

  subroutine ringstors(n,i,j,k,l,c,s,rings)
    !***********************************
    !* Smallest ring containing torsion i-j-k-l; 0 if none.
    !***********************************
    implicit none
    integer n,i,j,k,l,rings
    integer c(10,20,n),s(20,n)
    integer itest,rings1,rings2,rings3,rings4,m,a

    if (s(20,i) .eq. 0.or.s(20,j) .eq. 0.or.s(20,k) .eq. 0.or.s(20,l) .eq. 0) then
      rings = 0
      return
    end if

    rings1 = 99
    rings2 = 99
    rings3 = 99
    rings4 = 99

    do m = 1,s(20,i)    ! all rings of atom i
      itest = 0
      do a = 1,s(m,i)  ! all atoms of ring m
        if (c(a,m,i) .eq. j.or.c(a,m,i) .eq. k.or.c(a,m,i) .eq. l) itest = itest+1
      end do
      if (itest .eq. 3.and.s(m,i) .lt. rings1) then
        rings1 = s(m,i)
      end if
    end do
    do m = 1,s(20,j)    ! all rings of atom j
      itest = 0
      do a = 1,s(m,j)  ! all atoms of ring m
        if (c(a,m,j) .eq. i.or.c(a,m,j) .eq. k.or.c(a,m,j) .eq. l) itest = itest+1
      end do
      if (itest .eq. 3.and.s(m,j) .lt. rings2) then
        rings2 = s(m,j)
      end if
    end do
    do m = 1,s(20,k)    ! all rings of atom k
      itest = 0
      do a = 1,s(m,k)  ! all atoms of ring m
        if (c(a,m,k) .eq. i.or.c(a,m,k) .eq. j.or.c(a,m,k) .eq. l) itest = itest+1
      end do
      if (itest .eq. 3.and.s(m,k) .lt. rings3) then
        rings3 = s(m,k)
      end if
    end do
    do m = 1,s(20,l)    ! all rings of atom k
      itest = 0
      do a = 1,s(m,l)  ! all atoms of ring m
        if (c(a,m,l) .eq. i.or.c(a,m,l) .eq. k.or.c(a,m,l) .eq. j) itest = itest+1
      end do
      if (itest .eq. 3.and.s(m,l) .lt. rings4) then
        rings4 = s(m,l)
      end if
    end do

    rings = min(rings1,rings2,rings3,rings4)
    if (rings .eq. 99) then
      rings = 0
    end if

  end subroutine ringstors

  subroutine ringstorl(n,i,j,k,l,c,s,ringl)
    !***********************************
    !* Largest ring containing torsion i-j-k-l; 0 if none.
    !***********************************
    implicit none
    integer n,i,j,k,l,ringl
    integer c(10,20,n),s(20,n)
    integer itest,rings1,rings2,rings3,rings4,m,a

    if (s(20,i) .eq. 0.or.s(20,j) .eq. 0.or.s(20,k) .eq. 0.or.s(20,l) .eq. 0) then
      ringl = 0
      return
    end if

    rings1 = -99
    rings2 = -99
    rings3 = -99
    rings4 = -99

    do m = 1,s(20,i)    ! all rings of atom i
      itest = 0
      do a = 1,s(m,i)  ! all atoms of ring m
        if (c(a,m,i) .eq. j.or.c(a,m,i) .eq. k.or.c(a,m,i) .eq. l) itest = itest+1
      end do
      if (itest .eq. 3.and.s(m,i) .gt. rings1) then
        rings1 = s(m,i)
      end if
    end do
    do m = 1,s(20,j)    ! all rings of atom j
      itest = 0
      do a = 1,s(m,j)  ! all atoms of ring m
        if (c(a,m,j) .eq. i.or.c(a,m,j) .eq. k.or.c(a,m,j) .eq. l) itest = itest+1
      end do
      if (itest .eq. 3.and.s(m,j) .gt. rings2) then
        rings2 = s(m,j)
      end if
    end do
    do m = 1,s(20,k)    ! all rings of atom k
      itest = 0
      do a = 1,s(m,k)  ! all atoms of ring m
        if (c(a,m,k) .eq. i.or.c(a,m,k) .eq. j.or.c(a,m,k) .eq. l) itest = itest+1
      end do
      if (itest .eq. 3.and.s(m,k) .gt. rings3) then
        rings3 = s(m,k)
      end if
    end do
    do m = 1,s(20,l)    ! all rings of atom k
      itest = 0
      do a = 1,s(m,l)  ! all atoms of ring m
        if (c(a,m,l) .eq. i.or.c(a,m,l) .eq. k.or.c(a,m,l) .eq. j) itest = itest+1
      end do
      if (itest .eq. 3.and.s(m,l) .gt. rings4) then
        rings4 = s(m,l)
      end if
    end do

    ringl = max(rings1,rings2,rings3,rings4)
    if (ringl .eq. -99) then
      ringl = 0
    end if

  end subroutine ringstorl

  logical function chktors(n,xyz,i,j,k,l,iTrj,iTrk,iTrl,neigh)
    !***********************************
    !* True if angle j-i-k or i-j-l is within 10 deg of linear (torsion ill-defined).
    !***********************************
    implicit none
    type(TNeigh),intent(in) :: neigh
    integer n,i,j,k,l,iTrj,iTrk,iTrl
    real(wp) xyz(3,n),phi

    chktors = .true.

    call banglPBC(1,xyz,j,i,k,iTrj,iTrk,neigh%transVec,phi)
    if (phi*180./3.1415926d0 .gt. 170.0d0) return
    call banglPBC(2,xyz,i,j,l,iTrj,iTrl,neigh%transVec,phi)
    if (phi*180./3.1415926d0 .gt. 170.0d0) return

    chktors = .false.

  end function chktors

  logical function chkrng(nn,n,c)
    !***********************************
    !* True if c(1:n) holds n mutually distinct indices in [1,nn].
    !***********************************
    implicit none
    integer n,idum(nn),nn,c(10),i,j
    chkrng = .true.
    idum = 0
    do i = 1,n
      idum(c(i)) = idum(c(i))+1
    end do
    j = 0
    do i = 1,nn
      if (idum(i) .eq. 1) j = j+1
    end do
    if (j .ne. n) chkrng = .false.
  end function chkrng

  subroutine getring36(n,at,numnb,numctr,nbin,a0_in,cout,irout)
    !***********************************
    !* Enumerate rings of size 3-6 through atom a0_in by exhaustive walks
    !* over the neighbour list nbin, then drop duplicates. cout(1:size,m)
    !* is the atom list of ring m, irout(m) its size, irout(20) the count.
    !***********************************
    implicit none
    integer,intent(out) :: cout(10,20),irout(20)
    integer,intent(in) :: n,at(n),numnb,numctr,nbin(numnb,n),a0_in
    integer :: i,nb(numnb,n),a0
    integer i1,i2,i3,i4,i5,i6
    integer n0,n1,n2,n3,n4,n5,n6
    integer a1,a2,a3,a4,a5,a6
    integer maxr
    parameter(maxr=500)
    integer list(n),m,mm,nn,c(10),cdum(10,maxr),iring
    integer adum1(0:n),adum2(0:n),kk,j,idum(maxr),same(maxr)
    real(wp) w(n)

    cout = 0
    irout = 0
    if (n .le. 2..or.nbin(numnb-1,a0_in) .eq. 1) return

    nn = nbin(numnb,a0_in)

    cdum = 0
    kk = 0
    do m = 1,nn
      nb = nbin
      !>-- rings are searched only within the unit cell; nb already carries
      !>   neighbours from adjacent cells for that purpose
      if (nb(m,a0_in) .eq. 1) cycle
      do i = 1,n
        if (nb(numnb,i) .eq. 1) nb(numnb,i) = 0
      end do

      do mm = 1,nn
        w(mm) = dble(mm)
        list(mm) = mm
      end do
      w(m) = 0.0d0
      call ssort(nn,w,list)
      do mm = 1,nn
        nb(mm,a0_in) = nbin(list(mm),a0_in)
      end do

      iring = 0
      c = 0

      a0 = a0_in
      n0 = nb(numnb,a0)

      do i1 = 1,n0
        a1 = nb(i1,a0)
        if (a1 .eq. a0) cycle
        n1 = nb(numnb,a1)
        do i2 = 1,n1
          a2 = nb(i2,a1)
          if (a2 .eq. a1) cycle
          n2 = nb(numnb,a2)
          do i3 = 1,n2
            a3 = nb(i3,a2)
            n3 = nb(numnb,a3)
            if (a3 .eq. a2) cycle
            c(1) = a1
            c(2) = a2
            c(3) = a3
            if (a3 .eq. a0.and.chkrng(n,3,c)) then
              iring = 3
              if (kk .eq. maxr) goto 99
              kk = kk+1
              cdum(1:iring,kk) = c(1:iring)
              idum(kk) = iring
            end if
            do i4 = 1,n3
              a4 = nb(i4,a3)
              n4 = nb(numnb,a4)
              if (a4 .eq. a3) cycle
              c(4) = a4
              if (a4 .eq. a0.and.chkrng(n,4,c)) then
                iring = 4
                if (kk .eq. maxr) goto 99
                kk = kk+1
                cdum(1:iring,kk) = c(1:iring)
                idum(kk) = iring
              end if
              do i5 = 1,n4
                a5 = nb(i5,a4)
                n5 = nb(numnb,a5)
                if (a5 .eq. a4) cycle
                c(5) = a5
                if (a5 .eq. a0.and.chkrng(n,5,c)) then
                  iring = 5
                  if (kk .eq. maxr) goto 99
                  kk = kk+1
                  cdum(1:iring,kk) = c(1:iring)
                  idum(kk) = iring
                end if
                do i6 = 1,n5
                  a6 = nb(i6,a5)
                  n6 = nb(numnb,a6)
                  if (a6 .eq. a5) cycle
                  c(6) = a6
                  if (a6 .eq. a0.and.chkrng(n,6,c)) then
                    iring = 6
                    if (kk .eq. maxr) goto 99
                    kk = kk+1
                    cdum(1:iring,kk) = c(1:iring)
                    idum(kk) = iring
                  end if
                end do
              end do
            end do
          end do
        end do
      end do

99    continue

    end do

    !>-- flag rings that duplicate an earlier one of the same size and atom set
    same = 0
    do i = 1,kk
      do j = i+1,kk
        if (idum(i) .ne. idum(j)) cycle ! different ring size
        if (same(j) .eq. 1) cycle ! already double
        adum1 = 0
        adum2 = 0
        do m = 1,10
          i1 = cdum(m,i)
          i2 = cdum(m,j)
          adum1(i1) = 1
          adum2(i2) = 1
        end do
        if (sum(abs(adum1-adum2)) .ne. 0) then
          same(j) = 0
        else
          same(j) = 1
        end if
      end do
    end do

    m = 0
    do i = 1,kk
      if (same(i) .eq. 0) then
        m = m+1
        irout(m) = idum(i)     ! number of atoms in ring m
        nn = idum(i)
        cout(1:nn,m) = cdum(1:nn,i)
        if (m .gt. 19) then
          m = 19
          goto 999
        end if
      end if
    end do
999 irout(20) = m  ! number of rings for this atom

    return
  end subroutine getring36

  subroutine ssort(n,edum,ind)
    !***********************************
    !* Sort edum(1:n) ascending by selection sort, permuting the
    !* companion index array ind alongside it.
    !***********************************
    implicit none

    integer,intent(in) :: n
    real(wp),intent(inout) :: edum(n)
    integer,intent(inout) :: ind(n)

    integer :: i,k
    real(wp) :: temp_val
    integer :: temp_ind

    do i = 1,n-1
      k = minloc(edum(i:n),dim=1)+i-1

      if (k /= i) then
        temp_val = edum(i)
        edum(i) = edum(k)
        edum(k) = temp_val

        temp_ind = ind(i)
        ind(i) = ind(k)
        ind(k) = temp_ind
      end if
    end do
  end subroutine ssort

end module gfnff_topo_rings
