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

!> Wigner-Seitz image search per atom pair. The real space electrostatics
!> reaches only two cells and converges only if centred on the nearest image.
module gfnff_wsc
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  public :: wsc_images,wsc_maxcells

  !> largest number of cells the image search ever scans (3x3x3)
  integer,parameter :: wsc_maxcells = 27

  !> distance tolerance for treating two images as equidistant
  real(wp),parameter :: wsc_tol = 0.01_wp

contains  !> MODULE PROCEDURES START HERE

  pure subroutine wsc_images(nat,xyz,iat,jat,lattice,pbc,lattr,wc)
    !***********************************************************************
    !* Lattice translations lattr(:,1:wc) that bring atom jat closest to
    !* atom iat, all of them if several are equidistant within wsc_tol.
    !* lattice holds the direct lattice vectors as columns. wc is zero only
    !* for iat = jat in a cell without a periodic direction.
    !***********************************************************************
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in) :: iat,jat
    real(wp),intent(in) :: lattice(3,3)
    logical,intent(in) :: pbc(3)
    integer,intent(out) :: lattr(3,wsc_maxcells)
    integer,intent(out) :: wc

    integer :: rep(3),aa,bb,cc,c,img,minpos,nminpos
    integer :: cand(3,wsc_maxcells)
    real(wp) :: t(3),rw(3),dist(wsc_maxcells),mindist,nmindist
    logical :: avail(wsc_maxcells)

    where (pbc)
      rep = 1
    elsewhere
      rep = 0
    end where

    c = 0
    do aa = -rep(1),rep(1),1
      do bb = -rep(2),rep(2),1
        do cc = -rep(3),rep(3),1
          if ((aa .eq. 0.and.bb .eq. 0.and.cc .eq. 0).and.iat .eq. jat) cycle
          t = [aa,bb,cc]
          c = c+1
          cand(:,c) = [aa,bb,cc]
          rw = xyz(:,jat)+matmul(lattice,t)
          dist(c) = sqrt(sum((xyz(:,iat)-rw)**2))
        end do
      end do
    end do

    lattr(:,:) = 0
    wc = 0
    if (c .eq. 0) return

    !>-- nearest image first, then every further one within the tolerance
    avail(1:c) = .true.
    minpos = minloc(dist(1:c),dim=1)
    mindist = dist(minpos)
    avail(minpos) = .false.
    wc = 1
    lattr(:,1) = cand(:,minpos)
    !>-- bounded by c, since every candidate may fall inside the tolerance
    do img = 2,c
      nminpos = minloc(dist(1:c),dim=1,mask=avail(1:c))
      nmindist = dist(nminpos)
      if (abs(mindist-nmindist) .ge. wsc_tol) exit
      avail(nminpos) = .false.
      wc = wc+1
      lattr(:,wc) = cand(:,nminpos)
    end do

  end subroutine wsc_images

end module gfnff_wsc
