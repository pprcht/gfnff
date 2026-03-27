! ──────────────────────────────────────────────────────────────────────────────
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
! ──────────────────────────────────────────────────────────────────────────────
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert, Sebastian Spicher, Stefan Grimme
!> at https://github.com/grimme-lab/xtb
! ──────────────────────────────────────────────────────────────────────────────

!> defines the Wigner--Seitz-cell (wsc) data type
!
!  The Wigner--Seitz cell is used to define periodic boundary conditions
!  by the cyclic cluster model (CCM). This type is usually bound to
!  the molecule class but can in principle be used independently.
module gfnff_type_wsc
  use iso_fortran_env,only:wp => real64
  implicit none

  public :: gfnff_wsc,generate_wsc

  private

!> definition of the the Wigner--Seitz-cell (wsc) data type
  type :: gfnff_wsc
    integer  :: n = 0       !< number of atoms in the WSC
    integer  :: cells = 0   !< number of cells used to generate WSC
    integer  :: rep(3) = 0  !< translations defining the number of cells
    real(wp) :: lattice(3,3) = 0.0_wp    !< lattice parameters
    integer,allocatable :: at(:,:)      !< define species
    integer,allocatable :: lattr(:,:,:,:) !< lattice translation
    real(wp),allocatable :: w(:,:)       !< define weights
    integer,allocatable :: itbl(:,:)    !< define index table
  contains
    procedure :: allocate => allocate_wsc
    procedure :: deallocate => deallocate_wsc
    procedure :: write => write_wsc
  end type gfnff_wsc
  
! ══════════════════════════════════════════════════════════════════════════════
contains  !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════

!> constructor for Wigner--Seitz cell
  subroutine allocate_wsc(self,n,rep)
    implicit none
    class(gfnff_wsc),intent(inout) :: self
    integer,intent(in) :: n      !< number of atoms
    integer,intent(in) :: rep(3) !< translations
    integer :: cells
    cells = product(2*rep+1)
    self%n = n
    self%rep = rep
    self%cells = cells
    call self%deallocate
    allocate (self%at(n,n),source=0)
    allocate (self%lattr(3,cells,n,n),source=0)
    allocate (self%w(n,n),source=0.0_wp)
    allocate (self%itbl(n,n),source=0)
  end subroutine allocate_wsc

!> @brief deconstructor for Wigner--Seitz cell
  subroutine deallocate_wsc(self)
    implicit none
    class(gfnff_wsc),intent(inout) :: self
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%lattr)) deallocate (self%lattr)
    if (allocated(self%w)) deallocate (self%w)
    if (allocated(self%itbl)) deallocate (self%itbl)
  end subroutine deallocate_wsc

  subroutine write_wsc(self,iunit,comment)
    implicit none
    class(gfnff_wsc),intent(in) :: self
    integer,intent(in) :: iunit
    character(len=*),intent(in) :: comment
    character(len=*),parameter :: dfmt = '(1x,a,1x,"=",1x,g0)'

    write (iunit,'(72(">"))')
    write (iunit,'(1x,"*",1x,a)') "Writing 'gfnff_wsc' class"
    write (iunit,'(  "->",1x,a)') comment
    write (iunit,'(72("-"))')
    write (iunit,'(1x,"*",1x,a)') "status of the fields"
    write (iunit,dfmt) "integer :: n           ",self%n
    write (iunit,dfmt) "integer :: cells       ",self%cells
    write (iunit,dfmt) "integer :: rep(1)      ",self%rep(1)
    write (iunit,dfmt) "        &  rep(2)      ",self%rep(2)
    write (iunit,dfmt) "        &  rep(3)      ",self%rep(3)
    write (iunit,'(72("-"))')
    write (iunit,'(1x,"*",1x,a)') "allocation status"
    write (iunit,dfmt) "allocated? at(:)       ",allocated(self%at)
    write (iunit,dfmt) "allocated? lattr(:,:,:,:)",allocated(self%lattr)
    write (iunit,dfmt) "allocated? w(:,:)      ",allocated(self%w)
    write (iunit,dfmt) "allocated? itbl(:,:)   ",allocated(self%itbl)
    write (iunit,'(72("-"))')
    write (iunit,'(1x,"*",1x,a)') "size of memory allocation"
    if (allocated(self%at)) then
      write (iunit,dfmt) "size(1) :: at(*,:)     ",size(self%at,1)
      write (iunit,dfmt) "size(2) :: at(:,*)     ",size(self%at,2)
    end if
    if (allocated(self%lattr)) then
      write (iunit,dfmt) "size(1) :: lattr(*,:,:,:)",size(self%lattr,1)
      write (iunit,dfmt) "size(2) :: lattr(:,*,:,:)",size(self%lattr,2)
      write (iunit,dfmt) "size(3) :: lattr(:,:,*,:)",size(self%lattr,3)
      write (iunit,dfmt) "size(4) :: lattr(:,:,:,*)",size(self%lattr,4)
    end if
    if (allocated(self%w)) then
      write (iunit,dfmt) "size(1) :: w(*,:)      ",size(self%w,1)
      write (iunit,dfmt) "size(2) :: w(:,*)      ",size(self%w,2)
    end if
    if (allocated(self%w)) then
      write (iunit,dfmt) "size(1) :: itbl(*,:)   ",size(self%itbl,1)
      write (iunit,dfmt) "size(2) :: itbl(:,*)   ",size(self%itbl,2)
    end if
    write (iunit,'(72("<"))')

  end subroutine write_wsc

! ══════════════════════════════════════════════════════════════════════════════

  subroutine generate_wsc(nat,at,xyz,lattice,pbc,wsc)
    implicit none
    !> molecular structure informtion
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz(3,nat),lattice(3,3)
    logical,intent(in) :: pbc(3)
    !> Wigner--Seitz cell data type (might be contained in mol)
    type(gfnff_wsc),intent(inout) :: wsc
! ------------------------------------------------------------------------
!  Variables
! ------------------------------------------------------------------------
    integer  :: rep(3)
    integer  :: ii,jj
    integer  :: aa,bb,cc
    integer  :: c,wc
    integer  :: minpos
    integer  :: nminpos
    real(wp) :: t(3),rw(3)
    real(wp) :: mindist
    real(wp) :: nmindist
    !> overall WSC tolerance to consider atoms as WSC-images
    real(wp),parameter :: tol = 0.01_wp
    integer,allocatable,dimension(:,:)   :: lattr
    real(wp),allocatable,dimension(:)     :: dist
    logical,allocatable,dimension(:)     :: trans

    where (pbc)
      rep = 1
    elsewhere
      rep = 0
    end where

! ------------------------------------------------------------------------
!  allocate space for the WSC first
    call wsc%allocate(nat,rep)

! ------------------------------------------------------------------------
! Create the Wigner-Seitz Cell (WSC)
! ------------------------------------------------------------------------
    wsc%at = 0
    wsc%itbl = 0

    !$omp parallel default(none) &
    !$omp private(ii,jj,aa,bb,cc,wc,c,dist,trans,t,lattr,rw) &
    !$omp shared(nat,at,xyz,lattice,wsc,rep) &
    !$omp private(mindist,minpos,nmindist,nminpos)

    allocate (lattr(3,wsc%cells),source=0)
    allocate (dist(wsc%cells),source=0.0_wp)
    allocate (trans(wsc%cells),source=.true.)

    ! Each WSC of one atom consists of n atoms
    !$omp do collapse(2) schedule(dynamic,32)
    do ii = 1,nat
      do jj = 1,nat
        !if (ii.eq.jj) cycle
        ! find according neighbours
        c = 0
        dist = 0.0_wp
        lattr = 0
        do aa = -rep(1),rep(1),1
          do bb = -rep(2),rep(2),1
            do cc = -rep(3),rep(3),1
              if ((aa .eq. 0.and.bb .eq. 0.and.cc .eq. 0).and.ii .eq. jj) cycle
              t = [aa,bb,cc]
              c = c+1
              lattr(:,c) = [aa,bb,cc]
              rw = xyz(:,jj)+matmul(lattice,t)
              dist(c) = sqrt(sum((xyz(:,ii)-rw)**2))
            end do
          end do
        end do
        ! sanity check; otherwise code below crashes sometimes
        if (c .eq. 0) cycle
        ! get first image with same dist
        ! find minimum in dist-array and assign it to minpos = minimum position
        trans = .true.
        minpos = minloc(dist(:c),dim=1)
        ! minimum distance saved in mindist
        mindist = dist(minpos)
        trans(minpos) = .false.
        wc = 1
        wsc%lattr(:,wc,jj,ii) = lattr(:,minpos)
        ! get other images with same distance
        if (c > 1) then
          find_images: do
            nminpos = minloc(dist(:c),dim=1,mask=trans(:c))
            nmindist = dist(nminpos)
            if (abs(mindist-nmindist) .lt. tol) then
              trans(nminpos) = .false.
              wc = wc+1
              wsc%lattr(:,wc,jj,ii) = lattr(:,nminpos)
            else
              wsc%w(jj,ii) = 1.0_wp/real(wc,wp)
              wsc%itbl(jj,ii) = wc
              wsc%at(jj,ii) = jj
              exit find_images
            end if
          end do find_images
        else
          wsc%w(jj,ii) = 1.0_wp
          wsc%itbl(jj,ii) = 1
          wsc%at(jj,ii) = jj
        end if
      end do
    end do
    !$omp end do
    deallocate (lattr,dist,trans)
    !$omp end parallel

  end subroutine generate_wsc
! ══════════════════════════════════════════════════════════════════════════════
end module gfnff_type_wsc
