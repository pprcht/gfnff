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

!> Generator for lattice points
module gfnff_latticepoint
  use iso_fortran_env,only: wp => real64,stdout => output_unit
  use gfnff_helpers,only:crossproduct,bisectSearch,indexHeapSort
  implicit none
  private

  !> Tolerance for atomic distances
  real(wp),parameter :: tolSameDist = 1.0e-5_wp
  !> Tolerance for atomic square distances
  real(wp),parameter :: tolSameDist2 = tolSameDist**2

  public :: TLatticePoint,init_l

  !> Lattice point generator
  type :: TLatticePoint

    !> Boundary conditions for this generator
    integer :: boundaryCondition

    !> Ranges for generating translations
    integer :: ranges(2,3)

    !> Lattice vectors
    real(wp) :: lattice(3,3)

    !> Real space cutoff for the generation of lattice points
    real(wp) :: cutoff

    !> Exclude inversion symmetry in the lattice points
    logical :: excludeInversion

    !> Number of generated translations
    integer :: nTrans

    !> Lattice translations
    integer,allocatable :: trans(:,:)

    !> Distance of lattice points from origin
    real(wp),allocatable :: dist2(:)

  contains

    !> Update lattice points
    procedure :: update

    !> Returns all lattice points within a given cutoff
    procedure :: getLatticePoints

    !> Generate lattice points
    procedure :: generate

  end type TLatticePoint

  !> Initializes lattice point generator
  interface init_l
    module procedure :: initLatticePoint
    module procedure :: initLatticePointMolecule
  end interface init_l

contains

!> Initializes lattice point generator
  subroutine initLatticePointMolecule(self,nat,at,xyz,cutoff,excludeInversion)

    !> Instance of the lattice point generator
    type(TLatticePoint),intent(out) :: self

    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp) :: lattice(3,3) = 0.0_wp
    integer  :: boundaryCond = 0

    !> Real space cutoff for the generation of lattice points
    real(wp),intent(in) :: cutoff

    !> Exclude inversion symmetry in the lattice points
    logical,intent(in),optional :: excludeInversion

    call init_l(self,nat,at,xyz,lattice,boundaryCond,cutoff, &
       & excludeInversion)

  end subroutine initLatticePointMolecule

!> Initializes lattice point generator
  subroutine initLatticePoint(self,nat,at,xyz,lattice,boundaryCond,cutoff, &
        & excludeInversion)

    !> Instance of the lattice point generator
    type(TLatticePoint),intent(out) :: self

    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: lattice(3,3)
    integer,intent(in)  :: boundaryCond

    !> Real space cutoff for the generation of lattice points
    real(wp),intent(in) :: cutoff

    !> Exclude inversion symmetry in the lattice points
    logical,intent(in),optional :: excludeInversion

    self%boundaryCondition = boundaryCond
    self%lattice(:,:) = lattice
    self%cutoff = cutoff
    self%nTrans = 0

    if (present(excludeInversion)) then
      self%excludeInversion = excludeInversion
    else
      self%excludeInversion = .false.
    end if

    select case (boundaryCond)
    case (0)
      self%ranges(:,:) = 0

    case (3)
      call getRangesPBC3D(lattice,cutoff,self%ranges)

    end select

    call self%update(lattice)

  end subroutine initLatticePoint

!> Update lattice point generator
  subroutine update(self,lattice,updated)

    !> Source of the generated error
    character(len=*),parameter :: source = 'type_latticepoint_update'

    !> Instance of the lattice point generator
    class(TLatticePoint),intent(inout) :: self

    !> Lattice parameters
    real(wp),intent(in) :: lattice(:,:)

    !> New lattice points generated
    logical,intent(out),optional :: updated

    logical :: exitRun
    integer :: ranges(2,3),iostat

    self%lattice(:,:) = lattice

    select case (self%boundaryCondition)
    case (3)
      call getRangesPBC3D(lattice,self%cutoff,ranges)
      if (any(self%ranges /= ranges)) then
        self%ranges(:,:) = ranges
        self%nTrans = 0
      end if

    end select

    if (present(updated)) updated = self%nTrans == 0

    if (self%nTrans == 0) then
      call self%generate(iostat)
      if (iostat /= 0) then
        write (stdout,'("**ERROR** ",a,1x,a)') "Could not generate lattice points",source
        return
      end if
    end if

  end subroutine update

!> Returns all lattice points within a given cutoff
  subroutine getLatticePoints(self,latticePoint,cutoff)

    !> Source of the generated error
    character(len=*),parameter :: source = 'type_latticepoint_getLatticePoints'

    !> Instance of the lattice point generator
    class(TLatticePoint),intent(in) :: self

    !> Translation vectors for all lattice points
    real(wp),allocatable,intent(out) :: latticePoint(:,:)

    !> Cutoff for lattice points to be returned
    real(wp),intent(in),optional :: cutoff

    real(wp) :: cutoff2
    integer :: iTr,nTrans

    if (self%nTrans == 0) then
      return
    end if

    if (present(cutoff)) then
      cutoff2 = min(cutoff**2,self%cutoff**2)
    else
      cutoff2 = self%cutoff**2
    end if

    select case (self%boundaryCondition)
    case (0)
      allocate (latticePoint(3,1))
      latticePoint(:,:) = 0.0_wp

    case (3)
      call bisectSearch(nTrans,self%dist2(1:self%nTrans),cutoff2,tolSameDist2)
      allocate (latticePoint(3,nTrans))
      do iTr = 1,nTrans
        latticePoint(:,iTr) = self%lattice(:,1)*self%trans(1,iTr) &
              &              +self%lattice(:,2)*self%trans(2,iTr) &
              &              +self%lattice(:,3)*self%trans(3,iTr)
      end do

    end select

  end subroutine getLatticePoints

!> Generate lattice points
  subroutine generate(self,iostat)

    !> Source of the generated error
    character(len=*),parameter :: source = 'type_latticepoint_generate'

    !> Instance of the lattice point generator
    class(TLatticePoint),intent(inout) :: self

    integer,intent(out) :: iostat 

    iostat = 0

    select case (self%boundaryCondition)
    case (0)
      self%nTrans = 1
      if (allocated(self%trans)) then
        if (size(self%trans,dim=2) < self%nTrans) then
          deallocate (self%trans)
          allocate (self%trans(3,self%nTrans))
        end if
      else
        allocate (self%trans(3,self%nTrans))
      end if
      self%trans(:,1) = 0

    case (3)
      call generatePBC3D(self)

    case default
      write (stdout,'("**ERROR** ",a,1x,a)') "Boundary condition is not supported",source
      iostat = -1
    end select

  end subroutine generate

!> Generate lattice points under 3D infinite periodic boundary conditions
  subroutine generatePBC3D(self)

    !> Instance of the lattice point generator
    class(TLatticePoint),intent(inout) :: self

    integer :: mTrans,nTrans
    integer :: iTr,iTr1,iTr2,iTr3
    real(wp) :: cutoff2,r2,point(3)
    integer,allocatable :: indx(:)
    integer,allocatable :: trans(:,:)

    mTrans = product(self%ranges(2,:)-self%ranges(1,:)+1)
    if (allocated(self%trans)) then
      if (size(self%trans,dim=2) < mTrans) then
        deallocate (self%trans)
        allocate (self%trans(3,mTrans))
      end if
    else
      allocate (self%trans(3,mTrans))
    end if

    if (allocated(self%dist2)) then
      if (size(self%dist2) < mTrans) then
        deallocate (self%dist2)
        allocate (self%dist2(mTrans))
      end if
    else
      allocate (self%dist2(mTrans))
    end if

    cutoff2 = self%cutoff**2

    nTrans = 0
    do iTr1 = self%ranges(1,1),self%ranges(2,1)
      do iTr2 = self%ranges(1,2),self%ranges(2,2)
        do iTr3 = self%ranges(1,3),self%ranges(2,3)
          if (self%excludeInversion) then
            if (iTr1 < 0) cycle
            if (iTr2 < 0.and.iTr1 == 0) cycle
            if (iTr3 < 0.and.iTr2 == 0.and.iTr1 == 0) cycle
          end if
          point(:) = self%lattice(:,1)*iTr1 &
             &     +self%lattice(:,2)*iTr2 &
             &     +self%lattice(:,3)*iTr3
          r2 = sum(point**2)
          if (r2 <= cutoff2) then
            nTrans = nTrans+1
            self%dist2(nTrans) = r2
            self%trans(:,nTrans) = [iTr1,iTr2,iTr3]
          end if
        end do
      end do
    end do

    allocate (indx(nTrans))
    call indexHeapSort(indx(:nTrans),self%dist2(:nTrans),tolSameDist2)
    self%dist2(:nTrans) = self%dist2(indx(:nTrans))

    self%nTrans = nTrans
    call move_alloc(self%trans,trans)
    allocate (self%trans(3,nTrans))
    do iTr = 1,nTrans
      self%trans(:,iTr) = trans(:,indx(iTr))
    end do

  end subroutine generatePBC3D

!> Calculate the range of images of the central cell that interact
  subroutine getRangesPBC3D(lattice,cutoff,ranges)

    !> Lattice vectors
    real(wp),intent(in) :: lattice(:,:)

    !> Real space cutoff
    real(wp),intent(in) :: cutoff

    !> Array of the two extremal points
    integer,intent(out) :: ranges(2,3)

    integer :: tMax1,tMax2,tMax3
    real(wp) :: cos1,cos2,cos3
    real(wp) :: normal1(3),normal2(3),normal3(3)

    !> Get normals to lattice vectors
    normal1 = crossproduct(lattice(:,2),lattice(:,3))
    normal2 = crossproduct(lattice(:,3),lattice(:,1))
    normal3 = crossproduct(lattice(:,1),lattice(:,2))

    !> Normalize to unit length
    normal1 = normal1/norm2(normal1)
    normal2 = normal2/norm2(normal2)
    normal3 = normal3/norm2(normal3)

    !> Get angle between lattice vectors and normals
    cos1 = dot_product(normal1,lattice(:,1))
    cos2 = dot_product(normal2,lattice(:,2))
    cos3 = dot_product(normal3,lattice(:,3))

    !> Determine maximal translations
    tMax1 = ceiling(abs(cutoff/cos1))
    tMax2 = ceiling(abs(cutoff/cos2))
    tMax3 = ceiling(abs(cutoff/cos3))

    !> Save maximal translations as ranges
    ranges(:,1) = [-tMax1,tMax1]
    ranges(:,2) = [-tMax2,tMax2]
    ranges(:,3) = [-tMax3,tMax3]

  end subroutine getRangesPBC3D

end module gfnff_latticepoint
