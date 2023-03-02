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
! along with gfnff.  If not, see <https://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------------!
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert
!> at https://github.com/grimme-lab/xtb
!================================================================================!

!> Abstract solvation model
module solvation_type_solvation
   use iso_fortran_env, only : wp => real64
   implicit none
   private

   public :: TSolvation


   type, abstract :: TSolvation
   contains

      !> Update coordinates and internal state
      procedure(update), deferred :: update

      !> Add potential shift
      procedure(addShift), deferred :: addShift

      !> Calculate solvation energy
      procedure(getEnergy), deferred :: getEnergy

      !> Calculate derivatives of solvation energy
      procedure(addGradient), deferred :: addGradient

   end type TSolvation


   abstract interface
      !> Update coordinates and internal state
      subroutine update(self, num, xyz)
         import :: TSolvation, wp

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Atomic numbers
         integer, intent(in) :: num(:)

         !> Cartesian coordinates
         real(wp), intent(in) :: xyz(:, :)

      end subroutine update

      !> Add potential shift
      subroutine addShift(self, qat, qsh, atomicShift, shellShift)
         import :: TSolvation, wp

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Atomic partial charges
         real(wp), intent(in) :: qat(:)

         !> Shell-resolved partial charges
         real(wp), intent(in) :: qsh(:)

         !> Atomic potential shift
         real(wp), intent(inout) :: atomicShift(:)

         !> Shell-resolved potential shift
         real(wp), intent(inout) :: shellShift(:)

      end subroutine addShift


      !> Calculate solvation energy
      subroutine getEnergy(self, qat, qsh, energy)
         import :: TSolvation, wp

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Atomic partial charges
         real(wp), intent(in) :: qat(:)

         !> Shell-resolved partial charges
         real(wp), intent(in) :: qsh(:)

         !> Total solvation energy
         real(wp), intent(out) :: energy

      end subroutine getEnergy


      !> Calculate derivatives of solvation energy
      subroutine addGradient(self, num, xyz, qat, qsh, gradient)
         import :: TSolvation, wp

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Atomic numbers
         integer, intent(in) :: num(:)

         !> Cartesian coordinates
         real(wp), intent(in) :: xyz(:, :)

         !> Atomic partial charges
         real(wp), intent(in) :: qat(:)

         !> Shell-resolved partial charges
         real(wp), intent(in) :: qsh(:)

         !> Molecular gradient
         real(wp), intent(inout) :: gradient(:, :)

      end subroutine addGradient

   end interface


end module solvation_type_solvation
