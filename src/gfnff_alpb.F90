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
! along with gfnff. If not, see <https://www.gnu.org/licenses/>.
!==============================================================================!

!> module gfnff_gbsa
!> An interface to GBSA and ALPB calculations for gfnff
!>
!> NOTE: The GBSA for GFN-FF technically does not exist.
!>       We are using the improved ALPB model instead.
module gfnff_gbsa
  use iso_fortran_env,only:wp => real64,stdout => output_unit
#ifdef WITH_GBSA
  use solvation_solv_input, only : TSolvInput
  use solvation_solv_state, only : solutionState
  use solvation_solv_kernel, only : gbKernel
  use solvation_solv_model, only: TSolvModel, init, newBornModel, info
  use solvation_solv_gbsa, only: TBorn
#endif
  implicit none
  private

#ifndef WITH_GBSA
  !> TBorn placeholder if compiled without GBSA support
  type :: TBorn
    integer :: dummy
  end type TBorn
#endif

  public :: TBorn, gfnff_gbsa_init, gfnff_solvation, gfnff_gbsa_print

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

subroutine gfnff_gbsa_init(nat,at,solv,gbsa)
!> set up the GBSA (i.e., ALPB) parametrization
   implicit none
   !> INPUT
   integer,intent(in)          :: nat      !> number of atoms
   integer,intent(in)          :: at(nat)  !> atom types
   character(len=*),intent(in) :: solv     !> which solvent to use
   !> OUTPUT
   type(TBorn),intent(out)     :: gbsa     !> gbsa type returned 
   !> LOCAL
#ifdef WITH_GBSA
   type(TSolvInput) :: input
   type(TSolvModel) :: model

   input = TSolvInput(solvent=trim(solv), alpb=.true., kernel=gbKernel%p16)
   call init(model, input, 0)
   !call info(model,stdout)
   call newBornModel(model, gbsa, at)
#else
   gbsa%dummy = 0
#endif

   return
end subroutine gfnff_gbsa_init

subroutine gfnff_gbsa_print(gbsa,iunit)
   implicit none
   integer :: iunit
   type(TBorn),intent(in)     :: gbsa     !> gbsa type
#ifdef WITH_GBSA
   call gbsa%info(iunit)
#endif
end subroutine gfnff_gbsa_print

!========================================================================================!

subroutine gfnff_solvation(nat,at,xyz,qat,gbsa,gsolv,gradient)
!> get energy and gradient contribution from SASA term
   implicit none
   !> INPUT
   integer,intent(in)          :: nat        !> number of atoms
   integer,intent(in)          :: at(nat)    !> atom types 
   real(wp),intent(in)         :: xyz(3,nat) !> coordinats (bohr)
   real(wp),intent(in)         :: qat(nat)   !> charges
   type(TBorn),intent(inout)   :: gbsa       !> gbsa type returned
   !> OUTPUT
   real(wp),intent(out)        :: gsolv      !> solv contribution to E
   real(wp),intent(inout)      :: gradient(3,nat) !> atomic gradient (added to)
   !> LOCAL
   real(wp) :: gborn,ghb,gsasa,gshift
#ifdef WITH_GBSA
   call gbsa%update(at, xyz)
   call gbsa%addGradient(at,xyz,qat,qat,gradient)
   call gbsa%getEnergyParts(qat,qat,gborn,ghb,gsasa, &
   & gshift)
   gsolv = gsasa+gborn+ghb+gshift
#else
   gsolv = 0.0_wp
#endif
end subroutine gfnff_solvation
!========================================================================================!
end module gfnff_gbsa
