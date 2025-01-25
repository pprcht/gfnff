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
!================================================================================!
program gfnff_main_tester
    use iso_fortran_env, only: wp=>real64,stdout=>output_unit
    use testmol
    use gfnff_interface
    implicit none
    
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    integer :: chrg
    integer :: uhf
    integer :: i,j,k,l

!========================================================================================!

   real(wp) :: energy
   real(wp),allocatable :: gradient(:,:)
   real(wp) :: gnorm

   logical :: fail,pr
   integer :: io
   type(gfnff_data) :: dat

!========================================================================================!
    fail = .false.
    pr = .true.

    nat = testnat
    allocate(at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz
    chrg = 0

    energy = 0.0_wp
    gnorm  = 0.0_wp
    allocate(gradient(3,nat),source=0.0_wp)

    write(*,*) nat
    write(*,*)
    do i=1,nat
      write(*,'(2x,i3,3x,3f16.6)') at(i),xyz(1:3,i)
    enddo
    call writetestcoord()

!=======================================================================================!
!=======================================================================================!
!> STANDARD USAGE
!=======================================================================================!
!=======================================================================================!

    write(*,*)
    write(*,*) '================================================================='
    write(*,*) '==================== GFN-FF SINGLEPOINT ========================='
    write(*,*) '================================================================='
  
!>--- First, setup of parametrisation and topology into DAT
    write(*,*)
    write(*,*) 'Initializing GFN-FF'
    call gfnff_initialize(nat,at,xyz,DAT,verbose=.false.,ichrg=chrg,iostat=io)
    write(*,*) 'Setup exit status:',io

!>--- Then call to the singlepoint routine
    write(*,*)
    write(*,*) 'Performing GFN-FF singlepoint calculation'
    call  gfnff_singlepoint(nat,at,xyz,DAT, energy,gradient, pr,iostat=io)
    write(*,*) 'Singlepoint exit status:',io

!>--- And printout of the results from DAT%res
    write(*,*) 
    write(*,*) 'GFN-FF results'
    call print_gfnff_results(stdout, DAT%res, allocated(DAT%solvation))

    call DAT%deallocate()

    write(*,*)
    write(*,*) '=========================== END ================================='
    write(*,*) '==================== GFN-FF SINGLEPOINT ========================='
    write(*,*) '=========================== END ================================='


!=======================================================================================!
!=======================================================================================!
!> USAGE WITH ALPB
!=======================================================================================!
!=======================================================================================!

    write(*,*)
    write(*,*) '================================================================='
    write(*,*) '========= GFN-FF + ALPB implicit solvation SINGLEPOINT =========='
    write(*,*) '================================================================='

!>--- First, setup of parametrisation and topology into DAT
    write(*,*)
    write(*,*) 'Initializing GFN-FF + ALPB(water)'
    DAT%solvent = 'h2o'
    call gfnff_initialize(nat,at,xyz,DAT,print=.true.,ichrg=chrg,iostat=io)
    write(*,*) 'Setup exit status:',io

!>--- Then call to the singlepoint routine
    write(*,*)
    write(*,*) 'Performing GFN-FF singlepoint calculation'
    call  gfnff_singlepoint(nat,at,xyz,DAT, energy,gradient, pr,iostat=io)
    write(*,*) 'Singlepoint exit status:',io

!>--- And printout of the results from DAT%res
    write(*,*)
    write(*,*) 'GFN-FF results'
    call print_gfnff_results(stdout, DAT%res, allocated(DAT%solvation))

    write(*,*)
    write(*,*) '========================== END =================================='
    write(*,*) '========= GFN-FF + ALPB implicit solvation SINGLEPOINT =========='
    write(*,*) '========================== END =================================='


!=======================================================================================!
   deallocate(gradient)
   deallocate(xyz,at)
!=======================================================================================!
end program gfnff_main_tester
