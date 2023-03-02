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

!>--- First, setup of parametrisation and topology into DAT
    write(*,*)
    write(*,*) 'Initializing GFN-FF'
    call gfnff_initialize(nat,at,xyz,DAT,'no file',.false.,.false.,ichrg=chrg,iostat=io)
    write(*,*) 'Setup exit status:',io

!>--- Then call to the singlepoint routine
    write(*,*)
    write(*,*) 'Performing GFN-FF singlepoint calculation'
    call  gfnff_singlepoint(nat,at,xyz,DAT, energy,gradient, pr,iostat=io)
    write(*,*) 'Singlepoint exit status:',io

!>--- And printout of the results from DAT%res
    write(*,*) 
    write(*,*) 'GFN-FF results'
    call print_gfnff_results(stdout, DAT%res, pr, allocated(DAT%solvation))

!=======================================================================================!
   deallocate(gradient)
   deallocate(xyz,at)
!=======================================================================================!
end program gfnff_main_tester
