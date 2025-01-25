!================================================================================!
! This file is part of gfnff.
!
! Copyright (C) 2025 Philipp Pracht
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
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use omp_lib
  use gfnff_interface
  use gfnff_type_timer
  use xyzreader
  implicit none

  integer :: nat
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:)
  integer :: i,j,k,l
  character(len=50) :: atmp
!========================================================================================!
  integer :: ichrg
  character(len=:),allocatable :: alpbsolvent
  real(wp) :: energy
  real(wp),allocatable :: gradient(:,:)
  real(wp) :: gnorm
  logical :: fail,pr
  integer :: io
  type(gfnff_data) :: ffdata
  type(gfnff_timer) :: timer

  character(len=1028) :: inputfile
  integer :: threads
  real(wp),parameter :: autoaa = 0.529177249_wp
!========================================================================================!
  interface
    subroutine ParseCommandLineArgs(threads,inputfile,ichrg,alpbsolvent)
      implicit none
      integer,intent(inout) :: threads
      integer,intent(inout) :: ichrg
      character(len=*),intent(inout) :: inputfile
      character(len=:),allocatable,intent(inout) :: alpbsolvent
    end subroutine ParseCommandLineArgs
  end interface
!========================================================================================!
  call print_header()

  !> DEFAULTS
  pr = .true.

  ichrg = 0 !> molecular charge
  energy = 0.0_wp
  gnorm = 0.0_wp

!=======================================================================================!
!> parse optional command line args and their defaults
  threads = 1
  call ParseCommandLineArgs(threads,inputfile,ichrg,alpbsolvent)

!> update setting based on input
#ifdef WITH_OpenMP
  call OMP_Set_Num_Threads(threads)
#ifdef WITH_MKL
  call MKL_Set_Num_Threads(threads)
#endif
  call ompprint_intern(atmp)
#endif

!> read xyz file
  if (.not.exists(inputfile)) then
    error stop 'Input file not found!'
  end if
  call readxyz(inputfile,nat,at,xyz)
  allocate (gradient(3,nat),source=0.0_wp)

  write (*,*) 'Input coords:'
  call printxyz(stdout,nat,at,xyz)
  write (*,*)

  if (allocated(alpbsolvent)) then
    ffdata%solvent = alpbsolvent
  end if

!=======================================================================================!
!=======================================================================================!
!> STANDARD USAGE
!=======================================================================================!
!=======================================================================================!

  call timer%new(2,.true.)
!>--- First, setup of parametrisation and topology into ffdata
  call timer%measure(1,'GFN-FF topology setup')
  call gfnff_initialize(nat,at,xyz,ffdata,print=.true.,ichrg=ichrg,iostat=io)
  write (*,*)
  if (io == 0) then
    write (*,*) 'Topology setup successful!'
  else
    write (*,*) 'Topology setup exited with errors!'
    error stop
  end if
  call timer%measure(1)
  call timer%write_timing(stdout,1,verbose=.false.)

!>--- Then call to the singlepoint routine
  write (*,*)
  write (stdout,'(1x,a)',advance='no') 'Performing GFN-FF singlepoint calculation ... '
  call timer%measure(2,'GFN-FF energy evaluation')
  call gfnff_singlepoint(nat,at,xyz,ffdata,energy,gradient,pr,iostat=io)
  if (io == 0) then
    write (*,*) 'success!'
  else
    write (*,*) 'FAILED!'
    error stop
  end if
  call timer%measure(2)
  call timer%write_timing(stdout,2,verbose=.false.)

!>--- And printout of the results from ffdata%res
  write (*,*)
  write (*,*) 'GFN-FF results'
  call ffdata%resultprint()

  write (*,*)
  write (*,'(A)') "Gradient (x, y, z in Eh/a0) for each atom:"
  do i = 1,nat
    write (*,'("  Atom ",I3," : ",3(ES12.5," "))') i,gradient(1,i),gradient(2,i),gradient(3,i)
  end do

  call timer%write(stdout,'GFN-FF timings ')

!=======================================================================================!
  deallocate (gradient)
  deallocate (xyz,at)
!=======================================================================================!
end program gfnff_main_tester

!=======================================================================================!
subroutine ompprint_intern(str)
  use omp_lib
  implicit none
  integer :: nproc,TID
  character(len=*) :: str
!$OMP PARALLEL PRIVATE(TID)
  TID = OMP_GET_THREAD_NUM()
  IF (TID .EQ. 0) THEN
    nproc = OMP_GET_NUM_THREADS()
    write (str,'(a,i0,a)') 'Total runtime (',nproc,' threads)'
  END IF
!$OMP END PARALLEL
end subroutine ompprint_intern

!=======================================================================================!
subroutine ParseCommandLineArgs(threads,inputfile,ichrg,alpbsolvent)
  use iso_fortran_env,only:wp => real64
  use xyzreader
  implicit none
  character(len=256) :: arg,arg2 ! Buffer to hold each argument
  integer :: numArgs,i,j     ! Variables to store argument count and loop index
  integer :: io,dumi
  real(wp) :: dum
  !> IN/OUTPUTS
  integer,intent(inout) :: threads
  integer,intent(inout) :: ichrg
  character(len=*),intent(inout) :: inputfile
  character(len=:),allocatable,intent(inout) :: alpbsolvent

  ! Get the number of command-line arguments
  numArgs = COMMAND_ARGUMENT_COUNT()

  ! Check if there are any arguments passed
  if (numArgs == 0) then
    return
  end if

  ! Loop through the command-line arguments
  do i = 1,numArgs
    ! Fetch each argument and store it in 'arg'
    call GET_COMMAND_ARGUMENT(i,arg)

    !> First argument can be input file
    if (i == 1) then
      if (exists(arg)) then
        inputfile = trim(arg)
      end if
    end if

    select case (trim(arg))
    case ('-T','--threads')
      !> parallelization
      call GET_COMMAND_ARGUMENT(i+1,arg2)
      read (arg2,*,iostat=io) dumi
      if (io == 0) threads = dumi

    case ('-i','--input')
      !> input file
      call GET_COMMAND_ARGUMENT(i+1,arg2)
      if (exists(arg2)) inputfile = trim(arg2)

    case ('-c','-chrg','--charge')
      !> Probe radius
      call GET_COMMAND_ARGUMENT(i+1,arg2)
      read (arg2,*,iostat=io) dum
      if (io == 0) ichrg = nint(dum)

    case ('-alpb','--alpb')
      !> ALPB implicit solvation
      call GET_COMMAND_ARGUMENT(i+1,arg2)
      alpbsolvent = trim(arg2)

    case ('-h','--help')
      !> display info about these flags and program usage
      call printhelp()

    end select
  end do

end subroutine ParseCommandLineArgs

subroutine printhelp()
  !> Print a brief description of the program
  write (*,"(1x,a)") "This program performs a GFN-FF single point energy calculation."
  write (*,"(1x,a)") "Note, the topology setup (an expensive pre-processing step) is"
  write (*,"(1x,a)") "automatically performed in each program execution."
  write (*,"(1x,a)") "Usage: gfnff <inputfile>.xyz [options]"
  write (*,*)
  write (*,"(1x,a)") "Options (all optional):"
  !> Describing each command-line argument
  write (*,"(1x,a)") " -T, --threads <int>      Number of threads for parallelization."
  write (*,"(1x,a)") " -i, --input <str>        Specify input file (xyz format)."
  write (*,"(1x,a)") "                          Alternatively, the first argument can be the input file."
  write (*,"(1x,a)") " -c, --charge <int>       Set the molecular charge (default 0)"
  write (*,"(1x,a)") " -alpb <str>              Specify a solvent for ALPB solvation"
  write (*,*)
  write (*,"(1x,a)") " -h, --help               Display this help message and exit."

  !> Examples of common usage
  write (*,*)
  write (*,"(1x,a)") "Examples:"
  write (*,"(1x,a)") " xhcff inputfile.xyz -p 2.0"
  write (*,"(1x,a)") " xhcff -i inputfile.xyz -T 4"
  write (*,"(1x,a)") " xhcff --pressure 1.0 --proberad 1.2"

  stop
end subroutine printhelp

subroutine print_header()
  implicit none
  write (*,*)
  write (*,"(10x,a)") repeat('+',42)
  write (*,"(10x,a)") "+"//repeat(' ',12)//"GFN-FF App v0.0.1"//repeat(' ',11)//'+'
  write (*,"(10x,a)") repeat('+',42)
  write (*,"(10x,a)") "Author: P.Pracht"
  write (*,"(10x,a)") "Based on the method by Spicher and Grimme"
  write (*,"(10x,a)") "(see https://doi.org/10.1002/anie.202004239)"
  write (*,*)

end subroutine print_header
