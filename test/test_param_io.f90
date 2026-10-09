! This file is part of gfnff.
! SPDX-Identifier: LGPL-3.0-or-later
!> Parametrisation I/O and the force field version toggle.
module test_param_io
  use iso_fortran_env,only:wp => real64
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use gfnff_interface
  use gfnff_param_io,only:param_write_toml,param_read_toml,param_io_available
  use gfnff_param,only:gffVersion,gfnff_load_param,gfnff_set_param
  use gfnff_data_types,only:TGFFData,TGFFGenerator
  use coffeine,only:testnat,testat,testxyz
  implicit none
  private

  public :: collect_param_io

  !> Registered to run serially (see run_parallel in main.f90): toml-f must
  !> not be called from inside an OpenMP parallel region.

  real(wp),parameter :: thr = 1.0e-12_wp

contains

  subroutine collect_param_io(testsuite)
    !***********************************************************************
    !* Registers the parametrisation I/O tests.
    !***********************************************************************
    type(unittest_type),allocatable,intent(out) :: testsuite(:)
    testsuite = [ &
       & new_unittest("TOML round trip preserves every parameter",test_roundtrip), &
       & new_unittest("written set reproduces the energy",test_energy_identical), &
       & new_unittest("edited set changes the energy",test_edit_changes_energy), &
       & new_unittest("partial file overlays the internal set",test_partial_overlay), &
       & new_unittest("missing file is reported, not ignored",test_missing_file), &
       & new_unittest("unregistered version is rejected",test_unknown_version), &
       & new_unittest("mcGFN-FF is refused without a lattice",test_mc_needs_lattice), &
       & new_unittest("mcGFN-FF is accepted with a lattice",test_mc_with_lattice) &
       & ]
  end subroutine collect_param_io

  subroutine reference_param(param,gen)
    !***********************************************************************
    !* Fills a parametrisation the way the library does for a real system.
    !***********************************************************************
    type(TGFFData),intent(out) :: param
    type(TGFFGenerator),intent(out) :: gen
    logical :: ex
    call gfnff_load_param(gffVersion%angewChem2020_2,param,ex)
    call gfnff_set_param(testnat,gen,param)
  end subroutine reference_param

  subroutine test_roundtrip(error)
    !***********************************************************************
    !* Round-trips every parameter family through TOML into an identically
    !* seeded state; a match proves the file carries values, not reader defaults.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(TGFFData) :: p1,p2
    type(TGFFGenerator) :: g1,g2
    character(len=:),allocatable :: name,msg
    character(len=*),parameter :: fname = 'test_roundtrip.toml'
    integer :: io,ver,u


    if (.not.param_io_available()) return   !> built without toml-f

    call reference_param(p1,g1)
    call reference_param(p2,g2)

    !>-- perturb every family so a field that fails to round trip can't pass by luck
    p1%chi(1:86) = p1%chi(1:86)+0.125_wp
    p1%gam(1:86) = p1%gam(1:86)-0.0625_wp
    p1%bond(1:86) = p1%bond(1:86)*1.5_wp
    p1%tors2(1:86) = p1%tors2(1:86)+0.25_wp
    g1%linthr = 155.5_wp
    g1%fringbo = 0.0123_wp
    g1%bstren(3) = 1.75_wp
    g1%hdiag(7) = -2.5_wp
    g1%bsmat(1,2) = 0.375_wp
    g1%tdist_thr = 11.5

    ver = gffVersion%angewChem2020_2
    call param_write_toml(fname,p1,g1,ver,'roundtrip probe',io,msg)
    call check(error,io,0,"write failed: "//msg)
    if (allocated(error)) return

    call param_read_toml(fname,p2,g2,ver,name,io,msg)
    call check(error,io,0,"read failed: "//msg)
    if (allocated(error)) return

    call check(error,name,'roundtrip probe')
    if (allocated(error)) return

    call check(error,maxval(abs(p1%chi(1:86)-p2%chi(1:86))),0.0_wp,thr=thr)
    if (allocated(error)) return
    call check(error,maxval(abs(p1%gam(1:86)-p2%gam(1:86))),0.0_wp,thr=thr)
    if (allocated(error)) return
    call check(error,maxval(abs(p1%bond(1:86)-p2%bond(1:86))),0.0_wp,thr=thr)
    if (allocated(error)) return
    call check(error,maxval(abs(p1%tors2(1:86)-p2%tors2(1:86))),0.0_wp,thr=thr)
    if (allocated(error)) return

    call check(error,g1%linthr,g2%linthr,thr=thr)
    if (allocated(error)) return
    call check(error,g1%fringbo,g2%fringbo,thr=thr)
    if (allocated(error)) return
    call check(error,g1%bstren(3),g2%bstren(3),thr=thr)
    if (allocated(error)) return
    call check(error,g1%hdiag(7),g2%hdiag(7),thr=thr)
    if (allocated(error)) return
    call check(error,g1%bsmat(1,2),g2%bsmat(1,2),thr=thr)  !> off-diagonal catches a row/column mix-up
    if (allocated(error)) return
    call check(error,real(g1%tdist_thr,wp),real(g2%tdist_thr,wp),thr=1.0e-6_wp)
    if (allocated(error)) return

    open (newunit=u,file=fname); close (u,status='delete')
  end subroutine test_roundtrip

  subroutine test_energy_identical(error)
    !***********************************************************************
    !* Runs a second calculator from the parametrisation the first one wrote;
    !* matching energy means what toml-f emits is what was in use.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: c1,c2
    character(len=:),allocatable :: msg
    character(len=*),parameter :: fname = 'test_energy_identical.toml'
    real(wp) :: e1,e2,g1(3,testnat),g2(3,testnat)
    integer :: io,u


    if (.not.param_io_available()) return

    call gfnff_initialize(testnat,testat,testxyz,c1,iostat=io,printlevel=0)
    call check(error,io,0,"first init failed")
    if (allocated(error)) return
    call c1%singlepoint(testnat,testat,testxyz,e1,g1,printlevel=0)

    call gfnff_write_parametrisation(c1,fname,'in-use set',io,msg)
    call check(error,io,0,"write failed: "//msg)
    if (allocated(error)) return

    c2%parametrisation = fname
    call gfnff_initialize(testnat,testat,testxyz,c2,iostat=io,printlevel=0)
    call check(error,io,0,"init from written set failed")
    if (allocated(error)) return
    call c2%singlepoint(testnat,testat,testxyz,e2,g2,printlevel=0)

    !>-- Not bit-exact: toml-f writes floats below 1e3 with f24.16 (16 digits
    !>   after the point, not 16 significant digits), losing precision in the
    !>   round trip. The resulting shift is one ULP (relative 2e-16); an
    !>   actual misread parameter would move the energy far more, as shown below.
    call check(error,e1,e2,thr=1.0e-13_wp)
    if (allocated(error)) then
      call test_failed(error,"energy moved by more than round-off after a "// &
         & "parametrisation round trip")
      return
    end if
    call check(error,maxval(abs(g1-g2)),0.0_wp,thr=1.0e-12_wp)
    if (allocated(error)) return

    call c1%deallocate(); call c2%deallocate()
    open (newunit=u,file=fname); close (u,status='delete')
  end subroutine test_energy_identical

  subroutine test_edit_changes_energy(error)
    !***********************************************************************
    !* Counter-test to test_energy_identical: editing a parameter must move
    !* the energy, or the reader is silently ignoring the file.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: c1,c2
    type(TGFFData) :: param
    type(TGFFGenerator) :: gen
    character(len=:),allocatable :: msg
    character(len=*),parameter :: fname = 'test_edit.toml'
    real(wp) :: e1,e2,g1(3,testnat),g2(3,testnat)
    integer :: io,ver,u


    if (.not.param_io_available()) return

    call gfnff_initialize(testnat,testat,testxyz,c1,iostat=io,printlevel=0)
    call c1%singlepoint(testnat,testat,testxyz,e1,g1,printlevel=0)

    call reference_param(param,gen)
    param%bond(1:86) = param%bond(1:86)*1.1_wp   !> 10 % stiffer bonds
    ver = gffVersion%angewChem2020_2
    call param_write_toml(fname,param,gen,ver,'stiffer bonds',io,msg)
    call check(error,io,0,"write failed: "//msg)
    if (allocated(error)) return

    c2%parametrisation = fname
    call gfnff_initialize(testnat,testat,testxyz,c2,iostat=io,printlevel=0)
    call check(error,io,0,"init from edited set failed")
    if (allocated(error)) return
    call c2%singlepoint(testnat,testat,testxyz,e2,g2,printlevel=0)

    if (abs(e1-e2) < 1.0e-8_wp) then
      call test_failed(error,"edited parameters did not change the energy; "// &
         & "the file is being ignored")
      return
    end if

    call c1%deallocate(); call c2%deallocate()
    open (newunit=u,file=fname); close (u,status='delete')
  end subroutine test_edit_changes_energy

  subroutine test_partial_overlay(error)
    !***********************************************************************
    !* A file naming only a few keys must leave the rest at internal values.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(TGFFData) :: param
    type(TGFFGenerator) :: gen
    character(len=:),allocatable :: name,msg
    character(len=*),parameter :: fname = 'test_partial.toml'
    integer :: io,ver,u
    real(wp) :: chi_ref


    if (.not.param_io_available()) return

    call reference_param(param,gen)
    chi_ref = param%chi(6)

    open (newunit=u,file=fname,action='write',status='replace')
    write (u,'(a)') '[generator]'
    write (u,'(a)') 'fringbo = 0.5'
    close (u)

    ver = gffVersion%angewChem2020_2
    call param_read_toml(fname,param,gen,ver,name,io,msg)
    call check(error,io,0,"read failed: "//msg)
    if (allocated(error)) return

    call check(error,gen%fringbo,0.5_wp,thr=thr)
    if (allocated(error)) return
    call check(error,param%chi(6),chi_ref,thr=thr)  !> untouched by the file
    if (allocated(error)) return

    open (newunit=u,file=fname); close (u,status='delete')
  end subroutine test_partial_overlay

  subroutine test_missing_file(error)
    !***********************************************************************
    !* Rejects a nonexistent parameter file instead of silently proceeding.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(TGFFData) :: param
    type(TGFFGenerator) :: gen
    character(len=:),allocatable :: name,msg
    integer :: io,ver

    if (.not.param_io_available()) return

    call reference_param(param,gen)
    ver = gffVersion%angewChem2020_2
    call param_read_toml('definitely_not_here.toml',param,gen,ver,name,io,msg)
    if (io == 0) then
      call test_failed(error,"a missing parameter file was accepted")
      return
    end if
  end subroutine test_missing_file

  subroutine test_unknown_version(error)
    !***********************************************************************
    !* An unregistered force field version must fail at setup, not fall back.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: io

    calc%version = 987654        !> not in gffVersion
    call gfnff_initialize(testnat,testat,testxyz,calc,iostat=io,printlevel=0)
    if (io == 0) then
      call test_failed(error,"an unregistered force field version was accepted")
      return
    end if
    call calc%deallocate()
  end subroutine test_unknown_version

  subroutine test_mc_needs_lattice(error)
    !***********************************************************************
    !* mcGFN-FF needs a lattice: one of its four scaling factors is
    !* periodic-only, so a molecular request must be refused, not under-parametrized.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    integer :: io
    call gfnff_initialize(testnat,testat,testxyz,calc, &
       & version=gffVersion%mcgfnff2023,iostat=io,printlevel=0)
    if (io == 0) then
      call test_failed(error,"mcGFN-FF was accepted for a molecular system")
      return
    end if
    call calc%deallocate()
  end subroutine test_mc_needs_lattice

  subroutine test_mc_with_lattice(error)
    !***********************************************************************
    !* Counter-test: mcGFN-FF must be accepted once a lattice is supplied.
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(gfnff_data) :: calc
    real(wp) :: lattice(3,3),energy,grad(3,testnat)
    integer :: io
    lattice = 0.0_wp
    lattice(1,1) = 30.0_wp
    lattice(2,2) = 30.0_wp
    lattice(3,3) = 30.0_wp
    call gfnff_initialize(testnat,testat,testxyz,calc,lattice=lattice,npbc=3, &
       & version=gffVersion%mcgfnff2023,iostat=io,printlevel=0)
    call check(error,io,0,"mcGFN-FF was refused for a periodic system")
    if (allocated(error)) return
    call calc%singlepoint(testnat,testat,testxyz,energy,grad, &
       & lattice=lattice,printlevel=0)
    call check(error,energy < 0.0_wp)  !> a real number came out: parametrisation is in use
    if (allocated(error)) return
    call calc%deallocate()
  end subroutine test_mc_with_lattice

end module test_param_io
