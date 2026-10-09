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

!> Topology-only EEQ charges (Goedecker model): computed from topology
!> distances rather than real ones, used only to scale the bond, angle and
!> torsion parameters. goedeckera solves for them, qheavy condenses H charges
!> onto heavy neighbors.
module gfnff_topo_eeq
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_data_types,only:TGFFTopology
  use gfnff_math_wrapper,only:sytrf_wrap,sytrs_wrap
  implicit none
  private

  public :: goedeckera,qheavy

contains  !> MODULE PROCEDURES START HERE

  subroutine goedeckera(n,at,pair,q,es,topo,printlevel,printunit,failed)
    !***********************************************************************
    !* EEQ charge solver (molecular). failed is set if the linear solve fails;
    !* q and es are then zero.
    !***********************************************************************
    implicit none
    character(len=*),parameter :: source = 'gfnff_topo_eeq_goedeckera'
    type(TGFFTopology),intent(in) :: topo
    integer,intent(in)  :: n          ! number of atoms
    integer,intent(in)  :: at(n)      ! ordinal numbers
    real(wp),intent(in)  :: pair(n*(n+1)/2)
    real(wp),intent(out) :: q(n)       ! output charges
    real(wp),intent(out) :: es         ! ES energy
    integer,intent(in),optional :: printlevel  !< verbosity (0=silent,1=errors,2=info,3=verbose)
    integer,intent(in),optional :: printunit   !< output unit (default: stdout)
    logical,intent(out),optional :: failed

    integer :: mylevel,myunit
    logical :: exitRun
    integer  :: m,i,j,ii
    integer  :: ij
    integer,allocatable :: ipiv(:)

    real(wp) :: gammij
    real(wp) :: r2
    real(wp) :: rij
    real(wp) :: tsqrt2pi
    real(wp) :: tmp
    real(wp),allocatable :: A(:,:)
    real(wp),allocatable :: x(:)

    integer :: io1,io2
    parameter(tsqrt2pi=0.797884560802866_wp)

    mylevel = 0
    if (present(printlevel)) mylevel = printlevel
    if (present(printunit)) then
      myunit = printunit
    else
      myunit = stdout
    end if

    m = n+topo%nfrag ! # atoms frag constrain
    allocate (A(m,m),x(m),ipiv(m))

    A = 0

    !>-- setup RHS
    do i = 1,n
      x(i) = topo%chieeq(i) ! EN of atom
      A(i,i) = topo%gameeq(i)+tsqrt2pi/sqrt(topo%alpeeq(i))
    end do

    !>-- setup A matrix
    do i = 1,n
      do j = 1,i-1
        ij = i*(i-1)/2+j
        rij = pair(ij)
        r2 = rij*rij
        gammij = 1.d0/sqrt(topo%alpeeq(i)+topo%alpeeq(j)) ! squared above
        tmp = erf(gammij*rij)/rij  ! apart from diagonal(=0), if ij non-bonded rij=1.0d+12
        A(j,i) = tmp
        A(i,j) = tmp
      end do
    end do

    !>-- fragment charge constrain
    do i = 1,topo%nfrag
      x(n+i) = topo%qfrag(i)
      do j = 1,n
        if (topo%fraglist(j) .eq. i) then
          A(n+i,j) = 1
          A(j,n+i) = 1
        end if
      end do
    end do

    call sytrf_wrap(a,ipiv,io1)
    call sytrs_wrap(a,x,ipiv,io2)

    exitRun = (io1 /= 0).or.(io2 /= 0)
    if (present(failed)) failed = exitRun
    if (exitRun) then
      if (mylevel >= 1) write (myunit,'("**ERROR**",a,1x,a)') 'Solving linear equations failed',source
      q = 0.0_wp
      es = 0.0_wp
      return
    end if

    q(1:n) = x(1:n)

    if (n .eq. 1) q(1) = topo%qfrag(1)

    !>-- energy
    es = 0.0_wp
    do i = 1,n
      ii = i*(i-1)/2
      do j = 1,i-1
        ij = ii+j
        rij = pair(ij)
        gammij = 1.d0/sqrt(topo%alpeeq(i)+topo%alpeeq(j)) ! squared above
        tmp = erf(gammij*rij)/rij
        es = es+q(i)*q(j)*tmp/rij
      end do
      es = es-q(i)*topo%chieeq(i) &
     &        +q(i)*q(i)*0.5d0*(topo%gameeq(i)+tsqrt2pi/sqrt(topo%alpeeq(i)))
    end do

  end subroutine goedeckera

  subroutine qheavy(n,at,numnb,numctr,nb,q)
    !***********************************
    !* Condenses hydrogen charges onto their heavy-neighbor(s), based on
    !* the topological neighbor list.
    !***********************************
    implicit none
    integer,intent(in)  :: numnb,numctr
    integer,intent(in)   ::  n,nb(numnb,n,numctr),at(n)
    real(wp),intent(inout) ::  q(n)

    integer i,j,k,iTr
    real(wp) qtmp(n)
    qtmp = q
    do i = 1,n
      if (at(i) .ne. 1) cycle
      qtmp(i) = 0
      do iTr = 1,numctr
        do j = 1,nb(numnb,i,iTr)
          k = nb(j,i,iTr)
          qtmp(k) = qtmp(k)+q(i)/dble(sum(nb(numnb,i,:)))  ! could be a bridging H
        end do
      end do
    end do

    q = qtmp

  end subroutine qheavy

end module gfnff_topo_eeq
