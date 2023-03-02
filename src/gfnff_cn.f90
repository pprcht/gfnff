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
!--------------------------------------------------------------------------------!
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert, Sebastian Spicher, Stefan Grimme
!> at https://github.com/grimme-lab/xtb
!================================================================================!
module gfnff_cn
  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData
  implicit none
  private
  public :: gfnff_dlogcoord

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CN routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> logCN derivative saved in dlogCN array
  subroutine gfnff_dlogcoord(n,at,xyz,rab,logCN,dlogCN,thr2,param)
    implicit none
    type(TGFFData),intent(in) :: param
    integer,intent(in)  :: n
    integer,intent(in)  :: at(n)
    real(wp),intent(in)  :: xyz(3,n)
    real(wp),intent(in)  :: rab(n*(n+1)/2)
    real(wp),intent(out) :: logCN(n)
    real(wp),intent(out) :: dlogCN(3,n,n)
    real(wp),intent(in)  :: thr2
    real(wp)             :: cn(n)
    !> counter
    integer  :: i,j,ij,ii
    !> distances
    real(wp) :: r
    real(wp) :: r0
    real(wp) :: dr
    real(wp),dimension(3) :: rij
    !> treshold
    real(wp) :: thr
    !> cn parameter
    real(wp) :: erfCN
    real(wp) :: dlogdcni
    real(wp) :: dlogdcnj
    real(wp) :: derivative
    !> local parameter
    real(wp),parameter :: kn = -7.5_wp
    !> initialize vaiiables
    cn = 0.0_wp
    logCN = 0.0_wp
    dlogCN = 0.0_wp
    dlogdcni = 0.0_wp
    dlogdcnj = 0.0_wp
    derivative = 0.0_wp
    r = 0.0_wp
    r0 = 0.0_wp
    rij = 0.0_wp
    thr = sqrt(thr2)
    !> create error function CN
    do i = 2,n
      ii = i*(i-1)/2
      do j = 1,i-1
        ij = ii+j
        r = rab(ij)
        if (r .gt. thr) cycle
        r0 = (param%rcov(at(i))+param%rcov(at(j)))
        dr = (r-r0)/r0
        !> hier kommt die CN funktion hin
!           erfCN = create_expCN(16.0d0,r,r0)
        erfCN = 0.5_wp*(1.0_wp+erf(kn*dr))
        cn(i) = cn(i)+erfCN
        cn(j) = cn(j)+erfCN
      end do
    end do
    !> create cutted logarithm CN + derivatives
    do i = 1,n
      ii = i*(i-1)/2
      logCN(i) = create_logCN(cn(i),param)
      !> get dlogCN/dCNi
      dlogdcni = create_dlogCN(cn(i),param)
      do j = 1,i-1
        ij = ii+j
        !> get dlogCN/dCNj
        dlogdcnj = create_dlogCN(cn(j),param)
        r = rab(ij)
        if (r .gt. thr) cycle
        r0 = (param%rcov(at(i))+param%rcov(at(j)))
        !> get derfCN/dRij
        derivative = create_derfCN(kn,r,r0)
!           derivative = create_dexpCN(16.0d0,r,r0)
        rij = derivative*(xyz(:,j)-xyz(:,i))/r
        !> project rij gradient onto cartesians
        dlogCN(:,j,j) = dlogdcnj*rij+dlogCN(:,j,j)
        dlogCN(:,i,j) = -dlogdcnj*rij
        dlogCN(:,j,i) = dlogdcni*rij
        dlogCN(:,i,i) = -dlogdcni*rij+dlogCN(:,i,i)
      end do
    end do

  contains

    pure elemental function create_logCN(cn,param) result(count)
      type(TGFFData),intent(in) :: param
      real(wp),intent(in) :: cn
      real(wp) :: count
      count = log(1+exp(param%cnmax))-log(1+exp(param%cnmax-cn))
    end function create_logCN

    pure elemental function create_dlogCN(cn,param) result(count)
      type(TGFFData),intent(in) :: param
      real(wp),intent(in) :: cn
      real(wp) :: count
      count = exp(param%cnmax)/(exp(param%cnmax)+exp(cn))
    end function create_dlogCN

    pure elemental function create_erfCN(k,r,r0) result(count)
      real(wp),intent(in) :: k
      real(wp),intent(in) :: r
      real(wp),intent(in) :: r0
      real(wp) :: count
      real(wp) :: dr
      dr = (r-r0)/r0
      count = 0.5_wp*(1.0_wp+erf(kn*dr))
    end function create_erfCN

    pure elemental function create_derfCN(k,r,r0) result(count)
      real(wp),intent(in) :: k
      real(wp),intent(in) :: r
      real(wp),intent(in) :: r0
      real(wp),parameter :: sqrtpi = 1.77245385091_wp
      real(wp) :: count
      real(wp) :: dr
      dr = (r-r0)/r0
      count = k/sqrtpi*exp(-k**2*dr*dr)/r0
    end function create_derfCN

    pure elemental function create_expCN(k,r,r0) result(count)
      real(wp),intent(in) :: k
      real(wp),intent(in) :: r
      real(wp),intent(in) :: r0
      real(wp) :: count
      count = 1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))
    end function create_expCN

    pure elemental function create_dexpCN(k,r,r0) result(count)
      real(wp),intent(in) :: k
      real(wp),intent(in) :: r
      real(wp),intent(in) :: r0
      real(wp) :: count
      real(wp) :: expterm
      expterm = exp(-k*(r0/r-1._wp))
      count = (-k*r0*expterm)/(r**2*((expterm+1._wp)**2))
    end function create_dexpCN

  end subroutine gfnff_dlogcoord

end module gfnff_cn

