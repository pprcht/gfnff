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


!> D3 van-der-Waals radii
module solvation_param_vdwradd3
   use iso_fortran_env, only : wp => real64
   implicit none
   private

   public :: getVanDerWaalsRadD3, vanDerWaalsRadD3

   !>  convert bohr (a.u.) to Ångström and back
   real(wp),public,parameter :: autoaa = 0.52917726_wp
   real(wp),public,parameter :: aatoau = 1.0_wp/autoaa

   !> Get van-der-Waals Rad for a species
   interface getVanDerWaalsRadD3
      module procedure :: getVanDerWaalsRadD3Number
   end interface getVanDerWaalsRadD3


   !> D3 pairwise van-der-Waals radii (only homoatomic pairs present here)
   real(wp), parameter :: vanDerWaalsRadD3(1:94) = aatoau * [&
      & 1.09155_wp, 0.86735_wp, 1.74780_wp, 1.54910_wp, &
      & 1.60800_wp, 1.45515_wp, 1.31125_wp, 1.24085_wp, &
      & 1.14980_wp, 1.06870_wp, 1.85410_wp, 1.74195_wp, &
      & 2.00530_wp, 1.89585_wp, 1.75085_wp, 1.65535_wp, &
      & 1.55230_wp, 1.45740_wp, 2.12055_wp, 2.05175_wp, &
      & 1.94515_wp, 1.88210_wp, 1.86055_wp, 1.72070_wp, &
      & 1.77310_wp, 1.72105_wp, 1.71635_wp, 1.67310_wp, &
      & 1.65040_wp, 1.61545_wp, 1.97895_wp, 1.93095_wp, &
      & 1.83125_wp, 1.76340_wp, 1.68310_wp, 1.60480_wp, &
      & 2.30880_wp, 2.23820_wp, 2.10980_wp, 2.02985_wp, &
      & 1.92980_wp, 1.87715_wp, 1.78450_wp, 1.73115_wp, &
      & 1.69875_wp, 1.67625_wp, 1.66540_wp, 1.73100_wp, &
      & 2.13115_wp, 2.09370_wp, 2.00750_wp, 1.94505_wp, &
      & 1.86900_wp, 1.79445_wp, 2.52835_wp, 2.59070_wp, &
      & 2.31305_wp, 2.31005_wp, 2.28510_wp, 2.26355_wp, &
      & 2.24480_wp, 2.22575_wp, 2.21170_wp, 2.06215_wp, &
      & 2.12135_wp, 2.07705_wp, 2.13970_wp, 2.12250_wp, &
      & 2.11040_wp, 2.09930_wp, 2.00650_wp, 2.12250_wp, &
      & 2.04900_wp, 1.99275_wp, 1.94775_wp, 1.87450_wp, &
      & 1.72280_wp, 1.67625_wp, 1.62820_wp, 1.67995_wp, &
      & 2.15635_wp, 2.13820_wp, 2.05875_wp, 2.00270_wp, &
      & 1.93220_wp, 1.86080_wp, 2.53980_wp, 2.46470_wp, &
      & 2.35215_wp, 2.21260_wp, 2.22970_wp, 2.19785_wp, &
      & 2.17695_wp, 2.21705_wp]


contains

!> Get van-der-Waals radius for species with a given atomic number
elemental function getVanDerWaalsRadD3Number(number) result(rad)

   !> Atomic number
   integer, intent(in) :: number

   !> van-der-Waals radius
   real(wp) :: rad

   if (number > 0 .and. number <= size(vanDerWaalsRadD3, dim=1)) then
      rad = vanDerWaalsRadD3(number)
   else
      rad = -1.0_wp
   end if

end function getVanDerWaalsRadD3Number


end module solvation_param_vdwradd3

! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Cosmo van-der-Waals radii

module solvation_param_vdwradcosmo
   use iso_fortran_env, only : wp => real64
   implicit none
   private

   public :: getVanDerWaalsRadCosmo, vanDerWaalsRadCosmo

   !>  convert bohr (a.u.) to Ångström and back
   real(wp),public,parameter :: autoaa = 0.52917726_wp
   real(wp),public,parameter :: aatoau = 1.0_wp/autoaa

   !> Get van-der-Waals Rad for a species
   interface getVanDerWaalsRadCosmo
      module procedure :: getVanDerWaalsRadCosmoNumber
   end interface getVanDerWaalsRadCosmo


   !> Default value for unoptimized van-der-Waals radii
   real(wp), parameter :: cosmostub = 2.223_wp

   !> COSMO optimized van-der-Waals radii
   real(wp), parameter :: vanDerWaalsRadCosmo(94) = aatoau * [ &
       & 1.3000_wp, 1.6380_wp, 1.5700_wp, 1.0530_wp, &   ! h-be
       & 2.0480_wp, 2.0000_wp, 1.8300_wp, 1.7200_wp, &   ! B-O
       & 1.7200_wp, 1.8018_wp, 1.8000_wp, 1.6380_wp, &   ! F-Mg
       & 2.1530_wp, 2.2000_wp, 2.1060_wp, 2.1600_wp, &   ! Al-S
       & 2.0500_wp, 2.2000_wp, 2.2230_wp, cosmostub, &   ! Cl-Ca
       & cosmostub, 2.2930_wp, cosmostub, cosmostub, &   ! Sc-Cr
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Mn-Ni
       & cosmostub, 1.6260_wp, cosmostub, 2.7000_wp, &   ! Cu-Ge
       & 2.3500_wp, 2.2000_wp, 2.1600_wp, 2.3630_wp, &   ! As-Kr
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Rb-Zr
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Nb-Ru
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Rh-Cd
       & 2.2580_wp, 2.5500_wp, 2.4100_wp, 2.4100_wp, &   ! In-Te
       & 2.3200_wp, 2.5270_wp, cosmostub, cosmostub, &   ! I-Ba
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! La-Nd
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Pm-Gd
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Tb-Er
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Tm-Hf
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ta-Os
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ir-Hg
       & cosmostub, 2.3600_wp, 2.4220_wp, 2.3050_wp, &   ! Tl-Po
       & 2.3630_wp, 2.5740_wp, cosmostub, cosmostub, &   ! At-Ra
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ac-U
       & cosmostub, cosmostub]                           ! Np-Pu


contains

!> Get van-der-Waals radius for species with a given atomic number
elemental function getVanDerWaalsRadCosmoNumber(number) result(rad)

   !> Atomic number
   integer, intent(in) :: number

   !> van-der-Waals radius
   real(wp) :: rad

   if (number > 0 .and. number <= size(vanDerWaalsRadCosmo, dim=1)) then
      rad = vanDerWaalsRadCosmo(number)
   else
      rad = -1.0_wp
   end if

end function getVanDerWaalsRadCosmoNumber


end module solvation_param_vdwradcosmo
