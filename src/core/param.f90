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

!> Turns the fixed tables in gfnff_param_tables into a populated TGFFData,
!> reads a parameter file when supplied, and owns the accuracy-to-cutoff mapping.
module gfnff_param
  use iso_fortran_env,only:wp => real64
  use gfnff_data_types,only:TGFFData,TGFFGenerator,init
  !> No only-list: this module turns essentially every table into a TGFFData.
  use gfnff_param_tables
  implicit none
  private

  public :: gfnff_set_param,gfnff_load_param,gfnff_read_param
  public :: gfnff_thresholds
  !> Re-exported from gfnff_param_tables so a consumer needs only one use statement.
  public :: gffVersion,pse,sqrtZr4r2,covalentRadD3,paulingEN

contains  !> MODULE PROCEDURES START HERE

  subroutine gfnff_set_param(n,gen,param)
    !***********************************************************************
    !* Populates param (TGFFData) with the fitted/derived scalars used by
    !* the energy and gradient terms: angle/torsion damping, HB/XB cutoffs,
    !* acidities and basicities, and the 3-body/D3 prefactors built from gen.
    !***********************************************************************
    implicit none
    integer,intent(in)  :: n
    type(TGFFGenerator),intent(out) :: gen
    type(TGFFData),intent(inout) :: param
    integer   :: i,j,k
    real(wp)  :: dum

    if (.false.) write (*,*) n  ! silences -Wunused-dummy-argument

    call newGFNFFGenerator(gen)

    param%cnmax = 4.4         ! max. CN considered ie all larger values smoothly set to this val
    param%atcuta = 0.595_wp     ! angle damping
    param%atcutt = 0.505_wp     ! torsion angle damping
    param%atcuta_nci = 0.395_wp ! nci angle damping in HB term
    param%atcutt_nci = 0.305_wp ! nci torsion angle damping in HB term
    param%repscalb = 1.7583      ! bonded rep. scaling
    param%repscaln = 0.4270      ! non-bonded rep. scaling
    param%hbacut = 49.0_wp       ! HB angle cut-off
    param%hbscut = 22.0_wp       ! HB SR     "   "
    param%xbacut = 70.0_wp       ! same for XB
    param%xbscut = 5.0_wp        !
    param%hbsf = 1.0_wp          ! charge dep.
    param%hbst = 15.0_wp         ! 10 is better for S22, 20 better for HCN2 and S30L
    param%xbsf = 0.03_wp         !
    param%xbst = 15.0_wp         !
    param%hbalp = 6.0_wp         ! damp
    param%hblongcut = 85.0_wp    ! values larger than 85 yield large RMSDs for P26
    param%hblongcut_xb = 70.0_wp ! values larger than 70 yield large MAD for HAL28
    param%hbabmix = 0.80         !
    param%hbnbcut = 11.20        !
    param%tors_hb = 0.94         ! torsion potential shift in HB term
    param%bend_hb = 0.20         ! bending potential shift in HB term
    param%vbond_scale = 0.9      ! vbond(2) scaling for CN(H) = 1
    param%xhaci_globabh = 0.268  ! A-H...B gen. scaling
    param%xhaci_coh = 0.350      ! A-H...O=C gen. scaling
    param%xhaci_glob = 1.50      ! acidity
    param%xhbas(:) = 0.0_wp
    param%xhbas(6) = 0.80_wp     ! basicities (XB and HB), i.e., B...X-A or B...H..A
    param%xhbas(7) = 1.68_wp
    param%xhbas(8) = 0.67_wp
    param%xhbas(9) = 0.52_wp
    param%xhbas(14) = 4.0_wp
    param%xhbas(15) = 3.5_wp
    param%xhbas(16) = 2.0_wp
    param%xhbas(17) = 1.5_wp
    param%xhbas(35) = 1.5_wp
    param%xhbas(53) = 1.9_wp
    param%xhbas(33) = param%xhbas(15)
    param%xhbas(34) = param%xhbas(16)
    param%xhbas(51) = param%xhbas(15)
    param%xhbas(52) = param%xhbas(16)
    param%xhaci(:) = 0.0_wp
    param%xhaci(6) = 0.75               ! HB acidities, a bit weaker for CH
    param%xhaci(7) = param%xhaci_glob+0.1
    param%xhaci(8) = param%xhaci_glob
    param%xhaci(9) = param%xhaci_glob
    param%xhaci(15) = param%xhaci_glob
    param%xhaci(16) = param%xhaci_glob
    param%xhaci(17) = param%xhaci_glob+1.0
    param%xhaci(35) = param%xhaci_glob+1.0
    param%xhaci(53) = param%xhaci_glob+1.0
    param%xbaci(:) = 0.0_wp
    param%xbaci(15) = 1.0_wp              ! XB acidities
    param%xbaci(16) = 1.0_wp
    param%xbaci(17) = 0.5_wp
    param%xbaci(33) = 1.2_wp
    param%xbaci(34) = 1.2_wp
    param%xbaci(35) = 0.9_wp
    param%xbaci(51) = 1.2_wp
    param%xbaci(52) = 1.2_wp
    param%xbaci(53) = 1.2_wp

    !>-- 3-body bond prefactors and D3 R0^2 lookup
    k = 0
    do i = 1,86
      dum = dble(i)
      param%zb3atm(i) = -dum*gen%batmscal**(1.0_wp/3.0_wp)  ! inlcude pre-factor
      do j = 1,i
        k = k+1
        dum = sqrtZr4r2(i)*sqrtZr4r2(j)*3.0_wp
        param%d3r0(k) = (gen%d3a1*dsqrt(dum)+gen%d3a2)**2   ! save R0^2 for efficiency reasons
      end do
    end do
    param%zb3atm(1) = -0.25_wp*gen%batmscal**(1.0_wp/3.0_wp) ! slightly better than 1.0

  end subroutine gfnff_set_param

  subroutine gfnff_thresholds(accuracy,dispthr,cnthr,repthr,hbthr1,hbthr2)
    !***********************************************************************
    !* Maps the accuracy dial to pair cutoffs (dispersion, CN, repulsion,
    !* HB) used elsewhere to skip negligible interactions.
    !***********************************************************************
    real(wp),intent(in) :: accuracy
    real(wp),intent(out) :: dispthr
    real(wp),intent(out) :: cnthr
    real(wp),intent(out) :: repthr
    real(wp),intent(out) :: hbthr1
    real(wp),intent(out) :: hbthr2
    dispthr = 1500.0_wp-log10(accuracy)*1000._wp
    cnthr = 100.0_wp-log10(accuracy)*50.0_wp
    repthr = 400.0_wp-log10(accuracy)*100.0_wp
    hbthr1 = 200.0_wp-log10(accuracy)*50.0_wp
    hbthr2 = 400.0_wp-log10(accuracy)*50.0_wp
  end subroutine gfnff_thresholds

  subroutine gfnff_load_param(version,param,exist)
    !***********************************************************************
    !* Fills param with the fixed elemental tables and the per-element set
    !* matching version; exist is false if version has no parameter set.
    !***********************************************************************
    implicit none
    integer,intent(in) :: version
    type(TGFFData),intent(out) :: param
    logical,intent(out) :: exist

    exist = .false.

    call init(param,103)

    param%en(:) = en
    param%rad(:) = rad
    param%rcov(:) = covalentRadD3(1:103)
    param%metal(:) = metal
    param%group(:) = group
    param%normcn(:) = normcn
    param%repz(:) = repz

    select case (version)
    case (gffVersion%angewChem2020,gffVersion%angewChem2020_1, &
    &     gffVersion%angewChem2020_2,gffVersion%harmonic2020,gffVersion%mcgfnff2023, &
    &     gffVersion%conformer2020)
      call loadGFNFFAngewChem2020(param)
      exist = .true.
    end select

  end subroutine gfnff_load_param

  subroutine loadGFNFFAngewChem2020(param)
    !***********************************************************************
    !* Copies the AngewChem2020 per-element parameter tables into param.
    !***********************************************************************
    type(TGFFData),intent(inout) :: param
    param%chi(:) = chi_angewChem2020
    param%gam(:) = gam_angewChem2020
    param%cnf(:) = cnf_angewChem2020
    param%alp(:) = alp_angewChem2020
    param%bond(:) = bond_angewChem2020
    param%repa(:) = repa_angewChem2020
    param%repan(:) = repan_angewChem2020
    param%angl(:) = angl_angewChem2020
    param%angl2(:) = angl2_angewChem2020
    param%tors(:) = tors_angewChem2020
    param%tors2(:) = tors2_angewChem2020
  end subroutine loadGFNFFAngewChem2020

  subroutine gfnff_read_param(iunit,param)
    !***********************************************************************
    !* Fills param with the fixed elemental tables, then reads the 86
    !* per-element parameters (columns 2-12) from a text file on iunit.
    !***********************************************************************
    use gfnff_strings,only:readl
    implicit none
    integer,intent(in)  :: iunit
    type(TGFFData),intent(out) :: param
    integer  :: i,nn
    real(wp) :: xx(20)
    character(len=256) :: atmp

    call init(param,103)

    param%en(:) = en
    param%rad(:) = rad
    param%rcov(:) = covalentRadD3(1:103)
    param%metal(:) = metal
    param%group(:) = group
    param%normcn(:) = normcn
    param%repz(:) = repz

    do i = 1,86
      read (iunit,'(a)') atmp
      call readl(atmp,xx,nn)
      param%chi(i) = xx(2)
      param%gam(i) = xx(3)
      param%cnf(i) = xx(4)
      param%alp(i) = xx(5)
      param%bond(i) = xx(6)
      param%repa(i) = xx(7)
      param%repan(i) = xx(8)
      param%angl(i) = xx(9)
      param%angl2(i) = xx(10)
      param%tors(i) = xx(11)
      param%tors2(i) = xx(12)
    end do

  end subroutine gfnff_read_param

  subroutine newGFNFFGenerator(gen)
    !***********************************************************************
    !* Populates gen with the fitted GFN-FF generator constants: angle and
    !* torsion damping, HB/XB terms, Hueckel iteration, the bond-strength
    !* matrix bsmat, and related topology thresholds.
    !***********************************************************************
    type(TGFFGenerator),intent(out) :: gen

    gen%cnmax = 4.4         ! max. CN considered ie all larger values smoothly set to this val
    gen%linthr = 160.       ! angle considered linear above this; closer to 170 better for metals but unclear for Sc, kept at 160
    gen%fcthr = 1.d-3       ! skip torsion and bending if potential is small
    gen%tdist_thr = 12.     ! R threshold in Angstroem for cov distance estimated used in apprx EEQ
    gen%rthr = 1.25         ! bond determination threshold, critical for topo setup; large values yield more 1.23
    gen%rthr2 = 1.00        ! decrease if a metal is present, larger values yield smaller CN
    gen%rqshrink = 0.23     ! change of R0 for topo with charge qa, larger values yield smaller CN for metals in particular
    gen%hqabthr = 0.01      ! H charge (qa) threshold for H in HB list 18
    !>-- larger values are better for S30L but worse in PubChem RMSD checks
    gen%qabthr = 0.10       ! AB charge (qa) threshold for AB in HB list, avoids HBs with positive atoms
    gen%srb1 = 0.3731      ! bond params
    gen%srb2 = 0.3171      !
    gen%srb3 = 0.2538      !
    gen%qrepscal = 0.3480     ! change of non-bonded rep. with q(topo)
    gen%nrepscal = -0.1270    !   "    "      "       "   CN
    gen%hhfac = 0.6290        ! HH repulsion
    gen%hh13rep = 1.4580      !
    gen%hh14rep = 0.7080      !
    gen%bstren(1) = 1.00_wp      ! single bond
    gen%bstren(2) = 1.24_wp      ! double bond
    gen%bstren(3) = 1.98_wp      ! triple bond
    gen%bstren(4) = 1.22_wp      ! hyperval bond
    gen%bstren(5) = 1.00_wp      ! M-X
    gen%bstren(6) = 0.78_wp      ! M eta
    gen%bstren(7) = 3.40_wp      ! M-M
    gen%bstren(8) = 3.40_wp      ! M-M
    gen%qfacBEN = -0.54_wp     ! bend FC change with polarity
    gen%qfacTOR = 12.0_wp      ! torsion FC change with polarity
    gen%fr3 = 0.3          ! tors FC 3-ring
    gen%fr4 = 1.0          ! tors FC 4-ring
    gen%fr5 = 1.5          ! tors FC 5-ring
    gen%fr6 = 5.7          ! tors FC 6-ring
    gen%torsf(1) = 1.00        ! single bond
    gen%torsf(2) = 1.18        ! pi bond
    gen%torsf(3) = 1.05        ! improper
    gen%torsf(5) = 0.50        ! pi part improper
    gen%torsf(6) = -0.90       ! extra sp3 C
    gen%torsf(7) = 0.70        ! extra sp3 N
    gen%torsf(8) = -2.00       ! extra sp3 O
    gen%fbs1 = 0.50            ! small bend corr.
    gen%batmscal = 0.30_wp     ! bonded ATM scal
    gen%mchishift = -0.09_wp
    gen%rabshift = -0.110      ! gen shift
    gen%rabshifth = -0.050     ! XH
    gen%hyper_shift = 0.03     ! hypervalent
    gen%hshift3 = -0.11        ! heavy
    gen%hshift4 = -0.11        !
    gen%hshift5 = -0.06        !
    gen%metal1_shift = 0.2     ! group 1+2 metals
    gen%metal2_shift = 0.15    ! TM
    gen%metal3_shift = 0.05    ! main group metals
    gen%eta_shift = 0.040      ! eta bonded
    gen%qfacbm(0) = 1.0_wp     ! bond charge dep.

    gen%qfacbm(1:2) = -0.2_wp  !
    gen%qfacbm(3) = 0.70_wp    !
    gen%qfacbm(4) = 0.50_wp    !
    gen%qfacbm0 = 0.047        !
    gen%rfgoed1 = 1.175        ! topo dist scaling
    gen%htriple = 1.45_wp      ! decrease Hueckel off-diag for triple bonds because they are less well conjugated 1.4
    gen%hueckelp2 = 1.00_wp    ! increase pot depth depending on P
    gen%hueckelp3 = -0.24_wp   ! diagonal element change with qa
    gen%hdiag(5) = -0.5_wp     ! diagonal element relative to C
    gen%hdiag(6) = 0.00_wp     !
    gen%hdiag(7) = 0.14_wp     !
    gen%hdiag(8) = -0.38_wp    !
    gen%hdiag(9) = -0.29_wp    !
    gen%hdiag(16) = -0.30_wp   !
    gen%hdiag(17) = -0.30_wp   !
    gen%hoffdiag(5) = 0.5_wp   ! Hückel off-diag constants
    gen%hoffdiag(6) = 1.00_wp  !
    gen%hoffdiag(7) = 0.66_wp  !
    gen%hoffdiag(8) = 1.10_wp  !
    gen%hoffdiag(9) = 0.23_wp  !
    gen%hoffdiag(16) = 0.60_wp !
    gen%hoffdiag(17) = 1.00_wp !
    gen%hiter = 0.700_wp       ! iteration mixing
    gen%hueckelp = 0.340_wp    ! diagonal qa dep.
    gen%bzref = 0.370_wp       ! ref P value R shift
    gen%bzref2 = 0.315_wp      !  "  "  "    k stretch
    gen%pilpf = 0.530_wp       ! 2el diag shift
    gen%maxhiter = 5           ! the Hückel iterations can diverge so take only a few steps
    gen%d3a1 = 0.58_wp         ! D3, s8 fixed = 2
    gen%d3a2 = 4.80_wp
    gen%split0 = 0.670_wp      ! mixing of sp^n with sp^n-1
    gen%fringbo = 0.020_wp     ! str ring size dep.
    gen%aheavy3 = 89.          ! three coord. heavy eq. angle
    gen%aheavy4 = 100.         ! four   "       "    "    "
    gen%split1 = 1.0_wp-gen%split0
    gen%bsmat = -999.
    gen%bsmat(0,0) = gen%bstren(1)
    gen%bsmat(3,0) = gen%bstren(1)
    gen%bsmat(3,3) = gen%bstren(1)
    gen%bsmat(2,2) = gen%bstren(2)
    gen%bsmat(1,1) = gen%bstren(3)
    gen%bsmat(1,0) = gen%split0*gen%bstren(1)+gen%split1*gen%bstren(3)
    gen%bsmat(3,1) = gen%split0*gen%bstren(1)+gen%split1*gen%bstren(3)
    gen%bsmat(2,1) = gen%split0*gen%bstren(2)+gen%split1*gen%bstren(3)
    gen%bsmat(2,0) = gen%split0*gen%bstren(1)+gen%split1*gen%bstren(2)
    gen%bsmat(3,2) = gen%split0*gen%bstren(1)+gen%split1*gen%bstren(2)
    gen%bstren(9) = 0.5*(gen%bstren(7)+gen%bstren(8))

  end subroutine newGFNFFGenerator

end module gfnff_param
