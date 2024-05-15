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
!> Topological data for force field type calculations and neighbor lists
module gfnff_data_types
  use iso_fortran_env, only: wp => real64, sp=>real32
  implicit none
  private

  !> public types and routines
  public :: TGFFTopology
  public :: TGFFNeighbourList,new
  public :: TGFFData,init
  public :: TDispersionData,initgffdispersion
  public :: TGFFGenerator

!========================================================================================!

  !> Data for the dispersion contribution
  type :: TDispersionData

    !> Damping parameters
    real(wp) :: s6 = 1.0_wp
    real(wp) :: s8 = 0.0_wp
    real(wp) :: s10 = 0.0_wp
    real(wp) :: a1 = 0.0_wp
    real(wp) :: a2 = 0.0_wp
    real(wp) :: s9 = 0.0_wp
    integer  :: alp = 16

    !> Weighting factor for Gaussian interpolation
    real(wp) :: wf

    !> Charge steepness
    real(wp) :: g_a

    !> Charge height
    real(wp) :: g_c

    !> Reference data for the dispersion
    integer,allocatable :: atoms(:)
    integer,allocatable :: nref(:)
    integer,allocatable :: ncount(:,:)
    real(wp),allocatable :: cn(:,:)
    real(wp),allocatable :: q(:,:)
    real(wp),allocatable :: alpha(:,:,:)
    real(wp),allocatable :: c6(:,:,:,:)

  end type TDispersionData
!========================================================================================!

  !> Topology information for a given system
  type :: TGFFTopology

    !> some reference files
    character(len=:),allocatable :: filename
    character(len=:),allocatable :: refcharges

    !> number of terms
    integer  :: nbond
    integer  :: nangl
    integer  :: ntors
    integer  :: nathbH
    integer  :: nathbAB
    integer  :: natxbAB
    integer  :: nbatm
    integer  :: nfrag
    integer  :: maxsystem   !> max. number of fragmentsfor hessian
    integer  :: bond_hb_nr  !> number of unique AH...B HB/bond terms
    integer  :: b_max       !> number of B atoms per unique AH bond

    !> numbers that are rewritten, so must be stored for allocation
    integer  :: nbond_blist
    integer  :: nbond_vbond
    integer  :: nangl_alloc
    integer  :: ntors_alloc

    !> file type read
    integer  :: read_file_type

    !> lists
    integer,allocatable ::     nb(:,:)   ! neighbors nb(20,i) is the # neigbors
    integer,allocatable ::      hyb(:)   ! hybridization of every atom
    integer,allocatable ::    bpair(:)   ! # of cov. between atoms
    integer,allocatable ::  blist(:,:)   ! bonded atoms
    integer,allocatable ::  alist(:,:)   ! angles
    integer,allocatable ::  tlist(:,:)   ! torsions
    integer,allocatable :: b3list(:,:)   ! bond atm
    integer,allocatable :: sTorsl(:,:)   ! triple bonded carbon potential
    real(wp),allocatable::      pbo(:)   ! bond order from Hückel
    !-----------------------------------------------
    integer,allocatable :: nr_hb(:)        ! Nr. of H bonds per O-H or N-H bond
    integer,allocatable :: bond_hb_AH(:,:) ! A, H atoms in bonds that are also part of HBs
    integer,allocatable :: bond_hb_B(:,:)  ! B atoms in bonds that are also part of HBs
    integer,allocatable :: bond_hb_Bn(:)   ! Nr. of B atoms for one AH bond pair
    !-----------------------------------------------
    integer,allocatable :: hbatABl(:,:) ! AB atoms for HB
    integer,allocatable :: xbatABl(:,:) ! AB atoms for XB
    integer,allocatable :: hbatHl(:)    ! H  atoms for HB
    integer,allocatable :: fraglist(:)  ! atoms in molecular fragments (for EEQ)
    integer,allocatable :: qpdb(:)      ! atomic charge in residues from PDB file

    !potential parameters used in energy-gradient routine
    real(wp),allocatable :: vbond(:,:) ! bonds
    real(wp),allocatable :: vangl(:,:) ! angles
    real(wp),allocatable :: vtors(:,:) ! torsions
    real(wp),allocatable :: chieeq(:)  ! atomic ENs for EEQ
    real(wp),allocatable :: gameeq(:)  ! atomic gamma for EEQ
    real(wp),allocatable :: alpeeq(:)  ! atomic alpha for EEQ, squared
    real(wp),allocatable :: alphanb(:) ! non-bonded exponent for atom pairs
    real(wp),allocatable :: qa(:)      ! estimated atomic charges (fixed and obtained from topology EEQ)
    real(wp),allocatable :: xyze0(:,:) ! atom xyz, starting geom. (for Efield energy)
    real(wp),allocatable :: zetac6(:)  ! D4 scaling factor product
    real(wp),allocatable :: qfrag(:)   ! fragment charge (for EEQ)
    real(wp),allocatable :: hbbas(:)   ! HB donor atom basicity
    real(wp),allocatable :: hbaci(:)   ! HB acceptor atom acidity

    integer,allocatable  :: ispinsyst(:,:)
    integer,allocatable  :: nspinsyst(:)
    integer              :: nsystem

    type(TDispersionData) :: dispm

  contains

    procedure :: zero

  end type TGFFTopology
!========================================================================================!

  !> Neighbourlist storage
  type :: TGFFNeighbourList
    logical :: initialized = .false.
    logical :: force_hbond_update = .false.
    integer :: nhb1
    integer :: nhb2
    integer :: nxb
    !> atomic charges (obtained from EEQ)
    real(wp),allocatable :: q(:)
    !> atom xyz, used to check for HB list update
    real(wp),allocatable :: hbrefgeo(:,:)
    !> HBs loose
    integer,allocatable :: hblist1(:,:)
    !> HBs bonded
    integer,allocatable :: hblist2(:,:)
    !> XBs
    integer,allocatable :: hblist3(:,:)
  end type TGFFNeighbourList

  interface new
    module procedure :: newGFFNeighbourList
  end interface
!========================================================================================!

  !> Parametrisation data for the force field
  type :: TGFFData

    !> repulsion scaling
    real(wp) :: repscaln
    real(wp) :: repscalb

    !> bend/tors angle damping
    real(wp) :: atcuta
    real(wp) :: atcutt

    !> bend/tors nci angle damping for HB term
    real(wp) :: atcuta_nci
    real(wp) :: atcutt_nci

    !> damping HB
    real(wp) :: hbacut
    real(wp) :: hbscut

    !> damping XB
    real(wp) :: xbacut
    real(wp) :: xbscut

    !> damping HB/XB
    real(wp) :: hbalp

    !> damping HB/XB
    real(wp) :: hblongcut
    real(wp) :: hblongcut_xb

    !> charge scaling HB/XB
    real(wp) :: hbst
    real(wp) :: hbsf
    real(wp) :: xbst
    real(wp) :: xbsf

    !> HB AH-B
    real(wp) :: xhaci_globabh

    !> HB AH-O=C
    real(wp) :: xhaci_coh

    !> acidity
    real(wp) :: xhaci_glob

    !> HB AH-B
    real(wp) :: hbabmix

    !> new parameter for neighbour angle
    real(wp) :: hbnbcut

    !> new parameter for HB NCI angle term
    real(wp) :: tors_hb

    !> new parameter for HB NCI torsion term
    real(wp) :: bend_hb

    !> new parameter for FC scaling of bonds in HB
    real(wp) :: vbond_scale

    !> max CN cut-off
    real(wp) :: cnmax

    !> D3 scaling
    real(wp) :: dispscale

    !> Constant data
    real(wp),allocatable :: en(:)
    real(wp),allocatable :: rad(:)
    real(wp),allocatable :: rcov(:)
    integer,allocatable :: metal(:)
    integer,allocatable :: group(:)
    integer,allocatable :: normcn(:)

    !> rep alpha bond
    real(wp),allocatable :: repa(:)
    real(wp),allocatable :: repan(:)

    !> prefactor (Zval), 3atm bond term
    real(wp),allocatable :: repz(:)
    real(wp),allocatable :: zb3atm(:)

    !> HB/XB
    real(wp),allocatable :: xhaci(:)
    real(wp),allocatable :: xhbas(:)
    real(wp),allocatable :: xbaci(:)

    !> EN dep. in EEQ.
    real(wp),allocatable :: chi(:)
    real(wp),allocatable :: gam(:)
    real(wp),allocatable :: cnf(:)
    real(wp),allocatable :: alp(:)

    !> Elem. bond param.
    real(wp),allocatable :: bond(:)

    !> Elem. angular param.
    real(wp),allocatable :: angl(:)

    !> Elem. angular param.
    real(wp),allocatable :: angl2(:)

    !> Elem. torsion param_alloc.
    real(wp),allocatable :: tors(:)

    !> Elem. torsion param.
    real(wp),allocatable :: tors2(:)

    !> BJ radii set in gnff_ini()
    real(wp),allocatable :: d3r0(:)

  end type TGFFData

  !> Initialize a new instance of the parametrisation data
  interface init
    module procedure :: initGFFData
  end interface init

!========================================================================================!


   !> Generator for the force field topology
   type TGFFGenerator

      !> when is an angle close to linear ? (GEODEP)
      !  for metals values closer to 170 (than to 160) are better
      real(wp) :: linthr

      !> skip torsion and bending if potential is small
      real(wp) :: fcthr

      !> R threshold in Angstroem for cov distance estimated used in apprx EEQ
      real(sp) :: tdist_thr

      !> important bond determination threshold, large values yield more 1.23
      real(wp) :: rthr

      !> decrease if a metal is present, larger values yield smaller CN
      real(wp) :: rthr2

      !> change of R0 for topo with charge qa
      !  larger values yield smaller CN for metals in particular
      real(wp) :: rqshrink

      !> H charge (qa) threshold for H in HB list 18
      real(wp) :: hqabthr

      !> AB charge (qa) threshold for AB in HB list
      !  - avoids HBs with positive atoms,
      !  - larger val. better for S30L but worse in PubChem RMSD checks
      real(wp) :: qabthr

      !> Parameter
      real(wp) :: srb1
      real(wp) :: srb2
      real(wp) :: srb3

      !> change of non-bonded rep. with q(topo)
      real(wp) :: qrepscal

      !> change of non-bonded rep. with CN
      real(wp) :: nrepscal

      !> HH repulsion
      real(wp) :: hhfac
      real(wp) :: hh13rep
      real(wp) :: hh14rep
      real(wp) :: bstren(9)

      !> bend FC change with polarity
      real(wp) :: qfacBEN

      !> torsion FC change with polarity
      real(wp) :: qfacTOR

      !> tors FC 3-ring
      real(wp) :: fr3

      !> tors FC 4-ring
      real(wp) :: fr4

      !> tors FC 5-ring
      real(wp) :: fr5

      !> tors FC 6-ring
      real(wp) :: fr6

      !> bonds
      real(wp) :: torsf(8)

      !> small bend corr.
      real(wp) :: fbs1

      !> bonded ATM scal
      real(wp) :: batmscal

      !> Shifts
      real(wp) :: mchishift

      !> gen shift
      real(wp) :: rabshift

      !> XH
      real(wp) :: rabshifth

      !> hypervalent
      real(wp) :: hyper_shift

      !> heavy
      real(wp) :: hshift3
      real(wp) :: hshift4
      real(wp) :: hshift5

      !> group 1+2 metals
      real(wp) :: metal1_shift

      !> TM
      real(wp) :: metal2_shift

      !> main group metals
      real(wp) :: metal3_shift

      !> eta bonded
      real(wp) :: eta_shift

      !> Charge Param
      real(wp) :: qfacbm(0:4)

      !> bond charge dependent
      real(wp) :: qfacbm0

      !> topo dist scaling
      real(wp) :: rfgoed1

      !> Hückel Param
      !  decrease Hueckel off-diag for triple bonds because they are less well conjugated 1.4
      real(wp) :: htriple

      !> increase pot depth depending on P
      real(wp) :: hueckelp2

      !> diagonal element change with qa
      real(wp) :: hueckelp3

      !> diagonal element relative to C
      real(wp) :: hdiag(17)

      !> Huckel off-diag constants
      real(wp) :: hoffdiag(17)

      !> iteration mixing
      real(wp) :: hiter

      !> diagonal qa dep.
      real(wp) :: hueckelp

      !> ref P value R shift
      real(wp) :: bzref

      !> ref P value k stretch
      real(wp) :: bzref2

      !> 2el diag shift
      real(wp) :: pilpf

      !> the Hückel iterations can diverge so take only a few steps
      real(wp) :: maxhiter

      !> D3 Param
      real(wp) :: d3a1

      !> D3
      real(wp) :: d3a2

      !> mixing of sp^n with sp^n-1
      real(wp) :: split0

      !> mixing of sp^n with sp^n-1
      real(wp) :: split1

      !> str ring size dep.
      real(wp) :: fringbo

      !> three coord. heavy eq. angle
      real(wp) :: aheavy3

      !> four coord. heavy eq. angle
      real(wp) :: aheavy4
      real(wp) :: bsmat(0:3,0:3)

      !> max CN cut-off
      real(wp) :: cnmax

   end type TGFFGenerator


!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine zero(self)
    class(TGFFTopology),intent(inout) :: self

    self%nbond = 0
    self%nangl = 0
    self%ntors = 0
    self%nathbH = 0
    self%nathbAB = 0
    self%natxbAB = 0
    self%nbatm = 0
    self%nfrag = 0
    self%maxsystem = 0
    self%bond_hb_nr = 0
    self%b_max = 0

    self%nbond_blist = 0
    self%nbond_vbond = 0
    self%nangl_alloc = 0
    self%ntors_alloc = 0

    self%read_file_type = 0
  end subroutine zero
!========================================================================================!

  !> Initialize new instance for the neighbourlist
  subroutine newGFFNeighbourList(self,n,nhb1,nhb2,nxb)
    type(TGFFNeighbourList),intent(out) :: self
    integer,intent(in) :: n
    integer,intent(in) :: nhb1
    integer,intent(in) :: nhb2
    integer,intent(in) :: nxb
    self%initialized = .true.
    self%nhb1 = nhb1
    self%nhb2 = nhb2
    self%nxb = nxb
    allocate (self%q(n),source=0.0_wp)
    allocate (self%hbrefgeo(3,n),source=0.0_wp)
    allocate (self%hblist1(3,self%nhb1),source=0)
    allocate (self%hblist2(3,self%nhb2),source=0)
    allocate (self%hblist3(3,self%nxb),source=0)
  end subroutine newGFFNeighbourList
!========================================================================================!

  !> Initialize a new instance of the parametrisation data
  subroutine initGFFData(self,ndim)

    !> Instance of the parametrisation data
    type(TGFFData),intent(out) :: self

    !> Dimension for allocating space
    integer,intent(in) :: ndim

    allocate (self%en(ndim))
    allocate (self%rad(ndim))
    allocate (self%metal(ndim))
    allocate (self%group(ndim))
    allocate (self%normcn(ndim))
    allocate (self%rcov(ndim))

    allocate (self%repa(ndim))
    allocate (self%repan(ndim))

    allocate (self%repz(ndim))
    allocate (self%zb3atm(ndim))

    allocate (self%xhaci(ndim))
    allocate (self%xhbas(ndim))
    allocate (self%xbaci(ndim))

    allocate (self%chi(ndim))
    allocate (self%gam(ndim))
    allocate (self%cnf(ndim))
    allocate (self%alp(ndim))

    allocate (self%bond(ndim))

    allocate (self%angl(ndim))

    allocate (self%angl2(ndim))

    allocate (self%tors(ndim))

    allocate (self%tors2(ndim))

    allocate (self%d3r0(ndim*(1+ndim)/2))

  end subroutine initGFFData

!========================================================================================!
  subroutine initGFFDispersion(self)
    type(TDispersionData),intent(out) :: self
    integer :: elem,ref,freq

    elem = 118
    ref = 7
    freq = 23

    allocate (self%atoms(elem),source=0)
    allocate (self%nref(elem),source=0)
    allocate (self%ncount(ref,elem),source=0)
    allocate (self%cn(ref,elem),source=0.0_wp)
    allocate (self%q(ref,elem),source=0.0_wp)
    allocate (self%alpha(freq,ref,elem),source=0.0_wp)
    allocate (self%c6(ref,ref,elem,elem),source=0.0_wp)

    self%g_a = 3.0_wp
    self%g_c = 2.0_wp
    self%wf = 6.0_wp

    self%a1 = 0.8000000_wp
    self%a2 = 4.6000000_wp
    self%s8 = 2.8500000_wp

  end subroutine initGFFDispersion

!========================================================================================!
end module gfnff_data_types
