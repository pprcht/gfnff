!> SiO2 (alpha-quartz) unit cell.
!>
!> Coordinates from $coord block (Bohr), periodic 3D cell with hexagonal lattice:
!>   a = b = 9.28422449595511046 Bohr
!>   c = 10.21434769907115 Bohr
!>   alpha = beta = 90 deg, gamma = 120 deg
module sio2
  use iso_fortran_env, only: wp => real64
  implicit none
  private

  public :: sio2nat, sio2at, sio2xyz, sio2lattice

  integer, parameter :: sio2nat = 9

  integer, parameter :: sio2at(sio2nat) = [8,8,8,8,8,8,14,14,14]

  real(wp), parameter :: sio2xyz(3, sio2nat) = reshape([ &
  &   2.82781861325240_wp,  2.96439280874170_wp,  3.12827803849279_wp, &
  &   7.19124230791576_wp,  0.98723342603994_wp,  4.89004701836746_wp, &
  &   4.95491880597601_wp,  4.82830910314898_wp,  8.74847811174740_wp, &
  &   0.19290883043307_wp,  2.30645007856310_wp,  8.72969832061507_wp, &
  &  -2.01592208020090_wp,  6.16478744235115_wp,  4.87273962147340_wp, &
  &   0.66183062221384_wp,  7.07392578563696_wp,  0.27767968372345_wp, &
  &   4.55701736204879_wp,  0.06291337111965_wp,  3.31745840478609_wp, &
  &  -2.10064209975148_wp,  3.63969476409878_wp,  6.81014625000326_wp, &
  &   2.31009832827224_wp,  4.12572862149043_wp,  0.08842485276656_wp  &
  & ], shape(sio2xyz))

  !> Lattice matrix (column-major: columns are lattice vectors).
  !> a1 = [a,0,0], a2 = [a*cos(120),a*sin(120),0], a3 = [0,0,c]
  real(wp), parameter :: sio2_a   = 9.28422449595511046_wp
  real(wp), parameter :: sio2_c   = 10.21434769907115_wp
  real(wp), parameter :: cos120   = -0.5_wp
  real(wp), parameter :: sin120   = 0.86602540378443865_wp  !> sqrt(3)/2

  real(wp), parameter :: sio2lattice(3,3) = reshape([ &
  &  sio2_a,           0.0_wp,  0.0_wp, &  ! a1
  &  sio2_a * cos120,  sio2_a * sin120,  0.0_wp, &  ! a2
  &  0.0_wp,           0.0_wp,  sio2_c  &  ! a3
  & ], shape(sio2lattice))

end module sio2
