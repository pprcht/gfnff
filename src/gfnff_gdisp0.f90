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
module gfnff_gdisp0
   use iso_fortran_env, only: wp=>real64,stdout=>output_unit
   use gfnff_data_types, only: TDispersionData,initgffdispersion
   use gfnff_helpers, only: lin
   implicit none
   public :: d3_gradient,newD3Model
   private


!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

subroutine newD3Model(dispm,nat,at)
   use disp_dftd3param
   use dftd4param
   implicit none
   type(TDispersionData), intent(out) :: dispm

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)

   integer  :: i,ia,is,icn,j,ja,ii,jj,ij
   integer  :: cncount(0:18)
   real(wp) :: alpha(23),c6

   intrinsic :: nint

!   call init(dispm, maxElem=maxval(at))
   call initGFFDispersion(dispm)

   dispm%atoms = 0
   dispm%nref = 0.0_wp

   !> set up ncount und alpha, also obtain the dimension of the dispmat
   do i = 1, nat
      ia = at(i)
      if (dispm%atoms(ia).eq.0) then
         dispm%nref(ia) = refn(ia)
         do j = 1, dispm%nref(ia)
            is = refsys(j,ia)
            alpha = sscale(is)*secaiw(:,is)
            dispm%cn(j,ia) = refcn(j,ia)
            dispm%alpha(:,j,ia) = max(ascale(j,ia)*(alphaiw(:,j,ia) &
               &                - hcount(j,ia)*alpha), 0.0_wp)
         enddo
      endif
      dispm%atoms(ia) = dispm%atoms(ia)+1
   enddo

   ! integrate C6 coefficients
   do i = 1, maxval(at)
      do j = 1, i
         if (dispm%atoms(i) > 0 .and. dispm%atoms(j) > 0) then
            do ii = 1, dispm%nref(i)
               do jj = 1, dispm%nref(j)
                  alpha = dispm%alpha(:,ii,i)*dispm%alpha(:,jj,j)
                  c6 = thopi * trapzd(alpha)
                  dispm%c6(ii,jj,i,j) = c6
                  dispm%c6(jj,ii,j,i) = c6
               enddo
            enddo
         endif
      enddo
   enddo

end subroutine newD3Model

!> Calculate the weights of the reference system and the derivatives w.r.t.
!  coordination number for later use.
subroutine weight_references(nat, atoms, wf, cn, gwvec, gwdcn)
   use disp_dftd3param
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> Atomic numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> weighting for the atomic reference systems
   real(wp), intent(out) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(out) :: gwdcn(:, :)

   integer :: iat, ati, iref, icount
   real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk

   gwvec = 0.0_wp
   gwdcn = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)
      norm = 0.0_wp
      dnorm = 0.0_wp
      do iref = 1, number_of_references(ati)
         gw = weight_cn(wf, cn(iat), reference_cn(iref, ati))
         norm = norm + gw
         dnorm = dnorm + 2*wf*(reference_cn(iref, ati) - cn(iat)) * gw
      end do
      norm = 1.0_wp / norm
      do iref = 1, number_of_references(ati)
         expw = weight_cn(wf, cn(iat), reference_cn(iref, ati))
         expd = 2*wf*(reference_cn(iref, ati) - cn(iat)) * expw

         gwk = expw * norm
         if (gwk /= gwk) then
            if (maxval(reference_cn(:number_of_references(ati), ati)) &
               & == reference_cn(iref, ati)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         gwvec(iref, iat) = gwk

         dgwk = expd*norm-expw*dnorm*norm**2
         if (dgwk /= dgwk) then
            dgwk = 0.0_wp
         endif
         gwdcn(iref, iat) = dgwk

      end do
   end do

end subroutine weight_references
!========================================================================================!

!> Calculate the weights of the reference system and the derivatives w.r.t.
!  coordination number for later use.
subroutine weight_references_d4(dispm, nat, atoms, wf, cn, gwvec, gwdcn)
   type(TDispersionData), intent(in) :: dispm
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> Atomic numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> weighting for the atomic reference systems
   real(wp), intent(out) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(out) :: gwdcn(:, :)

   integer :: iat, ati, iref, icount
   real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk

   gwvec = 0.0_wp
   gwdcn = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)
      norm = 0.0_wp
      dnorm = 0.0_wp
      do iref = 1, dispm%nref(ati)
         gw = weight_cn(wf, cn(iat), dispm%cn(iref, ati))
         norm = norm + gw
         dnorm = dnorm + 2*wf*(dispm%cn(iref, ati) - cn(iat)) * gw
      end do
      norm = 1.0_wp / norm
      do iref = 1, dispm%nref(ati)
         expw = weight_cn(wf, cn(iat), dispm%cn(iref, ati))
         expd = 2*wf*(dispm%cn(iref, ati) - cn(iat)) * expw

         gwk = expw * norm
         if (gwk /= gwk) then
            if (maxval(dispm%cn(:dispm%nref(ati), ati)) &
               & == dispm%cn(iref, ati)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         gwvec(iref, iat) = gwk

         dgwk = expd*norm-expw*dnorm*norm**2
         if (dgwk /= dgwk) then
            dgwk = 0.0_wp
         endif
         gwdcn(iref, iat) = dgwk

      end do
   end do

end subroutine weight_references_d4
!========================================================================================!

!> calculate atomic dispersion coefficients and their derivatives w.r.t.
!  the coordination number.
subroutine get_atomic_c6(nat, atoms, gwvec, gwdcn, c6, dc6dcn)
   use disp_dftd3param
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in) :: gwdcn(:, :)
   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)
   !> derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out) :: dc6dcn(:, :)

   integer :: iat, jat, ati, atj, iref, jref
   real(wp) :: refc6, dc6, dc6dcni, dc6dcnj

   c6 = 0.0_wp
   dc6dcn = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)
      do jat = 1, iat
         atj = atoms(jat)
         dc6 = 0.0_wp
         dc6dcni = 0.0_wp
         dc6dcnj = 0.0_wp
         do iref = 1, number_of_references(ati)
            do jref = 1, number_of_references(atj)
               refc6 = get_c6(iref, jref, ati, atj)
               dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcni = dc6dcni + gwdcn(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcnj = dc6dcnj + gwvec(iref, iat) * gwdcn(jref, jat) * refc6
            end do
         end do
         c6(iat, jat) = dc6
         c6(jat, iat) = dc6
         dc6dcn(iat, jat) = dc6dcni
         dc6dcn(jat, iat) = dc6dcnj
      end do
   end do
end subroutine get_atomic_c6
!========================================================================================!

!> calculate atomic dispersion coefficients and their derivatives w.r.t.
!  the coordination number.
subroutine get_atomic_c6_d4(dispm, nat, atoms, gwvec, gwdcn, c6, dc6dcn)
   use disp_dftd3param, only : get_c6
   type(TDispersionData), intent(in) :: dispm
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in) :: gwdcn(:, :)
   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)
   !> derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out) :: dc6dcn(:, :)

   integer :: iat, jat, ati, atj, iref, jref
   real(wp) :: refc6, dc6, dc6dcni, dc6dcnj

   c6 = 0.0_wp
   dc6dcn = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)
      do jat = 1, iat
         atj = atoms(jat)
         dc6 = 0.0_wp
         dc6dcni = 0.0_wp
         dc6dcnj = 0.0_wp
         do iref = 1, dispm%nref(ati)
            do jref = 1, dispm%nref(atj)
               refc6 = dispm%c6(iref, jref, ati, atj)
               dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcni = dc6dcni + gwdcn(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcnj = dc6dcnj + gwvec(iref, iat) * gwdcn(jref, jat) * refc6
            end do
         end do
         c6(iat, jat) = dc6
         c6(jat, iat) = dc6
         dc6dcn(iat, jat) = dc6dcni
         dc6dcn(jat, iat) = dc6dcnj
      end do
   end do
end subroutine get_atomic_c6_d4
!========================================================================================!

subroutine d3_gradient(dispm, nat, at, xyz, npair, pairlist, zeta_scale, radii, &
      & r4r2, weighting_factor, dispscale, cn, dcndr, energy, gradient)
   use disp_dftd3param
   type(TDispersionData), intent(in) :: dispm
   integer, intent(in) :: nat
   integer, intent(in) :: at(:)
   real(wp), intent(in) :: xyz(:, :)
   integer, intent(in) :: npair
   integer, intent(in) :: pairlist(:, :)
   real(wp), intent(in) :: zeta_scale(:)
   real(wp), intent(in) :: r4r2(:)
   real(wp), intent(in) :: radii(:)

   real(wp), intent(in) :: dispscale
   real(wp), intent(in) :: weighting_factor
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)

   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: gradient(:, :)
   real(wp) :: sigma(3, 3)

   integer :: max_ref
   integer :: iat, jat, ati, atj, ij, img

   real(wp) :: r4r2ij, r0, rij(3), r2, t6, t8, t10, d6, d8, d10
   real(wp) :: dE, dG(3), dS(3, 3), disp, ddisp

   real(wp), allocatable :: gw(:, :), dgwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: energies(:), dEdcn(:)

   max_ref = maxval(dispm%nref(at))
   allocate(gw(max_ref, nat), dgwdcn(max_ref, nat), c6(nat, nat), &
      &     dc6dcn(nat, nat), energies(nat), dEdcn(nat), source=0.0_wp)

   call weight_references_d4(dispm, nat, at, weighting_factor, cn, gw, dgwdcn)

   !gw = gw*zeta_scale
   !dgwdcn = dgwdcn*zeta_scale
   call get_atomic_c6_d4(dispm, nat, at, gw, dgwdcn, c6, dc6dcn)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:energies, gradient, sigma, dEdcn) &
   !$omp shared(at, xyz, npair, pairlist, zeta_scale, dispscale, radii, r4r2, c6, &
   !$omp&       dc6dcn) &
   !$omp private(ij, img, iat, jat, ati, atj, r2, rij, r4r2ij, r0, t6, t8, t10, &
   !$omp&        d6, d8, d10, disp, ddisp, dE, dG, dS)
   do img = 1, npair
      iat = pairlist(1, img)
      jat = pairlist(2, img)
      ij = jat + iat*(iat-1)/2
      ati = at(iat)
      atj = at(jat)
      rij = xyz(:, iat) - xyz(:, jat)
      r2 = sum(rij**2)

      r4r2ij = 3*r4r2(ati)*r4r2(atj)
      r0 = radii(lin(ati, atj))

      t6 = 1._wp/(r2**3+r0**3)
      t8 = 1._wp/(r2**4+r0**4)

      d6 = -6*r2**2*t6**2
      d8 = -8*r2**3*t8**2

      disp = (t6 + 2*r4r2ij*t8) * zeta_scale(ij)*dispscale
      ddisp= (d6 + 2*r4r2ij*d8) * zeta_scale(ij)*dispscale

      dE = -c6(iat, jat)*disp * 0.5_wp
      dG = -c6(iat, jat)*ddisp*rij
      dS = spread(dG, 1, 3) * spread(rij, 2, 3) * 0.5_wp

      energies(iat) = energies(iat) + dE
      dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * disp
      sigma = sigma + dS
      if (iat /= jat) then
         energies(jat) = energies(jat) + dE
         dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * disp
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
         sigma = sigma + dS
      endif

   enddo
   !$omp end parallel do

   call dgemv('n', 3*nat, nat, 1.0_wp, dcndr, 3*nat, dEdcn, 1, 1.0_wp, gradient, 1)
   !call dgemv('n', 9, nat, 1.0_wp, dcndL, 9, dEdcn, 1, 1.0_wp, sigma, 1)

   energy = sum(energies)

end subroutine d3_gradient

real(wp) pure elemental function weight_cn(wf,cn,cnref) result(cngw)
   real(wp),intent(in) :: wf, cn, cnref
   intrinsic :: exp
   cngw = exp ( -wf * ( cn - cnref )**2 )
end function weight_cn

real(wp) pure elemental function pair_scale(iat, jat) result(scale)
   integer, intent(in) :: iat, jat
   if (iat == jat) then
      scale = 0.5_wp
   else
      scale = 1.0_wp
   endif
end function pair_scale

pure function trapzd(pol)
   real(wp),intent(in) :: pol(23)
   real(wp)            :: trapzd

   real(wp)            :: tmp1, tmp2
   real(wp),parameter  :: freq(23) = (/ &
&   0.000001_wp,0.050000_wp,0.100000_wp, &
&   0.200000_wp,0.300000_wp,0.400000_wp, &
&   0.500000_wp,0.600000_wp,0.700000_wp, &
&   0.800000_wp,0.900000_wp,1.000000_wp, &
&   1.200000_wp,1.400000_wp,1.600000_wp, &
&   1.800000_wp,2.000000_wp,2.500000_wp, &
&   3.000000_wp,4.000000_wp,5.000000_wp, &
&   7.500000_wp,10.00000_wp /)
!  just precalculate all weights and get the job done
   real(wp),parameter :: weights(23) = 0.5_wp * (/ &
&  ( freq (2) - freq (1) ),  &
&  ( freq (2) - freq (1) ) + ( freq (3) - freq (2) ),  &
&  ( freq (3) - freq (2) ) + ( freq (4) - freq (3) ),  &
&  ( freq (4) - freq (3) ) + ( freq (5) - freq (4) ),  &
&  ( freq (5) - freq (4) ) + ( freq (6) - freq (5) ),  &
&  ( freq (6) - freq (5) ) + ( freq (7) - freq (6) ),  &
&  ( freq (7) - freq (6) ) + ( freq (8) - freq (7) ),  &
&  ( freq (8) - freq (7) ) + ( freq (9) - freq (8) ),  &
&  ( freq (9) - freq (8) ) + ( freq(10) - freq (9) ),  &
&  ( freq(10) - freq (9) ) + ( freq(11) - freq(10) ),  &
&  ( freq(11) - freq(10) ) + ( freq(12) - freq(11) ),  &
&  ( freq(12) - freq(11) ) + ( freq(13) - freq(12) ),  &
&  ( freq(13) - freq(12) ) + ( freq(14) - freq(13) ),  &
&  ( freq(14) - freq(13) ) + ( freq(15) - freq(14) ),  &
&  ( freq(15) - freq(14) ) + ( freq(16) - freq(15) ),  &
&  ( freq(16) - freq(15) ) + ( freq(17) - freq(16) ),  &
&  ( freq(17) - freq(16) ) + ( freq(18) - freq(17) ),  &
&  ( freq(18) - freq(17) ) + ( freq(19) - freq(18) ),  &
&  ( freq(19) - freq(18) ) + ( freq(20) - freq(19) ),  &
&  ( freq(20) - freq(19) ) + ( freq(21) - freq(20) ),  &
&  ( freq(21) - freq(20) ) + ( freq(22) - freq(21) ),  &
&  ( freq(22) - freq(21) ) + ( freq(23) - freq(22) ),  &
&  ( freq(23) - freq(22) ) /)

   trapzd = sum(pol*weights)

end function trapzd



!========================================================================================!
end module gfnff_gdisp0
