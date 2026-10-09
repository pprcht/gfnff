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

!> Implementation of the Born radii integrator
module gfnff_solv_born
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  public :: compute_bornr
  public :: split_range

  !> van der Waals to Lee-Richard's surface correction (GBOBCII parameter)
  real(wp),parameter :: alp = 1._wp
  real(wp),parameter :: bet = 0.8_wp
  real(wp),parameter :: gam = 4.85_wp

contains

  subroutine compute_bornr(nat,nnrad,nnlistr,ddpair,vdwr,rho,svdw,c1, &
        & brad,brdr,psiout)
    !***********************************************************************
    !* Born radii and their Cartesian derivatives from the pairwise
    !* descreening integral with GBOBCII rescaling. Inputs: rho descreened
    !* vdW radii, svdw vdW radii with offset, c1 GBMV2-like radius scaling.
    !* psiout optionally returns the integral before brad overwrites it;
    !* the closed-form Hessian needs it.
    !***********************************************************************

    integer,intent(in) :: nat

    integer,intent(in) :: nnrad

    integer,intent(in) :: nnlistr(:,:)

    real(wp),intent(in) :: ddpair(:,:)

    real(wp),intent(in) :: vdwr(:)

    real(wp),intent(in) :: rho(:)

    real(wp),intent(in) :: svdw(:)

    real(wp),intent(in) :: c1

    real(wp),intent(out) :: brad(:)

    real(wp),intent(out) :: brdr(:,:,:)

    real(wp),intent(out),optional :: psiout(:)

    integer :: iat
    real(wp) :: br,dpsi,svdwi,vdwri,s1,v1,s2,arg,arg2
    real(wp) :: th,ch

    call compute_psi(nat,nnrad,nnlistr,ddpair,vdwr,rho,brad,brdr)
    if (present(psiout)) psiout(1:nat) = brad(1:nat)

    !>-- rescale psi to a Born radius, in parallel since brdr is 3*nat*nat
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(nat,brad,brdr,svdw,vdwr,c1) &
    !$omp private(iat,br,dpsi,svdwi,vdwri,s1,v1,s2,arg,arg2,th,ch)
    do iat = 1,nat

      br = brad(iat)

      svdwi = svdw(iat)
      vdwri = vdwr(iat)
      s1 = 1.0_wp/svdwi
      v1 = 1.0_wp/vdwri
      s2 = 0.5_wp*svdwi

      br = br*s2

      arg2 = br*(gam*br-bet)
      arg = br*(alp+arg2)
      arg2 = 2.0_wp*arg2+alp+gam*br*br

      th = tanh(arg)
      ch = cosh(arg)

      br = 1.0_wp/(s1-v1*th)
      br = c1*br

      dpsi = ch*(s1-v1*th)
      dpsi = s2*v1*arg2/(dpsi*dpsi)
      dpsi = c1*dpsi

      brad(iat) = br
      brdr(:,:,iat) = brdr(:,:,iat)*dpsi

    end do
    !$omp end parallel do

  end subroutine compute_bornr

  subroutine compute_psi(nat,nnrad,nnlistr,ddpair,vdwr,rho,psi,dpsidr)
    !***********************************************************************
    !* Descreening integral psi(i) and its derivatives, d psi(i)/d x_j in
    !* dpsidr(:,j,i). Own-row parallel loop: the thread owning atom a writes
    !* only row a and sums its pairs in serial order, so no atomics or thread
    !* copies are needed and the result matches the serial loop bit for bit.
    !***********************************************************************
!$  use omp_lib,only:omp_get_max_threads

    integer,intent(in) :: nat

    integer,intent(in) :: nnrad

    integer,intent(in) :: nnlistr(:,:)

    real(wp),intent(in) :: ddpair(:,:)

    real(wp),intent(in) :: vdwr(:)

    real(wp),intent(in) :: rho(:)

    real(wp),intent(out) :: psi(:)

    real(wp),intent(out) :: dpsidr(:,:,:)

    integer  :: kk,idx,a,t,nproc,lo,hi,run
    integer  :: ii,jj,nn
    real(wp) :: dr(3),r,rhoi,rhoj
    real(wp) :: gi,gj,dgi,dgj
    real(wp) :: drjj(3)
    real(wp) :: rvdwi,rvdwj
    logical  :: havei,havej
    integer,allocatable :: pstart(:),pcount(:),plist(:)
    integer,allocatable :: cnt(:,:),off(:,:)
    integer :: ppos(nat)

    !>-- index the pair list by atom, preserving pair order: per-thread counts
    !>   over contiguous ranges plus a prefix sum (a serial sort would dominate)
    nproc = 1
!$  nproc = omp_get_max_threads()
    nproc = max(1,min(nproc,nnrad))
    allocate (pstart(nat),pcount(nat),plist(2*nnrad))
    allocate (cnt(nat,nproc),off(nat,nproc))

    !$omp parallel do default(none) schedule(static) &
    !$omp shared(nproc,nnrad,nnlistr,cnt,nat) private(t,lo,hi,kk)
    do t = 1,nproc
      call split_range(nnrad,nproc,t,lo,hi)
      cnt(:,t) = 0
      do kk = lo,hi
        cnt(nnlistr(1,kk),t) = cnt(nnlistr(1,kk),t)+1
        cnt(nnlistr(2,kk),t) = cnt(nnlistr(2,kk),t)+1
      end do
    end do
    !$omp end parallel do

    run = 1
    do a = 1,nat
      pstart(a) = run
      do t = 1,nproc
        off(a,t) = run
        run = run+cnt(a,t)
      end do
      pcount(a) = run-pstart(a)
    end do

    !$omp parallel do default(none) schedule(static) &
    !$omp shared(nproc,nnrad,nnlistr,off,plist,nat) &
    !$omp private(t,lo,hi,kk,ii,jj,ppos)
    do t = 1,nproc
      call split_range(nnrad,nproc,t,lo,hi)
      ppos = off(:,t)
      do kk = lo,hi
        ii = nnlistr(1,kk)
        plist(ppos(ii)) = kk
        ppos(ii) = ppos(ii)+1
        jj = nnlistr(2,kk)
        plist(ppos(jj)) = kk
        ppos(jj) = ppos(jj)+1
      end do
    end do
    !$omp end parallel do
    deallocate (cnt,off)

    !$omp parallel do default(none) schedule(dynamic,8) &
    !$omp shared(nat,nnlistr,ddpair,vdwr,rho,psi,dpsidr,pstart,pcount,plist) &
    !$omp private(a,idx,kk,ii,jj,nn,dr,r,rhoi,rhoj,rvdwi,rvdwj, &
    !$omp&        gi,dgi,gj,dgj,havei,havej,drjj)
    do a = 1,nat

      psi(a) = 0.0_wp
      dpsidr(:,:,a) = 0.0_wp

      do idx = pstart(a),pstart(a)+pcount(a)-1
        kk = plist(idx)

        ii = nnlistr(1,kk)
        jj = nnlistr(2,kk)
        nn = nnlistr(3,kk)

        r = ddpair(1,nn)
        dr(:) = ddpair(2:4,nn)

        rhoi = rho(ii)
        rhoj = rho(jj)
        rvdwi = vdwr(ii)
        rvdwj = vdwr(jj)

        call psi_pair(r,rhoi,rhoj,rvdwi,rvdwj,a == ii,a == jj, &
           & gi,dgi,havei,gj,dgj,havej)

        if (a == ii) then
          if (havei) then
            psi(ii) = psi(ii)+gi
            drjj(:) = dgi*dr(:)
            dpsidr(:,ii,ii) = dpsidr(:,ii,ii)+drjj(:)
            dpsidr(:,jj,ii) = dpsidr(:,jj,ii)-drjj(:)
          end if
        else
          if (havej) then
            psi(jj) = psi(jj)+gj
            drjj(:) = dgj*dr(:)
            dpsidr(:,jj,jj) = dpsidr(:,jj,jj)-drjj(:)
            dpsidr(:,ii,jj) = dpsidr(:,ii,jj)+drjj(:)
          end if
        end if

      end do
    end do
    !$omp end parallel do

  end subroutine compute_psi

  pure subroutine split_range(n,nchunk,ichunk,lo,hi)
    !***********************************************************************
    !* Bounds of the ichunk-th of nchunk contiguous, near-equal slices of
    !* 1..n, in order, so both counting-sort passes agree on pair ownership.
    !***********************************************************************
    integer,intent(in) :: n,nchunk,ichunk
    integer,intent(out) :: lo,hi
    integer :: base,rest
    base = n/nchunk
    rest = n-base*nchunk
    lo = (ichunk-1)*base+min(ichunk-1,rest)+1
    hi = lo+base-1
    if (ichunk .le. rest) hi = hi+1
  end subroutine split_range

  pure subroutine psi_pair(r,rhoi,rhoj,rvdwi,rvdwj,wanti,wantj, &
        & gi,dgi,havei,gj,dgj,havej)
    !***********************************************************************
    !* Descreening contribution of one pair to psi(i) and psi(j): value g,
    !* d/dr dg, and a have flag each; wanti/wantj select what is evaluated.
    !* Sphere overlap is tested per direction. Equal reduced radii reuse one
    !* evaluation from rhoj for both atoms as in the original code; that
    !* test has a tolerance, so the shortcut must be reproduced exactly.
    !***********************************************************************
    real(wp),intent(in) :: r,rhoi,rhoj,rvdwi,rvdwj
    logical,intent(in) :: wanti,wantj
    real(wp),intent(out) :: gi,dgi,gj,dgj
    logical,intent(out) :: havei,havej

    real(wp) :: ap,am,lnab,rhab,ab
    real(wp) :: rh1,rhr1,r24,r1,aprh1,r12
    logical  :: ovij,ovji

    gi = 0.0_wp; dgi = 0.0_wp; havei = .false.
    gj = 0.0_wp; dgj = 0.0_wp; havej = .false.

    ovij = r .lt. (rvdwi+rhoj)
    ovji = r .lt. (rhoi+rvdwj)
    r1 = 1.0_wp/r

    !>-- neither sphere overlaps: the closed descreening integral applies
    if (.not.ovij.and..not.ovji) then

      if (abs(rhoi-rhoj) .lt. 1.d-8) then
        ap = r+rhoj
        am = r-rhoj
        ab = ap*am
        rhab = rhoj/ab
        lnab = 0.5_wp*log(am/ap)*r1
        gi = rhab+lnab
        dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
        gj = gi
        dgj = dgi
        havei = .true.
        havej = .true.
        return
      end if

      if (wanti) then
        ap = r+rhoj
        am = r-rhoj
        ab = ap*am
        rhab = rhoj/ab
        lnab = 0.5_wp*log(am/ap)*r1
        gi = rhab+lnab
        dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
        havei = .true.
      end if
      if (wantj) then
        ap = r+rhoi
        am = r-rhoi
        ab = ap*am
        rhab = rhoi/ab
        lnab = 0.5_wp*log(am/ap)*r1
        gj = rhab+lnab
        dgj = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
        havej = .true.
      end if
      return
    end if

    if (wanti) then
      if (.not.ovij) then
        ap = r+rhoj
        am = r-rhoj
        ab = ap*am
        rhab = rhoj/ab
        lnab = 0.5_wp*log(am/ap)*r1
        gi = rhab+lnab
        dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
        havei = .true.
      else if ((r+rhoj) .gt. rvdwi) then
        r12 = 0.5_wp*r1
        r24 = r12*r12

        ap = r+rhoj
        am = r-rhoj
        rh1 = 1.0_wp/rvdwi
        rhr1 = 1.0_wp/ap
        aprh1 = ap*rh1
        lnab = log(aprh1)

        gi = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

        dgi = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
           &         rhoj*r24*(rhr1-rh1*aprh1)+ &
           &         r12*(r1*lnab-rhr1)
        dgi = dgi*r1
        havei = .true.
      end if
    end if

    if (wantj) then
      if (.not.ovji) then
        ap = r+rhoi
        am = r-rhoi
        ab = ap*am
        rhab = rhoi/ab
        lnab = 0.5_wp*log(am/ap)*r1
        gj = rhab+lnab
        dgj = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
        havej = .true.
      else if ((r+rhoi) .gt. rvdwj) then
        r12 = 0.5_wp*r1
        r24 = r12*r12

        ap = r+rhoi
        am = r-rhoi
        rh1 = 1.0_wp/rvdwj
        rhr1 = 1.0_wp/ap
        aprh1 = ap*rh1
        lnab = log(aprh1)

        gj = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

        dgj = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
           &         rhoi*r24*(rhr1-rh1*aprh1)+ &
           &         r12*(r1*lnab-rhr1)
        dgj = dgj*r1
        havej = .true.
      end if
    end if

  end subroutine psi_pair

end module gfnff_solv_born
