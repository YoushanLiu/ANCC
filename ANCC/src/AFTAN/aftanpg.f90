! This file is part of ANCC.
!
! AFTAN is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! AFTAN is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
!
!======================================================================
! aftanpg function. Provides regular ftan analysis, jumps correction,
! and amplitude map output for input periods.
!
!
! Autor: M.Barmine, CIEI, CU. Date: Jun 15,2006. Version: 2.0
!
module aftanpg_m

implicit none


contains


subroutine aftanpg(PIover4, n, seis, t0, dt, delta, vmin, vmax, tmin, tmax, &
                 tresh, ffact, perc, npoints, taperl, nfin, nphpr, phprper, &
                 phprvel, nfout1, arr1, nfout2, arr2, tamp, nrow, ncol, amp, ierr)

!======================================================================
! Parameters for aftanpg function:
! Input parameters:
! PIover4   - phase shift = PI/4*PIover4, for cross-correlation
!             PIover4 should be -1.0 (real(8))
! n         - number of input samples, (integer(4))
! seis      - input waveform length of n, (real(8))
! t0        - time shift of SAC file in seconds, (real(8))
! dt        - sampling step, in seconds, (real(8))
! delta     - distance, km (real(8))
! vmin      - minimal group velocity, km/s (real(8))
! vmax      - maximal value of the group velocity, km/s (real(8))
! tmin      - minimal period, s (real(8))
! tmax      - maximal period, s (real(8))
! tresh     - treshold, usualy = 10, (real(8))
! ffact     - factor to automatic filter parameter, usualy =1, (real(8))
! perc      - minimal length of of output segment vs freq. range, % (real(8))
! npoints   - max number points in jump, (integer(4))
! taperl    - factor for the left end seismogram tapering,
!             taper = taperl*tmax,    (real(8))
! nfin      - starting number of frequencies, nfin <= 100,(integer(4))
! nphpr     - length of phprper and phprvel arrays
! phprper   - predicted phase velocity periods, s
! phprvel   - predicted phase velocity for corresponding periods, s
! NOTE.       If nphpr = 0, aftanpg subroutine does not output phase
!                velocity branch, instead of it outputs phase.
! ==========================================================
! Output parameters are placed in 2-D arrays arr1 and arr2,
! arr1 contains preliminary results and arr2 - final.
! ==========================================================
! nfout1    - output number of frequencies for arr1, (integer(4))
! arr1      - preliminary results.
!             Description: real(8) arr1(8,n), n >= nfin)
! arr1(1,:) - central periods, s
! arr1(2,:) - observed periods, s
! arr1(3,:) - group velocities, km/s or phase if nphpr=0, rad
! arr1(4,:) - phase velocities, km/s
! arr1(5,:) - amplitudes, Db
! arr1(6,:) - discrimination function
! arr1(7,:) - signal/noise ratio, Db
! arr1(8,:) - maximum half width, s
! arr2      - final results
! nfout2    - output number of frequencies for arr2, (integer(4))
!             Description: real(8) arr2(7,n), n >= nfin)
!             If nfout2 == 0, no final result.
! arr2(1,:) - central periods, s
! arr2(2,:) - observed periods, s
! arr2(3,:) - group velocities, km/sor phase if nphpr=0, rad
! arr2(4,:) - phase velocities, km/s
! arr2(5,:) - amplitudes, Db
! arr2(6,:) - signal/noise ratio, Db
! arr2(7,:) - maximum half width, s
! tamp      - time to the beginning of ampo table, s (real(8))
! nrow      - number of rows in array ampo, (integer(4))
! ncol      - number of columns in array ampo, (integer(4))
! amp       - FTAN amplitude array, Db, (real(8))
! ierr      - completion status, =0 - O.K.,           (integer(4))
!                                =1 - some problems occures
!                                =2 - no final results
!======================================================================

use taper_m

implicit none

include 'fftw3.f'

integer(4), intent(in) :: n, npoints
integer(4), intent(in) :: nfin, nphpr

real(8), intent(in) :: tresh, ffact, perc
real(8), intent(in) :: t0, dt, delta, taperl
real(8), intent(in) :: vmin, vmax, tmin, tmax

real(8), intent(in) :: PIover4

real(8), intent(in) :: seis(n)

real(8), intent(in) :: phprper(nphpr), phprvel(nphpr)

integer(4), intent(out) :: nrow, ncol
integer(4), intent(out) :: nfout1, nfout2, ierr

real(8), intent(out) :: tamp

real(8), dimension(:,:), allocatable, intent(out) :: amp, arr1, arr2


integer(4) i, j, k, m, nf, nq
integer(4) ici, icj, iciflag, ia, ier
integer(4) mm, mi, iml, imr, indl, indr
integer(4) ip, iflag, nindx, imax, iimax
integer(4) ki, kk, istrt, ibeg, iend, ima
integer(4) ipos, ist, ibe, ne, nb, ns, ntall
integer(4) njump, nijmp, nii, ntapb, ntape, ntime

integer(8) planf, planb

real(8) t, dom, omb, ome, step, dph, tm, ph
real(8) time, tim, dmaxt, wor, lm, rm, alpha, amax

real(8), parameter :: PI = 4.d0*datan(1.d0)

complex(8), parameter :: czero = dcmplx(0.d0, 0.d0)

integer(4), dimension(:), allocatable :: ii, ijmp, indx

integer(4), dimension(:,:), allocatable :: ind

real(8), dimension(:), allocatable :: phgrc, trig1
real(8), dimension(:), allocatable :: om, per, om1, per1, ftrig, ftrig1
real(8), dimension(:), allocatable :: grvel, tvis, ampgr, phgr, snr, wdth
real(8), dimension(:), allocatable :: grvel1, tvis1, ampgr1, phgr1, snr1, wdth1
real(8), dimension(:), allocatable :: grvelt, tvist, ampgrt, phgrt, snrt, wdtht

real(8), dimension(:,:), allocatable :: ipar, pha, ampo

complex(8), dimension(:), allocatable :: s, sf, fils, tmp



ierr = 0


ip = 1
iml = 0
imr = 0
lm = 0.d0
rm = 0.d0


! number of FTAN filters
nf = nfin

! automatic width of filters * factor ffact
alpha = ffact*20.d0*sqrt(delta*1.d-3)

! number of samples for tapering, left end
ntapb = nint(taperl*tmax/dt)

! number of samples for tapering, right end
ntape = nint(tmax/dt)

! [omb,ome] - frequency range
omb = 2.d0*PI/tmax
ome = 2.d0*PI/tmin

! seismgram tapering
nb = max(2, nint((delta/vmax-t0)/dt))
tamp = (nb-1)*dt+t0
ne = min(n, nint((delta/vmin-t0)/dt))



nrow = nfin
ncol = ne-nb+1
!! times for FTAN map
!allocate(v(nb:ne), stat=ier)
!do k = nb, ne, 1
!   time = (k-1)*dt
!   ! velocity for FTAN map
!   v(k) = delta/(t0 + time)
!enddo



ntime = ne-nb+1
! tapering both ends of seismogram
call taper(max(nb,ntapb+1), min(ne,n-ntape), n, seis, ntapb, ntape, ns, s)

! prepare FTAN filters
dom = 2.d0*PI/(ns*dt)
step =(log(omb) - log(ome))/(nf - 1)

! log scaling for frequency
allocate(om(1:nf), per(1:nf), stat=ier)
do k = 1, nf, 1
   om(k) = exp(log(ome) + (k-1)*step)
   per(k) = 2.d0*PI/om(k)
enddo

! make backward FFT for seismogram: s ==> sf
allocate(sf(1:ns), stat=ier)
call dfftw_plan_dft_1d(planf, ns, s, sf, FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute_dft(planf, s, sf)
call dfftw_destroy_plan(planf)

! filtering and FTAN amplitude diagram construction
allocate(fils(1:ns), tmp(1:ns), stat=ier)
call dfftw_plan_dft_1d(planb, ns, fils, tmp, FFTW_BACKWARD, FFTW_ESTIMATE)


! main loop by frequency
ntall = ntime + 2
allocate(pha(1:ntall,1:nf), amp(1:ntall,1:nf), ampo(1:ntall,1:nf), stat=ier)
nq = ns/2 + 1
do k = 1, nf, 1

   ! filtering
   call ftfilt(alpha, om(k), dom, ns, sf, fils)

   ! fill with zeros half spectra for Hilbert transformation and
   ! spectra ends ajustment
   fils(nq+1:ns) = czero

   fils(1) = cmplx(real(0.5d0*fils(1)), 0.d0)
   fils(nq) = 0.50*fils(nq)
   ! forward FFT: fils ==> tmp
   call dfftw_execute_dft(planb, fils, tmp)

   tmp(1:ns) = tmp(1:ns) / real(ns)

   ! extraction from FTAN map area of investigation
   j = 1
   do m = nb-1, ne+1, 1
      pha(j,k) = atan2(imag(tmp(m)), real(tmp(m)))
      wor = abs(tmp(m))
      ampo(j,k) = wor
      amp(j,k) = 20.d0*log10(wor)
      j = j + 1
   enddo

enddo
call dfftw_destroy_plan(planb)
deallocate(s, sf, fils, tmp)


! normalization amp diagram to 100 Db with three decade cutting
amax = maxval(amp(1:ntall,1:nf))
amp(1:ntall,1:nf) = max(amp(1:ntall,1:nf) - amax + 100.d0, 40.d0)


! construction reference indices table ind. It points to local maxima.
! table ipar contains three parameter for each local maximum:
! tim - group time tvis - observed period ampgr - amplitude values in Db
nfout1 = 0
allocate(tvis(1:nf), ampgr(1:nf), stat=ier)
allocate(grvel(1:nf), snr(1:nf), wdth(1:nf), phgr(1:nf), stat=ier)
! tmp arrays
allocate(ind(1:2,1:ntall*nf), stat=ier)
allocate(ipar(1:6,1:ntall*nf), stat=ier)
ici = 0
do k = 1, nf, 1

   ! find local maxima on the FTAN amplitude diagram map
   iciflag = 0
   do j = 2, ntall-1, 1

      if ((amp(j,k) > amp(j-1,k)) .and. (amp(j,k) > amp(j+1,k))) then
         iciflag = iciflag + 1
         ici = ici + 1
         ind(1,ici) = k
         ind(2,ici) = j
      endif

   enddo

   if (0 == iciflag) then
      ici = ici + 1
      ind(1,ici) = k
      ind(2,ici) = ntall-1
      if (ampo(2,k) > ampo(ntall-1,k)) then
         ind(2,ici) = 2
      end if
      iciflag = 1
   endif

   ! compute parameters for each maximum
   ia = 1
   icj = ici-iciflag+1
   amax = -huge(amax)
   do j = icj, ici, 1

      m = ind(2,j)

      call fmax(amp(m-1,k), amp(m,k), amp(m+1,k), pha(m-1,k), pha(m,k), &
                          pha(m+1,k), PIover4, om(k), dt, t, dph, tm, ph)

      ipar(1,j) = (nb+m-3+t)*dt
      ipar(2,j) = 2*PI*dt/dph
      ipar(3,j) = tm
      ipar(6,j) = ph

      if (tm > amax) then
         amax = tm
         ia = j
      endif

   enddo

   ! Compute signal to noise ratio (SNR) ------------
   mm = 0
   do j = icj, ici, 1

      mm = mm + 1

      m = ind(2,j)

      ! find boundaries around local maximum
      if (1 == mm) then
         iml = 1
         imr = ntall
         if (iciflag > 1) then
            imr = ind(2,j+1)
         end if
      elseif (iciflag == mm) then
         iml = 1
         imr = ntall
         if (iciflag > 1) then
            iml = ind(2,j-1)
         end if
      else
         iml = ind(2,j-1)
         imr = ind(2,j+1)
      endif

      ! compute left minimum -------
      indl = 1
      lm = ampo(m,k)
      do mi = iml, m, 1
         if (ampo(mi,k) <= lm) then
            lm = ampo(mi,k)
            indl = mi
         endif
      enddo

      ! compute right minimum -------
      indr = 1
      rm = ampo(m,k)
      do mi = m, imr, 1
         if (ampo(mi,k) <= rm) then
            rm = ampo(mi,k)
            indr = mi
         endif
      enddo

      ipar(4,j) = 20.d0*log10(ampo(m,k)/sqrt(lm*rm))
      if ((1 == indl) .and. (indr == ntall)) then
         ipar(4,j) = ipar(4,j) + 100.d0
      end if
      ipar(5,j) = 0.5d0*dt*(abs(dble(m-indl)) + abs(dble(m-indr)))

   enddo

   ! End of SNR computations
   tim      = ipar(1,ia)
   tvis(k)  = ipar(2,ia)
   ampgr(k) = ipar(3,ia)
   grvel(k) = delta / (tim + t0)
   snr(k)   = ipar(4,ia)
   wdth(k)  = ipar(5,ia)
   phgr(k)  = ipar(6,ia)

enddo
deallocate(pha, ampo)
nfout1 = nf


! ==========================================
! Check dispersion curve for jumps
! ==========================================
nfout2 = 0
allocate(trig1(1:nf), ftrig(1:nf), stat=ier)
call trigger(grvel, om, nf, tresh, trig1, ftrig, ierr)


allocate(per1(1:nf), om1(1:nf), ftrig1(1:nf), stat=ier)
allocate(ii(1:nf), ijmp(1:nf), indx(1:nf), stat=ier)
allocate(tvist(1:nf), ampgrt(1:nf), grvelt(1:nf), &
         snrt(1:nf), wdtht(1:nf), phgrt(1:nf), stat=ier)
allocate(tvis1(1:nf), ampgr1(1:nf), grvel1(1:nf), &
         snr1(1:nf), wdth1(1:nf), phgr1(1:nf), stat=ier)



if (0 /= ierr) then

   do k = 1, nf, 1
      grvelt(k) = grvel(k)
      tvist(k)  = tvis(k)
      ampgrt(k) = ampgr(k)
      phgrt(k)  = phgr(k)
      snrt(k)   = snr(k)
      wdtht(k)  = wdth(k)
   enddo

   ! find all jumps
   njump = 0
   nijmp = 0
   do i = 1, nf-1, 1
      if (abs(trig1(i+1)-trig1(i)) > 1.5d0) then
         nijmp = nijmp + 1
         ijmp(nijmp) = i
      endif
   enddo

   nii = 0
   do i = 1, nijmp-1, 1
      if (ijmp(i+1)-ijmp(i) <= npoints) then
         nii = nii + 1
         ii(nii) = i
      endif
   enddo

   ! main loop by jumps
   if (0 /= nii) then

      do ki = 1, nii, 1

         kk = ii(ki)

         do i = 1, nf, 1
            tvis1(i)  = tvist(i)
            ampgr1(i) = ampgrt(i)
            grvel1(i) = grvelt(i)
            snr1(i)   = snrt(i)
            wdth1(i)  = wdth(i)
            phgr1(i)  = phgrt(i)
         enddo

         istrt = ijmp(kk)
         ibeg = istrt+1
         iend = ijmp(kk+1)

         ima = 0
         do k = ibeg, iend, 1

            dmaxt = huge(dmaxt)

            do j = 1, ici, 1

               if (k == ind(1,j)) then

                  wor = abs(delta/(ipar(1,j)+t0) - grvel1(k-1))
                  if (wor < dmaxt) then
                     ima = j
                     dmaxt = wor
                  endif

               endif

            enddo

            grvel1(k) = delta / (ipar(1,ima) + t0)
            tvis1(k)  = ipar(2,ima)
            ampgr1(k) = ipar(3,ima)
            phgr1(k)  = ipar(6,ima)
            snr1(k)   = ipar(4,ima)
            wdth1(k)  = ipar(5,ima)

         enddo

         call trigger(grvel1, om, nf, tresh, trig1, ftrig1, ier)

         iflag = 0
         do k = istrt, iend+1, 1
            if (abs(trig1(k)) >= 0.5d0) then
               iflag = 1
               exit
            end if
         enddo

         if (0 == iflag) then

            do i = 1, nf, 1
               tvist(i)  = tvis1(i)
               ampgrt(i) = ampgr1(i)
               grvelt(i) = grvel1(i)
               snrt(i)   = snr1(i)
               wdtht(i)  = wdth1(i)
               phgrt(i)  = phgr1(i)
               njump     = njump + 1
            enddo

         endif

      enddo

   endif ! if (0 /= nii) then

   do i = 1, nf, 1
      tvis1(i)  = tvist(i)
      ampgr1(i) = ampgrt(i)
      grvel1(i) = grvelt(i)
      snr1(i)   = snrt(i)
      wdth1(i)  = wdtht(i)
      phgr1(i)  = phgrt(i)
   enddo


   ! ===============================================================
   ! after removing possible jumps, we cut frequency range to single
   ! segment with max length
   ! ===============================================================
   call trigger(grvel1, om, nf, tresh, trig1, ftrig1, ier)

   if (0 /= ier) then

      nindx = 1
      indx(1) = 1

      do i = 1, nf, 1

         if (abs(trig1(i)) >= 0.5d0) then
            nindx = nindx + 1
            indx(nindx) = i
         endif

      enddo

      nindx = nindx+1
      indx(nindx) = nf

      imax = 0
      ipos = 0
      do i =1, nindx-1, 1

         iimax = indx(i+1)-indx(i)
         if (iimax > imax) then
            ipos = i
            imax = iimax
         endif

      enddo

      ist = max(indx(ipos), 1)
      ibe = min(indx(ipos+1), nf)
      nfout2 = ibe-ist+1

      do i = ist, ibe, 1
         per1(i-ist+1)   = per(i)
         tvis1(i-ist+1)  = tvis1(i)
         ampgr1(i-ist+1) = ampgr1(i)
         grvel1(i-ist+1) = grvel1(i)
         snr1(i-ist+1)   = snr1(i)
         wdth1(i-ist+1)  = wdth1(i)
         phgr1(i-ist+1)  = phgr1(i)
         om1(i-ist+1)    = om(i)
      enddo


      call trigger(grvel1(1:nfout1), om1(1:nfout1), nfout2, tresh, &
                             trig1(1:nfout1), ftrig1(1:nfout1), ier)
      if (nfout2 < 0.01d0*perc*nf) then
         ier = 1
         nfout2 = 0
      endif

   else

      nfout2 = nf
      per1(1:nf) = per(1:nf)
      !do i = 1, nf, 1
      !   per1(i) = per(i)
      !enddo

   endif ! if (0 /= ier) then

else

   ierr = 0
   nfout2 = nf
   do i = 1, nf, 1
      per1(i)   = per(i)
      tvis1(i)  = tvis(i)
      ampgr1(i) = ampgr(i)
      grvel1(i) = grvel(i)
      snr1(i)   = snr(i)
      wdth1(i)  = wdth(i)
      phgr1(i)  = phgr(i)
   enddo

endif
deallocate(ii, ijmp, indx)
deallocate(om, om1, trig1, ind, ipar)
deallocate(grvelt, tvist, snrt, wdtht, phgrt)
! ==========================================
! fill out output data arrays
! ==========================================
if ((0 /= nfout1) .and. (0 /= nphpr)) then
   allocate(phgrc(1:nfout1), stat=ier)
   call phtovel(delta, ip, nf, tvis, grvel, phgr, nphpr, phprper, phprvel, nfout1, phgrc)
   phgr(1:nfout1) = phgrc(1:nfout1)
   deallocate(phgrc)
endif


allocate(arr1(1:8,1:nfout1), stat=ier)
do i = 1, nfout1, 1
   arr1(1,i) = per(i)
   arr1(2,i) = tvis(i)
   arr1(3,i) = grvel(i)
   arr1(4,i) = phgr(i)
   arr1(5,i) = ampgr(i)
   arr1(6,i) = ftrig(i)
   arr1(7,i) = snr(i)
   arr1(8,i) = wdth(i)
enddo
deallocate(per, tvis, grvel, phgr)
deallocate(ampgr, ftrig, snr, wdth)




if (0 /= nfout2) then

   if (0 /= nphpr) then

      allocate(phgrc(1:nfout2), stat=ier)
      call phtovel(delta, ip, nf, tvis1, grvel1, phgr1, nphpr, phprper, phprvel, nfout2, phgrc)
      phgr1(1:nfout2) = phgrc(1:nfout2)
      deallocate(phgrc)

   endif

   allocate(arr2(1:7,1:nfout2), stat=ier)
   do i = 1, nfout2, 1
      arr2(1,i) = per1(i)
      arr2(2,i) = tvis1(i)
      arr2(3,i) = grvel1(i)
      arr2(4,i) = phgr1(i)
      arr2(5,i) = ampgr1(i)
      arr2(6,i) = snr1(i)
      arr2(7,i) = wdth1(i)
   enddo

else

   ierr = 2

endif
deallocate(per1, tvis1, grvel1, phgr1)
deallocate(ampgr1, ftrig1, snr1, wdth1)


return


end subroutine aftanpg


end module aftanpg_m
