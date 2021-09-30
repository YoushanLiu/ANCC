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
! aftanipg function. Provides ftan analysis with phase match filter,
! jumps correction, phase velocity computation and amplitude map
! output for input periods.
!
!
! Autor: M. Barmine, CIEI, CU. Date: Jun 15, 2006. Version: 2.00
! Modified by Youshan Liu, at Institute of Geology and Geophysics, Chinese Academy of Sciences.
!
module aftanipg_m

implicit none


contains


subroutine aftanipg(PIover4, n, seis, t0, dt, delta, vmin, vmax, &
            tmin, tmax, tresh, ffact, perc, npoints, taperl, nfin, &
            fsnr, fmatch, npred, pred, nphpr, phprper, phprvel, &
            nfout1, arr1, nfout2, arr2, tamp, nrow, ncol, amp, ierr)


!======================================================================
! Parameters for aftanipg function:
! Input parameters:
! PIover4   - phase shift = PI/4*PIover4, for cross -correlation
!             PIover4 should be -1.0 (real(8))
! n         - number of input samples, (integer(4))
! seis      - input array length of n, (real(4))
! t0        - time shift of SAC file in seconds, (real(4))
! dt        - sampling rate in seconds, (real(4))
! delta     - distance, km (real(4))
! vmin      - minimal group velocity, km/s (real(4))
! vmax      - maximal value of the group velocity, km/s (real(4))
! tmin      - minimal period, s (real(4))
! tmax      - maximal period, s (real(4))
! tresh     - treshold, usualy = 10, (real(4))
! ffact     - factor to automatic filter parameter, usualy =1, (real(4))
! perc      - minimal length of of output segment vs freq. range, % (real(4))
! npoints   - max number points in jump, (integer(4))
! taperl    - factor for the left end seismogram tapering,
!             taper = taperl*tmax,    (real(4))
! nfin      - starting number of frequencies, nfin <= 100,(integer(4))
! fsnr      - phase match filter parameter, spectra ratio to
!             determine cutting point   (real(4))
! fmatch    - factor to length of phase matching window (real(4))
! npred     - length of the group velocity prediction table
! pred      - group velocity prediction table:    (real(4))
! pred(0,:) - periods of group velocity prediction table, s
! pred(1,:) - pedicted group velocity, km/s
! nphpr     - length of phprper and phprvel arrays
! phprper   - predicted phase velocity periods, s
! phprvel   - predicted phase velocity for corresponding periods, s
! ==========================================================
! Output parameters are placed in 2-D arrays arr1 and arr2,
! arr1 contains preliminary results and arr2 - final.
! ==========================================================
! nfout1    - output number of frequencies for arr1, (integer(4))
! arr1      - preliminary results.
!             Description: real(4) arr1(8,n), n >= nfin)
! arr1(1,:) - central periods, s (real(4))
! arr1(2,:) - apparent periods, s (real(4))
! arr1(3,:) - group velocities, km/s (real(4))
! arr1(4,:) - phase velocities, km/s (real(4))
! arr1(5,:) - amplitudes, Db (real(4))
! arr1(6,:) - discrimination function, (real(4))
! arr1(7,:) - signal/noise ratio, Db (real(4))
! arr1(8,:) - maximum half width, s (real(4))
! arr2      - final results
! nfout2    - output number of frequencies for arr2, (integer(4))
!             Description: real(4) arr2(7,n), n >= nfin)
!             If nfout2 == 0, no final results.
! arr2(1,:) - central periods, s (real(4))
! arr2(2,:) - apparent periods, s (real(4))
! arr2(3,:) - group velocities, km/s (real(4))
! arr1(4,:) - phase velocities, km/s (real(4))
! arr2(5,:) - amplitudes, Db (real(4))
! arr2(6,:) - signal/noise ratio, Db (real(4))
! arr2(7,:) - maximum half width, s (real(4))
! tamp      - time to the beginning of ampo table, s (real(4))
! nrow      - number of rows in array ampo, (integer(4))
! ncol      - number of columns in array ampo, (integer(4))
! amp       - FTAN amplitude array, Db, (real(4))
! ierr      - completion status, =0 - O.K.,           (integer(4))
!                                =1 - some problems occures
!                                =2 - no final results
!=====================================================================

use taper_m
use spline_m


implicit none

include 'fftw3.f'

integer(4), intent(in) :: n, npoints
integer(4), intent(in) :: nfin, nphpr, npred

real(4), intent(in) :: seis(n)

real(4), intent(in) :: t0, dt, delta
real(4), intent(in) :: tresh, ffact, perc
real(4), intent(in) :: fsnr, fmatch, taperl
real(4), intent(in) :: vmin, vmax, tmin, tmax

real(8), intent(in) :: PIover4

real(4), intent(in) :: phprper(npred), phprvel(npred)
real(4), intent(in) :: pred(npred,2)


integer(4), intent(out) :: nrow, ncol
integer(4), intent(out) :: nfout1, nfout2, ierr

real(4), intent(out) :: tamp

real(4), dimension(:,:), allocatable, intent(out) :: amp, arr1, arr2


integer(4) i, j, k, m, nf, ier
integer(4) njump, nijmp, nii, nq
integer(4) ici, icj, iciflag, ia
integer(4) iflag, nindx, imax, iimax
integer(4) mm, mi, iml, imr, indl, indr
integer(4) ki, kk, istrt, ibeg, iend, ima
integer(4) ipos, ist, ibe, ip, inds, inde
integer(4) ntapb, ntape, ne, nb, ns, ntime, ntall

integer(8) plan1, plan2, plan3, plan4

real(4) alpha, omb, ome, dom, step
real(4) lm, rm, om0, tg0, omstart, dw
real(4) t, dmaxt, wor, amax, dph, tm, ph
real(4) pha_corr, ome1, omb1, minTpr, maxTpr

real(8), parameter :: PI = 4.d0*datan(1.d0)

complex czero, dci, dc2

integer(4), dimension(:), allocatable :: indx, ii, ijmp
integer(4), dimension(:,:), allocatable :: ind

real(4), dimension(:), allocatable :: tim, phgrc, trig1, omdom, ampdom
real(4), dimension(:), allocatable :: om, per, om1, per1, ftrig, ftrig1
real(4), dimension(:), allocatable :: grvel, tvis, ampgr, phgr, snr, wdth
real(4), dimension(:), allocatable :: grvel1, tvis1, ampgr1, phgr1, snr1, wdth1
real(4), dimension(:), allocatable :: grvelt, tvist, ampgrt, phgrt, snrt, wdtht

real(4), dimension(:,:), allocatable :: ipar, pha, ampo

complex, dimension(:), allocatable :: pha_cor, env, spref
complex, dimension(:), allocatable :: s, sf, fils, tmp



ierr = 0

czero = cmplx(0.0, 0.0)
dci   = cmplx(0.0, 1.0)
dc2   = cmplx(2.0, 0.0)


iml = 0
imr = 0
lm = 0.0
rm = 0.0


! number of FTAN filters
nf = nfin

! automatic width of filters * factor ffact
alpha = ffact*20.0*sqrt(delta*1.e-3)

! number of samples for tapering, left end
ntapb = nint(taperl*tmax/dt)

! number of samples for tapering, right end
ntape = nint(tmax/dt)

! [omb,ome] - frequency range
omb = 2.0*PI/tmax
ome = 2.0*PI/tmin

! find min/max of prediction period
maxTpr = maxval(pred(1:npred,1))
minTpr = minval(pred(1:npred,2))


! evaluation of spline polynomial forms for phase match filter
ip = 1
call pred_cur(ip, delta, sqrt(omb*ome), npred, pred, om0, tg0)


! seismgram tapering
nb = max(2, nint((delta/vmax-t0)/dt))
ne = min(n, nint((delta/vmin-t0)/dt))
tamp = (nb-1)*dt + t0


nrow = nfin
ncol = ne-nb+1
!! times for FTAN map
!allocate(time(nb:ne), v(nb:ne), stat=ier)
!do k = nb, ne, 1
!   time(k) = (k-1)*dt
!   ! velocity for FTAN map
!   v(k) = delta / (time(k) + t0)
!enddo
ntime = ne-nb+1


! tapering both ends of seismogram
call taper(max(nb,ntapb+1), min(ne,n-ntape), n, seis, ntapb, ntape, ns, s)


! prepare FTAN filters
dom = 2.0*PI/(ns*dt)
step =(log(omb) - log(ome))/real(nf - 1)

! log scaling for frequency
allocate(om(1:nf), per(1:nf), stat=ier)
do k = 1, nf, 1
   om(k) = exp(log(ome) + (k-1)*step)
   per(k) = 2.0*PI/om(k)
enddo

!==================================================================
! Phase match filtering
! make backward FFT for seismogram: s ==> sf
allocate(sf(1:ns), stat=ier)
call sfftw_plan_dft_1d(plan1, ns, s, sf, FFTW_FORWARD, FFTW_ESTIMATE)
call sfftw_execute(plan1)
call sfftw_destroy_plan(plan1)


! filtering and FTAN amplitude diagram construction
allocate(fils(1:ns), tmp(1:ns), stat=ier)
call sfftw_plan_dft_1d(plan2, ns, fils, tmp, FFTW_BACKWARD, FFTW_ESTIMATE)

nq = ns/2 + 1
sf(1) = 0.50*sf(1)
sf(nq) = cmplx(real(sf(nq)), 0.0)


! spectra tapering
ome1 = min(ome, 2.0*PI/minTpr)
omb1 = max(omb, 2.0*PI/maxTpr)
allocate(omdom(1:ns), ampdom(1:ns), stat=ier)
call tapers(omb1, ome1, dom, alpha, ns, omstart, inds, inde, omdom, ampdom)
omstart = real(nint(omstart/dom))*dom


inde = min(inde, nq+1)
allocate(pha_cor(1:ns), stat=ier)
do i = 1, ns, 1

   if (i < inds) then
      pha_cor(i) = czero
   elseif (i > inde) then
      pha_cor(i) = czero
   else
      call msplint(ip+1, npred, om0, omdom(i), pha_corr, ierr)
      pha_cor(i) = cmplx(pha_corr, 0.0)
   endif

   pha_cor(i) = exp(dci*pha_cor(i))
   sf(i) = sf(i)*pha_cor(i)*ampdom(i)

enddo
deallocate(omdom, ampdom)
call free_mspline()


!  forward FFT to get signal envelope
allocate(env(1:ns), stat=ier)
call sfftw_plan_dft_1d(plan3, ns, sf, env, FFTW_BACKWARD, FFTW_ESTIMATE)
call sfftw_execute(plan3)
call sfftw_destroy_plan(plan3)


env(1:ns) = env(1:ns) * dc2 / real(ns)

! cutting impulse response in time
dw = om(1) - om(nf)

allocate(spref(1:ns), stat=ier)
call tgauss(fsnr, tg0, dw, dt, ns, fmatch, env, spref)
deallocate(env)

! back to spectra after filtering
call sfftw_plan_dft_1d(plan4, ns, spref, sf, FFTW_FORWARD, FFTW_ESTIMATE)
call sfftw_execute(plan4)
call sfftw_destroy_plan(plan4)
deallocate(spref)


! apply back phase correction for spectra
sf(1:ns) = sf(1:ns) / pha_cor(1:ns)
deallocate(pha_cor)
!==================================================================
ntall = ntime + 2
allocate(pha(1:ntall,1:nf), amp(1:ntall,1:nf), ampo(1:ntall,1:nf), stat=ier)
! main loop by frequency
do k = 1, nf, 1

   ! filtering
   call ftfilt(alpha, om(k), dom, ns, sf, fils)

   ! fill with zeros half spectra for Hilbert transformation and spectra ends ajastment
   fils(nq+1:ns) = czero

   fils(1) = cmplx(0.50*real(fils(1)), 0.0)
   fils(nq) = dcmplx(real(fils(nq)), 0.0)

   ! forward FFT: fils ==> tmp
   call sfftw_execute(plan2)

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
call sfftw_destroy_plan(plan2)
deallocate(s, sf, fils, tmp)

! normalization amp diagram to 100 Db with three decade cutting
amax = maxval(amp(1:ntall,1:nf))
amp(1:ntall,1:nf) = max(amp(1:ntall,1:nf) + 100.0 - amax, 40.0)


! construction reference indices table ind. It points to local maxima.
! table ipar contains three parameter for each local maximum:
! tim - group time tvis - observed period ampgr - amplitude values in Db
nfout1 = 0
allocate(tim(1:nf), tvis(1:nf), ampgr(1:nf), stat=ier)
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

      ipar(4,j) = 20.0*log10(ampo(m,k)/sqrt(lm*rm))
      if((1 == indl) .and. (indr == ntall)) then
         ipar(4,j) = ipar(4,j) + 100.0
      end if

      ipar(5,j) = 0.50*dt*(abs(dble(m-indl)) + abs(dble(m-indr)))

   enddo
   ! End of SNR computations

   tim(k)   = ipar(1,ia)
   tvis(k)  = ipar(2,ia)
   ampgr(k) = ipar(3,ia)
   grvel(k) = delta / (tim(k) + t0)
   snr(k)   = ipar(4,ia)
   wdth(k)  = ipar(5,ia)
   phgr(k)  = ipar(6,ia)

enddo
deallocate(ind, pha, ampo, tim)
nfout1 = nf


! ==========================================
!       Check dispersion curve for jumps
! ==========================================
nfout2 = 0
allocate(trig1(1:nf), ftrig(1:nf), stat=ier)
call trigger(grvel, om, nf, tresh, trig1, ftrig, ierr)


allocate(per1(1:nf), om1(1:nf), ftrig1(1:nf), stat=ier)
allocate(ii(1:nf), ijmp(1:nf), indx(1:nf), stat=ier)
allocate(grvelt(1:nf), tvist(1:nf), ampgrt(1:nf), &
         phgrt(1:nf), snrt(1:nf), wdtht(1:nf), stat=ier)
allocate(grvel1(1:nf), tvis1(1:nf), ampgr1(1:nf), &
         phgr1(1:nf), snr1(1:nf), wdth1(1:nf), stat=ier)



if (0 /= ierr) then

   do k = 1, nf, 1
      grvelt(k) = grvel(k)
      tvist(k)  = tvis(k)
      ampgrt(k) = ampgr(k)
      ampgrt(k) = ampgr(k)
      phgrt(k)  = phgr(k)
      wdtht(k)  = wdth(k)
   enddo

   ! find all jumps
   njump = 0
   nijmp = 0
   do i = 1, nf-1, 1
      if (abs(trig1(i+1)-trig1(i)) > 1.50) then
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
            grvel1(i) = grvelt(i)
            tvis1(i)  = tvist(i)
            ampgr1(i) = ampgrt(i)
            phgr1(i)  = phgrt(i)
            snr1(i)   = snrt(i)
            wdth1(i)  = wdth(i)
         enddo

         istrt = ijmp(kk)
         ibeg = istrt+1
         iend = ijmp(kk+1)

         ima = 0
         do k = ibeg, iend, 1

            dmaxt = huge(dmaxt)

            do j = 1, ici, 1

               if (ind(1,j) == k) then

                  wor = abs(delta/(ipar(1,j)+t0) - grvel1(k-1))
                  if (wor < dmaxt) then
                     ima = j
                     dmaxt = wor
                  endif

               endif

            enddo

            grvel1(k) = delta / (ipar(1,ima)+t0)
            tvis1(k)  = ipar(2,ima)
            ampgr1(k) = ipar(3,ima)
            phgr1(k)  = ipar(6,ima)
            snr1(k)   = ipar(4,ima)
            wdth1(k)  = ipar(5,ima)

         enddo

         call trigger(grvel1, om, nf, tresh, trig1, ftrig1, ier)

         iflag = 0
         do k = istrt, iend+1, 1
            if (abs(trig1(k)) >= 0.50) then
               iflag = 1
               exit
            end if
         enddo

         if (0 == iflag) then

            do i = 1, nf, 1
               grvelt(i) = grvel1(i)
               tvist(i)  = tvis1(i)
               ampgrt(i) = ampgr1(i)
               phgrt(i)  = phgr1(i)
               snrt(i)   = snr1(i)
               wdtht(i)  = wdth1(i)
               njump = njump + 1
            enddo

         endif

      enddo

   endif

   do i = 1, nf, 1
      grvel1(i) = grvelt(i)
      tvis1(i)  = tvist(i)
      ampgr1(i) = ampgrt(i)
      phgr1(i)  = phgrt(i)
      snr1(i)   = snrt(i)
      wdth1(i)  = wdtht(i)
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

         if (abs(trig1(i)) >= 0.50) then
            nindx = nindx+1
            indx(nindx) = i
         endif

      enddo

      nindx = nindx+1
      indx(nindx) = nf

      imax = 0
      ipos = 0
      do i = 1, nindx-1, 1

         iimax = indx(i+1)-indx(i)
         if(iimax > imax) then
            ipos = i
            imax = iimax
         endif

      enddo

      ist = max(indx(ipos), 1)
      ibe = min(indx(ipos+1), nf)
      nfout2 = ibe -ist+1

      do i = ist, ibe, 1
         per1(i-ist+1)   = per(i)
         grvel1(i-ist+1) = grvel1(i)
         tvis1(i-ist+1)  = tvis1(i)
         ampgr1(i-ist+1) = ampgr1(i)
         phgr1(i-ist+1)  = phgr1(i)
         snr1(i-ist+1)   = snr1(i)
         wdth1(i-ist+1)  = wdth1(i)
         om1(i-ist+1)    = om(i)
      enddo

      call trigger(grvel1, om1, nfout2, tresh, trig1, ftrig1, ier)
      if(nfout2 < 0.010*nf*perc) then
         ier = 1
         nfout2 = 0
      endif

   else

      nfout2 = nf
      !do i = 1, nf, 1
      !   per1(i) = per(i)
      !enddo
      per1(1:nf) = per(1:nf)

   endif

else

   ierr = 0
   nfout2 = nf
   do i = 1, nf, 1
      per1(i)   = per(i)
      tvis1(i)  = tvis(i)
      ampgr1(i) = ampgr(i)
      phgr1(i)  = phgr(i)
      grvel1(i) = grvel(i)
      snr1(i)   = snr(i)
      wdth1(i)  = wdth(i)
   enddo

endif
deallocate(ii, ijmp, indx)
deallocate(om, om1, trig1, ipar)
deallocate(grvelt, tvist, snrt, wdtht, phgrt)
! ===========================================================
! fill out output data arrays
! ===========================================================
allocate(phgrc(1:nfout1), phgr(1:nfout1), stat=ier)
if ((0 /= nfout1) .and. (0 /= nphpr)) then
   call phtovel(delta, ip, nfout1, tvis, grvel, phgr, nphpr, phprper, phprvel, phgrc)
   phgr(1:nfout1) = phgrc(1:nfout1)
endif
deallocate(phgrc)


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

      allocate(phgrc(1:nfout2), phgr1(1:nfout2), stat=ier)
      call phtovel(delta, ip, nfout2, tvis1, grvel1, phgr1, nphpr, phprper, phprvel, phgrc)
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
   deallocate(per1, tvis1, grvel1, phgr1)
   deallocate(ampgr1, ftrig1, snr1, wdth1)

else

   ierr = 2

endif


return


end subroutine aftanipg


end module aftanipg_m
