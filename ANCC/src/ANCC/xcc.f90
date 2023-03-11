! This file is part of ANCC.
!
! ANCC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! ANCC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
module xcc_m

use, intrinsic :: iso_c_binding     ! Allow to define the equivalents of C data types (e.g. c_ptr, C_INT)

use db_m                 ! imported module containing data type and variable definitions.
use math_m
use sac_io_m
use string_m
use date_time_m
use bindata_io_m


implicit none


include 'fftw3.f'                   ! fftw3.f: contains the Fortran constant definitions and interface definitions (e.g. FFTW_BACKWARD)


real(DBL), parameter :: PI = 4.d0*datan(1.d0)



contains


! =======================================================================================
! Process the sac file for one particular is_save_record and fill in the sdb.rec info.
! evpath: event path containing the sac data [input]
! sdb: sac_db struct [input and output]
! iev: event iterator [input]
! ist: station iterator [input]
! channel: channel name [input]
! is_verbose: verbose indicator [input]
! =======================================================================================
subroutine mk_one_rec(evpath, iev, ist, channel, sdb, sdb_tmp)

implicit none

integer, intent(in) :: iev, ist

character(len=*), intent(in) :: evpath
character(len=*), intent(in) :: channel


type(sac_db), intent(in) :: sdb
type(sac_db), intent(inout) :: sdb_tmp


integer npts
integer nstrArray, ier

logical is_existed

real(DBL) frac, dt, t0

character(len=512) sacinfile, sacoutfile, sacname, str_date

character(len=128), allocatable, dimension(:) :: strArray



! ***************************************************************
! Copy one station sac file from evpath to target folder if it exists.
! ***************************************************************

sacname = trim(adjustl(sdb%st(ist)%n_name))// &
            '.'//trim(adjustl(channel))//'.SAC'
sacinfile = trim(adjustl(evpath))//'/'//trim(adjustl(sacname))
if (is_overwrite_data) then
   sacoutfile = trim(adjustl(sacinfile))
else
   sacoutfile = trim(adjustl(sdb_tmp%ev(iev)%name))//'/'//trim(adjustl(sacname))
end if


! Return if it doesn't exist.
inquire(file=sacinfile, exist=is_existed)
if (.not.(is_existed)) return



! ***************************************************************
if (is_verbose) then
   call split_string(evpath, '/', strArray, nstrArray)
   str_date = strArray(nstrArray)
   write(*,"('Event: ',A,'   Station: ',A)") trim(adjustl(str_date)), &
                                      trim(adjustl(sdb%st(ist)%n_name))
   call flush(6)
   if (allocated(strArray)) then
      deallocate(strArray)
   end if
end if



! ***************************************************************
! Correct the fraction time and set the reference time to be the beginning time (b=0)).
! ***************************************************************
call correct_sac_file(sacinfile, sacoutfile, frac, npts, dt, t0, ier)
if (0 /= ier) then
   write(*,"(A)") 'Error: correct_sac_file failed ! '//trim(adjustl(sacname))
   call flush(6)
   return
end if


! ***************************************************************
! Fill the elements in this is_save_record.
! ***************************************************************
if (0 /= npts) then
   sdb_tmp%rec(ist,iev)%npts = npts
   sdb_tmp%rec(ist,iev)%frac = frac
   sdb_tmp%rec(ist,iev)%dt = dt
   sdb_tmp%rec(ist,iev)%t0 = t0
   sdb_tmp%rec(ist,iev)%name = ''
   if (is_overwrite_data) then
      sdb_tmp%rec(ist,iev)%name = trim(adjustl(sacinfile))
   else
      sdb_tmp%rec(ist,iev)%name = trim(adjustl(sacoutfile))
   end if
   sdb_tmp%rec(ist,iev)%channel = ''
   sdb_tmp%rec(ist,iev)%channel = trim(adjustl(channel))
else
   sdb_tmp%rec(ist,iev)%npts = 0
   sdb_tmp%rec(ist,iev)%frac = 0.0
end if


end subroutine mk_one_rec



! =======================================================================================
! Correct the fraction time and set the reference time to be the beginning time (b=0)).
! if t=10.6s, then t=11s, frac=-0.4s; if t=10.4s, then t=10s, frac=0.4s
! fname: sac filename [input]
! frac: fraction time of this sac file  [output]
! npts: npts of this sac header [output]
! dt: time sampling interval [output]
! t: beginning time of the first data point [output]
! ier: status indicator [output]
! =======================================================================================
subroutine correct_sac_file(finname, foutname, frac, npts, dt, t, ier)

implicit none

integer, intent(out) :: npts, ier

real(DBL), intent(out) :: frac, dt, t

character(len=*), intent(in) :: finname, foutname


integer k, nf

real(SGL) coeff
real(DBL) tf, sec

type(sachead) shd

real(SGL), allocatable, dimension(:) :: seis_data



if (0 == len_trim(adjustl(finname))) return


! ***************************************************************
! read the sac file
call sacio_readsac(finname, shd, seis_data, ier)
if (0 /= ier) then
   write(*,*) "Error: Cannot read: "//trim(adjustl(finname))
   call flush(6)
   return
end if


! npts and dt
npts = shd%npts
! remove round error
dt = nint(shd%delta*1e6)*1.d-6


! Get the initial beginning time of the first data point.
t = datetime2timestamp(shd%nzyear, shd%nzjday, shd%nzhour, shd%nzmin, shd%nzsec+1.d-3*shd%nzmsec+shd%b)


! Apply taper
do k = 1, ntaper, 1
   coeff = 0.50*(1.0 + cos(PI*(ntaper-k+1)/dble(ntaper)))
   seis_data(k) = seis_data(k)*coeff
   seis_data(npts-k+1) = seis_data(npts-k+1)*coeff
end do


! ***************************************************************
! Make the time fraction correction
! ***************************************************************
!tf = floor(t)
!nf = nint(sngl(t - tf)/dt))
!frac = t - (tf + nf*dt)
!t = tf + nf*dt
!!if (frac > 0.5*dt) then
!!   t = t + dt
!!   frac = frac - dt
!!end if
!tf = floor(t)
!t = tf + nint(sngl((t - tf)*1.d6))*1.d-6
!frac = nint(sngl(frac*1.d6))*1.d-6



tf = nint(t)
nf = nint((t - tf)/dt)
tf = tf + nf*dt
frac = t - tf
t = tf




! ***************************************************************
! Change the sac header to make sure b=0
! ***************************************************************
call timestamp2datetime(t, shd%nzyear, shd%nzjday, shd%nzhour, shd%nzmin, sec)
shd%nzsec = int(sec)
shd%nzmsec = int((sec-int(sec))*1000.d0)
shd%b = 0.0
!shd%user1 = real(frac)


! Overwrite the sac file
call sacio_writesac(foutname, shd, seis_data, ier)
if (0 /= ier) then
   write(*,*) "Error: Cannot overwrite: "//trim(adjustl(foutname))
   call flush(6)
   return
end if


deallocate(seis_data)



end subroutine correct_sac_file



! =======================================================================================
! Write the info database into a ascii file if is_save_record is true.
! sdb: sac_db struct [input]
! filename: name of the output ascii file [input]
! =======================================================================================
subroutine sacdb_to_asc(sdb, filename)

implicit none

character(len=*), intent(in) :: filename

type(sac_db), intent(in) :: sdb


integer iev, ist
integer nstrArray, ier

type(sachead) shd

character(len=512) str_date, sacname

character(len=128), allocatable, dimension(:) :: strArray




if (0 == len_trim(adjustl(filename))) return



! ***************************************************************
open(unit=17, file=filename, status='replace', action='write', iostat=ier)

   ! ***************************************************************
   ! Write the number of stations and events in the first line.
   ! ***************************************************************
   write(17,"(A,I8,5X,A,I6)") 'Number of events:', sdb%nev, 'Number of stations:',sdb%nst
   write(17,"(A)") '===================================================================='
   call flush(17)

   ! ***************************************************************
   ! Write the data is_save_record.
   ! ***************************************************************
   do iev = 1, sdb%nev, 1
      do ist = 1, sdb%nst, 1

         call split_string(sdb%ev(iev)%name, '/', strArray, nstrArray)
         str_date = strArray(nstrArray)

         write(17,"(A20,$)") trim(adjustl(str_date))

         if (0 == sdb%rec(ist,iev)%npts) then    ! Write "NO DATA" if rec[ie][is] == 0

            write(17,"(A)") 'NO DATA at '//trim(adjustl(sdb%st(ist)%n_name))

         else

            ! read the sac file to retrive its header information into shd struct.
            call sacio_readhead(sdb%rec(ist,iev)%name, shd, ier)

            ! Write sac file name, t0 (reference time), frac(time fraction),
            ! data length (npts*delta)
            call split_string(sdb%rec(ist,iev)%name, '/', strArray, nstrArray)
            sacname = strArray(nstrArray)

            if (allocated(strArray)) then
               deallocate(strArray)
            end if

            !write(17,"(A30,3X,'t0: ',I4,'/',I3.3,'/',I2.2,':',I2.2,':',A,4X,'Frac:', &
            !     &F10.5,'s',4X,'Record Length:',F10.2,'s')") trim(adjustl(sacname)), &
            !                          shd%nzyear, shd%nzjday, shd%nzhour, shd%nzmin, &
            !                 trim(adjustl(padzero(shd%nzsec+0.001*shd%nzmsec,2,3))), &
            !                                 sdb%rec(ist,iev)%frac, shd%delta*shd%npts
            write(17,"(A30,3X,'t0: ',I4,'/',I3.3,'/',I2.2,':',I2.2,':',F10.5,4X,'Frac:', &
                     &F10.5,'s',4X,'Record Length:',F10.2,'s')") trim(adjustl(sacname)), &
                              shd%nzyear, shd%nzjday, shd%nzhour, shd%nzmin, shd%nzsec + &
                           0.001d0*shd%nzmsec, sdb%rec(ist,iev)%frac, shd%delta*shd%npts

         end if

         call flush(17)

      end do
   end do

   if (allocated(strArray)) then
      deallocate(strArray)
   end if

close(unit=17)


end subroutine sacdb_to_asc



! =======================================================================================
! Remove the instrument response
! sdb: sac_db struct [input]
! iev: event iterator [input]
! ist: station iterator [input]
! myrank: process id number [input]
! f1, f2, f3, f4: freqency limits [input]
! is_verbose: verbose indicator [input]
! =======================================================================================
subroutine remove_RESP(sdb, iev, ist, f1, f2, f3, f4, pzfolder)

implicit none

integer, intent(in) :: iev, ist

real(SGL), intent(in) :: f1, f2, f3, f4

character(len=*), intent(in) :: pzfolder

type(sac_db), intent(in) :: sdb



logical is_existed

character(len=128) str_myrank

character(len=512) pzfile


write(pzfile,"(A)") trim(adjustl(pzfolder))//'/'//trim(adjustl(sdb%st(ist)%n_name))// &
                                   '..'//trim(adjustl(sdb%rec(ist,iev)%channel))//'.PZ'

inquire(file=trim(adjustl(pzfile)), exist=is_existed)
if (.not.(is_existed)) return


! ***************************************************************
if ((f1 > 0.0) .and. (f2 > f1) .and. (f3 > f2) .and. (f4 > f3)) then

   if (sdb%rec(ist,iev)%npts > 0) then

      ! Each process has its own sac script

      write(str_myrank, "(A, I6.6)") './tmp/', myrank
      open(unit=18, file=trim(adjustl(str_myrank))//'.sh', status='replace', action='write')
         write(18, "(A)") 'sac<<EOF'
         write(18, "(A)") 'r '//trim(adjustl(sdb%rec(ist,iev)%name))
         write(18, "(A)") 'rmean'
         write(18, "(A)") 'rtrend'
         write(18, "(A)") 'taper'
         write(18, "(A,F10.4,F10.4,F10.4,F10.4)") 'transfer from polezero subtype '// &
                                trim(adjustl(pzfile))//' to vel freq ', f1, f2, f3, f4
         !write(18, "(A)") 'mul 1.e9'
         write(18, "(A)") 'w over'
         write(18, "(A)") 'quit'
         write(18, "(A)") 'EOF'
         call flush(18)
      close(unit=18)

      ! ***************************************************************
      ! Remove the instrument response to obtain the velocity mesurement
      ! using transfer command in SAC with frequency limits [f1,f2,f3,f4].
      ! ***************************************************************
      if (is_verbose) then
         call system('sh '//trim(adjustl(str_myrank))//'.sh')
      else
         call system('sh '//trim(adjustl(str_myrank))//'.sh > /dev/null')
      end if

      call system('rm -rf '//trim(adjustl(str_myrank))//'.sh')

   end if

else

   write(*,"(A)") 'The corner periods should satisfy f1 < f2 < f3 < f4 !'
   call flush(6)

end if


end subroutine remove_RESP



! =======================================================================================
! Apply fractional time correction, temporal domain normalization, spectra whitening,
! band-rejection filtering [optional], band-pass filtering, cut data, and forward FFT
! sdb: sac_db struct [input]
! iev: event iterator [input]
! ist: station iterator [input]
! myrank: process id number [input]
! is_running_time_average: time normalization indicator [input]
! is_onebit: is_onebit normalization indicator [input]
! is_suppress_notch: is_suppress_notch indicator [input]
! f1, f2, f3, f4: frequency limits [input]
! is_bandpass_earthquake: if earthquake band_pass filtering at [fr1 fr2] [input]
! fr1, fr2: nperiod limits for earthquake band-pass filtering in time normalization [input]
! npow_costaper: power of cosine tapering function [input]
! nwt, nwf: half-window length time normalization and spectral whitening [input]
! freqmin: retaining factor for is_suppress_notch repressing [input]
! t0: starting time [input]
! tlen: data length [input]
! is_verbose: verbose indicator [input]
! =======================================================================================
subroutine preprocess(sdb, iev, ist, npow_costaper, nwt, nwf, &
                   f1, f2, f3, f4, fr1, fr2, freqmin, t0, tlen)

implicit none


include 'fftw3.f'


integer, intent(in) :: iev, ist
integer, intent(in) :: nwt, nwf, npow_costaper

real(SGL), intent(in) :: f1, f2, f3, f4
real(SGL), intent(in) :: fr1, fr2, freqmin

real(DBL), intent(in) :: t0, tlen

type(sac_db), intent(in) :: sdb



integer n, nfft, nq
integer k, ngap, norder
integer nstrArray, n1, n2, ier

logical is_existed

real(SGL) df, dfft, coeff, wtr

real(DBL) dt, frac, tend
real(DBL) time, trb, tre, sec

type(c_ptr) planf, planb

type(sachead) shd

character(len=128) str_myrank

character(len=512) sacname

real(SGL), allocatable, dimension(:) :: seis_data, abs_data, wgt_data

real(DBL), allocatable, dimension(:) :: a, b, tr

character(len=128), allocatable, dimension(:) :: strArray

complex(SGL), allocatable, dimension(:) :: s, sf



! ***************************************************************
!if (0 == sdb%rec(ist,iev)%npts) return
inquire(file=trim(adjustl(sdb%rec(ist,iev)%name)), exist=is_existed)
if (.not.(is_existed)) return


! Each process has its own sac script
write(str_myrank, '(A, I6.6)') './tmp/', myrank


! read the sac file
call sacio_readsac(trim(adjustl(sdb%rec(ist,iev)%name)), shd, seis_data, ier)


call split_string(sdb%rec(ist,iev)%name, '/', strArray, nstrArray)
sacname = trim(adjustl(strArray(nstrArray-1)))//'/'//trim(adjustl(strArray(nstrArray)))




n = shd%npts
! remove round error
dt = nint(shd%delta*1e6)*1.d-6



! =================================================================================================
! ======================================= Fractional correction ===================================
! =================================================================================================


! Obtain the useful header information.
!frac = shd%user1
frac = sdb%rec(ist,iev)%frac


if (abs(frac) > 0.01*dt) then


   ! Determine the power for FFT
   nfft = 2**ceiling(log(dble(n))/log(2.d0))   ! nfft: number of points for FFT
   nq = nfft/2 + 1
   dfft = 1.0 / dble(nfft)

   df = dfft / dt     ! df: frequency interval


   ! Allocate memory for s and sf.
   allocate(s(nfft), sf(nfft), stat=ier)


   call sfftw_plan_dft_1d(planf, nfft, s, sf, FFTW_FORWARD, FFTW_ESTIMATE)
   call sfftw_plan_dft_1d(planb, nfft, sf, s, FFTW_BACKWARD, FFTW_ESTIMATE)



   ! Initialize s with complex zero.
   s = czero


   ! Apply taper
   do k = 1, ntaper, 1
      coeff = 0.50*(1.0 + cos(PI*(ntaper-k+1)/dble(ntaper)))
      seis_data(k) = seis_data(k)*coeff
      seis_data(n-k+1) = seis_data(n-k+1)*coeff
   end do


   ! Fill s with real data.
   s(1:n) = cmplx(seis_data(1:n), 0.0)


   ! Make forward FFT for the seismogram: s => sf
   call sfftw_execute_dft(planf, s, sf)


   ! Make time fraction correction
   do k = 1, nq, 1
      sf(k) = sf(k) * exp(-ci*(2.0*PI*(k-1)*df)*frac)
   end do


   ! Kill half spectra
   sf(nq+1:nfft) = czero


   ! Correct the ends
   sf(1) = 0.50*sf(1)
   !sf(nq) = cmplx(real(sf(nq)), 0.0)


   ! ***************************************************************
   ! Bandpass filtering.
   ! ***************************************************************
   call bandpass_filter(f1, f2, f3, f4, df, nq, npow_costaper, sf)


   ! Make backward FFT for the seismogram: sf => s
   call sfftw_execute_dft(planb, sf, s)


   ! Get the final result.
   seis_data(1:n) = 2.0*real(s(1:n))*dfft  ! 2 is introduced because half of the spectra is set as complex zero.



   if (is_verbose) then
      write(*,"(A,' bandpass filtering and fractional time correction is done ... ')") trim(adjustl(sacname))
      call flush(6)
   end if



   call sfftw_destroy_plan(planb)
   call sfftw_destroy_plan(planf)


   if (allocated(s)) then
      deallocate(s)
   end if
   if (allocated(sf)) then
      deallocate(sf)
   end if

end if


! =================================================================================================
! ======================================= Temporal Normalization ==================================
! =================================================================================================

! ***********************************************************************
! Perform one-bit normalization if is_onebit == .true.
! ***********************************************************************
if (is_onebit) then

   ! Apply one-bit normalization.
   !where (seis_data > 0.0)
   !   seis_data = +1.0
   !else where (seis_data < 0.0)
   !   seis_data = -1.0
   !end where
   seis_data = sign(1.0, seis_data)

   if (is_verbose) then
      write(*,*) trim(adjustl(sacname))//' one-bit normalization is done ... '
      call flush(6)
   end if

else

   ! ***********************************************************************
   ! Perform time domain running average if is_running_time_average is TRUE and is_onebit is FALSE.
   ! ***********************************************************************
   if (is_running_time_average) then

      allocate(abs_data(1:n), wgt_data(1:n), stat=ier)

      abs_data(1:n) = abs(seis_data(1:n))

      if (is_bandpass_earthquake) then
         norder = 2
         norder = 2*norder
         allocate(a(0:norder), b(0:norder), tr(1:n), stat=ier)
         call buttbp(norder, dt, fr1, fr2, a, b)
         tr(1:n) = seis_data(1:n)
         call filtfilt(norder, n-1, a, b, tr)
         abs_data(1:n) = abs(tr(1:n))
         deallocate(a, b, tr)
      end if

      do k = 1, n, 1
         n1 = max(1, k-nwt)
         n2 = min(n, k+nwt)
         wgt_data(k) = sum(abs_data(n1:n2)) / real(n2-n1+1)
      end do
      !wtr = 1.e-6*maxval(abs(wgt_data))
      !seis_data(1:n) = seis_data(1:n) / max(wgt_data(1:n), wtr)
      !seis_data(1:n) = min(seis_data(1:n)/wgt_data(1:n), hugeval)
      seis_data(1:n) = tan(atan2(seis_data(1:n), wgt_data(1:n)))

      deallocate(abs_data, wgt_data)

      if (is_verbose) then
         write(*,"(A)") trim(adjustl(sacname))//' time-domain running absolute average is done ... '
         call flush(6)
      end if

   end if

end if



! =================================================================================================
! ======================================= Cut data ================================================
! =================================================================================================

! Record info.
!dt = sdb%rec(ist,iev)%dt
!N = sdb%rec(ist,iev)%npts         ! N:    number of real data points
!Nlen = nint(sngl(tlen/dt))        ! Nlen: number of intercepted data points


! tend: desired ending time of the signal
tend = t0 + tlen


! trb: real data beginning time relative to the event time (1970_4_1_0_0_0) in seconds
trb = sdb%rec(ist,iev)%t0 - sdb%ev(iev)%t0


! tre: real data ending time relative to the event time (1970_4_1_0_0_0) in seconds
tre = trb + (n-1)*dt


! ***************************************************************
! If the real data beginning time is larger than t1 or the real data ending time
! is smaller than t2, the sac file will not be processed.
! ***************************************************************
if ((trb > t0) .or. (tre < tend)) then
   if (is_verbose) then
      write(*,"(A,A,A,F10.2,A,F10.2,A)") 'Short length file: ', trim(adjustl(sacname)), &
                                '  Beginning time:', trb, 's  Bad length:', (N-1)*dt, 's'
      call flush(6)
   end if
   call system('rm -rf '//trim(adjustl(sdb%rec(ist,iev)%name)))
   return
end if


!ngap = nint(sngl((t0 - trb)/dt))
ngap = nint((t0 - trb)/dt)


seis_data(1:Nlen) = seis_data(ngap+1:ngap+Nlen)
seis_data(Nlen+1:n) = 0.0


if (is_verbose) then
   write(*,"(A)") trim(adjustl(sacname))//' cutting data is done ... '
   call flush(6)
end if


if (maxval(abs(seis_data)) < tinyval) then
   call system('rm -rf '//trim(adjustl(sdb%rec(ist,iev)%name)))
   return
end if



! =================================================================================================
! ====================================== Forward FFT===============================================
! =================================================================================================

nfft = 2**ceiling(log(dble(2*Nlen-1))/log(2.d0))   ! nfft: number of points for FFT
nq = nfft/2 + 1
dfft = 1.0 / dble(nfft)


df = dfft / dt     ! df: frequency interval


! Allocate memory for s and sf.
allocate(s(nfft), sf(nfft), stat=ier)


call sfftw_plan_dft_1d(planf, nfft, s, sf, FFTW_FORWARD, FFTW_ESTIMATE)


! Initialize s with complex zero.
s = czero


! Apply taper
do k = 1, ntaper, 1
   coeff = 0.50*(1.0 + cos(PI*(ntaper-k+1)/dble(ntaper)))
   seis_data(k) = seis_data(k)*coeff
   seis_data(Nlen-k+1) = seis_data(Nlen-k+1)*coeff
end do


! Fill s with real data.
s(1:Nlen) = cmplx(seis_data(1:Nlen), 0.0)


! Make forward FFT for the seismogram: s => sf
call sfftw_execute_dft(planf, s, sf)


call sfftw_destroy_plan(planf)



! Kill half spectra
sf(nq+1:nfft) = czero


! Correct the ends
sf(1) = 0.50*sf(1)
!sf(nq) = cmplx(real(sf(nq)), 0.0)



! =================================================================================================
! ======================================= Spectral Whitening ======================================
! =================================================================================================

if (is_specwhitenning) then

   ! ***************************************************************
   ! Apply spectra whitening.
   ! ***************************************************************
   call whiten_spectra(f1, f4, df, nq, sf, nwf)


   ! ***************************************************************
   ! Reject the spike at the nperiod band [25s 27s].
   ! ***************************************************************
   if (is_suppress_notch) then
      call bandstop_filter(0.0350, 0.0360, 0.0390, 0.0400, df, nq, sf, npow_costaper, freqmin)
   end if


   if (is_verbose) then
      write(*,"(A)") trim(adjustl(sacname))//' spectral whitenning is done ... '
      call flush(6)
   end if

end if




! ***************************************************************
! Bandpass filtering.
! ***************************************************************
call bandpass_filter(f1, f2, f3, f4, df, nq, npow_costaper, sf)



! Destroy the original SAC file.
call system('rm -rf '//trim(adjustl(sdb%rec(ist,iev)%name)))


! ***********************************************************************
! Write the complex value of the FFT results to a local file, using the
! original sac name.
! ***********************************************************************
call write_bindata(trim(adjustl(sdb%rec(ist,iev)%name)), nq, sf(1:nq), ier)



if (is_verbose) then
   write(*,"(A)") trim(adjustl(sacname))//' fourier spectrum is done ... '
   call flush(6)
end if



! =================================================================================================
! =================================================================================================
! =================================================================================================



if (allocated(s)) then
   deallocate(s)
end if
if (allocated(sf)) then
   deallocate(sf)
end if
if (allocated(seis_data)) then
   deallocate(seis_data)
end if
if (allocated(strArray)) then
   deallocate(strArray)
end if


end subroutine preprocess



! =======================================================================================
! Spectra whitening algorithm. It works the same as running average amplitude in the
! time domain, and it is equivalent to do 'smooth mean h nwt' and 'divf avg.amp' in SAC.
! f1, f4: frequency band to do spectral whitening [input]
! df: frequency interval [input]
! nk: half-length of the data points in the frequency domain [input]
! sf: FFT values in complex form [input and output]
! nwf: half-window length in spectral whitening [input]
! =======================================================================================
subroutine whiten_spectra(f1, f4, df, nq, sf, nwf)

implicit none

integer, intent(in) :: nq, nwf

real(SGL), intent(in) :: f1, f4, df

complex(SGL), dimension(:), intent(inout) :: sf ! sf: assumed-shape dummy array


integer iw1, iw2
integer k, k1, k2

real(SGL) f, rsum, dw, wtr

real(SGL), dimension(:), allocatable :: sf_amp, sf_weight ! temporary arrays



! Return if 0 == nwf
if (0 == nwf) then
   write(*,"(A)") 'Error: nwf should be a positive integer !'
   call flush(6)
   return
end if


allocate(sf_amp(nq), sf_weight(nq))
! ***************************************************************
! compute the amplitude of the spectra and water level
! ***************************************************************
sf_amp(1:nq) = abs(sf(1:nq))
!wtr = 1.e-6*maxval(sf_amp(1:nq))


! ***************************************************************
! Loop on each frequency point
! ***************************************************************
sf_weight = 0.0
k1 = max(1, floor(f1/df)+1)
k2 = min(nq, ceiling(f4/df)+1)
do k = k1, k2, 1

   f = (k-1) * df      ! f: frequency value

   ! Only compute the weight at frequency band [f1 f4] with half-window length nwf.
   ! Set the weight at frequency band <f1 and >f4 to be zero.
   if ((f >= f1) .and. (f <= f4)) then
      iw1 = max(1 , k-nwf)
      iw2 = min(nq, k+nwf)
      dw = real(iw2 - iw1 + 1)
      rsum = sum(sf_amp(iw1:iw2))
      !sf_weight(k) = dw / max(rsum, wtr)
      !sf_weight(k) = min(dw/rsum, hugeval)
      sf_weight(k) = tan(atan2(dw, rsum))
   end if

end do


! ***************************************************************
! Obtain the whitened spectra (running averaged amplitude) at frequency band [f1 f4].
! Set the spectra at frequency band <f1 and >f4 to be zero.
! ***************************************************************
sf(1:nq) = sf(1:nq) * sf_weight(1:nq)


deallocate(sf_amp, sf_weight)


end subroutine whiten_spectra



! =======================================================================================
! Band-pass filtering computed in the frequency domain
! f1, f2, f3, f4: frequency limits [input]
! df: frequency interval [input]
! nk: half-length of the data points in the frequency domain [input]
! npow_costaper: power of the cosine tapering function [input]
! sf: FFT values in complex form [input and output]
! =======================================================================================
subroutine bandpass_filter(f1, f2, f3, f4, df, nq, npow_costaper, sf)

implicit none

integer, intent(in) :: nq, npow_costaper

real(SGL), intent(in) :: f1, f2, f3, f4, df

complex(SGL), dimension(:), intent(inout) :: sf


integer j, k, k1, k2

real(SGL) tmp, f

real(SGL), dimension(:), allocatable :: alpha



allocate(alpha(nq))


! Initialize alpha with zero
alpha = 0.0
k1 = max(1, floor(f1/df))
k2 = min(nq, ceiling(f4/df))
do k = k1, k2, 1

   f = (k-1)*df

   ! Keep alpha to be zero if f <= f1

   ! alpha = 0.50*(1.0 + cos(PI*(f2-f)/(f2-f1))) if f1 < f <= f2
   if ((f > f1) .and. (f <= f2)) then

      tmp = 1.0
      do j = 1, npow_costaper, 1
         tmp = tmp*0.50*(1.0 + cos(PI*(f2-f)/(f2-f1)))
      end do

      alpha(k) = tmp

   ! alpha = 1 if f2 < f <= f3
   else if (f <= f3) then

      alpha(k) = 1.0

   else if (f <= f4) then

      tmp = 1.0
      do j = 1, npow_costaper, 1
         tmp = tmp*0.50*(1.0 + cos(PI*(f-f3)/(f4-f3)))
      end do

      alpha(k) = tmp

   end if

   ! Keep alpha to be zero if f > f4

end do


! Apply the cosine tapering.
sf(1:nq) = sf(1:nq) * alpha(1:nq)


deallocate(alpha)


end subroutine bandpass_filter



! =======================================================================================
! Band-rejection filtering
! This function works just like the opposite of the band-pass filtering with
! two fliped consine taper function acting at [f1 f2] and [f3 f4], respectively.
! f1, f2, f3, f4: frequency limits [input]
! df: frequency interval [input]
! nk: half-length of the data points in the frequency domain [input]
! sf: FFT values in complex form [input and output]
! npow_costaper: power of the cosine tapering function [input]
! freqmin: retaining factor for the spectral whitening
! freqmin is the percentage (0.5 means 50%) of amplitude we try to retain
! =======================================================================================
subroutine bandstop_filter(f1, f2, f3, f4, df, nq, sf, npow_costaper, freqmin)


implicit none

integer, intent(in) :: nq, npow_costaper

real(SGL), intent(in) :: f1, f2, f3, f4
real(SGL), intent(in) :: freqmin

real(SGL), intent(in) :: df

complex(SGL), dimension(:), intent(inout) :: sf ! sf: assumed-shape dummy array


integer j, k, k1, k2

real(SGL) tmp, f

real(SGL), dimension(:), allocatable :: alpha



allocate(alpha(nq))


! Initialize alpha with 1
alpha = 1.0
k1 = max(1, floor(f1/df))
k2 = min(nq, ceiling(f4/df))
do k = 1, nq, 1

   f = (k-1)*df

   ! Keep alpha to be 1 if f <= f1.

   ! alpha = 0.50*(1.0 + cos(PI*(f-f1)/(f2-f1)))*(1.0 - freqmin) + freqmin if f1 < f <= f2
   if ((f > f1) .and. (f <= f2)) then

      tmp = 1.0
      do j = 1, npow_costaper, 1
         tmp = tmp*(0.50*(1.0 + cos(PI*(f-f1)/(f2-f1)))*(1.0 - freqmin) + freqmin)
      end do

      alpha(k) = tmp

   ! alpha = 1 if f2 < f <= f3
   else if (f <= f3) then

      alpha(k) = freqmin

   ! alpha = 0.50*(1.0 + cos(PI*(f-f1)/(f2-f1)))*(1.0 - freqmin) + freqmin if f3 < f <= f4
   else if (f <= f4) then

      tmp = 1.0
      do j = 1, npow_costaper, 1
         tmp = tmp*(0.50*(1.0 + cos(PI*(f4-f)/(f4-f3)))*(1.0 - freqmin) + freqmin)
      end do

      alpha(k) = tmp

  end if

  ! Keep alpha to be 1 if f >= f4

end do

! Apply the cosine tapering.
sf(1:nq) = sf(1:nq) * alpha(1:nq)


deallocate(alpha)


end subroutine bandstop_filter



! =======================================================================================
! Do the cross-correlation computation
! sdb: sac_db struct [input]
! nlag: lat time of the cross-correlation function [input]
! tarfolder: target folder to store the cross-correlation functions [input]
! N_bs: number of repeating times of the BOOTSTRAP method (e.g., 500)
! bs_type: which type does the BOOTSTRAP method apply to (e.g., 2_2)
! ist1, ist2: station indicies [input]
! myrank: process id number [input]
! is_verbose: verbose indicator [input]
! is_save_record: if output cross-correlation is_save_records [input]
! =======================================================================================
subroutine cc_and_aftan(sdb, ist1, ist2, nlag, N_bs, is_pws, str_pws, &
                   str_weight, str_per1, str_per2, bs_type, tarfolder)

implicit none

integer, intent(in) :: ist1, ist2, N_bs, nlag
logical, intent(in) :: is_pws

type(sac_db), intent(in) :: sdb

character(len=*), intent(in) :: tarfolder, bs_type
character(len=*), intent(in) :: str_per1, str_per2
character(len=*), intent(in) :: str_pws, str_weight



integer i, k, ii, jj, ier
integer iev, nev, nstack, nv
integer nrow1, nrow2, nout, nlen
integer nperiod, nperiod1, nperiod2

logical is_existed

real(SGL) dt, delta

real(SGL) groupV, phaseV

real(SGL) u_mean, u_std, c_mean, c_std

type(sachead) shd

character(len=128) str_myrank, str_stack

character(len=512) staname1, staname2
character(len=512) path, path_ls, path_pws
character(len=512) disp_name, bootstrap_name
character(len=512) stapair_path, stapair_name
character(len=512) sacname, listname, filename
character(len=512) sacfile_prefix, str_bootstrap

integer, allocatable, dimension(:) :: rand_array

real(SGL), allocatable, dimension(:) :: dataout

real(SGL), allocatable, dimension(:) :: tmpcorr, xcorr_bs

real(SGL), allocatable, dimension(:) :: grv_mean, grv_std, phv_mean, phv_std

real(SGL), allocatable, dimension(:,:) :: grv_2darr, phv_2darr

real(SGL), allocatable, dimension(:,:) :: matrix1, matrix2

real(DBL), allocatable, dimension(:) :: rand_tmp

complex(SGL), allocatable, dimension(:) :: fftdata1, fftdata2




! Initiate stacking number and cross-correlation function.
nev = sdb%nev
if (nev <= 0) return



! ***************************************************************
! Return if the corresponding dispersion file already exists.
! ***************************************************************
stapair_name = trim(adjustl(sdb%st(ist1)%n_name))//'_'//trim(adjustl(sdb%st(ist2)%n_name))
stapair_path = trim(adjustl(sdb%st(ist1)%n_name))//'/'//trim(adjustl(stapair_name))
path_ls = trim(adjustl(tarfolder))//'/FINAL/LINEAR/'//trim(adjustl(stapair_path))


if (is_pws) then
   path_pws = trim(adjustl(tarfolder))//'/FINAL/PWS/'//trim(adjustl(stapair_path))
   listname = trim(adjustl(path_pws))//'.dat'
   inquire(file=listname, exist=is_existed)
   if (is_existed) then
      write(*,"(A)") trim(adjustl(listname))//' exist, skip !'
      call flush(6)
      return
   end if
else
   listname = trim(adjustl(path_ls))//'.dat'
   inquire(file=listname, exist=is_existed)
   if (is_existed) then
      write(*,"(A)") trim(adjustl(listname))//' exist, skip !'
      call flush(6)
      return
   end if
end if


!! ======================================================================
!write(str_pws,'(I6)') ipws


if (is_stack) then
   ! Each process has its own process id
   write(str_myrank, '(I6.6)') myrank

   ! Create tmp directory to save single cross-correlation data.
   path = './tmp/'//trim(adjustl(str_myrank))
   call system('rm -rf '//trim(adjustl(path)))
   call system('mkdir '//trim(adjustl(path)))
else
   path = trim(adjustl(tarfolder))//'/CC_AFTAN/'//trim(adjustl(stapair_path))//'/prestack'
   !path = trim(adjustl(tarfolder))//'/CC_AFTAN/'//trim(adjustl(stapair_path))
   !call system('rm -rf '//trim(adjustl(path)))
   call system('mkdir -p '//trim(adjustl(path)))
end if
sacfile_prefix = trim(adjustl(path))//'/'//trim(adjustl(stapair_name))//'_'



! construct SAC header
call sacio_nullhead(shd)
shd%iztype = 11                 ! IO=11
shd%iftype = 1                  ! ITIME=1
shd%leven = 1                   ! TRUE=1
shd%npts = 2*nlag + 1
shd%o = 0.0
shd%evla = sdb%st(ist1)%lat
shd%evlo = sdb%st(ist1)%lon
shd%stla = sdb%st(ist2)%lat
shd%stlo = sdb%st(ist2)%lon
shd%kevnm = trim(adjustl(sdb%st(ist1)%name))
shd%kstnm = trim(adjustl(sdb%st(ist2)%name))
shd%kuser1 = trim(adjustl(sdb%st(ist1)%n_name))
shd%kuser2 = trim(adjustl(sdb%st(ist2)%n_name))
call geodist(shd)



allocate(tmpcorr(1:2*nlag+1))

nstack = 0
! Loop on the events.
do iev = 1, nev, 1

   staname1 = ''
   staname1 = trim(adjustl(sdb%rec(ist1,iev)%name))
   staname2 = ''
   staname2 = trim(adjustl(sdb%rec(ist2,iev)%name))
   ! ***************************************************************
   ! Check if there are FFT data for this station pair
   ! at this event. check_data is an internal procedure.
   ! ***************************************************************
   if (check_data(staname1, staname2)) then

      ! ***************************************************************
      ! read in the FFT data for the two stations.
      ! ***************************************************************

      call read_bindata(staname1, nlen, fftdata1, ier)
      call read_bindata(staname2, nlen, fftdata2, ier)

      ! ***************************************************************
      ! Compute the cross-correlation in the frequency domain.
      ! ***************************************************************
      call xcorr(fftdata1, fftdata2, nlen, dataout, nout)

      ! ***************************************************************
      ! Assign the cross-correlation resulted from frequency domain
      ! computation to the time domain cross-correlation series.
      ! ***************************************************************
      tmpcorr(nlag+1) = dataout(1)
      do k = 2, nlag+1, 1
         tmpcorr(nlag+2-k) = dataout(nout+2-k)
         tmpcorr(nlag+k) = dataout(k)
      end do

      ! **************************************************************
      ! Count the cross-correlation times.
      ! **************************************************************
      nstack = nstack + 1

      ! ***************************************************************
      ! Save single cross-correlation function into tmpfolder
      ! ***************************************************************
      write(str_stack,"(I6)") nstack
      sacname = trim(adjustl(sacfile_prefix))//trim(adjustl(str_stack))//'.SAC'


      if (1 == nstack) then
         ! Get the time interval.
         dt = nint(sdb%rec(ist1,iev)%dt*1e6)*1.d-6
		 shd%delta = dt
		 shd%b = -nlag*dt
		 shd%e =  nlag*dt

         !! Get the time interval.
         !dt = nint(sdb%rec(ist1,iev)%dt*1e6)*1.d-6
         !call sacio_newhead(shd, dt, 2*nlag+1, -nlag*dt)
         !shd%evla = sdb%st(ist1)%lat
         !shd%evlo = sdb%st(ist1)%lon
         !shd%stla = sdb%st(ist2)%lat
         !shd%stlo = sdb%st(ist2)%lon
         !shd%kevnm = trim(adjustl(sdb%st(ist1)%name))
         !shd%kstnm = trim(adjustl(sdb%st(ist2)%name))
         !shd%kuser1 = trim(adjustl(sdb%st(ist1)%n_name))
         !shd%kuser2 = trim(adjustl(sdb%st(ist2)%n_name))
         !call geodist(shd)
      end if


      shd%nzyear = sdb%ev(iev)%yy
	  shd%nzjday = date2jday(shd%ev(iev)%yy, shd%ev(iev)%mm, shd%ev(iev)%dd)
	  shd%nzhour = shd%ev(iev)%h
	  shd%nzmin  = shd%ev(iev)%m
	  shd%nzsec  = int(shd%ev(iev)%s)
	  shd%nzmsec = int((shd%ev(iev)%s - shd%nzsec)*1000)



      ! Write the single cross-correlation function.
      call sacio_writesac(sacname, shd, tmpcorr, ier)

   end if

end do




! ***************************************************************
! Write cross-correlation log if is_save_record is true.
! ***************************************************************
if (is_save_record) then
   write(str_stack,"(I6)") nstack
   call system('echo "'//trim(adjustl(stapair_name))//' '//trim(adjustl(str_stack))//'" | column -t >> CCRecord.lst')
end if



if ((0 == nstack) .or. (.not.(is_stack))) then
   if (allocated(fftdata1)) then
      deallocate(fftdata1)
   end if
   if (allocated(fftdata2)) then
       deallocate(fftdata2)
   end if
   if (allocated(tmpcorr)) then
       deallocate(tmpcorr)
   end if
   !call system('rm -rf '//trim(adjustl(tarfolder))//'/'//trim(adjustl(str_myrank)))
   if (is_verbose) then
      write(*,"(A)") 'Cross-correlation between '//trim(adjustl(sdb%st(ist1)%n_name))// &
                             ' and '//trim(adjustl(sdb%st(ist2)%n_name))//' is done ... '
      call flush(6)
   end if
   return
end if




! ***************************************************************
! After the event iteration, write the final cross-correlation to
! a local binary SAC file, measure the dispersion curves and do
! the bootstrap measurements [optional].
! ***************************************************************


if ((nstack > 0) .and. is_stack) then

   path = trim(adjustl(tarfolder))//'/CC_AFTAN/'//trim(adjustl(stapair_path))
   call system('mkdir -p '//trim(adjustl(path)))

   ! ***************************************************************
   ! Apply phase weighted stacking procedure, outputing both linear
   ! and phase weighted stacking final cross-correlations
   ! ***************************************************************
   sacfile_prefix = trim(adjustl(path))//'/'//trim(adjustl(stapair_name))
   call system('ls ./tmp/'//trim(adjustl(str_myrank))//'/*.SAC'// &
               ' | TF_PWS -B '//trim(adjustl(str_per2))//' -E '//trim(adjustl(str_per1))// &
               ' -W '//trim(adjustl(str_weight))//' -O '//trim(adjustl(sacfile_prefix))// &
               ' -P '//trim(adjustl(str_pws)))

   ! ***************************************************************
   ! Convert cross-correlation from sac format to ascii format
   ! in order to generate input files for Huajian Yao' code.
   ! ***************************************************************
   sacname = trim(adjustl(sacfile_prefix))//'_ls.SAC'
   call sac2asc(trim(adjustl(sacname)))


   sacname = trim(adjustl(sacfile_prefix))//'_pws.SAC'
   call sac2asc(trim(adjustl(sacname)))



   if (is_verbose) then
      write(*,"(A)") 'Stacking cross-correlation between '//trim(adjustl(sdb%st(ist1)%n_name))// &
                                      ' and '//trim(adjustl(sdb%st(ist2)%n_name))//' is done ... '
      call flush(6)
   end if



   if (.not.(is_onlycc)) then

      ! ***************************************************************
      ! Do the AFTAN for linear result.
      ! ***************************************************************
      sacname = trim(adjustl(sacfile_prefix))//'_ls.SAC'
      !call system('printf "r '//trim(adjustl(sacname))//'\nwh over\nq\n" | sac')


      ! Retrive distance header.
      !call sacio_readhead(sacname, shd, ier)
      delta = shd%dist

      call system('echo '//trim(adjustl(sacname))//' > '//trim(adjustl(str_myrank))//'.lst')
      call system('AFTAN '//trim(adjustl(str_myrank))//'.lst')


      if (is_pws) then
         ! ***************************************************************
         ! Do the AFTAN for PWS result.
         ! ***************************************************************
         sacname = trim(adjustl(sacfile_prefix))//'_pws.SAC'
         !call system('printf "r '//trim(adjustl(sacname))//'\nwh over\nq\n" | sac')

         ! Retrieve distance header.
         !call sacio_readhead(sacname, shd, ier)
         !delta = shd%dist

         call system('echo '//trim(adjustl(sacname))//' > '//trim(adjustl(str_myrank))//'.lst')
         call system('AFTAN '//trim(adjustl(str_myrank))//'.lst')
      end if

      call system('rm -rf '//trim(adjustl(str_myrank))//'.lst')


      if (is_pws) then
         ! ***************************************************************
         ! Write final dispersion data based on pws cross-correlation.
         ! ***************************************************************
         disp_name = trim(adjustl(sacfile_prefix))//'_pws_'//trim(adjustl(bs_type))//'.dat'
         inquire(file=disp_name, exist=is_existed)

         if (.not.(is_existed)) then

            if (is_verbose) then
               write(*,"(A)") '  NO final disperion data for pws cross-correlation !'
               call flush(6)
            end if

         else

            call system('mkdir -p '//trim(adjustl(tarfolder))//'/FINAL/PWS/'//trim(adjustl(sdb%st(ist1)%n_name)))

            listname = trim(adjustl(path_pws))//'.dat'
            open(unit=29, file=listname, status='replace', action='write', iostat=ier)

               write(29, "(A,2X,A)") trim(adjustl(sdb%st(ist1)%n_name)), trim(adjustl(sdb%st(ist2)%n_name))
               write(29, "(4F10.4,F14.4)") sdb%st(ist1)%lon, sdb%st(ist1)%lat, &
                                           sdb%st(ist2)%lon, sdb%st(ist2)%lat, delta
               write(29, "(A)") " Period  GroupV    PhaseV       SNR"
               call flush(29)

            close(unit=29)

            call system('cat '//trim(adjustl(disp_name))//' >> '//trim(adjustl(listname)))

         end if

      end if ! if (is_pws) then


      ! ***************************************************************
      ! Write final dispersion data based on linear stacking cross-correlation.
      ! ***************************************************************
      disp_name = trim(adjustl(sacfile_prefix))//'_ls_'//trim(adjustl(bs_type))//'.dat'

      inquire(file=disp_name, exist=is_existed)
      if (.not.(is_existed)) then
         if (is_verbose) then
            write(*,"(A)") '  NO final disperion data for linear stacking cross-correlation !'
            call flush(6)
         end if
         call system('rm -rf '//trim(adjustl(tarfolder))//'/'//trim(adjustl(str_myrank)))
         return
      end if



      if (.not.(is_sbs)) then

         call system('mkdir -p '//trim(adjustl(tarfolder))//'/FINAL/LINEAR/'//trim(adjustl(sdb%st(ist1)%n_name)))

         listname = trim(adjustl(path_ls))//'.dat'

         open(unit=30, file=listname, status='replace', action='write', iostat=ier)

            write(30, "(A,2X,A)") trim(adjustl(sdb%st(ist1)%n_name)), trim(adjustl(sdb%st(ist2)%n_name))
            write(30, "(4F10.4,F14.4)") sdb%st(ist1)%lon, sdb%st(ist1)%lat, &
                                        sdb%st(ist2)%lon, sdb%st(ist2)%lat, delta
            write(30, "(A)") " Period  GroupV    PhaseV       SNR"
            call flush(30)

         close(unit=30)

         call system('cat '//trim(adjustl(disp_name))//' >> '//trim(adjustl(listname)))

      else

         ! Allocate memory.
         allocate(xcorr_bs(1:2*nlag+1))
         allocate(rand_tmp(nstack), rand_array(nstack))
         allocate(grv_2darr(N_bs,100), phv_2darr(N_bs,100))

         ! Initialize the BOOTSTRAP matrix with zero
         grv_2darr = 0.0
         phv_2darr = 0.0



         ! ***************************************************************
         ! do the BOOTSTRAP
         ! ***************************************************************
         do i = 1, N_bs, 1

            xcorr_bs = 0.0

            ! Generate random integer data in [1,nstack]
            call init_random_seed()
            call random_number(rand_tmp)
            rand_tmp = rand_tmp*(nstack-1)+1
            rand_array = nint(rand_tmp)


            ! Stack selected daily cross-collelations
            do k = 1, nstack, 1

               write(sacname,"(I5)") rand_array(k)
               sacname = trim(adjustl(tarfolder))//'/'//trim(adjustl(str_myrank))//'/'// &
                          trim(adjustl(stapair_name))//'_'//trim(adjustl(sacname))//'.SAC'

               call sacio_readsac(sacname, shd, dataout, ier)
               tmpcorr = dataout
               xcorr_bs = xcorr_bs + tmpcorr

            end do


            ! Fill in the sac header.
            !call sacio_newhead(shd, dt, 2*nlag+1, -nlag*dt)
            !shd%evla = sdb%st(ist1)%lat
            !shd%evlo = sdb%st(ist1)%lon
            !shd%stla = sdb%st(ist2)%lat
            !shd%stlo = sdb%st(ist2)%lon
            !shd%kevnm = trim(adjustl(sdb%st(ist1)%name))
            !shd%kstnm = trim(adjustl(sdb%st(ist2)%name))
            !shd%kuser1 = trim(adjustl(sdb%st(ist1)%n_name))
            !shd%kuser2 = trim(adjustl(sdb%st(ist2)%n_name))
            shd%user0 = nstack

            sacname = trim(adjustl(tarfolder))//'/'//trim(adjustl(str_myrank))//'/'// &
                                                    trim(adjustl(stapair_name))//'.SAC'

            ! Write the bootstrap cross-correlation.
            call sacio_writesac(sacname, shd, xcorr_bs, ier)

            ! Update the sac header (e.g., delta).
            call system('printf "r '//trim(adjustl(sacname))//'\nwh over\nq\n" | sac')

            ! Do the AFTAN
            call system('echo '//trim(adjustl(sacname))//' > '//trim(adjustl(str_myrank))//'.lst')
            call system('AFTAN '//trim(adjustl(str_myrank))//'.lst')
            call system('rm -rf '//trim(adjustl(str_myrank))//'.lst')

            ! read the dispersion data file.
            filename = trim(adjustl(sacname))//'_'//trim(adjustl(bs_type))
            open(unit=25, file=filename, status='old', action='read', iostat=ier)

               ! Fill the BOOTSTRAP matrix with dispersion data.
               do
                   read(25, *, iostat=ier) nperiod, groupV, phaseV
                   if (0 /= ier) exit
                   grv_2darr(i, nperiod) = groupV
                   phv_2darr(i, nperiod) = phaseV
               end do

            close(unit=25)

         end do


         ! ***************************************************************
         ! Calculate the mean and standard deviation of the BOOTSTRAP measurements
         ! ***************************************************************
         call matrix_mean_std(grv_2darr, grv_mean, grv_std, 0.0, nv)
         call matrix_mean_std(phv_2darr, phv_mean, phv_std, 0.0, nv)


         path = trim(adjustl(tarfolder))//'/BOOTSTRAP/'//trim(adjustl(sdb%st(ist1)%n_name))
         call system('mkdir -p '//trim(adjustl(path)))

         str_bootstrap = trim(adjustl(path))//'/'//trim(adjustl(stapair_name))
         bootstrap_name = trim(adjustl(str_bootstrap))//'.dat'

         open(unit=26, file=bootstrap_name, status='replace', action='write', iostat=ier)

            do k = 1, nv, 1
               if ((grv_mean(k) > 0.0) .or. (grv_std(k) > 0.0) .or. (phv_mean(k) > 0.0) .or. (phv_std(k) > 0.0)) then
                  write(26, "(I4,4F12.6)") k, grv_mean(k), grv_std(k), phv_mean(k), phv_std(k)
               end if
            end do

            call flush(26)

         close(unit=26)

         deallocate(xcorr_bs)
         deallocate(grv_mean, grv_std)
         deallocate(phv_mean, phv_std)
         deallocate(rand_tmp, rand_array)
         deallocate(grv_2darr, phv_2darr)


         ! ***************************************************************
         ! Merge the dispersion and bootstrap data together.
         ! ***************************************************************
         nrow1 = 0       ! nrow1: number of rows of the dispersion file
         nrow2 = 0       ! nrow2: number of rows of the bootstrap file

         ! ***************************************************************
         ! Count the rows of the dispersion file and load the data.
         ! ***************************************************************
         open(unit=27, file=disp_name, status='old', action='read', iostat=ier)

            do
               read(27, *, iostat=ier)
               if (0 /= ier) exit
               nrow1 = nrow1 + 1
            end do

            if (allocated(matrix1)) then
               deallocate(matrix1)
            end if
            allocate(matrix1(nrow1,4))

            rewind(unit=27)

            read(27,*) ((matrix1(ii,jj), jj=1,4), ii=1,nrow1)

         close(unit=27)

         ! ***************************************************************
         ! Count the rows of the bootstrap file and load the data.
         ! ***************************************************************
         !bootstrap_name = trim(adjustl(str_bootstrap))//'.dat'
         open(unit=28, file=bootstrap_name, status='old', action='read', iostat=ier)

            do
               read(28, *, iostat=ier)
               if (0 /= ier) exit
               nrow2 = nrow2 + 1
            end do

            if (allocated(matrix2)) then
               deallocate(matrix2)
            end if
            allocate(matrix2(nrow2,5))

            rewind(unit=28)

            read(28,*) ((matrix2(ii,jj), jj=1,5), ii=1,nrow2)

         close(unit=28)


         ! ***************************************************************
         ! Write the final result.
         ! ***************************************************************
         call system('mkdir -p '//trim(adjustl(tarfolder))//'/FINAL/LINEAR/'//trim(adjustl(sdb%st(ist1)%n_name)))

         listname = trim(adjustl(path_ls))//'.dat'

         open(unit=29, file=listname, status='replace', action='write', iostat=ier)

            write(29, "(A,2X,A)") trim(adjustl(sdb%st(ist1)%n_name)), trim(adjustl(sdb%st(ist2)%n_name))
            write(29, "(4F10.4,F14.4)") sdb%st(ist1)%lon,sdb%st(ist1)%lat, &
                                     sdb%st(ist2)%lon,sdb%st(ist2)%lat,delta
            write(29, "(A)") " Period  GroupV     gMean     gStd     PhaseV     pMean     pStd      SNR"

            do ii = 1, nrow1, 1

               u_mean = 0.0
               u_std = 0.0
               c_mean = 0.0
               c_std = 0.0

               nperiod1 = int(matrix1(ii, 1))

               do jj = 1, nrow2, 1
                  nperiod2 = int(matrix2(jj, 1))
                  if (nperiod1 == nperiod2) then
                     u_mean = matrix2(jj, 2)
                     u_std = matrix2(jj, 3)
                     c_mean = matrix2(jj, 4)
                     c_std = matrix2(jj, 5)
                  end if
               end do

               write(29, "(I5,7F10.4)") nperiod1, matrix1(ii,2), u_mean, u_std, matrix1(ii,3), c_mean, c_std, matrix1(ii,4)

            end do

            call flush(29)

         close(unit=29)

         deallocate(matrix1, matrix2)

      end if ! if (.not. is_sbs) then

   end if ! if (.not.(is_onlycc)) then

end if ! if (nstack > 0) then


! Remove tmp directory.
call system('rm -rf ./tmp/'//trim(adjustl(str_myrank)))



if (allocated(tmpcorr)) then
   deallocate(tmpcorr)
end if
if (allocated(fftdata1)) then
   deallocate(fftdata1)
end if
if (allocated(fftdata2)) then
    deallocate(fftdata2)
end if
if (allocated(dataout)) then
    deallocate(dataout)
end if
if (allocated(rand_tmp)) then
    deallocate(rand_tmp)
end if
if (allocated(rand_array)) then
    deallocate(rand_array)
end if
if (allocated(grv_2darr)) then
    deallocate(grv_2darr)
end if
if (allocated(phv_2darr)) then
    deallocate(phv_2darr)
end if
if (allocated(matrix1)) then
    deallocate(matrix1)
end if
if (allocated(matrix2)) then
    deallocate(matrix2)
end if



contains



! ***************************************************************
! Internal procedure to check if there are FFT data for one
! station pair at one particular event.
! staname1, staname2: station directions [input]
! ***************************************************************
logical function check_data(staname1, staname2)

implicit none

character(len=*), intent(in) :: staname1, staname2

logical is_existed



check_data = .false.


if ((0 == len_trim(adjustl(staname1))) .or. (0 == len_trim(adjustl(staname2)))) then
   return
end if


inquire(file=trim(adjustl(staname1)), exist=is_existed)
if (.not.(is_existed)) return


inquire(file=trim(adjustl(staname2)), exist=is_existed)
if (.not.is_existed) return

check_data = .true.


end function check_data


end subroutine cc_and_aftan



! =======================================================================================
! Compute the cross-correlation in the frequency domain
! sf1: half-length FFT data at station 1 [input]
! sf2: half-length FFT data at station 2 [input]
! nlen: number of data points [input]
! dataout: output time domain cross-correlation data [output]
! nout: number of output time domain cross-correlation data [output]
! =======================================================================================
subroutine xcorr(sf1, sf2, nlen, dataout, nout)

implicit none

integer, intent(in) :: nlen

complex(SGL), dimension(nlen), intent(in) :: sf1, sf2


integer, intent(out) :: nout

real(SGL), allocatable, dimension(:), intent(out) :: dataout


integer k, ier

real(SGL) denom, wtr

type(c_ptr) plan

complex(SGL), allocatable, dimension(:) :: scorr, sfcorr



nout = (nlen-1)*2


allocate(dataout(nout), scorr(nout), sfcorr(nout), stat=ier)
if (0 /= ier) then
   write(*,"(A)") "Error: Allocating memory for dataout, scorr and sfcorr failed!"
   call flush(6)
   deallocate(dataout, scorr, sfcorr)
   return
end if


!wtr = 1.e-6*sum(abs(sf1(1:nlen))*abs(sf2(1:nlen)))/dble(nlen)
do k = 1, nlen, 1
   sfcorr(k) = conjg(sf1(k))*sf2(k)
   !sfcorr(k) = conjg(sf1(k))*sf2(k) / (abs(sf1(k))*abs(sf2(k)) + wtr)
end do


sfcorr(nlen+1:nout) = czero


! Make forward FFT for the cross-correlation: sfconj => s
call sfftw_plan_dft_1d(plan, nout, sfcorr, scorr, FFTW_BACKWARD, FFTW_ESTIMATE)
call sfftw_execute_dft(plan, sfcorr, scorr)
call sfftw_destroy_plan(plan)


dataout = 2.0*real(scorr)/real(nout)


if (allocated(scorr)) then
   deallocate(scorr)
end if
if (allocated(sfcorr)) then
   deallocate(sfcorr)
end if


end subroutine xcorr




subroutine geodist(shd)

implicit none

type(sachead), intent(inout) :: shd

real(8), parameter :: deg2rad = PI / 180.d0
real(8), parameter :: R = 6371.0

real(8) stla, stlo, evla, evlo, c, theta


evla = shd%evla * deg2rad
evlo = shd%evlo * deg2rad
stla = shd%stla * deg2rad
stlo = shd%stlo * deg2rad

c = sin(stla)*sin(evla) + cos(stla)*cos(evla)*cos(stlo - evlo)

if (abs(c - 1.0) < tinyval) then
   theta = 0.0
else if (abs(c + 1.0) > tinyval) then
   theta = PI
else
   theta = acos(c)
end if

shd%dist = R * theta


end subroutine geodist



end module xcc_m
