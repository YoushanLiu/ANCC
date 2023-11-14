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
real(DBL), parameter :: twoPI = 2.d0*PI



contains


! =======================================================================================
! Process the SAC file for one particular is_save_record and fill in the sdb.rec info.
! evtpath: event path containing the sac data [input]
! sdb: sac_db struct [input and output]
! iev: event struct [input]
! ist: station struct [input]
! channel: channel name [input]
! is_verbose: verbose indicator [input]
! =======================================================================================
subroutine mk_one_rec(sdb, iev, ist, npow_costaper, f1, f2, f3, f4, channel, evtpath, npts, dt, sdb_loc)

implicit none

integer, intent(in) :: iev, ist
integer, intent(in) :: npow_costaper

real(SGL), intent(in) :: f1, f2, f3, f4

character(len=*), intent(in) :: channel
character(len=*), intent(in) :: evtpath

type(sac_db), intent(in) :: sdb


integer, intent(out) :: npts
real(DBL), intent(out) :: dt

type(sac_db), intent(inout) :: sdb_loc


integer nstr, ier

logical is_existing

real(DBL) t0

character(len=512) str_date, sacname
character(len=512) sacinfile, sacoutfile

character(len=128), allocatable, dimension(:) :: strs



! ***************************************************************
! Copy one station SAC file from evtpath to target folder if it exists.
! ***************************************************************

sacname = trim(adjustl(sdb%st(ist)%ns_name))// &
             '.'//trim(adjustl(channel))//'.SAC'
sacinfile = trim(adjustl(evtpath))//'/'//trim(adjustl(sacname))
if (is_overwrite_data) then
   sacoutfile = trim(adjustl(sacinfile))
else
   sacoutfile = trim(adjustl(sdb_loc%ev(iev)%evtpath))//'/'//trim(adjustl(sacname))
end if


! Return if it doesn't exist.
inquire(file=sacinfile, exist=is_existing)
if (.not.(is_existing)) return



! ***************************************************************
if (is_verbose) then
   call split_string(evtpath, '/', nstr, strs)
   str_date = strs(nstr)
   write(*,"('Event: ',A,'   Station: ',A)") trim(adjustl(str_date)), &
                                     trim(adjustl(sdb%st(ist)%ns_name))
   call flush(6)
   if (allocated(strs)) then
      deallocate(strs)
   end if
end if



! ***************************************************************
! Correct the fraction time and set the reference time to be the beginning time (b=0)).
! ***************************************************************
call correct_sac_file(sacinfile, sacoutfile, sacname, npow_costaper, f1, f2, f3, f4, npts, dt, t0, ier)
if (0 /= ier) then
   write(*,"(A)") 'Error: Failed to correct fractional time for '//trim(adjustl(sacname))
   call flush(6)
   return
end if


! ***************************************************************
! Fill the elements in this is_save_record.
! ***************************************************************
sdb_loc%rec(ist,iev)%t0 = t0
sdb_loc%rec(ist,iev)%sacfile = ''
if (is_overwrite_data) then
   sdb_loc%rec(ist,iev)%sacfile = trim(adjustl(sacinfile))
else
   sdb_loc%rec(ist,iev)%sacfile = trim(adjustl(sacoutfile))
end if


end subroutine mk_one_rec



! =======================================================================================
! Correct the fraction time and set the reference time to be the beginning time (b=0)).
! if t=10.6s, then t=11s, frac=-0.4s; if t=10.4s, then t=10s, frac=0.4s
! fname: SAC filename [input]
! frac: fraction time of this SAC file  [output]
! npts: npts of this SAC header [output]
! dt: time sampling interval [output]
! t: beginning time of the first data point [output]
! ier: status indicator [output]
! =======================================================================================
subroutine correct_sac_file(finname, foutname, sacname, npow_costaper, f1, f2, f3, f4, npts, dt, t, ier)

implicit none

integer, intent(in) :: npow_costaper

real(SGL), intent(in) :: f1, f2, f3, f4

character(len=*), intent(in) :: finname, foutname, sacname


integer, intent(out) :: npts, ier

real(DBL), intent(out) :: dt, t


integer k, nf, nzsec

real(SGL) coeff
real(DBL) frac, tf, sec

type(sachead) shd

real(SGL), allocatable, dimension(:) :: seis_data



ier = -1
if (0 == len_trim(adjustl(finname))) return


! ***************************************************************
! Read the SAC file
call sacio_readsac(finname, shd, seis_data, ier)
if (0 /= ier) then
   write(*,*) "Error: Cannot read: "//trim(adjustl(finname))
   call flush(6)
   return
end if


! npts and dt
npts = shd%npts
! Remove round error
dt = nint(shd%delta*1e6)*1.d-6


! Get the initial beginning time of the first data point.
t = datetime2timestamp(shd%nzyear, shd%nzjday, shd%nzhour, shd%nzmin, shd%nzsec+1.d-3*shd%nzmsec+shd%b)


!! Apply taper
!do k = 1, ntaper, 1
!   coeff = 0.50*(1.0 + cos(PI*(ntaper-k+1)/ntaper))
!   seis_data(k) = seis_data(k)*coeff
!   seis_data(npts-k+1) = seis_data(npts-k+1)*coeff
!end do


! ***************************************************************
! Make the time fractional correction
! ***************************************************************
!tf = floor(t)
!nf = nint((t - tf)/dt)
!frac = t - (tf + nf*dt)
!t = tf + nf*dt
!!if (frac > 0.5*dt) then
!!   t = t + dt
!!   frac = frac - dt
!!end if
!tf = floor(t)
!t = tf + nint((t - tf)*1.d6)*1.d-6
!frac = nint(frac*1.d6)*1.d-6



tf = nint(t)
nf = nint((t - tf)/dt)
tf = tf + nf*dt
frac = t - tf
t = tf




! ***************************************************************
! Change the SAC header to make sure b=0
! ***************************************************************
call timestamp2datetime(t, shd%nzyear, shd%nzjday, shd%nzhour, shd%nzmin, sec)
sec = nint(sec*1e3)
nzsec = int(sec*0.001)
shd%nzsec = nzsec
shd%nzmsec = int(sec-nzsec*1000.0)
shd%b = 0.0
shd%user9 = sngl(frac)


! Correct fractional time
if (abs(frac) > 0.01*dt) then
   call correct_fractional_time(npts, npow_costaper, f1, f2, f3, f4, dt, frac, seis_data)
   if (is_verbose) then
      write(*,"(A)") trim(adjustl(sacname))//' fractional time correction and bandpass filtering is done ... '
      call flush(6)
   end if
end if



! Overwrite the SAC file
call sacio_writesac(foutname, shd, seis_data, ier)
if (0 /= ier) then
   write(*,*) "Error: Cannot overwrite: "//trim(adjustl(foutname))
   call flush(6)
   return
end if


deallocate(seis_data)


end subroutine correct_sac_file


! =======================================================================================
! =======================================================================================
! Peform time fraction correction and apply the Bandpass filtering
! sdb: sac_db struct [input]
! iev: event struct [input]
! ist: station struct [input]
! f1, f2, f3, f4: frequency limits [input]
! npow: power of cosine tapering function [input]
! isverbose: verbose indicator [input]
! =======================================================================================
! =======================================================================================
subroutine correct_fractional_time(npts, npow_costaper, f1, f2, f3, f4, dt, frac, seis_data)

implicit none

include 'fftw3.f'


integer, intent(in) :: npts, npow_costaper

real(SGL), intent(in) :: f1, f2, f3, f4

real(DBL), intent(in) :: dt, frac

real(SGL), intent(inout) :: seis_data(1:npts)


integer k, nfft, nq, ier

real(SGL) df, dfft, coeff

integer(8) fwd, bwd

complex(SGL), allocatable, dimension(:) :: s, sf


! Determine the power for FFT
nfft = 2**ceiling(log(dble(npts) + 16.0/dt)/log(2.d0))   ! nfft: number of points for FFT
nq = nfft/2 + 1
dfft = 1.0 / dble(nfft)

df = dfft / dt     ! df: frequency interval


! Allocate memory for s and sf.
allocate(s(nfft), sf(nfft), stat=ier)


call sfftw_plan_dft_1d(fwd, nfft, s, sf, FFTW_FORWARD, FFTW_ESTIMATE)
call sfftw_plan_dft_1d(bwd, nfft, sf, s, FFTW_BACKWARD, FFTW_ESTIMATE)



! Initialize s with complex zero.
s = czero


! Apply taper
do k = 1, ntaper, 1
   coeff = 0.50*(1.0 + cos(PI*(ntaper-k+1)/ntaper))
   seis_data(k) = seis_data(k)*coeff
   seis_data(npts-k+1) = seis_data(npts-k+1)*coeff
end do


! Fill s with real data.
s(1:npts) = cmplx(seis_data(1:npts), 0.0)


! Make forward FFT for the seismogram: s => sf
call sfftw_execute_dft(fwd, s, sf)


! Make time fraction correction
do k = 1, nq, 1
   !         tau > 0, left-shifted
   !         tau < 0, right-shifted
   !                   - j      omega       tau
   sf(k) = sf(k) * exp(-ci*(twoPI*(k-1)*df)*frac)
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
call sfftw_execute_dft(bwd, sf, s)


! Get the final result.
seis_data(1:npts) = 2.0*real(s(1:npts))*dfft  ! 2 is introduced because half of the spectra is set as complex zero.


call sfftw_destroy_plan(bwd)
call sfftw_destroy_plan(fwd)


if (allocated(s)) then
   deallocate(s)
end if
if (allocated(sf)) then
   deallocate(sf)
end if


end subroutine correct_fractional_time


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
integer nstr, ier

type(sachead) shd

character(len=512) str_date, sacname, sacfile

character(len=128), allocatable, dimension(:) :: strs




if (0 == len_trim(adjustl(filename))) return



! ***************************************************************
open(unit=17, file=filename, status='replace', action='write', iostat=ier)

   ! ***************************************************************
   ! Write the number of stations and events in the first line.
   ! ***************************************************************
   write(17,"(A,I8,5X,A,I6)") 'Number of events: ', sdb%nev, 'Number of stations: ', sdb%nst
   write(17,"(A)") '===================================================================='
   call flush(17)

   ! ***************************************************************
   ! Write the data is_save_record.
   ! ***************************************************************
   do iev = 1, sdb%nev, 1
      do ist = 1, sdb%nst, 1

         call split_string(sdb%ev(iev)%evtpath, '/', nstr, strs)
         str_date = strs(nstr)

         write(17,"(A20,$)") trim(adjustl(str_date))

         !if (0 == sdb%rec(ist,iev)%npts) then    ! Write "NO DATA" if rec[ist][iev]%npts == 0
         if (sdb%rec(ist,iev)%t0 < 0.0) then    ! Write "NO DATA" if rec[ist][iev]%t0 < 0

            write(17,"(A)") 'NO DATA at '//trim(adjustl(sdb%st(ist)%ns_name))

         else

			sacfile = trim(adjustl(sdb%rec(ist,iev)%sacfile))
            ! Read the SAC file to retrive its header information into shd struct.
            call sacio_readhead(sacfile, shd, ier)

            ! Write SAC file name, t0 (reference time), frac(time fraction),
            ! data length (npts*delta)
            call split_string(sacfile, '/', nstr, strs)
            sacname = strs(nstr)

            if (allocated(strs)) then
               deallocate(strs)
            end if

            !write(17,"(A30,3X,'t0: ',I4,'/',I3.3,'/',I2.2,':',I2.2,':',A,4X,'Frac:', &
            !     &F10.5,'s',4X,'Record Length:',F10.2,'s')") trim(adjustl(sacname)), &
            !                          shd%nzyear, shd%nzjday, shd%nzhour, shd%nzmin, &
            !                 trim(adjustl(padzero(shd%nzsec+0.001*shd%nzmsec,2,3))), &
            !                                 sdb%rec(ist,iev)%frac, shd%delta*shd%npts
            !write(17,"(A30,3X,'t0: ',I4,'/',I3.3,'/',I2.2,':',I2.2,':',F10.5,4X,'Frac:', &
            !         &F10.5,'s',4X,'Record Length:',F10.2,'s')") trim(adjustl(sacname)), &
            !                  shd%nzyear, shd%nzjday, shd%nzhour, shd%nzmin, shd%nzsec + &
            !               0.001d0*shd%nzmsec, sdb%rec(ist,iev)%frac, shd%delta*shd%npts
            write(17,"(A30,3X,'t0: ',I4,'/',I3.3,'/',I2.2,':',I2.2,':',F10.5,4X,'Frac:', &
                     &F10.5,'s',4X,'Record Length:',F10.2,'s')") trim(adjustl(sacname)), &
                              shd%nzyear, shd%nzjday, shd%nzhour, shd%nzmin, shd%nzsec + &
                                         0.001d0*shd%nzmsec, shd%user9, shd%delta*shd%npts

         end if

         call flush(17)

      end do
   end do

   if (allocated(strs)) then
      deallocate(strs)
   end if

close(unit=17)

end subroutine sacdb_to_asc



! =======================================================================================
! Remove the instrument response
! sdb: sac_db struct [input]
! iev: event struct [input]
! ist: station struct [input]
! myrank: process ID number [input]
! f1, f2, f3, f4: freqency limits [input]
! is_verbose: verbose indicator [input]
! =======================================================================================
subroutine remove_RESP(sdb, iev, ist, f1, f2, f3, f4, channel, pzfolder)

implicit none

integer, intent(in) :: iev, ist

real(SGL), intent(in) :: f1, f2, f3, f4

character(len=*), intent(in) :: channel
character(len=*), intent(in) :: pzfolder

type(sac_db), intent(in) :: sdb



logical is_existing

character(len=8) str_myrank

character(len=512) pzfile


!write(pzfile,"(A)") trim(adjustl(pzfolder))//'/'//trim(adjustl(sdb%st(ist)%ns_name))// &
!                                   '..'//trim(adjustl(sdb%rec(ist,iev)%channel))//'.PZ'
write(pzfile,"(A)") trim(adjustl(pzfolder))//'/'//trim(adjustl(sdb%st(ist)%ns_name))// &
                                                     '..'//trim(adjustl(channel))//'.PZ'

inquire(file=pzfile, exist=is_existing)
if (.not.(is_existing)) return


! ***************************************************************
if ((f1 > 0.0) .and. (f2 > f1) .and. (f3 > f2) .and. (f4 > f3)) then

   !if (sdb%rec(ist,iev)%npts > 0) then
   if (sdb%rec(ist,iev)%t0 > 0.0) then

      ! Each process has its own SAC script

      str_myrank = ''
      write(str_myrank, "(A, I6.6)") './tmp/', myrank
      open(unit=18, file=trim(adjustl(str_myrank))//'.sh', status='replace', action='write')
         write(18, "(A)") 'sac<<EOF'
         write(18, "(A)") 'r '//trim(adjustl(sdb%rec(ist,iev)%sacfile))
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
      ! Remove the instrument response to obtain the velocity measurement
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
! band-rejection filtering [optional], Bandpass filtering, cut data, and forward FFT
! sdb: sac_db struct [input]
! iev: event struct [input]
! ist: station struct [input]
! myrank: process ID number [input]
! is_running_time_average: time normalization indicator [input]
! is_onebit: is_onebit normalization indicator [input]
! is_suppress_notch: is_suppress_notch indicator [input]
! f1, f2, f3, f4: frequency limits [input]
! is_bandpass_earthquake: if earthquake band_pass filtering at [fr1 fr2] [input]
! fr1, fr2: period limits for earthquake Bandpass filtering in time normalization [input]
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
integer, intent(in) :: nwt, nwf
integer, intent(in) :: npow_costaper

real(SGL), intent(in) :: f1, f2, f3, f4
real(SGL), intent(in) :: fr1, fr2, freqmin

real(DBL), intent(in) :: t0, tlen

type(sac_db), intent(in) :: sdb



integer n, nfft, nq
integer k, nskip, norder
integer nstr, n1, n2, ier

logical is_existing

real(SGL) df, dfft, coeff, wtr

real(DBL) dt, tend
!real(DBL) dt, frac, tend
real(DBL) time, trb, tre, sec

integer(8) fwd, bwd

type(sachead) shd

character(len=512) sacname, sacfile

real(SGL), allocatable, dimension(:) :: seis_data
real(SGL), allocatable, dimension(:) :: abs_data, wgt_data

real(DBL), allocatable, dimension(:) :: a, b, tr

character(len=128), allocatable, dimension(:) :: strs

complex(SGL), allocatable, dimension(:) :: s, sf



! ***************************************************************
!if (0 == sdb%rec(ist,iev)%npts) return
sacfile = trim(adjustl(sdb%rec(ist,iev)%sacfile))
inquire(file=sacfile, exist=is_existing)
if (.not.(is_existing)) return



! Read the SAC file
call sacio_readsac(sacfile, shd, seis_data, ier)


call split_string(sacfile, '/', nstr, strs)
sacname = trim(adjustl(strs(nstr-1)))//'/'//trim(adjustl(strs(nstr)))




n = shd%npts
! Remove round error
dt = nint(shd%delta*1e6)*1.d-6



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

      if (is_bandpass_earthquake) then
         norder = 2
         norder = 2*norder
         allocate(a(0:norder), b(0:norder), tr(1:n), stat=ier)
         call buttbp(norder, dt, fr1, fr2, a, b)
         tr(1:n) = seis_data(1:n)
         call filtfilt(norder, n-1, a, b, tr)
         abs_data(1:n) = abs(tr(1:n))
         deallocate(a, b, tr)
      else
         abs_data(1:n) = abs(seis_data(1:n))
      end if

      do k = 1, n, 1
         n1 = max(1, k-nwt)
         n2 = min(n, k+nwt)
         wgt_data(k) = sum(abs_data(n1:n2)) / real(n2-n1+1)
      end do
      !wtr = 1.e-6*maxval(abs(wgt_data))
      !seis_data(1:n) = seis_data(1:n) / max(wgt_data(1:n), wtr)
      !seis_data(1:n) = min(seis_data(1:n)/wgt_data(1:n), HUGEVAL)
      !seis_data(1:n) = tan(atan2(seis_data(1:n), wgt_data(1:n)))
      wgt_data(1:n) = abs(tan(atan2(1.0, wgt_data(1:n))))
      seis_data(1:n) = seis_data(1:n) * wgt_data(1:n)

      !call hilbert(n, seis_data, wgt_data)
      !wgt_data(1:n) = abs(tan(atan2(1.0, wgt_data(1:n))))
      !seis_data(1:n) = seis_data(1:n) * wgt_data(1:n)

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
!Nlen = nint(tlen/dt)              ! Nlen: number of intercepted data points


! tend: desired ending time of the signal
tend = t0 + tlen


! trb: real data beginning time relative to the event time (1970_4_1_0_0_0) in seconds
trb = sdb%rec(ist,iev)%t0 - sdb%ev(iev)%t0


! tre: real data ending time relative to the event time (1970_4_1_0_0_0) in seconds
tre = trb + (n-1)*dt


! ***************************************************************
! If the real data beginning time is larger than t1 or the real data ending time
! is smaller than t2, the SAC file will not be processed.
! ***************************************************************
if ((trb > t0) .or. (tre < tend)) then
   if (is_verbose) then
      write(*,"(A,A,A,F10.2,A,F10.2,A)") 'Short length file: ', trim(adjustl(sacname)), &
                                '  Beginning time:', trb, 's  Bad length:', (N-1)*dt, 's'
      call flush(6)
   end if
   call system('rm -rf '//trim(adjustl(sacfile)))
   return
end if


nskip = nint((t0 - trb)/dt)


seis_data(1:Nlen) = seis_data(nskip+1:nskip+Nlen)
seis_data(Nlen+1:n) = 0.0


if (is_verbose) then
   write(*,"(A)") trim(adjustl(sacname))//' cutting data is done ... '
   call flush(6)
end if


if (maxval(abs(seis_data)) < TINYVAL) then
   call system('rm -rf '//trim(adjustl(sacfile)))
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


call sfftw_plan_dft_1d(fwd, nfft, s, sf, FFTW_FORWARD, FFTW_ESTIMATE)


! Initialize s with complex zero.
s = czero


! Apply taper
do k = 1, ntaper, 1
   coeff = 0.50*(1.0 + cos(PI*(ntaper-k+1)/ntaper))
   seis_data(k) = seis_data(k)*coeff
   seis_data(Nlen-k+1) = seis_data(Nlen-k+1)*coeff
end do


! Fill s with real data.
s(1:Nlen) = cmplx(seis_data(1:Nlen), 0.0)


! Make forward FFT for the seismogram: s => sf
call sfftw_execute_dft(fwd, s, sf)


call sfftw_destroy_plan(fwd)



! Kill half spectra
sf(nq+1:nfft) = czero


! Correct the ends
sf(1) = 0.50*sf(1)
!sf(nq) = cmplx(real(sf(nq)), 0.0)



! =================================================================================================
! ======================================= Spectral whitening ======================================
! =================================================================================================

if (is_specwhitening) then

   ! ***************************************************************
   ! Apply spectral whitening.
   ! ***************************************************************
   call whiten_spectra(f1, f4, df, nq, nwf, sf)


   if (is_verbose) then
      write(*,"(A)") trim(adjustl(sacname))//' spectral whitening is done ... '
      call flush(6)
   end if

end if



! ***************************************************************
! Reject the spike at the period band [25s 27s].
! ***************************************************************
if (is_suppress_notch) then
   call bandstop_filter(0.0350, 0.0360, 0.0390, 0.0400, df, nq, npow_costaper, freqmin, sf)
end if


! ***************************************************************
! Bandpass filtering.
! ***************************************************************
call bandpass_filter(f1, f2, f3, f4, df, nq, npow_costaper, sf)



! Destroy the original SAC file.
call system('rm -rf '//trim(adjustl(sacfile)))


! ***********************************************************************
! Write the complex value of the FFT results to a local file, using the
! original sac name.
! ***********************************************************************
call write_bindata(sacfile, nq, sf(1:nq), ier)



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
if (allocated(strs)) then
   deallocate(strs)
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
subroutine whiten_spectra(f1, f4, df, nq, nwf, sf)

implicit none

integer, intent(in) :: nq, nwf

real(SGL), intent(in) :: df
real(SGL), intent(in) :: f1, f4

complex(SGL), dimension(:), intent(inout) :: sf ! sf: assumed-shape dummy array


integer iw1, iw2
integer k, k1, k2

real(SGL) f, rsum, dw, wtr

real(SGL), dimension(:), allocatable :: sf_amp, sf_wgt ! temporary arrays



! Return if 0 == nwf
if (nwf < 0) then
   write(*,"(A)") 'Error: nwf should be a nonnegative integer !'
   call flush(6)
   return
end if


allocate(sf_amp(nq), sf_wgt(nq))
! ***************************************************************
! Compute the amplitude of the spectra and water level
! ***************************************************************
sf_amp(1:nq) = abs(sf(1:nq))
!wtr = 1.e-6*maxval(sf_amp(1:nq))


! ***************************************************************
! Loop on each frequency point
! ***************************************************************
sf_wgt = 0.0
k1 = max(1, floor(f1/df)+1)
k2 = min(nq, ceiling(f4/df)+1)
do k = k1, k2, 1

   f = (k-1) * df      ! f: frequency value

   ! Only compute the weight at frequency band [f1 f4] with half-window length nwf.
   ! Set the weight at frequency band <f1 and >f4 to be zero.
   if ((f >= f1) .and. (f <= f4)) then
      iw1 = max(1 , k-nwf)
      iw2 = min(nq, k+nwf)
      dw = real(iw2-iw1+1)
      rsum = sum(sf_amp(iw1:iw2))
      !sf_wgt(k) = dw / max(rsum, wtr)
      !sf_wgt(k) = min(dw/rsum, HUGEVAL)
      sf_wgt(k) = abs(tan(atan2(dw, rsum)))
   end if

end do


! ***************************************************************
! Obtain the whitenned spectra (running averaged amplitude) at frequency band [f1 f4].
! Set the spectra at frequency band <f1 and >f4 to be zero.
! ***************************************************************
sf(1:nq) = sf(1:nq) * sf_wgt(1:nq)


deallocate(sf_amp, sf_wgt)


end subroutine whiten_spectra



! =======================================================================================
! Bandpass filtering computed in the frequency domain
! f1, f2, f3, f4: frequency limits [input]
! df: frequency interval [input]
! nk: half-length of the data points in the frequency domain [input]
! npow_costaper: power of the cosine tapering function [input]
! sf: FFT values in complex form [input and output]
! =======================================================================================
subroutine bandpass_filter(f1, f2, f3, f4, df, nq, npow_costaper, sf)

implicit none

integer, intent(in) :: nq, npow_costaper

real(SGL), intent(in) :: df
real(SGL), intent(in) :: f1, f2, f3, f4

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
! This function works just like the opposite of the Bandpass filtering with
! two fliped consine taper function acting at [f1 f2] and [f3 f4], respectively.
! f1, f2, f3, f4: frequency limits [input]
! df: frequency interval [input]
! nk: half-length of the data points in the frequency domain [input]
! sf: FFT values in complex form [input and output]
! npow_costaper: power of the cosine tapering function [input]
! freqmin: retaining factor for the spectral whitening
! freqmin is the percentage (0.5 means 50%) of amplitude we try to retain
! =======================================================================================
subroutine bandstop_filter(f1, f2, f3, f4, df, nq, npow_costaper, freqmin, sf)


implicit none

integer, intent(in) :: nq, npow_costaper

real(SGL), intent(in) :: f1, f2, f3, f4
real(SGL), intent(in) :: df, freqmin

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
subroutine hilbert(nt, x, y)

implicit none

include 'fftw3.f'

integer(4), intent(in) :: nt

real(4), intent(in) :: x(1:nt)

real(4), intent(out) :: y(1:nt)


integer(4) it
integer(4) n, n2

integer(8) :: fwd = 0, bwd = 0

real(8) PI, dn

complex, dimension(:), allocatable :: ctrf, ctrb


n = 2**(ceiling(log10(dble(nt))/log10(2.d0)))
n2 = int(n/2)
dn = dble(n)


allocate(ctrf(1:n))
allocate(ctrb(1:n))


call sfftw_plan_dft_1d(fwd, n, ctrf, ctrf, FFTW_FORWARD , FFTW_ESTIMATE)
call sfftw_plan_dft_1d(bwd, n, ctrb, ctrb, FFTW_BACKWARD, FFTW_ESTIMATE)

ctrf(1:n) = czero
ctrf(1:nt) = cmplx(x(1:nt), 0.0)
call sfftw_execute_dft(fwd, ctrf, ctrf)

ctrb(1) = ctrf(1)
do it = 2, n2+1, 1
   ctrb(it) = 2.0*ctrf(it)
end do
do it = n2+2, n, 1
   ctrb(it) = czero
end do

call sfftw_execute(bwd, ctrb, ctrb)

ctrb(1:n) = ctrb(1:n) / dn

y(1:nt) = abs(ctrb(1:nt))


call sfftw_destroy_plan(bwd)
call sfftw_destroy_plan(fwd)


deallocate(ctrf, ctrb)


end subroutine hilbert



! =======================================================================================
! Do the cross-correlation computation
! sdb: sac_db struct [input]
! nlag: lat time of the cross-correlation function [input]
! tarfolder: target folder to store the cross-correlation functions [input]
! num_bootstrap: number of repeating times of the BOOTSTRAP method (e.g., 500)
! bootstrap_type: which type does the BOOTSTRAP method apply to (e.g., 2_2)
! ist1, ist2: station indicies [input]
! myrank: process ID number [input]
! is_verbose: verbose indicator [input]
! is_save_record: if output cross-correlation is_save_records [input]
! =======================================================================================
subroutine cc_and_aftan(sdb, ist1, ist2, nlag, num_bootstrap, is_pws, str_pws, &
                    str_npow_pws, str_per1, str_per2, bootstrap_type, tarfolder)

implicit none

integer, intent(in) :: ist1, ist2
integer, intent(in) :: num_bootstrap, nlag

logical, intent(in) :: is_pws

type(sac_db), intent(in) :: sdb

character(len=*), intent(in) :: str_per1, str_per2
character(len=*), intent(in) :: str_pws, str_npow_pws
character(len=*), intent(in) :: bootstrap_type, tarfolder



integer i, j, k
integer iev, ier
integer nlen, nout
integer nperiod, iperiod
integer iperiod1, iperiod2
integer nperiod1, nperiod2
integer nev, nstack, nzsec

logical is_existing

real(SGL) dist

real(SGL) groupV, phaseV

real(SGL) u_mean, u_std, c_mean, c_std

real(DBL) dt, sec

type(sachead) shd

character(len=8) str_myrank, str_stack

character(len=512) binfile1, binfile2
character(len=512) path, path_ls, path_pws
character(len=512) disp_name, bootstrap_name
character(len=512) stapair_path, stapair_name
character(len=512) sacname, listname, filename
character(len=512) sacfile_prefix, path_bootstrap

integer, allocatable, dimension(:) :: rand_array

real(SGL), allocatable, dimension(:) :: dataout

real(SGL), allocatable, dimension(:) :: tmpcorr, xcorr_bootstrap

real(SGL), allocatable, dimension(:) :: grv_mean, grv_std, phv_mean, phv_std

real(SGL), allocatable, dimension(:,:) :: grv_2darr, phv_2darr

real(SGL), allocatable, dimension(:,:) :: matrix1, matrix2

real(DBL), allocatable, dimension(:) :: rand_tmp

complex(SGL), allocatable, dimension(:) :: fftdata1, fftdata2




! Initialize stacking number and cross-correlation function.
nev = sdb%nev
dt = sdb%dt
if (nev <= 0) return



! ***************************************************************
! Return if the corresponding dispersion file already exists.
! ***************************************************************
stapair_name = trim(adjustl(sdb%st(ist1)%ns_name))//'_'//trim(adjustl(sdb%st(ist2)%ns_name))
stapair_path = trim(adjustl(sdb%st(ist1)%ns_name))//'/'//trim(adjustl(stapair_name))
path_ls = trim(adjustl(tarfolder))//'/FINAL/LINEAR/'//trim(adjustl(stapair_path))


if (is_pws) then
   path_pws = trim(adjustl(tarfolder))//'/FINAL/PWS/'//trim(adjustl(stapair_path))
   listname = trim(adjustl(path_pws))//'.dat'
   inquire(file=listname, exist=is_existing)
   if (is_existing) then
      write(*,"(A)") trim(adjustl(listname))//' exist, skip !'
      call flush(6)
      return
   end if
else
   listname = trim(adjustl(path_ls))//'.dat'
   inquire(file=listname, exist=is_existing)
   if (is_existing) then
      write(*,"(A)") trim(adjustl(listname))//' exist, skip !'
      call flush(6)
      return
   end if
end if



!! ======================================================================
if (is_stack) then
   ! Each process has its own process ID
   str_myrank = ''
   write(str_myrank, '(I6.6)') myrank

   ! Create tmp directory to save single cross-correlation data.
   path = './tmp/'//trim(adjustl(str_myrank))
   !call system('rm -rf '//trim(adjustl(path)))
   call system("perl -e 'for(<"//trim(path)//"/*>){unlink}'")
   call system('mkdir '//trim(adjustl(path)))
else
   path = trim(adjustl(tarfolder))//'/CC_AFTAN/'//trim(adjustl(stapair_path))//'/prestack'
   !path = trim(adjustl(tarfolder))//'/CC_AFTAN/'//trim(adjustl(stapair_path))
   !call system('rm -rf '//trim(adjustl(path)))
   call system("perl -e 'for(<"//trim(path)//"/*>){unlink}'")
   call system('mkdir -p '//trim(adjustl(path)))
end if
sacfile_prefix = trim(adjustl(path))//'/'//trim(adjustl(stapair_name))//'_'




! Construct SAC header
call sacio_nullhead(shd)
shd%iztype = 11               ! IO=11
shd%iftype = 1                ! ITIME=1
shd%leven = 1                 ! TRUE=1
shd%npts = 2*nlag + 1
shd%delta = dt
shd%b = -nlag*dt
shd%e =  nlag*dt
shd%o = 0.0
shd%evla = sdb%st(ist1)%lat
shd%evlo = sdb%st(ist1)%lon
shd%stla = sdb%st(ist2)%lat
shd%stlo = sdb%st(ist2)%lon
shd%kevnm = trim(adjustl(sdb%st(ist1)%staname))
shd%kstnm = trim(adjustl(sdb%st(ist2)%staname))
!call geodist(shd)
shd%dist = geodist(shd%evla, shd%evlo, shd%stla, shd%stlo)



allocate(tmpcorr(1:2*nlag+1))

nstack = 0
! Loop on the events.
do iev = 1, nev, 1

   binfile1 = ''
   binfile1 = trim(adjustl(sdb%rec(ist1,iev)%sacfile))
   binfile2 = ''
   binfile2 = trim(adjustl(sdb%rec(ist2,iev)%sacfile))
   ! ***************************************************************
   ! Check if there are FFT data for this station pair
   ! at this event. check_data is an internal procedure.
   ! ***************************************************************
   if (check_data(binfile1, binfile2)) then

      ! ***************************************************************
      ! Read in the FFT data for the two stations.
      ! ***************************************************************

      call read_bindata(binfile1, nlen, fftdata1, ier)
      call read_bindata(binfile2, nlen, fftdata2, ier)

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
      str_stack = ''
      write(str_stack,"(I0)") nstack
      sacname = trim(adjustl(sacfile_prefix))//trim(adjustl(str_stack))//'.SAC'



      !if (1 == nstack) then
      !   ! Get the time interval.
      !   !dt = nint(sdb%rec(ist1,iev)%dt*1e6)*1.d-6
      !   shd%delta = dt
      !   shd%b = -nlag*dt
      !   shd%e =  nlag*dt
      !   !! Get the time interval.
      !   !dt = nint(sdb%rec(ist1,iev)%dt*1e6)*1.d-6
      !   !call sacio_newhead(shd, dt, 2*nlag+1, -nlag*dt)
      !   !shd%evla = sdb%st(ist1)%lat
      !   !shd%evlo = sdb%st(ist1)%lon
      !   !shd%stla = sdb%st(ist2)%lat
      !   !shd%stlo = sdb%st(ist2)%lon
      !   !shd%kevnm = trim(adjustl(sdb%st(ist1)%staname))
      !   !shd%kstnm = trim(adjustl(sdb%st(ist2)%staname))
      !   !call geodist(shd)
      !end if



      if (.not.(is_stack)) then
         shd%nzyear = sdb%ev(iev)%yy
         shd%nzjday = date2jday(sdb%ev(iev)%yy, sdb%ev(iev)%mm, sdb%ev(iev)%dd)
         shd%nzhour = sdb%ev(iev)%h
         shd%nzmin  = sdb%ev(iev)%m
         sec = nint(sdb%ev(iev)%s*1e3)
         nzsec = int(sec*0.001)
         shd%nzsec = nzsec
         shd%nzmsec = int(sec-nzsec*1000.0)
      end if



      ! Write the single cross-correlation function.
      call sacio_writesac(sacname, shd, tmpcorr, ier)

   end if

end do




! ***************************************************************
! Write cross-correlation log if is_save_record is true.
! ***************************************************************
if (is_save_record) then
   str_stack = ''
   write(str_stack,"(I0)") nstack
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
   if (is_verbose .and. (nstack > 0)) then
      if (is_ac) then
         write(*,"(A)") 'Auto-correlation between '//trim(adjustl(sdb%st(ist1)%ns_name))// &
                               ' and '//trim(adjustl(sdb%st(ist2)%ns_name))//' is done ... '
      else
         write(*,"(A)") 'Cross-correlation between '//trim(adjustl(sdb%st(ist1)%ns_name))// &
                                ' and '//trim(adjustl(sdb%st(ist2)%ns_name))//' is done ... '
      end if
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
               ' -P '//trim(adjustl(str_pws))//' -W '//trim(adjustl(str_npow_pws))// &
               ' -O '//trim(adjustl(sacfile_prefix)))

   ! ***************************************************************
   ! Convert cross-correlation from sac format to ascii format
   ! in order to generate input files for Huajian Yao' code.
   ! ***************************************************************
   sacname = trim(adjustl(sacfile_prefix))//'_ls.SAC'
   call sac2asc(trim(adjustl(sacname)))


   sacname = trim(adjustl(sacfile_prefix))//'_pws.SAC'
   call sac2asc(trim(adjustl(sacname)))



   if (is_verbose) then
      if (is_ac) then
         write(*,"(A)") 'Stacking auto-correlation between '//trim(adjustl(sdb%st(ist1)%ns_name))// &
                                        ' and '//trim(adjustl(sdb%st(ist2)%ns_name))//' is done ... '
      else
         write(*,"(A)") 'Stacking cross-correlation between '//trim(adjustl(sdb%st(ist1)%ns_name))// &
                                         ' and '//trim(adjustl(sdb%st(ist2)%ns_name))//' is done ... '
      end if
      call flush(6)
   end if



   if (.not.(is_only_cf)) then

      ! ***************************************************************
      ! Do the AFTAN for linear result.
      ! ***************************************************************
      sacname = trim(adjustl(sacfile_prefix))//'_ls.SAC'
      !call system('printf "r '//trim(adjustl(sacname))//'\nwh over\nq\n" | sac')


      ! Retrive distance header.
      !call sacio_readhead(sacname, shd, ier)
      dist = shd%dist

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
         !dist = shd%dist

         call system('echo '//trim(adjustl(sacname))//' > '//trim(adjustl(str_myrank))//'.lst')
         call system('AFTAN '//trim(adjustl(str_myrank))//'.lst')
      end if

      call system('rm -rf '//trim(adjustl(str_myrank))//'.lst')


      if (is_pws) then
         ! ***************************************************************
         ! Write final dispersion data based on pws cross-correlation.
         ! ***************************************************************
         disp_name = trim(adjustl(sacfile_prefix))//'_pws_'//trim(adjustl(bootstrap_type))//'.dat'
         inquire(file=disp_name, exist=is_existing)

         if (.not.(is_existing)) then

            if (is_verbose) then
               write(*,"(A)") '  NO final disperion data for pws cross-correlation !'
               call flush(6)
            end if

         else

            call system('mkdir -p '//trim(adjustl(tarfolder))//'/FINAL/PWS/'//trim(adjustl(sdb%st(ist1)%ns_name)))

            listname = trim(adjustl(path_pws))//'.dat'
            open(unit=29, file=listname, status='replace', action='write', iostat=ier)

               write(29, "(A,2X,A)") trim(adjustl(sdb%st(ist1)%ns_name)), trim(adjustl(sdb%st(ist2)%ns_name))
               write(29, "(4F10.4,F14.4)") sdb%st(ist1)%lon, sdb%st(ist1)%lat, &
                                           sdb%st(ist2)%lon, sdb%st(ist2)%lat, dist
               write(29, "(A)") " Period  GroupV    PhaseV       SNR"
               call flush(29)

            close(unit=29)

            call system('cat '//trim(adjustl(disp_name))//' >> '//trim(adjustl(listname)))

         end if

      end if ! if (is_pws) then


      ! ***************************************************************
      ! Write final dispersion data based on linear stacking cross-correlation.
      ! ***************************************************************
      disp_name = trim(adjustl(sacfile_prefix))//'_ls_'//trim(adjustl(bootstrap_type))//'.dat'

      inquire(file=disp_name, exist=is_existing)
      if (.not.(is_existing)) then
         if (is_verbose) then
            write(*,"(A)") '  NO final disperion data for linear stacking cross-correlation !'
            call flush(6)
         end if
         call system('rm -rf '//trim(adjustl(tarfolder))//'/'//trim(adjustl(str_myrank)))
         return
      end if



      if (.not.(is_bootstrap)) then

         call system('mkdir -p '//trim(adjustl(tarfolder))//'/FINAL/LINEAR/'//trim(adjustl(sdb%st(ist1)%ns_name)))

         listname = trim(adjustl(path_ls))//'.dat'

         open(unit=30, file=listname, status='replace', action='write', iostat=ier)

            write(30, "(A,2X,A)") trim(adjustl(sdb%st(ist1)%ns_name)), trim(adjustl(sdb%st(ist2)%ns_name))
            write(30, "(4F10.4,F14.4)") sdb%st(ist1)%lon, sdb%st(ist1)%lat, &
                                        sdb%st(ist2)%lon, sdb%st(ist2)%lat, dist
            write(30, "(A)") " Period  GroupV    PhaseV       SNR"
            call flush(30)

         close(unit=30)

         call system('cat '//trim(adjustl(disp_name))//' >> '//trim(adjustl(listname)))

      else

         ! Allocate memory.
         allocate(xcorr_bootstrap(1:2*nlag+1))
         allocate(rand_tmp(nstack), rand_array(nstack))
         allocate(grv_2darr(200,num_bootstrap), phv_2darr(200,num_bootstrap))

         ! Initialize the BOOTSTRAP matrix with zero
         grv_2darr = 0.0
         phv_2darr = 0.0



         ! ***************************************************************
         ! Do the BOOTSTRAP
         ! ***************************************************************
         do i = 1, num_bootstrap, 1

            xcorr_bootstrap = 0.0

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
               xcorr_bootstrap = xcorr_bootstrap + tmpcorr

            end do


            ! Fill in the SAC header.
            !call sacio_newhead(shd, dt, 2*nlag+1, -nlag*dt)
            !shd%evla = sdb%st(ist1)%lat
            !shd%evlo = sdb%st(ist1)%lon
            !shd%stla = sdb%st(ist2)%lat
            !shd%stlo = sdb%st(ist2)%lon
            !shd%kevnm = trim(adjustl(sdb%st(ist1)%staname))
            !shd%kstnm = trim(adjustl(sdb%st(ist2)%staname))
            !shd%kuser1 = trim(adjustl(sdb%st(ist1)%ns_name))
            !shd%kuser2 = trim(adjustl(sdb%st(ist2)%ns_name))
            shd%user0 = nstack

            sacname = trim(adjustl(tarfolder))//'/'//trim(adjustl(str_myrank))//'/'// &
                                                    trim(adjustl(stapair_name))//'.SAC'

            ! Write the bootstrap cross-correlation.
            call sacio_writesac(sacname, shd, xcorr_bootstrap, ier)

            ! Update the SAC header (e.g., dist).
            !call system('printf "r '//trim(adjustl(sacname))//'\nwh over\nq\n" | sac')

            ! Do the AFTAN
            call system('echo '//trim(adjustl(sacname))//' > '//trim(adjustl(str_myrank))//'.lst')
            call system('AFTAN '//trim(adjustl(str_myrank))//'.lst')
            call system('rm -rf '//trim(adjustl(str_myrank))//'.lst')

            ! Read the dispersion data file.
            filename = trim(adjustl(sacname))//'_'//trim(adjustl(bootstrap_type))
            open(unit=25, file=filename, status='old', action='read', iostat=ier)

               ! Fill the BOOTSTRAP matrix with dispersion data.
               do
                   read(25, *, iostat=ier) iperiod, groupV, phaseV
                   if (0 /= ier) exit
                   grv_2darr(iperiod,i) = groupV
                   phv_2darr(iperiod,i) = phaseV
               end do

            close(unit=25)

         end do ! end of do i = 1, num_bootstrap, 1


         ! ***************************************************************
         ! Calculate the mean and standard deviation of the BOOTSTRAP measurements
         ! ***************************************************************
         call matrix_mean_std(grv_2darr, 0.0, nperiod, grv_mean, grv_std)
         call matrix_mean_std(phv_2darr, 0.0, nperiod, phv_mean, phv_std)


         path = trim(adjustl(tarfolder))//'/BOOTSTRAP/'//trim(adjustl(sdb%st(ist1)%ns_name))
         call system('mkdir -p '//trim(adjustl(path)))

         path_bootstrap = trim(adjustl(path))//'/'//trim(adjustl(stapair_name))
         bootstrap_name = trim(adjustl(path_bootstrap))//'.dat'

         open(unit=26, file=bootstrap_name, status='replace', action='write', iostat=ier)

            do iperiod = 1, nperiod, 1
               if ((grv_mean(iperiod) > 0.0) .or. (grv_std(iperiod) > 0.0) .or. (phv_mean(iperiod) > 0.0) .or. (phv_std(iperiod) > 0.0)) then
                  write(26, "(I4,4F12.6)") iperiod, grv_mean(iperiod), grv_std(iperiod), phv_mean(iperiod), phv_std(iperiod)
               end if
            end do

            call flush(26)

         close(unit=26)

         deallocate(xcorr_bootstrap)
         deallocate(grv_mean, grv_std)
         deallocate(phv_mean, phv_std)
         deallocate(rand_tmp, rand_array)
         deallocate(grv_2darr, phv_2darr)


         ! ***************************************************************
         ! Merge the dispersion and bootstrap data together.
         ! ***************************************************************
         nperiod1 = 0       ! nperiod1: number of rows of the dispersion file
         nperiod2 = 0       ! nperiod2: number of rows of the bootstrap file

         ! ***************************************************************
         ! Count the rows of the dispersion file and load the data.
         ! ***************************************************************
         open(unit=27, file=disp_name, status='old', action='read', iostat=ier)

            do
               read(27, *, iostat=ier)
               if (0 /= ier) exit
               nperiod1 = nperiod1 + 1
            end do

            if (allocated(matrix1)) then
              deallocate(matrix1)
            end if
            allocate(matrix1(4,nperiod1))

            rewind(unit=27)

            read(27,*) ((matrix1(i,j), i=1,4), j=1,nperiod1)

         close(unit=27)

         ! ***************************************************************
         ! Count the rows of the bootstrap file and load the data.
         ! ***************************************************************
         !bootstrap_name = trim(adjustl(path_bootstrap))//'.dat'
         open(unit=28, file=bootstrap_name, status='old', action='read', iostat=ier)

            do
               read(28, *, iostat=ier)
               if (0 /= ier) exit
               nperiod2 = nperiod2 + 1
            end do

            if (allocated(matrix2)) then
              deallocate(matrix2)
            end if
            allocate(matrix2(5,nperiod2))

            rewind(unit=28)

            read(28,*) ((matrix2(i,j), i=1,5), j=1,nperiod2)

         close(unit=28)


         ! ***************************************************************
         ! Write the final result.
         ! ***************************************************************
         call system('mkdir -p '//trim(adjustl(tarfolder))//'/FINAL/LINEAR/'//trim(adjustl(sdb%st(ist1)%ns_name)))

         listname = trim(adjustl(path_ls))//'.dat'

         open(unit=29, file=listname, status='replace', action='write', iostat=ier)

            write(29, "(A,2X,A)") trim(adjustl(sdb%st(ist1)%ns_name)), trim(adjustl(sdb%st(ist2)%ns_name))
            write(29, "(4F10.4,F14.4)") sdb%st(ist1)%lon, sdb%st(ist1)%lat, &
                                        sdb%st(ist2)%lon, sdb%st(ist2)%lat, dist
            write(29, "(A)") " Period  GroupV     gMean     gStd     PhaseV     pMean     pStd      SNR"

            do i = 1, nperiod1, 1

               u_mean = 0.0
               u_std = 0.0
               c_mean = 0.0
               c_std = 0.0

               iperiod1 = int(matrix1(1,i))

               do j = 1, nperiod2, 1
                  iperiod2 = int(matrix2(1,j))
                  if (iperiod1 == iperiod2) then
                     u_mean = matrix2(2,j)
                     u_std = matrix2(3,j)
                     c_mean = matrix2(4,j)
                     c_std = matrix2(5,j)
                  end if
               end do

               write(29, "(I5,7F10.4)") iperiod1, matrix1(2,i), u_mean, u_std, matrix1(3,i), c_mean, c_std, matrix1(4,i)

            end do

            call flush(29)

         close(unit=29)

         deallocate(matrix1, matrix2)

      end if ! if (.not.(is_bootstrap)) then

   end if ! if (.not.(is_only_cf)) then

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
   deallocate(dataout)
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
! binfile1, binfile2: station files [input]
! ***************************************************************
logical function check_data(binfile1, binfile2)

implicit none

character(len=*), intent(in) :: binfile1, binfile2

logical is_existing



check_data = .false.


if ((0 == len_trim(adjustl(binfile1))) .or. (0 == len_trim(adjustl(binfile2)))) then
   return
end if


inquire(file=binfile1, exist=is_existing)
if (.not.(is_existing)) return


inquire(file=binfile2, exist=is_existing)
if (.not.(is_existing)) return

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

integer(8) plan

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




real function geodist(evlaf, evlof, stlaf, stlof)

real, intent(in) :: stlaf, stlof, evlaf, evlof


real(8), parameter :: deg2rad = PI / 180.d0
real(8), parameter :: R = 6371.0

real(8) stla, stlo, evla, evlo, c, theta


evla = evlaf * deg2rad
evlo = evlof * deg2rad
stla = stlaf * deg2rad
stlo = stlof * deg2rad

c = sin(stla)*sin(evla) + cos(stla)*cos(evla)*cos(stlo - evlo)

if (abs(c - 1.0) < TINYVAL) then
   theta = 0.0
else if (abs(c + 1.0) < TINYVAL) then
   theta = PI
else
   theta = acos(c)
end if


geodist = R * theta


return


end function geodist



!subroutine geodist(shd)

!implicit none

!type(sachead), intent(inout) :: shd


!real(8), parameter :: deg2rad = PI / 180.d0
!real(8), parameter :: R = 6371.0

!real(8) stla, stlo, evla, evlo, c, theta


!evla = shd%evla * deg2rad
!evlo = shd%evlo * deg2rad
!stla = shd%stla * deg2rad
!stlo = shd%stlo * deg2rad


!c = sin(stla)*sin(evla) + cos(stla)*cos(evla)*cos(stlo - evlo)

!if (abs(c - 1.0) < TINYVAL) then
!   theta = 0.0
!else if (abs(c + 1.0) < TINYVAL) then
!   theta = PI
!else
!   theta = acos(c)
!end if

!shd%dist = R * theta


!end subroutine geodist



end module xcc_m
