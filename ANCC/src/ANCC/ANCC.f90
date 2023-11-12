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
program ANCC

use mpi
use xcc_m         ! xcc_m: imported module which contains all the functions, subroutines
                  ! and also other imported modules (e.g. my_definition_m, string_m, math_m, date_time_m).
use db_m, only: myrank, nprocs


implicit none


! ***********************************************************************
! Variable declaration section.
! ***********************************************************************
character(64) :: version = ' (v6.6)'

integer i, j, k, n
integer jmax, nlag
integer iunit, ipws
integer nargin, myroot, ier
integer iev, ist, ist1, ist2
integer nstrArray, num_bootstrap
integer nev, nst, nev_loc, nev_gathered
integer npow_costaper, nwt, nwf, npow_pws
integer npts, npts_read, npts_min, npts_max
integer event_type, station_type, record_type

integer(2) errcode


integer(8) ndim, ndim1, ndim2, idim1, iproc
integer(4) nstxnev, nstxnev_loc, nstxnev_gathered


real(SGL) lat, lon, freqmin, tlag
real(SGL) f1, f2, f3, f4, fr1, fr2

real(DBL) dt, dt_read, dt_min, dt_max, t0, tlen


character(len=3) str_specwhitenning
character(len=3) str_overwrite_data
character(len=3) str_stack, str_pws
character(len=3) str_only_cf, str_ac
character(len=3) str_verbose, str_save_record
character(len=3) str_bootstrap, bootstrap_type
character(len=3) str_onebit, str_running_time_average
character(len=3) str_bandpass_earthquake, str_suppress_notch

character(len=8) netname, staname, channel

character(len=32) str_npow_pws
character(len=32) str_per1, str_per4

character(len=512) sacfile, msg
character(len=512) evtpath, str_date, path
character(len=512) sacfolder, pzfolder, tarfolder


type(event) evt_node
type(station) sta_node
type(record) rec_node
type(sac_db) sdb, sdb_loc


integer, allocatable, dimension(:) :: blocklen, types
integer(MPI_ADDRESS_KIND), allocatable, dimension(:) :: base, disp

integer, allocatable, dimension(:) :: recvns, displs

character(len=128), allocatable, dimension(:) :: strArray




! ***********************************************************************
! Initialize MPI.
! ***********************************************************************
call MPI_INIT(ier)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ier)
myroot = nprocs - 1


! ***********************************************************************
! Construct new event data type for later data sharing.
! ***********************************************************************
n = 4
allocate(base(n), disp(n), blocklen(n), types(n))
blocklen = (/ 128, 6, 1, 1 /)
types =(/ MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION /)
call MPI_GET_ADDRESS(evt_node%evtpath, base(1), ier)
call MPI_GET_ADDRESS(evt_node%yy, base(2), ier)
call MPI_GET_ADDRESS(evt_node%s, base(3), ier)
call MPI_GET_ADDRESS(evt_node%t0, base(4), ier)
disp(1) = 0
disp(2) = base(2) - base(1)
disp(3) = base(3) - base(1)
disp(4) = base(4) - base(1)
call MPI_TYPE_CREATE_STRUCT(n, blocklen, disp, types, event_type, ier)
call MPI_TYPE_COMMIT(event_type, ier)
deallocate(base, disp, blocklen, types)


! ***********************************************************************
! Construct new station data type for later data sharing.
! ***********************************************************************
n = 3
allocate(base(n), disp(n), blocklen(n), types(n))
blocklen = (/ 8, 17, 2 /)
types =(/ MPI_CHARACTER, MPI_CHARACTER, MPI_REAL /)
call MPI_GET_ADDRESS(sta_node%staname, base(1), ier)
call MPI_GET_ADDRESS(sta_node%ns_name, base(2), ier)
call MPI_GET_ADDRESS(sta_node%lat, base(3), ier)
disp(1) = 0
disp(2) = base(2) - base(1)
disp(3) = base(3) - base(1)
call MPI_TYPE_CREATE_STRUCT(n, blocklen, disp, types, station_type, ier)
call MPI_TYPE_COMMIT(station_type, ier)
deallocate(base, disp, blocklen, types)


! ***********************************************************************
! Construct new record data type for later data sharing.
! ***********************************************************************
n = 2
allocate(base(n), disp(n), blocklen(n), types(n))
blocklen = (/ 192, 1 /)
types =(/ MPI_CHARACTER, MPI_DOUBLE_PRECISION /)
call MPI_GET_ADDRESS(rec_node%sacfile, base(1), ier)
call MPI_GET_ADDRESS(rec_node%t0, base(2), ier)
disp(1) = 0
disp(2) = base(2) - base(1)
call MPI_TYPE_CREATE_STRUCT(n, blocklen, disp, types, record_type, ier)
call MPI_TYPE_COMMIT(record_type, ier)
deallocate(base, disp, blocklen, types)





!if (myrank == myroot) then

! Parse the number of command line argumet(s).
nargin = command_argument_count()

! Stop if the number of input argument(s) is wrong.
if (3 /= nargin) then
   write(*,"(A)") 'ANCC' // trim(version)
   write(*,"(A)") "Usage: ANCC sacfolder pzfolder tarfolder"
   call flush(6)
   call MPI_FINALIZE(ier)
   stop
end if


! =====================================================================================
! =============================== SECTION 1 BEGINS ====================================
! =====================================================================================
! This section parse the input parameters.
if (myrank == myroot) then
   write(*,"(A)")
   write(*,"(A)") 'This program computes cross/auto-correlation and/or does AFTAN' // trim(version)
   write(*,"(A)") 'Its efficiency has been improved significantly by removing any unneccessary '
   write(*,"(A)") 'MPI_SEND & MPI_RECV and paralleling all parts by Youshan-Liu'
   write(*,"(A)") 'All processors are used to compute instead of the master processor only for '
   write(*,"(A)") 'message passing just as those original version done'
   write(*,"(A)")
   write(*,"(A)") '***********************************************************************'
   write(*,"(A)") '                         SECTION 1 BEGINS'
   write(*,"(A)") '***********************************************************************'
   write(*,"(A)") 'Reading input parameters ...'
   write(*,"(A)") '***********************************************************************'
   call flush(6)
end if

call MPI_BARRIER(MPI_COMM_WORLD, ier)





! Parse the command line input argument.
call get_command_argument(1, sacfolder)
call get_command_argument(2, pzfolder)
call get_command_argument(3, tarfolder)



! ***********************************************************************
! Read the parameters from the 'input.dat' file.
! ***********************************************************************
iunit = myrank + 11
open(unit=iunit, file='input.dat', status='old', action='read', iostat=ier)

   if (0 /= ier) then
      write(*,"(A)") 'Error: Cannot open input.dat ! '
      call flush(6)
      close(unit=iunit)
      call MPI_FINALIZE(ier)
      stop
   end if

   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   read(iunit,*) channel
   read(iunit,*) f1, f2, f3, f4
   read(iunit,*) t0, tlen
   read(iunit,*) npow_costaper
   read(iunit,*) str_running_time_average, nwt, str_bandpass_earthquake, fr1, fr2
   read(iunit,*) str_onebit
   read(iunit,*) str_specwhitenning, nwf
   read(iunit,*) str_suppress_notch, freqmin
   read(iunit,*) tlag
   read(iunit,*) str_pws, npow_pws
   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   read(iunit,*) str_bootstrap, num_bootstrap
   read(iunit,*) bootstrap_type
   do i = 1, 20, 1
      read(iunit,*)
   end do
   read(iunit,*) str_verbose
   read(iunit,*) str_overwrite_data
   read(iunit,*) str_save_record
   read(iunit,*) str_only_cf
   read(iunit,*) str_stack
   read(iunit,*) str_ac

close(unit=iunit)



! =======================================================================
! Check parameter validation
if (nwt < 0) then
   call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
   call MPI_FINALIZE(ier)
   stop 'Error: The half window width [nwt] for temporal normalization must be nonnegative integer!'
end if

if (nwf < 0) then
   call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
   call MPI_FINALIZE(ier)
   stop 'Error: The half window width [nwf] for spectral whitenning must be nonnegative integer!'
end if

if (.not.((f1 > f2) .and. (f2 > f3) .and. (f3 > f4))) then
   !write(*,"(A)") 'Error: periods f1 > f2 > f3 > f4 [sec.] should be satisfied !'
   !call flush(6)
   call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
   call MPI_FINALIZE(ier)
   stop 'Error: periods f1 > f2 > f3 > f4 [sec.] should be satisfied !'
end if
! =======================================================================



! ***********************************************************************
! Convert the periods to frequencies.
! ***********************************************************************
write(str_per1,'(F8.2)') f1
write(str_per4,'(F8.2)') f4
f1 = 1.0/f1
f2 = 1.0/f2
f3 = 1.0/f3
f4 = 1.0/f4
fr1 = 1.0/fr1
fr2 = 1.0/fr2


! ***********************************************************************
! Initialize the logical variables.
! ***********************************************************************
is_running_time_average = .false.
is_bandpass_earthquake = .false.
is_onebit = .false.
is_specwhitenning = .false.
is_suppress_notch = .false.
ipws = 0
is_bootstrap = .false.
is_verbose = .false.
is_overwrite_data = .false.
is_save_record = .false.
is_only_cf = .false.
is_stack = .false.
is_ac = .false.
if ((str_running_time_average == 'Y') .or. (str_running_time_average == 'y')) is_running_time_average = .true.
if ((str_bandpass_earthquake == 'Y') .or. (str_bandpass_earthquake == 'y')) is_bandpass_earthquake = .true.
if ((str_onebit == 'Y') .or. (str_onebit == 'y')) is_onebit = .true.
if ((str_specwhitenning == 'Y') .or. (str_specwhitenning == 'y')) is_specwhitenning = .true.
if ((str_suppress_notch == 'Y') .or. (str_suppress_notch == 'y')) is_suppress_notch = .true.
if ((str_bootstrap == 'Y') .or. (str_bootstrap == 'y')) is_bootstrap = .true.
if ((str_pws == 'Y') .or. (str_pws == 'y')) ipws = 1
if ((str_verbose == 'Y') .or. (str_verbose == 'y')) is_verbose = .true.
if ((str_overwrite_data == 'Y') .or. (str_overwrite_data == 'y')) is_overwrite_data = .true.
if ((str_save_record == 'Y') .or. (str_save_record == 'y')) is_save_record = .true.
if ((str_only_cf == 'Y') .or. (str_only_cf == 'y')) is_only_cf = .true.
if ((str_stack == 'Y') .or. (str_stack == 'y')) is_stack = .true.
if ((str_ac == 'Y') .or. (str_ac == 'y')) is_ac = .true.
if (myrank == myroot) then
   write(*,"(A,/)") 'Reading input parameters is done ... '
   call flush(6)
end if



! =======================================================================
! Check parameter validation
if (.not.(is_only_cf) .and. (.not.(is_stack))) then
   write(*,"(A)") 'Error: Dispersion analysis cannot be done for prestack cross-correlation functions !'
   call flush(6)
   call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
   call MPI_FINALIZE(ier)
   stop
end if

if (is_bandpass_earthquake .and. (fr1 >= fr2)) then
   write(*,"(A)") 'Error: fr1 < fr2 should be satisfied !'
   call flush(6)
   call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
   call MPI_FINALIZE(ier)
   stop
end if
! =======================================================================



! ***********************************************************************
! Clear the [DATA] folder inside the target folder.
! ***********************************************************************
!!call system('rm -rf '//trim(adjustl(tarfolder)))
if (.not.(is_overwrite_data)) then
   call system('rm -rf '//'./tmp/DATA')
end if
! ***********************************************************************
! Clear the [tmp] folder in current folder.
! ***********************************************************************
call system('rm -rf '//'./tmp')
! ***********************************************************************
! Create the [tmp] folder in current folder.
! ***********************************************************************
call system('mkdir -p '//'./tmp')

call MPI_BARRIER(MPI_COMM_WORLD, ier)




! =====================================================================================
! =============================== SECTION 2 BEGINS ====================================
! =====================================================================================
! This section processes the SAC files and fill in the elements in the sdb struct.
if (myrank == myroot) then
   write(*,"(A)") '***********************************************************************'
   write(*,"(A)") '                         SECTION 2 BEGINS'
   write(*,"(A)") '***********************************************************************'
   write(*,"(A)") 'Constructing sdb struct ...'
   write(*,"(A)") '***********************************************************************'
   call flush(6)
end if

call MPI_BARRIER(MPI_COMM_WORLD, ier)




! ***********************************************************************
! Count the number of stations.
! ***********************************************************************
iunit = nprocs + myrank + 11
open(unit=iunit, file='stations.lst', status='old', action='read', iostat=ier)
   if (0 /= ier) then
      write(*,"(A)") 'Error: Cannot open stations.lst !'
      call flush(6)
      call MPI_FINALIZE(ier)
   close(unit=11)
   stop
   end if
   nst = 0
   do
      read(iunit, *, iostat=ier)
      if (0 /= ier) exit
      nst = nst + 1
   end do
close(unit=iunit)



! ***********************************************************************
! Count the number of events.
! ***********************************************************************
iunit = 2*nprocs + myrank + 11
open(unit=iunit, file='events.lst', status='old', action='read', iostat=ier)
   if (0 /= ier) then
      write(*,"(A)") 'Error: Cannot open events.lst !'
      call flush(6)
      call MPI_FINALIZE(ier)
      close(unit=iunit)
      stop
   end if
   nev = 0
   do
      read(iunit, *, iostat=ier)
      if (0 /= ier) exit
      nev = nev + 1
   end do
close(unit=iunit)

!end if  ! if (myrank == myroot) then
nstxnev = nst*nev



!! ***********************************************************************
!! Broadcast input parameters from the master processor to all other processors.
!! ***********************************************************************
!call MPI_BCAST(sacfolder, 512, MPI_CHARACTER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(pzfolder, 512, MPI_CHARACTER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(tarfolder, 512, MPI_CHARACTER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(nev, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(nst, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(f1, 1, MPI_REAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(f2, 1, MPI_REAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(f3, 1, MPI_REAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(f4, 1, MPI_REAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(t0, 1, MPI_DOUBLE_PRECISION, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(tlen, 1, MPI_DOUBLE_PRECISION, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(npow_costaper, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_running_time_average, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_bandpass_earthquake, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(fr1, 1, MPI_REAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(fr2, 1, MPI_REAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_onebit, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_specwhitenning, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(nwf, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_suppress_notch, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(freqmin, 1, MPI_REAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(tlag, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(ipws, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(npow_pws, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_bootstrap, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(num_bootstrap, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(bootstrap_type, 3, MPI_CHARACTER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_verbose, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_overwrite_data, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_save_record, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_only_cf, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_stack, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_ac, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
! ***********************************************************************
! Allocate memory for the station, event and record elements in sdb.
! ***********************************************************************
allocate(sdb%st(nst), stat=ier)
nev_loc = ceiling(real(nev)/real(nprocs))
allocate(sdb_loc%ev(nev_loc), sdb_loc%rec(nst,nev_loc))



! ***********************************************************************
! sdb elements are filled in the master processor.
! ***********************************************************************
!if (myrank == myroot) then

! ***********************************************************************
! Fill the station elements (lon, lat, sta, net.sta) into sdb.st.
! ***********************************************************************
iunit = 3*nprocs + myrank + 13
open(unit=iunit, file='stations.lst', status='old', action='read', iostat=ier)

   ist = 0

   do
      read(iunit, *, iostat=ier) netname, staname, lat, lon
      if (0 /= ier) exit
      ist = ist + 1
      sdb%st(ist)%staname = trim(adjustl(staname))
      sdb%st(ist)%ns_name = trim(adjustl(netname))//'.'//trim(adjustl(staname))
      sdb%st(ist)%lat = lat
      sdb%st(ist)%lon = lon
   end do

   ! Save the number of stations in sdb%nst.
   sdb%nst = nst

close(unit=iunit)



! ***********************************************************************
! Do the time correction and fill in the sdb.
! ***********************************************************************
errcode = 0

npts = 0
dt = 0.d0
iunit = 4*nprocs + myrank + 11
open(unit=iunit, file='events.lst', status='old', action='read', iostat=ier)

   iev = 0

   do j = myrank, nev-1, nprocs

      ! Skip head myrank lines, because they are processed by other processorsors
      do k = 1, myrank, 1
         read(iunit,*, iostat=ier)
         if (0 /= ier) exit
      end do

      read(iunit, "(A)", iostat=ier) evtpath


      if (0 /= ier) exit
      iev = iev + 1


      ! Split the input event path name.
      call split_string(evtpath, '/', nstrArray, strArray)
      str_date = strArray(nstrArray)


      ! ***********************************************************************
      ! Fill in the event time information into sdb.ev.
      ! ***********************************************************************
      read(str_date(1:4),*) sdb_loc%ev(iev)%yy
      read(str_date(5:6),*) sdb_loc%ev(iev)%mm
      read(str_date(7:8),*) sdb_loc%ev(iev)%dd
      read(str_date(10:11),*) sdb_loc%ev(iev)%h
      read(str_date(12:13),*) sdb_loc%ev(iev)%m
      read(str_date(14:15),*) sdb_loc%ev(iev)%s
      sdb_loc%ev(iev)%jday = date2jday(sdb_loc%ev(iev)%yy, sdb_loc%ev(iev)%mm, sdb_loc%ev(iev)%dd)
      sdb_loc%ev(iev)%t0 = datetime2timestamp(sdb_loc%ev(iev)%yy, sdb_loc%ev(iev)%jday, &
                                 sdb_loc%ev(iev)%h, sdb_loc%ev(iev)%m, sdb_loc%ev(iev)%s)


      ! ***********************************************************************
      ! Create result folder and tmp path to sdb%ev%name.
      ! ***********************************************************************
      call system('mkdir -p '//trim(adjustl(tarfolder)))
      if (is_overwrite_data) then
         sdb_loc%ev(iev)%evtpath = trim(adjustl(str_date))
      else
         path = './tmp/DATA/'//trim(adjustl(strArray(nstrArray-2)))//'/'// &
          trim(adjustl(strArray(nstrArray-1)))//'/'//trim(adjustl(str_date))
         sdb_loc%ev(iev)%evtpath = trim(adjustl(path))

         ! Create the target event folder.
         call system('mkdir -p '//trim(adjustl(path)))
      end if


      ! Loop the station to processor the SAC files and fill in the sdb elements.
      do ist = 1, nst, 1


         ! ***********************************************************************
         ! Initialize the sdb.rec.t0 and sdb.rec.sacfile elements.
         ! ***********************************************************************
         sdb_loc%rec(ist,iev)%t0 = -1.d0
         sdb_loc%rec(ist,iev)%sacfile = ''
   
         ! ***********************************************************************
         ! Process the SAC file for one record and fill in the sdb.rec info.
         ! ***********************************************************************
         call mk_one_rec(sdb, iev, ist, npow_costaper, f1, f2, f3, f4, channel, evtpath, npts_read, dt_read, sdb_loc)


         sacfile = trim(adjustl(sdb_loc%rec(ist,iev)%sacfile))
         ! Check consistence in sampling points
         if (0 /= npts_read) then
         	if (0 == npts) then
         	   npts = npts_read
         	else if (npts /= npts_read) then
         	   errcode = -1
               call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
               write(msg, "(A,I0,A)") 'Error: Inconsistence of sampling points (npts_read=', npts_read, ') in '//trim(sacfile)
               !stop trim(adjustl(msg))
               write(*,*) trim(adjustl(msg))
               call flush(6)
               stop
         	end if
         end if


         ! Check consistence in sampling interval
         if (dt_read > 0.0) then
         	if (abs(dt) < 1.e-15) then
         	   dt = dt_read
         	else if ((abs(dt - dt_read) > 1.e-15) .or. (dt_read < 0.0)) then
         	   errcode = -2
               call MPI_ABORT(MPI_COMM_WORLD, -2, ier)
               write(msg, "(A,I0,A)") 'Error: Inconsistence of sampling interval (dt_read=', dt_read, ') in '//trim(sacfile)
               !stop trim(adjustl(msg))
               write(*,*) trim(adjustl(msg))
               call flush(6)
               stop
         	end if
         end if

         ! Check the correctness of tlen
         if ((0 /= npts_read) .and. ((t0 + tlen) > (npts_read-1)*dt_read)) then
            write(*,"(A)") "Error: t0 + tlen > (npts-1)*dt !"
            write(*,"(A, F12.6)") "t0 + tlen   = ", t0 + tlen
            write(*,"(A, F12.6)") "(npts-1)*dt = ", (npts_read-1)*dt_read
            write(*,"(A)") "Error: Parameters t0 and tlen must be set wrongly, please reset !"
            call flush(6)
            errcode = -3
            call MPI_ABORT(MPI_COMM_WORLD, -3, ier)
            stop
         end if

         ! Check the correctness of tlen
         if ((0 /= npts_read) .and. (tlen < 0.80*(npts_read-1)*dt_read)) then
            write(*,"(A)") "Error: tlen < 0.80*(npts-1)*dt !"
            write(*,"(A, F12.6)") "tlen             = ", tlen
            write(*,"(A, F12.6)") "0.80*(npts-1)*dt = ", (npts_read-1)*dt_read
            write(*,"(A)") "Error: Parameters t0 and tlen must be set wrongly, please reset !"
            call flush(6)
            errcode = -4
            call MPI_ABORT(MPI_COMM_WORLD, -4, ier)
            stop
         end if

      end do

      ! Skip tail nprocs-(myrank+1) lines, because they are processed by other processorsors
      do k = myrank+2, nprocs, 1
         read(iunit,*, iostat=ier)
         if (0 /= ier) exit
      end do

   end do ! end of iev = 1, nev, 1

   ! Save the number of events in sdb%nev.
   !sdb%nev = nev

close(unit=iunit)

call MPI_BARRIER(MPI_COMM_WORLD, ier)


! =======================================================================
if (0 /= errcode) then
   select case(errcode)
   case(-1)
      write(msg,"(A)") "Error: Inconsistence of sampling points, please check your data !"
   case(-2)
      write(msg,"(A)") "Error: Inconsistence of sampling interval, please chech your data !"
   case(-3,-4)
      write(msg,"(A)") "Error: Parameters t0 and tlen must be set wrongly, please reset !"
   end select
   !call flush(6)
   call MPI_ABORT(MPI_COMM_WORLD, -100, ier)
   call MPI_FINALIZE(ier)
   !stop trim(adjustl(msg))
   write(*,*) trim(adjustl(msg))
   call flush(6)
   stop
end if
! =======================================================================




call MPI_ALLREDUCE(npts, npts_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)
call MPI_ALLREDUCE(npts, npts_min, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ier)
if (npts_min /= npts_max) then
   write(msg,"(A)") "Error: Inconsistence of sampling points across processors !"
   call MPI_ABORT(MPI_COMM_WORLD, -100, ier)
   call MPI_FINALIZE(ier)
   !stop trim(adjustl(msg))
   write(*,*) trim(adjustl(msg))
   call flush(6)
   stop
end if
call MPI_ALLREDUCE(dt, dt_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
call MPI_ALLREDUCE(dt, dt_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
if (abs(dt_max - dt_min) > 1.e-15) then
   write(msg,"(A)") "Error: Inconsistence of sampling intervals across processors !"
   call MPI_ABORT(MPI_COMM_WORLD, -100, ier)
   call MPI_FINALIZE(ier)
   !stop trim(adjustl(msg))
   write(*,*) trim(adjustl(msg))
   call flush(6)
   stop
end if


! =======================================================================
! Check parameter validation
if (dt < 0.0) then
   !write(*,"(A)") 'Error: dt < 0, please check your data !'
   !call flush(6)
   call MPI_ABORT(MPI_COMM_WORLD, -100, ier)
   call MPI_FINALIZE(ier)
   stop 'Error: dt < 0, please check your data !'
end if

if (f4 > 0.5/dt) then
   !write(*,"(A)") 'Error: 1/f4 < 0.5/dt should be satisfied !'
   !call flush(6)
   call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
   call MPI_FINALIZE(ier)
   stop 'Error: 1/f4 < 0.5/dt should be satisfied !'
end if
! =======================================================================
! Pass parameter validity test
! Save the sampling interval and sampling points in sdb%dt and sdb%npts, respectively.
sdb%dt = dt
sdb%npts = npts
! =======================================================================



! Compute the displ in GatherV
allocate(recvns(nprocs), displs(nprocs), stat=ier)
!call MPI_GATHER(iev, 1, MPI_INTEGER, recvns, 1, &
!        MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(recvns, nprocs, MPI_INTEGER, &
!                 myroot, MPI_COMM_WORLD, ier)
call MPI_ALLGATHER(iev, 1, MPI_INTEGER, recvns, 1, &
                   MPI_INTEGER, MPI_COMM_WORLD, ier)
nev_gathered = sum(recvns)
if ((nev /= nev_gathered) .and. (myrank == myroot)) then
   write(*,*) 'Error: nev is not equal to nev_gathered !'
   call flush(6)
end if



! Gather sdb%ev
displs(1) = 0
do iproc = 2, nprocs, 1
   displs(iproc) = displs(iproc-1) + recvns(iproc-1)
end do
allocate(sdb%ev(nev), sdb%rec(nst,nev), stat=ier)
!call MPI_GATHERV(sdb_loc%ev, iev, event_type, sdb%ev, recvns, &
!               displs, event_type, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%ev, nev, event_type, myroot, MPI_COMM_WORLD, ier)
call MPI_ALLGATHERV(sdb_loc%ev, iev, event_type, sdb%ev, &
          recvns, displs, event_type, MPI_COMM_WORLD, ier)
! Save the number of events in sdb%nev.
sdb%nev = nev



! Gather sdb%rec
recvns = recvns*nst
displs(1) = 0
do iproc = 2, nprocs, 1
   displs(iproc) = displs(iproc-1) + recvns(iproc-1)
end do
nstxnev_loc = nst*iev
nstxnev_gathered = sum(recvns)
if ((nstxnev /= nstxnev_gathered) .and. (myrank == myroot)) then
   write(*,*) 'Error: nstxnev is not equal to nstxnev_gathered !'
   call flush(6)
end if
!call MPI_GATHERV(sdb_loc%rec, nstxnev_loc, record_type, sdb%rec, &
!         recvns, displs, record_type, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%rec, nstxnev, record_type, myroot, MPI_COMM_WORLD, ier)
call MPI_ALLGATHERV(sdb_loc%rec, nstxnev_loc, record_type, sdb%rec, &
                    recvns, displs, record_type, MPI_COMM_WORLD, ier)


deallocate(sdb_loc%ev, sdb_loc%rec, recvns, displs)



! ***********************************************************************
! Write the info database into a ascii file if is_save_record is true.
! ***********************************************************************
if (is_save_record .and. (myrank == myroot)) then
   call sacdb_to_asc(sdb, 'DataRecord.lst')
   write(*,"(A)") 'SAC data records are written into DataRecord.lst.'
   call flush(6)
end if

if (myrank == myroot) then
   write(*,"(A,/)") 'Constructing sdb struct is done ... '
   call flush(6)
end if

!end if ! if (myrank == myroot) then

call MPI_BARRIER(MPI_COMM_WORLD, ier)
! =======================================================================



!nstxnev = nst*nev
! ***********************************************************************
! Broadcast all the elements in sdb to all other processors.
! ***********************************************************************
!call MPI_BCAST(sdb%ev, nev, event_type, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%st, nst, station_type, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%rec, nstxnev, record_type, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%nev, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%nst, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)


! ***********************************************************************
! Broadcast dt and nwt to all other processors.
! ***********************************************************************
!call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(nwt, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)




! ***********************************************************************
! Note: In the following version, all processors will be used for computation instead of
! the master processor used for message passing just as those done in the original version.
! ***********************************************************************

! ***********************************************************************
! Preprocess data, including remove instrument response, fractional time correction,
! temporal normalization, spectral whitening, cuting data, computing Fourier spectrum, etc.
! ***********************************************************************
if (myrank == myroot) then

   ! =====================================================================================
   ! =============================== SECTION 3 BEGINS ====================================
   ! =====================================================================================
   ! This section removes the instrument response, cut the data, do the band-pass filtering,
   ! correct the time fraction, do the time domain normalization and spectra whitening.
   ! All the tasks are done using the so called self-scheduling mode.

   write(*,"(A)") '***********************************************************************'
   write(*,"(A)") '                         SECTION 3 BEGINS'
   write(*,"(A)") '***********************************************************************'
   write(*,"(A)") 'Data preprocessing ...'
   write(*,"(A)") '***********************************************************************'
   call flush(6)

end if

call MPI_BARRIER(MPI_COMM_WORLD, ier)



! length of taper
ntaper = max(ceiling(max(t0, taper_min) / dt), ntaper_min)
!ntaper = max(ceiling(sngl(max(0.50*(npts*dt-tlen), taper_min) / dt)), ntaper_min)



! Nlen: number of intercepted data points
Nlen = nint(tlen/dt)



! Determine corresponding half-window length for time domain normalization
! Maximum allowed half-window length is 128 in SAC for smooth command
if (nwt > 0) then
   !nwt = min(int(nwt/dt), 128)
   nwt = min(nwt, Nlen/2)
end if



ndim1 = nev
ndim2 = nst
ndim = nstxnev - 1

do iproc = myrank, ndim, nprocs

   idim1 = int8(iproc/ndim2)
   iev = idim1 + 1
   ist = iproc - int8(idim1*ndim2) + 1

   if ((iev < 1) .or. (ist < 1) .or. (iev > ndim1) .or. (ist > ndim2)) cycle


   call remove_RESP(sdb, iev, ist, f1, f2, f3, f4, channel, pzfolder)
   call preprocess(sdb, iev, ist, npow_costaper, nwt, nwf, &
                f1, f2, f3, f4, fr1, fr2, freqmin, t0, tlen)

end do

call MPI_BARRIER(MPI_COMM_WORLD, ier)




! ***********************************************************************
! Compute cross-correlation and AFTAN etc.
! ***********************************************************************

nlag = nint(tlag/dt)
if (myrank == myroot) then

   ! =====================================================================================
   ! =============================== SECTION 4 BEGINS ====================================
   ! =====================================================================================
   ! This section computes the cross-correlation and measures the dispersion curves.

   write(*,"(A)") '***********************************************************************'
   write(*,"(A)") '                         SECTION 4 BEGINS'
   write(*,"(A)") '***********************************************************************'
   if (is_ac) then
      write(*,"(A)") 'Doing auto-correlation ...'
   else
      if (is_only_cf) then
         write(*,"(A)") 'Doing cross-correlation ...'
      else
         write(*,"(A)") 'Doing cross-correlation and AFTAN ...'
      end if
   end if
   write(*,"(A)") '***********************************************************************'
   call flush(6)

   ! Delete the CCRecord.lst file if is_save_record is true.
   if (is_save_record) then
      call system('rm -rf CCRecord.lst')
   end if

end if

call MPI_BARRIER(MPI_COMM_WORLD, ier)




if (0 == ipws) then
   is_pws = .false.
else
   is_pws = .true.
end if
write(str_pws,'(I3)') ipws
write(str_npow_pws,'(I6)') npow_pws


if (is_ac) then

   ndim1 = nst
   ndim2 = nst
   ndim = nst - 1

   do iproc = myrank, ndim, nprocs

      ist1 = iproc + 1
      ist2 = ist1

      if ((ist1 < 1) .or. (ist2 < 1) .or. (ist1 > ndim1) .or. (ist2 > ndim2)) cycle


      !! NOTE: You can comment here for some specific targets,
      !!       such as Z-N or Z-E cross-correlation
      !if (ist2 <= ist1) cycle


      call cc_and_aftan(sdb, ist1, ist2, nlag, num_bootstrap, is_pws, str_pws, &
                    str_npow_pws, str_per1, str_per4, bootstrap_type, tarfolder)

   end do

else

   ndim1 = nst
   ndim2 = nst
   ndim = nst*nst - 1

   do iproc = myrank, ndim, nprocs

      idim1 = int8(iproc/ndim2)
      ist1 = idim1 + 1
      ist2 = iproc - int8(idim1*ndim2) + 1

      if ((ist1 < 1) .or. (ist2 < 1) .or. (ist1 > ndim1) .or. (ist2 > ndim2)) cycle


      ! NOTE: You can comment here for some specific targets,
      !       such as Z-N or Z-E cross-correlation
      if (ist2 <= ist1) cycle


      call cc_and_aftan(sdb, ist1, ist2, nlag, num_bootstrap, is_pws, str_pws, &
                    str_npow_pws, str_per1, str_per4, bootstrap_type, tarfolder)

   end do

end if

call MPI_BARRIER(MPI_COMM_WORLD, ier)



! ***********************************************************************
! Deallocate memory for the elements in sdb and strArray.
! ***********************************************************************
deallocate(sdb%st, sdb%ev, sdb%rec)
if (allocated(strArray)) then
   deallocate(strArray)
end if


! ***********************************************************************
! Postprocess to remove temporary folders or files
! ***********************************************************************


if (myrank == myroot) then

   write(*,"(A)")
   write(*,"(A)") 'Postprocessing ...'
   write(*,"(A)")
   call flush(6)

   ! ***********************************************************************
   ! Remove the DATA folder.
   ! ***********************************************************************
   if (is_overwrite_data) then
      call system('rm -rf '//trim(adjustl(sacfolder)))
   end if

   ! ***********************************************************************
   ! Remove the tmp folder.
   ! ***********************************************************************
   call system('rm -rf ./tmp/')
   !call system("find ./tmp -depth -type 'd' -empty -exec rmdir {} \;")

   ! ***********************************************************************
   ! Remove possible empty folder(s) and file(s).
   ! ***********************************************************************
   !call system("find "//trim(adjustl(tarfolder))//" -depth -type 'd' -empty -exec rmdir {} \;")
   !call system("find "//trim(adjustl(tarfolder))//' -name "*" -type f -size 0c | xargs -n 1 rm -f')

   write(*,"(A)")
   if (is_ac) then
      write(*,"(A)") 'Auto-correlation is done ... '
   else
      if (is_only_cf) then
         write(*,"(A)") 'Cross-correlation is done ... '
      else
         write(*,"(A)") 'Cross-correlation and AFTAN is done ... '
      end if
   end if
   write(*,"(A)")
   call flush(6)

end if

call MPI_BARRIER(MPI_COMM_WORLD, ier)



call MPI_FINALIZE(ier)



end program ANCC
