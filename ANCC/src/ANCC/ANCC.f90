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
character(64) :: version = ' (v6.3)'

integer i, j, k, jmax, nlag
integer nargin, myroot, ier
integer iev, ist, ist1, ist2
integer nstrArray, N_bs, ipws, iunit
integer npow_costaper, nwt, nwf, nweight
integer event_type, station_type, record_type
integer nst, nev, nev_tmp, nev_gathered, npts, count

logical flag

integer(8) ndim, ndim1, ndim2, idim1, iproc
integer(4) nstxnev, nstxnev_tmp, nstxnev_gathered


real(SGL) lon, lat, freqmin, tlag
real(SGL) f1, f2, f3, f4, fr1, fr2

real(DBL) dt, t0, tlen


character(len=3) str_sbs, str_pws, str_save_record, str_verbose, bs_type
character(len=3) str_onlycc, str_overwrite_data, str_stack_cc, str_specwhitenning, str_ac
character(len=3) str_onebit, str_running_time_average, str_bandpass_earthquake, str_suppress_notch

character(len=8) netname, staname, channel

character(len=32) str_per1, str_per4, str_weight

character(len=512) evpath, str_tmp, path
character(len=512) sacfolder, pzfolder, tarfolder


type(event) evt_tmp
type(station) sta_tmp
type(record) rec_tmp
type(sac_db) sdb, sdb_tmp


integer, allocatable, dimension(:) :: blocklen, types
integer(MPI_ADDRESS_KIND), allocatable, dimension(:) :: base, disp

integer, allocatable, dimension(:) :: recvcounts, displs

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
count = 4
allocate(base(count), disp(count), blocklen(count), types(count))
blocklen = (/ 512, 6, 1, 1 /)
types =(/ MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION /)
call MPI_GET_ADDRESS(evt_tmp%name, base(1), ier)
call MPI_GET_ADDRESS(evt_tmp%yy, base(2), ier)
call MPI_GET_ADDRESS(evt_tmp%s, base(3), ier)
call MPI_GET_ADDRESS(evt_tmp%t0, base(4), ier)
disp(1) = 0
disp(2) = base(2) - base(1)
disp(3) = base(3) - base(1)
disp(4) = base(4) - base(1)
call MPI_TYPE_CREATE_STRUCT(count, blocklen, disp, types, event_type, ier)
call MPI_TYPE_COMMIT(event_type, ier)
deallocate(base, disp, blocklen, types)


! ***********************************************************************
! Construct new station data type for later data sharing.
! ***********************************************************************
count = 3
allocate(base(count), disp(count), blocklen(count), types(count))
blocklen = (/ 8, 16, 2 /)
types =(/ MPI_CHARACTER, MPI_CHARACTER, MPI_REAL /)
call MPI_GET_ADDRESS(sta_tmp%name, base(1), ier)
call MPI_GET_ADDRESS(sta_tmp%n_name, base(2), ier)
call MPI_GET_ADDRESS(sta_tmp%lat, base(3), ier)
disp(1) = 0
disp(2) = base(2) - base(1)
disp(3) = base(3) - base(1)
call MPI_TYPE_CREATE_STRUCT(count, blocklen, disp, types, station_type, ier)
call MPI_TYPE_COMMIT(station_type, ier)
deallocate(base, disp, blocklen, types)


! ***********************************************************************
! Construct new record data type for later data sharing.
! ***********************************************************************
count = 6
allocate(base(count), disp(count), blocklen(count), types(count))
blocklen = (/ 512, 8, 1, 1, 1, 1 /)
types =(/ MPI_CHARACTER, MPI_CHARACTER, MPI_DOUBLE_PRECISION, &
     MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_INTEGER /)
call MPI_GET_ADDRESS(rec_tmp%name, base(1), ier)
call MPI_GET_ADDRESS(rec_tmp%channel, base(2), ier)
call MPI_GET_ADDRESS(rec_tmp%t0, base(3), ier)
call MPI_GET_ADDRESS(rec_tmp%frac, base(4), ier)
call MPI_GET_ADDRESS(rec_tmp%dt, base(5), ier)
call MPI_GET_ADDRESS(rec_tmp%npts, base(6), ier)
disp(1) = 0
disp(2) = base(2) - base(1)
disp(3) = base(3) - base(1)
disp(4) = base(4) - base(1)
disp(5) = base(5) - base(1)
disp(6) = base(6) - base(1)
call MPI_TYPE_CREATE_STRUCT(count, blocklen, disp, types, record_type, ier)
call MPI_TYPE_COMMIT(record_type, ier)
deallocate(base, disp, blocklen, types)





!if (myrank == myroot) then

! Parse the number of command line argumet(s).
nargin = command_argument_count()

! Stop if the number of input argument(s) is wrong.
if (3 /= nargin) then
   write(*,"(A)") 'ANCC' // trim(version)
   write(*,"(A)") "Usage: AND sacfolder pzfolder tarfolder"
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
   write(*,"(A)") '************************************'
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
   read(iunit,*) str_pws, nweight
   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   read(iunit,*) str_sbs, N_bs
   read(iunit,*) bs_type
   do i = 1, 20, 1
      read(iunit,*)
   end do
   read(iunit,*) str_save_record
   read(iunit,*) str_verbose
   read(iunit,*) str_onlycc
   read(iunit,*) str_overwrite_data
   read(iunit,*) str_stack_cc

close(unit=iunit)


! ***********************************************************************
! Convert the periods to frequencies.
! ***********************************************************************
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
is_suppress_notch = .false.
is_sbs = .false.
is_save_record = .false.
is_verbose = .false.
ipws = 0
is_onlycc = .false.
is_overwrite_data = .false.
is_stack = .false.
is_specwhitenning = .false.
is_ac = .false.
if ((str_running_time_average == 'Y') .or. (str_running_time_average == 'y')) is_running_time_average = .true.
if ((str_bandpass_earthquake == 'Y') .or. (str_bandpass_earthquake == 'y')) is_bandpass_earthquake = .true.
if ((str_onebit == 'Y') .or. (str_onebit == 'y')) is_onebit = .true.
if ((str_suppress_notch == 'Y') .or. (str_suppress_notch == 'y')) is_suppress_notch = .true.
if ((str_sbs == 'Y') .or. (str_sbs == 'y')) is_sbs = .true.
if ((str_save_record == 'Y') .or. (str_save_record == 'y')) is_save_record = .true.
if ((str_verbose == 'Y') .or. (str_verbose == 'y')) is_verbose = .true.
if ((str_pws == 'Y') .or. (str_pws == 'y')) ipws = 1
if ((str_onlycc == 'Y') .or. (str_onlycc == 'y')) is_onlycc = .true.
if ((str_overwrite_data == 'Y') .or. (str_overwrite_data == 'y')) is_overwrite_data = .true.
if ((str_stack_cc == 'Y') .or. (str_stack_cc == 'y')) is_stack = .true.
if ((str_specwhitenning == 'Y') .or. (str_specwhitenning == 'y')) is_specwhitenning = .true.
if ((str_ac == 'Y') .or. (str_ac == 'y')) is_ac = .true.
if (myrank == myroot) then
   write(*,"(A,/)") 'Reading input parameters is done ... '
   call flush(6)
end if



! Check parameter validation
if (is_bandpass_earthquake .and. (fr1 >= fr2)) then
   write(*,"(A)") 'Error: fr1 < fr2 should be satisfied !'
   call flush(6)
   call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
   call MPI_FINALIZE(ier)
   stop
end if



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
! This section process the sac files and fill in the elements in the sdb struct.
if (myrank == myroot) then
   write(*,"(A)") '***********************************************************************'
   write(*,"(A)") '                         SECTION 2 BEGINS'
   write(*,"(A)") '***********************************************************************'
   write(*,"(A)") 'Constructing sdb struct ...'
   write(*,"(A)") '************************************'
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
!! Broadcast input parametes from the master proces to all other process.
!! ***********************************************************************
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
!call MPI_BCAST(is_onebit, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_suppress_notch, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(fr1, 1, MPI_REAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(fr2, 1, MPI_REAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(nwf, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(freqmin, 1, MPI_REAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_sbs, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(N_bs, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(bs_type, 3, MPI_CHARACTER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(tlag, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(ipws, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(nweight, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sacfolder, 512, MPI_CHARACTER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(pzfolder, 512, MPI_CHARACTER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(tarfolder, 512, MPI_CHARACTER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_save_record, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_verbose, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_onlycc, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_overwrite_data, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_stack, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(is_specwhitenning, 1, MPI_LOGICAL, myroot, MPI_COMM_WORLD, ier)




! ***********************************************************************
! Allocate memory for the station, event and record elements in sdb.
! ***********************************************************************
allocate(sdb%st(nst), stat=ier)
nev_tmp = ceiling(real(nev)/real(nprocs))
allocate(sdb_tmp%ev(nev_tmp), sdb_tmp%rec(nst,nev_tmp))



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
      read(iunit, *, iostat=ier) netname, staname, lon, lat
      if (0 /= ier) exit
      ist = ist + 1
      sdb%st(ist)%name = trim(adjustl(staname))
      sdb%st(ist)%n_name = trim(adjustl(netname))//'.'//trim(adjustl(staname))
      sdb%st(ist)%lon = lon
      sdb%st(ist)%lat = lat

   end do

   ! Save the number of stations in sdb%nst.
   sdb%nst = nst

close(unit=iunit)



! ***********************************************************************
! do the time correction and fill in the sdb.
! ***********************************************************************
flag = .false.

iunit = 4*nprocs + myrank + 11
open(unit=iunit, file='events.lst', status='old', action='read', iostat=ier)

   iev = 0

   do j = myrank, nev-1, nprocs

      ! Skip head myrank lines, because they are processed by other processors
      do k = 1, myrank, 1
         read(iunit,*, iostat=ier)
         if (0 /= ier) exit
      end do

      read(iunit, "(A)", iostat=ier) evpath


      if (0 /= ier) exit
      iev = iev + 1


      ! Split the input event path name.
      call split_string(evpath, '/', strArray, nstrArray)
      str_tmp = strArray(nstrArray)


      ! ***********************************************************************
      ! Fill in the event time information into sdb.ev.
      ! ***********************************************************************
      read(str_tmp(1:4),*) sdb_tmp%ev(iev)%yy
      read(str_tmp(5:6),*) sdb_tmp%ev(iev)%mm
      read(str_tmp(7:8),*) sdb_tmp%ev(iev)%dd
      read(str_tmp(10:11),*) sdb_tmp%ev(iev)%h
      read(str_tmp(12:13),*) sdb_tmp%ev(iev)%m
      read(str_tmp(14:15),*) sdb_tmp%ev(iev)%s
      sdb_tmp%ev(iev)%jday = date2jday(sdb_tmp%ev(iev)%yy, sdb_tmp%ev(iev)%mm, sdb_tmp%ev(iev)%dd)
      sdb_tmp%ev(iev)%t0 = datetime2timestamp(sdb_tmp%ev(iev)%yy, sdb_tmp%ev(iev)%jday, &
                                    sdb_tmp%ev(iev)%h, sdb_tmp%ev(iev)%m, dble(sdb_tmp%ev(iev)%s))


      ! ***********************************************************************
      ! Create result folder and tmp path to sdb%ev%name.
      ! ***********************************************************************
      call system('mkdir -p '//trim(adjustl(tarfolder)))
      if (is_overwrite_data) then
         path = trim(adjustl(tarfolder))//'/DATA/'//trim(adjustl(strArray(nstrArray-2)))// &
                      '/'//trim(adjustl(strArray(nstrArray-1)))//'/'//trim(adjustl(str_tmp))
         sdb_tmp%ev(iev)%name = trim(adjustl(path))
      else
         path = './tmp/DATA/'//trim(adjustl(strArray(nstrArray-2)))//'/'// &
           trim(adjustl(strArray(nstrArray-1)))//'/'//trim(adjustl(str_tmp))
         sdb_tmp%ev(iev)%name = trim(adjustl(path))

         ! Create the target event folder.
         call system('mkdir -p '//trim(adjustl(path)))
      end if


      ! Loop the station to processor the sac files and fill in the sdb elements.
      do ist = 1, nst, 1

         ! ***********************************************************************
         ! Initiate the sdb.rec.nrec and sdb.rec.frac elements.
         ! ***********************************************************************
         sdb_tmp%rec(ist,iev)%npts = 0
         sdb_tmp%rec(ist,iev)%frac = 0.0

         ! ***********************************************************************
         ! processor the sac file for one record and fill in the sdb.rec info.
         ! ***********************************************************************
         call mk_one_rec(evpath, iev, ist, channel, sdb, sdb_tmp)


         dt = sdb_tmp%rec(ist,iev)%dt
         npts = sdb_tmp%rec(ist,iev)%npts
         if ((0 /= npts ) .and. ((t0 + tlen) > (npts-1)*dt)) then
            write(*,"(A)") "Error: t0 + tlen > (npts-1)*dt !"
            write(*,"(A, F12.6)") "t0 + tlen   = ", t0 + tlen
            write(*,"(A, F12.6)") "(npts-1)*dt = ", (npts-1)*dt
            write(*,"(A)") "Error: parameters t0 and tlen must be set wrongly, please reset !"
            call flush(6)
            flag = .true.
            call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
            stop
         end if

         if ((0 /= npts ) .and. (tlen < 0.80*(npts-1)*dt)) then
            write(*,"(A)") "Error: tlen < 0.80*(npts-1)*dt !"
            write(*,"(A, F12.6)") "tlen             = ", tlen
            write(*,"(A, F12.6)") "0.80*(npts-1)*dt = ", (npts-1)*dt
            write(*,"(A)") "Error: parameters t0 and tlen must be set wrongly, please reset !"
            call flush(6)
            flag = .true.
            call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
            stop
         end if

      end do

      ! Skip tail nprocs-(myrank+1) lines, because they are processed by other processors
      do k = myrank+2, nprocs, 1
         read(iunit,*, iostat=ier)
         if (0 /= ier) exit
      end do

   end do ! end of iev = 1, nev, 1

   ! Save the number of events in sdb%nev.
   !sdb%nev = nev

close(unit=iunit)



if (flag) then
   write(*,"(A)") "Error: parameters t0 and tlen must be set wrongly, please reset !"
   call flush(6)
   call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
   call MPI_FINALIZE(ier)
   stop
end if




! compute the displ in GatherV
allocate(recvcounts(nprocs), displs(nprocs), stat=ier)
!call MPI_GATHER(iev, 1, MPI_INTEGER, recvcounts, 1, &
!            MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(recvcounts, nprocs, MPI_INTEGER, &
!                     myroot, MPI_COMM_WORLD, ier)
call MPI_ALLGATHER(iev, 1, MPI_INTEGER, recvcounts, 1, &
                        MPI_INTEGER, MPI_COMM_WORLD, ier)
nev_gathered = sum(recvcounts)
if ((nev /= nev_gathered) .and. (myrank == myroot)) then
   write(*,*) 'Error: nev is not equal to nev_gathered !'
   call flush(6)
end if



displs(1) = 0
do iproc = 2, nprocs, 1
   displs(iproc) = displs(iproc-1) + recvcounts(iproc-1)
end do
allocate(sdb%ev(nev), sdb%rec(nst,nev), stat=ier)
!call MPI_GATHERV(sdb_tmp%ev, iev, event_type, sdb%ev, &
!     recvcounts, displs, event_type, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%ev, nev, event_type, myroot, MPI_COMM_WORLD, ier)
call MPI_ALLGATHERV(sdb_tmp%ev, iev, event_type, sdb%ev, &
      recvcounts, displs, event_type, MPI_COMM_WORLD, ier)
! Save the number of events in sdb%nev.
sdb%nev = nev



recvcounts = recvcounts*nst
displs(1) = 0
do iproc = 2, nprocs, 1
   displs(iproc) = displs(iproc-1) + recvcounts(iproc-1)
end do
nstxnev_tmp = nst*iev
nstxnev_gathered = sum(recvcounts)
if ((nstxnev /= nstxnev_gathered) .and. (myrank == myroot)) then
   write(*,*) 'Error: nstxnev is not equal to nstxnev_gathered !'
   call flush(6)
end if
!call MPI_GATHERV(sdb_tmp%rec, nstxnev_tmp, record_type, sdb%rec, &
!     recvcounts, displs, record_type, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%rec, nstxnev, record_type, myroot, MPI_COMM_WORLD, ier)
call MPI_ALLGATHERV(sdb_tmp%rec, nstxnev_tmp, record_type, sdb%rec, &
                  recvcounts, displs, record_type, MPI_COMM_WORLD, ier)


deallocate(sdb_tmp%ev, sdb_tmp%rec, recvcounts, displs)



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



! Check consistence of dt
dt = -1.0
flag = .false.
do iev = 1, nev, 1
   if (flag) exit
   do ist = 1, nst, 1
      if ((sdb%rec(ist,iev)%npts) > 0) then
         dt = sdb%rec(ist,iev)%dt
         flag = .true.
         exit
      endif
   end do
end do


if (dt < 0.0) then
   write(*,"(A)") 'Error: input data is wrongly set, please check !'
   call flush(6)
   call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
   call MPI_FINALIZE(ier)
   stop
end if


! Check parameter validation
if (.not.( (f1 < f2) .and. (f2 < f3) .and. (f3 < f4))) then
   write(*,"(A)") 'Error: f1 > f2 > f3 > f4 [sec.] should be satisfied !'
   call flush(6)
   call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
   call MPI_FINALIZE(ier)
   stop
end if


if (f4 > 0.5/dt) then
   write(*,"(A)") 'Error: 1/f4 < 0.5/dt should be satisfied !'
   call flush(6)
   call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
   call MPI_FINALIZE(ier)
   stop
end if



! Determine corresponding half-window length for time domain normalization
! Maximum allowed half-window length is 128 in sac for smooth command
if (nwt > 0.0) then
   nwt = min(int(nwt/dt), 128)
end if



call MPI_BARRIER(MPI_COMM_WORLD, ier)



!nstxnev = nst*nev
! ***********************************************************************
! Broadcast all the elements in sdb to all other process.
! ***********************************************************************
!call MPI_BCAST(sdb%ev, nev, event_type, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%st, nst, station_type, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%rec, nstxnev, record_type, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%nev, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)
!call MPI_BCAST(sdb%nst, 1, MPI_INTEGER, myroot, MPI_COMM_WORLD, ier)


! ***********************************************************************
! Broadcast dt and nwt to all other process.
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
   write(*,"(A)") '************************************'
   call flush(6)

end if

call MPI_BARRIER(MPI_COMM_WORLD, ier)



! length of taper
ntaper = max(ceiling(max(t0, taper_min) / dt), ntaper_min)
!ntaper = max(ceiling(sngl(max(0.50*(npts*dt-tlen), taper_min) / dt)), ntaper_min)



! Nlen: number of intercepted data points
Nlen = nint(tlen/dt)




ndim1 = nev
ndim2 = nst
ndim = nstxnev - 1

do iproc = myrank, ndim, nprocs

   idim1 = int(int8(iproc/ndim2))
   iev = idim1 + 1
   ist = iproc - int8(idim1*ndim2) + 1

   if ((iev < 1) .or. (ist < 1) .or. (iev > ndim1) .or. (ist > ndim2)) cycle


   call remove_RESP(sdb, iev, ist, f1, f2, f3, f4, pzfolder)
   call preprocess(sdb, iev, ist, npow_costaper, nwt, nwf, &
                f1, f2, f3, f4, fr1, fr2, freqmin, t0, tlen)

end do



call MPI_BARRIER(MPI_COMM_WORLD, ier)




! ***********************************************************************
! compute cross-correlation and AFTAN etc.
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
      if (is_onlycc) then
         write(*,"(A)") 'Doing cross-correlation ...'
      else
         write(*,"(A)") 'Doing cross-correlation and AFTAN ...'
      end if
   end if
   write(*,"(A)") '****************************************'
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
write(str_per1,'(F8.2)') 1.0/f1
write(str_per4,'(F8.2)') 1.0/f4
write(str_weight,'(I6)') nweight


if (is_ac) then


   ndim1 = nst
   ndim2 = nst
   ndim = nst - 1

   do iproc = myrank, ndim, nprocs

      ist1 = iproc + 1
      ist2 = ist1

      if ((ist1 < 1) .or. (ist2 < 1) .or. (ist1 > ndim1) .or. (ist2 > ndim2)) cycle


      ! NOTE: You can comment here for some specific targets,
      !       such as Z-N or Z-E cross-correlation
      if (ist2 <= ist1) cycle


      call cc_and_aftan(sdb, ist1, ist2, nlag, N_bs, is_pws, str_pws, &
                  str_weight, str_per1, str_per4, bs_type, tarfolder)

   end do

else

   ndim1 = nst
   ndim2 = nst
   ndim = nst*nst - 1

   do iproc = myrank, ndim, nprocs

      idim1 = int(int8(iproc/ndim2))
      ist1 = idim1 + 1
      ist2 = iproc - int8(idim1*ndim2) + 1

      if ((ist1 < 1) .or. (ist2 < 1) .or. (ist1 > ndim1) .or. (ist2 > ndim2)) cycle


      ! NOTE: You can comment here for some specific targets,
      !       such as Z-N or Z-E cross-correlation
      if (ist2 <= ist1) cycle


      call cc_and_aftan(sdb, ist1, ist2, nlag, N_bs, is_pws, str_pws, &
                  str_weight, str_per1, str_per4, bs_type, tarfolder)

   end do

end if

call MPI_BARRIER(MPI_COMM_WORLD, ier)



! ***********************************************************************
! post-process to remove temporary folders or files
! ***********************************************************************

! ***********************************************************************
! Deallocate memory for the elements in sdb and strArray.
! ***********************************************************************
deallocate(sdb%st, sdb%ev, sdb%rec)
if (allocated(strArray)) then
   deallocate(strArray)
end if


write(*,"(A)") 'Postprocessing ...'

if (myrank == myroot) then

   ! ***********************************************************************
   ! Remove the DATA folder.
   ! ***********************************************************************
   if (.not.(is_overwrite_data)) then
      call system('rm -rf ./tmp/')
      !call system('rm -rf '//'./tmp/DATA')
   end if

   ! ***********************************************************************
   ! Remove possible empty folder(s) and file(s).
   ! ***********************************************************************
   !call system("find "//trim(adjustl(tarfolder))//" -depth -type 'd' -empty -exec rmdir {} \;")
   !call system("find "//trim(adjustl(tarfolder))//' -name "*" -type f -size 0c | xargs -n 1 rm -f')

   write(*,"(A)")
   if (is_ac) then
       write(*,"(A)") 'Auto-correlation is done ... '
   else
      if (is_onlycc) then
         write(*,"(A)") 'Cross-correlation is done ... '
      else
         write(*,"(A)") 'Cross-correlation and AFTAN is done ... '
      end if
   end if
   write(*,"(A)")
   call flush(6)

   ! ***********************************************************************
   ! Remove the tmp folder.
   ! ***********************************************************************
   !call system("find ./tmp -depth -type 'd' -empty -exec rmdir {} \;")

end if




call MPI_BARRIER(MPI_COMM_WORLD, ier)



call MPI_FINALIZE(ier)



end program ANCC
