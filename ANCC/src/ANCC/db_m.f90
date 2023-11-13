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
module db_m


implicit none


! ***********************************************************************************************
! ********************************* VARIABLE DEFINITION SECTION *********************************
! ***********************************************************************************************

! some constants
integer, parameter :: SGL = SELECTED_real_KIND(5,20)
integer, parameter :: DBL = SELECTED_real_KIND(13,100)
integer, parameter :: TIMEMAX = 10000000 ! TIMEMAX: maximum time length of the SAC file
real(SGL), parameter :: HUGEVAL = huge(0.0)
real(SGL), parameter :: TINYVAL = tiny(0.0)
complex(SGL), parameter :: ci = cmplx(0.0, 1.0)
complex(SGL), parameter :: czero = cmplx(0.0, 0.0)


! size and rank of MPI_COMM_WORLD
integer :: nprocs = 1
integer :: myrank = 0


! taper parameters
integer, parameter :: ntaper_min = 10
real, parameter :: taper_min = 1.0

integer :: Nlen = 0, ntaper = 0


! some logical variables
logical is_stack, is_pws
logical is_only_cf, is_ac
logical is_specwhitening
logical is_verbose, is_bootstrapp
logical is_overwrite_data, is_save_record
logical is_onebit, is_running_time_average
logical is_bandpass_earthquake, is_suppress_notch


! ***********************************************************************************************
! *********************************** type DEFINITION SECTION ***********************************
! ***********************************************************************************************
type station
   character(len=8)  :: staname = ''                       ! staname: station name (e.g. MONP)
   character(len=17) :: ns_name = ''                       ! ns_name: network.station name (e.g. AZ.MONP)
   real(SGL) :: lat, lon                                   ! lat, lon: latitude and longitude of the station
end type station

type event
   character(len=128) :: evtpath = ''                      ! evtpath: event path (e.g. path/20150401_000000)
   integer :: yy, mm, dd, h, m, jday                       ! yy, mm, dd, h, m, jday: year, month, day, hour, minute, julian day
   real(DBL) :: s = 0.d0                                   ! ss: second in decimal form
   real(DBL) :: t0 = -1.d0                                 ! t0: epoch time of the event (e.g. 123456789)
end type event

type record
   character(len=192) :: sacfile = ''                      ! sacname: SAC path and filename (e.g. path/20150401_000000/AZ.MONP.BHZ.SAC)
   real(DBL) :: t0 = -1.d0                                 ! t0: epoch time of the first data point
end type record

type sac_db
   type(event), allocatable, dimension(:) :: ev            ! ev: EVENT struct array
   type(station), allocatable, dimension(:) :: st          ! st: STATION struct array
   type(record), allocatable, dimension(:,:) :: rec        ! rec: record struct array
   integer :: nev, nst                                     ! nev, nst: number of events and stations
   integer :: npts = 0                                     ! npts: number of data points (should be default to 0)
   real(DBL) :: dt = 0.d0                                  ! dt: data sampling interval
end type sac_db


end module db_m
