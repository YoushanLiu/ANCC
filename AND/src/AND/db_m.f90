module db_m


implicit none


! ***********************************************************************************************
! ********************************* VARIABLE DEFINITION SECTION *********************************
! ***********************************************************************************************

! some constants
integer, parameter :: SGL = SELECTED_real_KIND(5,20)
integer, parameter :: DBL = SELECTED_real_KIND(13,100)
integer, parameter :: TIMEMAX = 10000000 ! TIMEMAX: maximum time length of the SAC file
real(SGL), parameter :: hugeval = huge(0.0)
real(SGL), parameter :: tinyval = tiny(0.0)
complex(SGL), parameter :: ci = cmplx(0.0, 1.0)
complex(SGL), parameter :: czero = cmplx(0.0, 0.0)


! size and rank of MPI_COMM_WORLD
integer :: nprocs = 1
integer :: myrank = 0


! taper parameters
real, parameter :: taper_min = 1.0
integer, parameter :: ntaper_min = 10
integer :: Nlen = 0, ntaper = 0


! some logical variables
logical :: is_onebit, is_running_time_average
logical :: is_suppress_notch, is_bandpass_earthquake
logical :: is_sbs, is_pws, is_save_record, is_overwrite_data
logical :: is_onlycc, is_verbose, is_stack, is_specwhitenning


! ***********************************************************************************************
! *********************************** type DEFINITION SECTION ***********************************
! ***********************************************************************************************
type :: station
   character(len=8)  :: name = ''       ! name: station name (e.g. MONP)
   character(len=16) :: n_name = ''     ! n_name: network.station name (e.g. AZ.MONP)
   real(SGL) :: lat, lon                ! lat, lon: latitude and longitude of the station
end type station

type :: event
   character(len=512) :: name = ''      ! name: event name (e.g. absolute_path/20150401_000000)
   integer :: yy, mm, dd, h, m, jday    ! yy,mm,dd,h,m,jday: year,month,day,hour,minute,julian day
   real(DBL) :: s = 0.0                 ! ss: second in decimal form
   real(DBL) :: t0 = 0.0                ! t0: epoch time of the event (e.g. 123456789)
end type event

type :: record
   character(len=512) :: name = ''      ! name: SAC file path  (e.g. absolute_path/20150401_000000/AZ.MONP.LHZ.SAC)
   character(len=8) :: channel = ''     ! channel: channel name (e.g. LHZ)
   real(DBL) :: t0 = 0.0                ! t0: epoch time of the first data point
   real(DBL) :: frac = 0.0              ! frac: time fraction of the SAC file (e.g. 0.0)
   real(DBL) :: dt = 0.0                ! dt: data sampling interval
   integer :: npts = 0                  ! npts: number of data points (should be default to 0)
end type record

type :: sac_db
   type(event), allocatable, dimension(:) :: ev            ! ev: EVENT struct array
   type(station), allocatable, dimension(:) :: st          ! st: STATION struct array
   type(record), allocatable, dimension(:,:) :: rec        ! rec: is_save_record struct array
   integer :: nev, nst                                     ! nev, nst: number of events and stations
end type sac_db


end module db_m
