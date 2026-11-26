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
module date_time_m

use db_m



implicit none



contains


! ====================================================================
! Convert the date to julian day
! year: year [input]
! month: month [input]
! day: day [input]
! date2jday: converted julian day [output]
! ====================================================================
integer function date2jday(year, month, day)

implicit none

integer, intent(in) :: year, month, day

integer i, leap_year


if (0 == mod(year,400)) then
   leap_year = 1
else if (0 == mod(year,100)) then
   leap_year = 0
else if (0 == mod(year,4)) then
   leap_year = 1
else
   leap_year = 0
end if


date2jday = day

do i = 1, month-1, 1

   select case (i)
      case (1,3,5,7,8,10,12)
         date2jday = date2jday + 31
      case (4,6,9,11)
         date2jday = date2jday + 30
      case (2)
         date2jday = date2jday + 28 + leap_year
   end select

end do


end function date2jday


! ====================================================================
! Convert datetime to timestamp since 1970-01-01T00:00:000000Z in seconds
! year: year [input]
! jday: julian day [input]
! hour: hour [input]
! min: minute [input]
! sec: second [input]
! datetime2timestamp: timestamp in seconds [output]
! ====================================================================
real(DBL) function datetime2timestamp(year, jday, hour, minute, sec)

use iso_fortran_env, only: int64

implicit none

integer, intent(in) :: year, jday, hour, minute
real(DBL), intent(in) :: sec

integer(int64) yeardiff, leap_days, days

yeardiff = int(year, int64) - 1970_int64

! leap days between 1970 and the beginning of year
leap_days = floor_div(yeardiff + 1_int64, 4_int64)   &
          + floor_div(yeardiff + 369_int64, 400_int64) &
          - floor_div(yeardiff + 69_int64, 100_int64)

days = 365_int64 * yeardiff + leap_days + int(jday, int64) - 1_int64

datetime2timestamp = real(days, DBL) * 86400.d0 &
                   + real(hour * 3600 + minute * 60, DBL) &
                   + sec

end function datetime2timestamp


! ====================================================================
! Convert timestamp to datetime
! t: timestamp since 1970-01-01T00:00:000.000Z in seconds [input]
! year: year [output]
! jday: julian day [output]
! hour: hour [output]
! min: minute [output]
! sec: second [output]
! ====================================================================
subroutine timestamp2datetime(t, year, jday, hour, minute, sec)

use iso_fortran_env, only: int64

implicit none

real(DBL), intent(in) :: t

integer, intent(out) :: year, jday, hour, minute
real(DBL), intent(out) :: sec

integer(int64) idate, iysupp, idsupp
integer(int64) irest, itime
real(DBL) fsec

! extract integer days since 1970-01-01.
! floor is required for negative timestamps.
idate = floor(t / 86400.d0, kind=int64)

! estimate year offset from 1970.
iysupp = floor_div(idate, 365_int64)

! correct the estimated year.
idsupp = iysupp * 365_int64                     &
       + floor_div(iysupp + 1_int64, 4_int64)   &
       + floor_div(iysupp + 369_int64, 400_int64) &
       - floor_div(iysupp + 69_int64, 100_int64)

do while (idate < idsupp)

    iysupp = iysupp - 1_int64

    idsupp = iysupp * 365_int64                    &
           + floor_div(iysupp + 1_int64, 4_int64) &
           + floor_div(iysupp + 369_int64, 400_int64) &
           - floor_div(iysupp + 69_int64, 100_int64)

end do

do

    itime = iysupp + 1_int64

    itime = itime * 365_int64                    &
          + floor_div(itime + 1_int64, 4_int64) &
          + floor_div(itime + 369_int64, 400_int64) &
          - floor_div(itime + 69_int64, 100_int64)

    if (idate < itime) exit

    iysupp = iysupp + 1_int64
    idsupp = itime

end do

! extract year and day of year.
year = int(1970_int64 + iysupp)
jday = int(idate - idsupp + 1_int64)


! extract time of day.
fsec = t - real(idate, DBL) * 86400.d0

irest = int(fsec, int64)
fsec = fsec - real(irest, DBL)

itime = irest

irest = itime / 60_int64
sec = real(itime - irest * 60_int64, DBL) + fsec

itime = irest

irest = itime / 60_int64
minute = int(itime - irest * 60_int64)

hour = int(irest)

end subroutine timestamp2datetime


! ====================================================================
pure function floor_div(a, b) result(q)

use iso_fortran_env, only: int64

implicit none

integer(int64), intent(in) :: a, b

integer(int64) q, r

q = a / b
r = mod(a, b)

if (r /= 0_int64 .and. ((r < 0_int64) .neqv. (b < 0_int64))) then
    q = q - 1_int64
end if

end function floor_div


end module date_time_m
