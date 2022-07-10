! This file is part of ANCC.
!
! AND is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! AND is distributed in the hope that it will be useful,
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

integer :: i, leap_year


if (mod(year, 400) == 0) then
   leap_year = 1
else if (mod(year, 100) == 0) then
   leap_year = 0
else if (mod(year, 4) == 0) then
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

implicit none

integer, intent(in) :: year, jday, hour, minute
real(DBL), intent(in) :: sec

integer(8) yeardiff, vis


yeardiff = year - 1970
vis = (yeardiff+1)/4 + (yeardiff+369)/400 - (yeardiff+69)/100
datetime2timestamp = (365*yeardiff + vis + jday - 1)*86400.d0 + (hour*60 + minute)*60.d0 + sec


end function datetime2timestamp


! ====================================================================
! Convert timestamp to datetime
! t: timestamp since 1970-01-01T00:00:000000Z in seconds [input]
! year: year [output]
! jday: julian day [output]
! hour: hour [output]
! min: minute [output]
! sec: second [output]
! ====================================================================
subroutine timestamp2datetime(t, year, jday, hour, minute, sec)

implicit none

real(DBL), intent(in) :: t

integer, intent(out) :: year, jday, hour, minute
real(DBL), intent(out) :: sec


integer(8) idate, itime, irest, iysupp, idsupp
real(DBL) fsec


idate = int(t/86400.d0)
iysupp = idate/365
idsupp = iysupp*365 + (iysupp+1)/4 + (iysupp+369)/400 - (iysupp+69)/100
if (idate < idsupp) iysupp = iysupp - 1
idsupp = iysupp*365 + (iysupp+1)/4 + (iysupp+369)/400 - (iysupp+69)/100


! extract year
year = int(1970 + iysupp, kind=4)
! extract jday
jday = int(idate - idsupp + 1, kind=4)
! time part
fsec = t - idate*86400.d0
irest = int(fsec)
fsec = fsec - irest
itime = irest
irest = itime/60
! extract second
sec = int(itime - irest*60 + fsec, kind=4)
itime = irest
irest = itime/60
! extract minute
minute = int(itime - irest*60, kind=4)
! extract hour
hour = int(irest, kind=4)


end subroutine timestamp2datetime


end module date_time_m
