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
module string_m


implicit none


contains


! =============================================================================
! Delete spaces within one string
! str: input and output string [input and output]
! =============================================================================
subroutine delspace(str)

implicit none

integer :: str_len

character(len=*), intent(inout) :: str
character(len=512) :: str_tmp, str_result



str_result = ''


do
    str = trim(adjustl(str))
    if (str == '') exit
    read(str, *) str_tmp
    str_len = len(trim(str_tmp))
    str(1:str_len) = ''
    str_result = trim(str_result) // trim(str_tmp)
end do

str = str_result


end subroutine delspace


! =============================================================================
! Split a string by the given delimiter(s)
! InStr: input string [input]
! delimiter: string containing the delimiter(s) (e.g. '/.') [input]
! strArray: output string array (allocatable) [output]
! nsize: number of output string array [output]
! =============================================================================
subroutine split_string(InStr, delimiter, strArray, nsize)

implicit none

character(len=*), intent(in) :: InStr, delimiter
character(len=*), allocatable, dimension(:), intent(out) :: strArray
integer, intent(out) :: nsize


integer i, j, istart
integer nstr, ndelim, ncount


nsize = 1
istart = 1

nstr = len(trim(InStr))
ndelim = len(trim(delimiter))

! Count how many segements the input string has
do i = 1, nstr, 1
   do j = 1, ndelim, 1
      if (InStr(i:i) == delimiter(j:j)) then
         nsize = nsize + 1
      end if
   end do
end do


allocate(strArray(nsize))

! Split the input string
ncount = 0
do i = 1, nstr, 1
   do j = 1, ndelim, 1

      if (InStr(i:i) == delimiter(j:j)) then
         if (i == 1) then        ! the first element of the string is the delimiter
            strArray(1) = ''
            ncount = ncount + 1
            istart = 2
         else
            ncount = ncount + 1
            strArray(ncount) = InStr(istart:i-1)
            istart = i+1
         end if
      end if

   end do
end do


if (nsize > 1) then
   if (istart <= nstr) then
      strArray(nsize) = InStr(istart:nstr)
   else        ! the last element of the string is the delimiter
      strArray(nsize) = ''
   end if
else            ! no spliting occured
    strArray(1) = InStr
end if


end subroutine split_string


! =============================================================================
! Convert a float point number to a string with leading zeros
! value: input floating point number [input]
! m: width of the integer part [input]
! n: width of the decimal part [input]
! padzero: output formatted string [output]
! =============================================================================
function padzero(numeric, m, n)

implicit none

real, intent(in) :: numeric
integer, intent(in) :: m, n

character(len=m+n+2) :: padzero


character(len=32) str_tmp1, str_tmp2
character(len=256) str_tmp



if (numeric >= 0.0) then
   write(str_tmp1, "(I0, A1, I0)") m, '.', m
   write(str_tmp2, "(I0, A1, I0)") 0, '.', n
   write(str_tmp,"('(A1,', 'I', A, ',F', A, ')')") trim(adjustl(str_tmp1)), trim(adjustl(str_tmp2))
   write(padzero,trim(adjustl(str_tmp))) '+', abs(int(numeric)), abs(numeric-int(numeric))
else
   write(str_tmp1, "(I0, A1, I0)") m, '.', m
   write(str_tmp2, "(I0, A1, I0)") 0, '.', n
   write(str_tmp,"('(A1,', 'I', A, ',F', A, ')')") trim(adjustl(str_tmp1)), trim(adjustl(str_tmp2))
   write(padzero,trim(adjustl(str_tmp))) '-', abs(int(numeric)), abs(numeric-int(numeric))
end if


return


end function padzero



end module string_m
