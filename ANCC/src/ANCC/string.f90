module string_m


implicit none


contains


! =============================================================================
! Delete spaces within one string
! str: input and output string [input and output]
! =============================================================================
subroutine delspace(str)

implicit none

character(len=*), intent(inout) :: str


integer nstr

character(len=512) :: str_tmp, str_result


str_result = ''

do
    str = trim(adjustl(str))
    if (str == '') exit
    read(str, *) str_tmp
    nstr = len_trim(str_tmp)
    str(1:nstr) = ''
    str_result = trim(str_result) // trim(str_tmp)
end do

str = str_result


end subroutine delspace


! =============================================================================
! Split a string by the given delimiter(s)
! InStr: input string [input]
! delimiter: string containing the delimiter(s) (e.g. '/.') [input]
! OutStr: output string array (allocatable) [output]
! n: number of output string array [output]
! =============================================================================
subroutine split_string(InStr, delimiter, nout, OutStr)

implicit none

character(len=*), intent(in) :: InStr, delimiter

integer, intent(out) :: nout
character(len=*), allocatable, dimension(:), intent(out) :: OutStr


integer i, j, istart, n
integer nstr, ndelim, ncount


n = 1
istart = 1

nstr = len_trim(InStr)
ndelim = len_trim(delimiter)

! Count how many segements the input string has
do i = 1, nstr, 1
   do j = 1, ndelim, 1
      if (InStr(i:i) == delimiter(j:j)) then
         n = n + 1
      end if
   end do
end do


allocate(OutStr(n))

! Split the input string
ncount = 0
do i = 1, nstr, 1
   do j = 1, ndelim, 1

      if (InStr(i:i) == delimiter(j:j)) then
         if (i == 1) then        ! the first element of the string is the delimiter
            OutStr(1) = ''
            ncount = ncount + 1
            istart = 2
         else
            ncount = ncount + 1
            OutStr(ncount) = trim(adjustl(InStr(istart:i-1)))
            istart = i+1
         end if
      end if

   end do
end do


nout = n
if (n > 1) then
   if (istart <= nstr) then
      OutStr(n) = InStr(istart:nstr)
   else        ! the last element of the string is the delimiter
      OutStr(n) = ''
   end if
else            ! no spliting occured
    OutStr(1) = InStr
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



function tolower(string) result(tolower_result)

implicit none

character(len=*), intent(in) :: string
character(len=len(string)) :: tolower_result


integer i, ii

do i = 1, len(string), 1

   ii = iachar(string(i:i))

   select case (ii)
      case (65:90)
         tolower_result(i:i) = achar(ii+32)
      case default
         tolower_result(i:i) = string(i:i)
   end select

end do

end function tolower




function toupper(string) result(toupper_result)

implicit none

character(len=*), intent(in) :: string
character(len=len(string)) :: toupper_result


integer i, ii

do i = 1, len(string), 1

   ii = iachar(string(i:i))

   select case (ii)
      case (97:122)
         toupper_result(i:i)=achar(ii-32)
      case default
         toupper_result(i:i)=string(i:i)
   end select

end do

end function toupper


end module string_m
