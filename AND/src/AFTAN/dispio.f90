! This file is part of ANCC.
!
! AFTAN is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! AFTAN is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
!
module dispio_m

implicit none

integer, parameter :: SGL = SELECTED_REAL_KIND(5,20)
integer, parameter :: DBL = SELECTED_REAL_KIND(13,100)


contains


! ==============================================================
! write the dispersion results after linear interpolation.
! array1: raw measurement matrix [input]
! nfout1: # of raw measurements (# of columns) [input]
! array2: clean measurement matrix [input]
! nfout2: # of clean measurements (# of columns) [input]
! filename: output file name [input]
! ==============================================================
subroutine write_data(array1, nfout1, array2, nfout2, filename)

implicit none

integer, intent(in) :: nfout1, nfout2

character(len=*), intent(in) :: filename

real(SGL), dimension(8,nfout1), intent(in) :: array1
real(SGL), dimension(7,nfout2), intent(in) :: array2


integer :: lowPeriod, highPeriod, n, i, ier

real(SGL), allocatable, dimension(:) :: x, y, z, snr



lowPeriod = nint(minval(array1(2,1:nfout1)))
highPeriod = nint(maxval(array1(2,1:nfout1)))
n = highPeriod - lowPeriod + 1


allocate(x(n), y(n), z(n), snr(n))

do i = 1, n, 1
   x(i) = lowPeriod + i - 1
end do

call linear_interpo(array1(2,1:nfout1), array1(3,1:nfout1), x, y, ier)
call linear_interpo(array1(2,1:nfout1), array1(4,1:nfout1), x, z, ier)
call linear_interpo(array1(2,1:nfout1), array1(7,1:nfout1), x, snr, ier)

open(unit=33, file=trim(adjustl(filename))//'_1', status='replace', action='write')
   do i = 1, n, 1
      write(33, "(I5,2F10.4,F12.4)") int(x(i)), y(i), z(i), snr(i)
   end do
close(unit=33)

deallocate(x, y, z, snr)

! ==============================================================
lowPeriod = nint(minval(array2(2,1:nfout2)))
highPeriod = nint(maxval(array2(2,1:nfout2)))
n = highPeriod - lowPeriod + 1

allocate(x(n), y(n), z(n), snr(n))

do i = 1, n
   x(i) = lowPeriod + i - 1
end do

call linear_interpo(array2(2,1:nfout2), array2(3,1:nfout2), x, y, ier)
call linear_interpo(array2(2,1:nfout2), array2(4,1:nfout2), x, z, ier)
call linear_interpo(array2(2,1:nfout2), array2(6,1:nfout2), x, snr, ier)

open(unit=33, file=trim(adjustl(filename))//'_2', status='replace', action='write')
   do i = 1, n, 1
      write(33, "(I5,2F10.4,F12.4)") int(x(i)), y(i), z(i), snr(i)
   end do
close(unit=33)


deallocate(x, y, z, snr)


end subroutine write_data


! ==============================================================
! Sort the array to make sure the array value increases
! mononically.
! array_x: x array (e.g., period) [input and output]
! array_y: y array (e.g., phase velocity) [input and output]
! ==============================================================
subroutine sort_array_xy(array_x, array_y)

implicit none

real(SGL), dimension(:), intent(inout) :: array_x, array_y


integer :: n, i, j, ip

real(SGL) :: temp, ip_x, ip_y


n = size(array_x)
do i = 1, n-1, 1

   ip = i
   ip_x = array_x(i)
   do j = i+1, n, 1
      if (array_x(j) < ip_x) then
         ip = j
         ip_x = array_x(j)
         ip_y = array_y(j)
      end if
   end do

   if (i /= ip) then
      temp = array_x(i)
      array_x(i) = ip_x
      array_x(ip) = temp
      temp = array_y(i)
      array_y(i) = ip_y
      array_y(ip) = temp
   end if

end do


end subroutine sort_array_xy


! ************************************************************
! Check if there are duplicates in the array
! array: input array [input]
! ************************************************************
logical function any_near_dupl(array)

implicit none

integer :: n, i

real(SGL), dimension(:), intent(in) :: array


any_near_dupl = .false.
n = size(array)
if (1 == n) return

do i = 1, n-1, 1

   if (array(i) == array(i+1)) then
      any_near_dupl = .true.
      return
   end if

end do


end function any_near_dupl

! ************************************************************
! Find the index of the value which is non-positive in the array.
! array: input array [input]
! index: obtained index [output]
! ************************************************************
subroutine find_max_nonpos(array, idx)

implicit none


real(SGL), dimension(:), intent(in) :: array

integer, intent(out) :: idx


integer :: n, i

real(SGL) :: max_value


idx = 0
max_value = -huge(max_value)

n = size(array)
do i = 1, n, 1
   if ((array(i) <= 0) .and. (array(i) > max_value)) then
      idx = i
      max_value = array(i)
   end if
end do
!idx = maxloc((array(1:n) <= 0), 1)


end subroutine find_max_nonpos

! ************************************************************
! Linear interpolation
! x0, y0: input data points [input]
! x: inquired x values [input]
! y: computed y values [output]
! ier: status indicator [output]
! ************************************************************
subroutine linear_interpo(x0, y0, x, y, ier)

implicit none

real(SGL), dimension(:), intent(in) :: x0, y0, x

integer, intent(out) :: ier

real(SGL), dimension(:), intent(out) :: y


integer :: n1, n2, k, idx

real(SGL) :: x1, x2, y1, y2

real(SGL), allocatable, dimension(:) :: x0_c, y0_c



ier = 1
n1 = size(x0)
n2 = size(x)

allocate(x0_c(n1), y0_c(n1))

x0_c = x0
y0_c = y0

call sort_array_xy(x0_c, y0_c)
if (any_near_dupl(x0_c)) then
   write(*,"(A)") 'NO duplicates allowed in x0 !'
   return
end if

do k = 1, n2, 1

   call find_max_nonpos(x0_c-x(k), idx)

   if (0 == idx) then
      x1 = x0_c(1)
      y1 = y0_c(1)
      x2 = x0_c(2)
      y2 = y0_c(2)
      y(k) = (y2-y1)/(x2-x1)*(x(k)-x1)+y1
   else if (n1 == idx) then
      x1 = x0_c(n1-1)
      y1 = y0_c(n1-1)
      x2 = x0_c(n1)
      y2 = y0_c(n1)
      y(k) = (y2-y1)/(x2-x1)*(x(k)-x1)+y1
   else
      x1 = x0_c(idx)
      y1 = y0_c(idx)
      x2 = x0_c(idx+1)
      y2 = y0_c(idx+1)
      y(k) = (y2-y1)/(x2-x1)*(x(k)-x1)+y1
   end if

end do

deallocate(x0_c, y0_c)


ier = 0


end subroutine linear_interpo


end module dispio_m
