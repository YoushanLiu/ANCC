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
module math_m

use db_m


implicit none



contains

! ***************************************************************
! This program is downloaded from http://fcode.cn/guide-96-1.html
! After invoke this subroutine, random_number can produce different
! random numbers in the range of [0 1]. Without using this subroutine,
! random_number will produce the same random numbers every time.
! ***************************************************************
subroutine init_random_seed()


#ifdef __INTEL_COMPILER
use IFPORT
#endif


implicit none


integer ised, i, pid

integer(DBL) t

integer, allocatable, dimension(:) :: sed


call random_seed(size = ised) ! Get the size of the seed


allocate(sed(ised)) ! Distribute the seed
call system_clock(t) ! Get the time
pid = getpid() ! Get the id of the processor
t = ieor(t, int(pid, kind(t))) ! XOR operation

do i = 1, ised, 1
   sed(i) = lcg(t) ! LCG operation
end do


call random_seed(put = sed) ! Set the seed value


end subroutine init_random_seed

! ***************************************************************
! Linear congruential generator
! ***************************************************************
function lcg(s)

integer(DBL), intent(inout) :: s


integer lcg


if (0 == s) then

   s = 104729

else

   s = mod(s, 4294967296_DBL)

end if

s = mod(s * 279470273_DBL, 4294967291_DBL)
lcg = int(mod(s, int(huge(0), DBL)), kind(0))


return


end function lcg


! ***************************************************************
! ***************************************************************
subroutine matrix_mean_std(A, skip_value, nout, mean, std)

implicit none

real(SGL), intent(in) :: skip_value

real(SGL), dimension(:,:), intent(in) :: A

integer, intent(out) :: nout

real(SGL), dimension(:), allocatable, intent(out) :: mean, std


integer i, j, ier
integer n, nrow, ncol

real(SGL) rn, rsum1, rsum2


nrow = size(A,1)
ncol = size(A,2)

nout = nrow

allocate(mean(1:nrow), std(1:nrow), stat=ier)
std = 0.0
mean = 0.0
do i = 1, nrow, 1

   n = 0

   rsum1 = 0.0
   rsum2 = 0.0
   do j = 1, ncol, 1
      if (abs(A(i,j)-skip_value) > epsilon(A(i,j))) then
         rsum1 = rsum1 + A(i,j)
         rsum2 = rsum2 + A(i,j)*A(i,j)
         n = n + 1
      end if
   end do

   rn = real(n)
   if (n >= 1) then
      mean(i) = rsum1 / rn
   end if

   if (n >= 2) then
      std(i) = sqrt(abs(rn*rsum2 - rsum1*rsum1)/(rn*(n-1)))
   end if

end do


end subroutine matrix_mean_std


end module math_m
