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
! - -----------------------------------------------------------
! test dispersion curve for jumps
! - -----------------------------------------------------------
subroutine trigger(grvel, om, nf, tresh, trig, ftrig, ier)

implicit none

integer(4), intent(in) :: nf

real(4), intent(in) :: tresh

real(4), intent(in) :: grvel(nf), om(nf)


integer(4), intent(out) :: ier

real(4), intent(out) :: trig(nf), ftrig(nf)


integer(4) i

real(4) hh1, hh2, hh3, r


ier = 0
ftrig(1) = 0.0
ftrig(nf) = 0.0
trig(1) = 0.0
trig(nf) = 0.0


do i = 1, nf-2, 1

   trig(i+1) = 0.0
   hh1 = om(i+1) - om(i)
   hh2 = om(i+2) - om(i+1)
   hh3 = hh1 + hh2

   r = 25.0*(grvel(i)/hh1 - (1.0/hh1 + 1.0/hh2)*grvel(i+1) + grvel(i+2)/hh2)*hh3

   ftrig(i+1) = r

   if (abs(r) > tresh) then
      trig(i+1) = sign(1.0, r)
      ier = 1
   endif

enddo


return


end subroutine trigger
