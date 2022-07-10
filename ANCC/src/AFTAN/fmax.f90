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
!----------------------------------------------------------------------
! parabolic interpolation of signal amplitude and phase,
! finding phase derivative
!----------------------------------------------------------------------
subroutine fmax(am1, am2, am3, ph1, ph2, ph3, PIover4, om, dt, t, dph, tm, ph)

implicit none

real(4), intent(in) :: om, dt
real(4), intent(in) :: am1, am2, am3
real(4), intent(in) :: ph1, ph2, ph3

real(8), intent(in) :: PIover4


real(4), intent(out) :: t, dph, tm, ph


integer(4) k

real(8), parameter :: twoPI = 8.d0*datan(1.d0)

real(4) a1, a2, a3, dd



dd = am1 + am3 - 2.0*am2

t = 0.0
if (0.0 /= dd) then
   t = 0.50*(am1 - am3)/dd
endif

!  phase derivative
a1 = ph1
a2 = ph2
a3 = ph3


!  check for 2*pi phase jump
k = nint((a2 - a1 - om*dt)/twoPI)
a2 = a2 - k*twoPI
k = nint((a3 - a2 - om*dt)/twoPI)
a3 = a3 - k*twoPI


! interpolation
dph = t*(a1 + a3 - 2.0*a2) + 0.50*(a3 - a1)
tm = 0.50*t*t*(am1 + am3 - 2.0*am2) + 0.50*t*(am3 - am1) + am2
ph = 0.50*t*t*(a1 + a3 - 2.0*a2) + 0.50*t*(a3 - a1) + a2 + 0.1250*twoPI*PIover4


return


end subroutine fmax
