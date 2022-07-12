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
!-----------------------------------------------------------------
! convert phase to phase velocity
!-----------------------------------------------------------------
subroutine phtovel(delta, ip, nf, per, U, pha, npr, prper, prvel, n, V)

!
! delta  - distance, km
! n      - # of points in per
! per    - apparent periods
! U      - group velocity
! pha    - phase
! npr    - # of points in prper
! prper  - pedicted periods
! prvel  - pedicted phase velocity
! V      - observed phase velocity

use spline_m


implicit none

integer(4), intent(in) :: ip, n, nf, npr

real(8), intent(in) :: delta

real(8), intent(in) :: per(nf), U(nf), pha(nf)
real(8), intent(in) :: prper(npr), prvel(npr)


real(8), intent(out) :: V(n)


integer(4) i, k, m, ier

real(8) Vpred, phpred, s, ss

real(8), parameter :: PI = 4.d0*datan(1.d0)

real(8), dimension(:), allocatable :: om, sU, t



allocate(om(1:n), sU(1:n), t(1:n), stat = ier)


Vpred = 0.d0
do i = 1, n, 1
   om(i) = 2.d0*PI/per(i)
   sU(i) = 1.d0/U(i)
   t(i) = delta*sU(i)
   V(i) = 0.d0
enddo


! find velocity for the largest period by spline interpolation
! with not a knot boundary conditions
call mspline(ip+2, npr, prper, prvel, 0, 0.d0, 0, 0.d0)
call msplder(ip+2, per(n), Vpred, s, ss, ier)
call free_mspline()


phpred = om(n)*(t(n) - delta/Vpred)
k = nint(0.5d0*(phpred - pha(n))/PI)
V(n) = delta/(t(n) - (pha(n) + 2.d0*k*PI)/om(n))

do m = n-1, 1, -1

   Vpred = 1.d0/((0.5d0*(sU(m) + sU(m+1))*(om(m) - om(m+1)) + om(m+1)/V(m+1))/om(m))
   phpred = om(m)*(t(m) - delta/Vpred)
   k = nint(0.5d0*(phpred - pha(m))/PI)
   V(m) = delta/(t(m) - (pha(m) + 2.d0*k*PI)/om(m))

enddo


deallocate(om, sU, t)


return


end subroutine phtovel
