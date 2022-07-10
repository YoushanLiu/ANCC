!-----------------------------------------------------------------
! convert phase to phase velocity
!-----------------------------------------------------------------
subroutine phtovel(delta, ip, n, per, U, pha, npr, prper, prvel, V)

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

integer(4), intent(in) :: ip, n, npr

real(4), intent(in) :: delta

real(4), intent(in) :: per(n), U(n), pha(n)
real(4), intent(in) :: prper(n), prvel(n)


real(4), intent(out) :: V(n)


integer(4) i, k, m, ier

real(4) Vpred, phpred, s, ss

real(8), parameter :: PI = 4.d0*datan(1.d0)

real(4), dimension(:), allocatable :: om, sU, t



allocate(om(1:n), sU(1:n), t(1:n), stat = ier)


Vpred = 0.0
do i = 1, n, 1
   om(i) = 2.0*PI/per(i)
   sU(i) = 1.0/U(i)
   t(i) = delta*sU(i)
   V(i) = 0.0
enddo


! find velocity for the largest period by spline interpolation
! with not a knot boundary conditions
call mspline(ip+2, npr, prper, prvel, 0, 0.0, 0, 0.0)
call msplder(ip+2, per(n), Vpred, s, ss, ier)


phpred = om(n)*(t(n) - delta/Vpred)
k = nint(0.50*(phpred - pha(n))/PI)
V(n) = delta/(t(n) - (pha(n) + 2.0*k*PI)/om(n))

do m = n-1, 1, -1

   Vpred = 1.0/((0.50*(sU(m) + sU(m+1))*(om(m) - om(m+1)) + om(m+1)/V(m+1))/om(m))
   phpred = om(m)*(t(m) - delta/Vpred)
   k = nint(0.50*(phpred - pha(m))/PI)
   V(m) = delta/(t(m) - (pha(m) + 2.0*k*PI)/om(m))

enddo


call free_mspline()


deallocate(om, sU, t)


return


end subroutine phtovel
