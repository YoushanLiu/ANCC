!--------------------------------------------------------------
! create phase prediction curve by group velocity
!--------------------------------------------------------------
subroutine pred_cur(ip, delta, om0, npred, pred, om1, Gt0)

use spline_m

implicit none

integer(4), intent(in) :: ip, npred

real(4), intent(in) :: delta, om0

real(4), intent(in) :: pred(npred,2)


real(4), intent(out) :: om1, Gt0


integer(4) i, ier

real(8), parameter :: PI = 4.d0*datan(1.d0)

real(4) s, ss

real(4), dimension(:), allocatable :: x, y


allocate(x(1:npred), y(1:npred), stat=ier)


! transfer prediction curve (T,vel) ==> (omrga,t)
do i = 1, npred, 1
   x(i) = 2.0*PI/pred(npred-i+1,1)
   y(i) = delta/pred(npred-i+1,2)
enddo

! get velocity for low integral boundary
call mspline(ip, npred, x, y, 0, 0.0, 0, 0.0)
call msplder(ip, om0, Gt0, s, ss, ier)

call free_mspline()


! construct spline for ph filter
do i = 1, npred, 1
   y(i) = y(i) - Gt0
enddo

call mspline(ip+1, npred, x, y, 0, 0.0, 0, 0.0)


om1 = om0


deallocate(x, y)


return


end subroutine pred_cur
