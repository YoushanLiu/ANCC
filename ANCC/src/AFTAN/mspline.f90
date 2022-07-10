! =======================================================================
! estimate peacewise cubic spline coefficients
! =======================================================================
module spline_m


implicit none


integer(4), dimension(:), allocatable :: nn, ndim

real(4), dimension(:,:), allocatable :: xx
real(4), dimension(:,:,:), allocatable :: sc, scc


contains


subroutine mspline(ip, n, x, y, ind1, d1, ind2, d2)

implicit none

integer(4), intent(in) :: ip, n
integer(4), intent(in) :: ind1, ind2

real(4), intent(in) :: d1, d2

real(4), intent(in) :: x(n), y(n)


integer(4) i, j, ier

real(4) s, z

real(4), dimension(:,:), allocatable :: c, cc

allocate(nn(1:n), ndim(1:n), stat=ier)

allocate(xx(1:4,1:n), sc(1:4,1:n,1:10), scc(1:5,1:n,1:10), stat=ier)



allocate(c(1:4,1:n), cc(1:5,1:n), stat=ier)


! validate input parameters
!if ((n < 2) .or. (n > 1001)) then
!   write(*, *) 'mspline: n=', n, ' is out of bounds. Should be in range 2-1001.'
!   stop
!endif
if ((ip < 1) .or. (ip > 10)) then
   write(*, *) 'mspline: ip=', ip, ' is out of bounds. Should be in range 1-10.'
   stop
endif
if ((ind1 < 0) .or. (ind1 > 2)) then
   write(*, *) 'mspline: ind1=', ind1, ' is out of bounds. Should be in range 0-2.'
   stop
endif
if ((ind2 < 0) .or. (ind2 > 2)) then
   write(*, *) 'mspline: ind2=', ind2, ' is out of bounds. Should be in range 0-2.'
   stop
endif


do i = 1, n, 1
   c(1,i) = y(i)
enddo
c(2,1) = d1
c(2,n) = d2
call cubspl(x, n, ind1, ind2, c)


! remove 2 and 6
do i = 1, n-1, 1
   c(3,i) = 0.50*c(3,i)
   c(4,i) = c(4,i)/6.0
enddo

! create integral coefficients
s = 0.0
cc = 0.0
do i = 1, n-1, 1
   do j = 1, 4, 1
      cc(j+1,i) = c(j,i)/dble(j)
   enddo
   cc(1,i) = s
   z = x(i+1) - x(i)
   s = s + z*(cc(2,i) + z*(cc(3,i) + z*(cc(4,i) + z*cc(5,i))))
enddo
! save coefficients into common block mspdat
nn(ip) = n
do i = 1, n, 1
   xx(i,ip) = x(i)
   do j = 1, 4, 1
      sc(j,i,ip) = c(j,i)
      scc(j,i,ip) = cc(j,i)
   enddo
   scc(5,i,ip) = cc(5,i)
enddo


deallocate(c, cc)


return


end subroutine mspline


! =======================================================================
! integral for spline
! =======================================================================
subroutine msplint(ip, sa, sb, sint, ier)

implicit none

integer(4), intent(in) :: ip

real(4), intent(in) :: sa, sb


integer(4), intent(out) :: ier

real(4), intent(out) :: sint


integer(4) i, j, ii, n

real(4) s1, s2, z

real(4), dimension(:), allocatable :: x
real(4), dimension(:,:), allocatable :: cc



ier = 0
! restore polinomial coefficients
n = nn(ip)
allocate(x(1:n), cc(1:5,1:n), stat=ier)
do i = 1, n, 1
   x(i) = xx(i,ip)
   do j = 1, 5, 1
      cc(j,i) = scc(j,i,ip)
   enddo
enddo
! compute integral for sa
ii = 0
do i = 1, n-1, 1
   if ((x(i) - sa)*(x(i+1) - sa) <= 0.0) then
      ii = i
      exit
   end if
enddo
if (0 == ii) then
   if (sa <= x(1)) then
      ii = 1
   end if
   if (sa >= x(n-1)) then
      ii = n - 1
   end if
endif
z = sa - x(ii)
s1 = cc(1,ii) + z*(cc(2,ii) + z*(cc(3,ii) + z*(cc(4,ii) + z*cc(5,ii))))
! compute integral for sb
ii = 0
do i = 1, n-1, 1
   if ((x(i) - sb)*(x(i+1) - sb) <= 0.0) then
      ii = i
      exit
   end if
enddo
if (0 == ii) then
   if (sb <= x(1)) then
      ii = 1
   end if
   if (sb >= x(n-1)) then
      ii = n - 1
   end if
endif
z = sb - x(ii)
s2 = cc(1,ii) + z*(cc(2,ii) + z*(cc(3,ii) + z*(cc(4,ii) + z*cc(5,ii))))
sint = s2 - s1


deallocate(cc, x)


return


end subroutine msplint
! =======================================================================
! spline interpolation of funcion and its 1st and 2nd derivatives
! =======================================================================
subroutine msplder(ip, xt, s, sd, sdd, ier)

implicit none

integer(4), intent(in) :: ip

real(4), intent(in) :: xt


integer(4), intent(out) :: ier

real(4), intent(out) :: s, sd, sdd


integer(4) i, j, ii, n

real(4) z

real(4), dimension(:), allocatable :: x
real(4), dimension(:,:), allocatable :: c




ier = 0
! restore polinomial coefficients
n = nn(ip)
allocate(x(1:n), c(1:4,1:n), stat=ier)
do i = 1, n, 1
   x(i) = xx(i,ip)
   do j = 1, 4, 1
      c(j,i) = sc(j,i,ip)
   enddo
enddo

! find interval for interpolation
ii = 0
do i = 1, n-1, 1
   if ((x(i) - xt)*(x(i+1) - xt) <= 0.0) then
      ii = i
      exit
   end if
enddo

if (0 == ii) then
   if (xt <= x(1)) then
      ii = 1
   end if
   if (xt >= x(n-1)) then
      ii = n - 1
   end if
endif

z = xt - x(ii)
s = c(1,ii) + z*(c(2,ii) + z*(c(3,ii) + z*c(4,ii)))
sd = c(2,ii) + z*(2.0*c(3,ii) + 3.0*z*c(4,ii))
sdd = 2.0*c(3,ii) + 6.0*z*c(4,ii)


deallocate(x, c)


return


end subroutine msplder
! =======================================================================
subroutine cubspl(tau, n, ibcbeg, ibcend, c)

!  from  * a practical guide to splines *  by c. de boor
!     ************************  input  ***************************
!     n = number of data points. assumed to be greater than 1.
!     (tau(i), c(1,i), i=1,...,n) = abscissae and ordinates of the
!        data points. tau is assumed to be strictly increasing.
!     ibcbeg, ibcend = boundary condition indicators, and
!     c(2,1), c(2,n) = boundary condition information. Specifically,
!        ibcbeg = 0  means no boundary condition at tau(1) is given.
!           in this case, the not-a-knot condition is used, i.e. the
!           jump in the third derivative across tau(2) is forced to
!           zero, thus the first and the second cubic polynomial pieces
!           are made to coincide.)
!        ibcbeg = 1  means that the slope at tau(1) is made to equal
!           c(2,1), supplied by input.
!        ibcbeg = 2  means that the second derivative at tau(1) is
!           made to equal c(2,1), supplied by input.
!        ibcend = 0, 1, or 2 has analogous meaning concerning the
!           boundary condition at tau(n), with the additional infor-
!           mation taken from c(2,n).
!     ***********************  output  **************************
!     c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
!        of the cubic interpolating spline with interior knots (or
!        joints) tau(2), ..., tau(n-1). precisely, in the interval
!        (tau(i), tau(i+1)), the spline f is given by
!           f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
!        where h = x - tau(i). the function program *ppvalu* may be
!        used to evaluate f or its derivatives from tau,c, l = n-1,
!        and k=4.

implicit none

integer(4), intent(in) :: n

real(4), intent(in) :: tau(n)

integer(4), intent(in) :: ibcbeg, ibcend


real(4), intent(out) :: c(4,n)


integer(4) i, j, l, m

real(4)  divdf1, divdf3, dtau, g


! a tridiagonal linear system for the unknown slopes s(i) of
!  f  at tau(i), i=1,...,n, is generated and then solved by gauss elimination,
!  with s(i) ending up in c(2,i), all i.
!     c(3,.) and c(4,.) are used initially for temporary storage.

l = n-1


! compute first differences of tau sequence and store in c(3,.). also,
! compute first divided difference of data and store in c(4,.).
do m = 2, n, 1
   c(3,m) = tau(m) - tau(m-1)
   c(4,m) = (c(1,m) - c(1,m-1))/c(3,m)
enddo


if (1 == ibcbeg) then

   c(4,1) = 1.0
   c(3,1) = 0.0

   if (2 == n) then

      if (ibcend > 1) then

         c(2,n) = 3.0*c(4,n) + 0.50*c(3,n)*c(2,n)
         c(4,n) = 2.0
         g = -1.0/c(4,n-1)
         c(4,n) = g*c(3,n-1) + c(4,n)
         c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)

      else if (ibcend < 1) then

         if (ibcbeg > 0) then
            c(2,n) = 2.0*c(4,n)
            c(4,n) = 1.0
            g = -1.0/c(4,n-1)
            c(4,n) = g*c(3,n-1) + c(4,n)
            c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
         else
            ! cannot arrive
            c(2,n) = c(4,n)
         end if

      end if

   else if (2 /= n) then

      do m = 2, l, 1
         g = -c(3,m+1)/c(4, m-1)
         c(2,m) = g*c(2,m-1) + 3.0*(c(3,m)*c(4,m+1) + c(3,m+1)*c(4,m))
         c(4,m) = g*c(3,m-1) + 2.0*(c(3,m) + c(3,m+1))
      enddo

      if (ibcend > 1) then
         c(2,n) = 3.0*c(4,n) + 0.50*c(3,n)*c(2,n)
         c(4,n) = 2.0
         g = -1.0/c(4,n-1)
      else if (ibcend < 1) then
         if ((3 == n) .and. (0 == ibcbeg)) then
            ! cannot arrive
            c(2,n) = 2.0*c(4,n)
            c(4,n) = 1.0
            g = -1.0/c(4,n-1)
         else
            g = c(3,n-1) + c(3,n)
            c(2,n) = ((c(3,n) + 2.0*g)*c(4,n)*c(3,n-1) + c(3,n)*c(3,n)*(c(1,n-1) - c(1,n-2))/c(3,n-1))/g
            g = -g/c(4,n-1)
            c(4,n) = c(3,n-1)
         end if
      end if
      c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)

   end if

elseif (ibcbeg > 1) then

   c(4,1) = 2.0
   c(3,1) = 1.0
   c(2,1) = 3.0*c(4,2) - 0.50*c(3,2)*c(2,1)

   if (2 == n) then

      if (ibcend > 1) then
         c(2,n) = 3.0*c(4,n) + 0.50*c(3,n)*c(2,n)
         c(4,n) = 2.0
         g = -1.0/c(4,n-1)
         c(4,n) = g*c(3,n-1) + c(4,n)
         c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
      else if (ibcend < 1) then
         if (ibcbeg > 0) then
            c(2,n) = 2.0*c(4,n)
            c(4,n) = 1.0
            g = -1.0/c(4,n-1)
            c(4,n) = g*c(3,n-1) + c(4,n)
            c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
         else
            ! cannot arrive
            c(2,n) = c(4,n)
         end if
      end if

   else if (2 /= n) then

      do m = 2, l, 1
         g = -c(3,m+1)/c(4,m-1)
         c(2,m) = g*c(2,m-1) + 3.0*(c(3,m)*c(4,m+1) + c(3,m+1)*c(4,m))
         c(4,m) = g*c(3,m-1) + 2.0*(c(3,m) + c(3,m+1))
      enddo

      if (ibcend > 1) then
         c(2,n) = 3.0*c(4,n) + 0.50*c(3,n)*c(2,n)
         c(4,n) = 2.0
         g = -1.0/c(4,n-1)
      else if (ibcend < 1) then
         if ((3 == n) .and. (0 == ibcbeg)) then
            ! cannot arrive
            c(2,n) = 2.0*c(4,n)
            c(4,n) = 2.0
            g = -1.0/c(4,n-1)
         else
            g = c(3,n-1) + c(3,n)
            c(2,n) = ((c(3,n) + 2.0*g)*c(4,n)*c(3,n-1) + c(3,n)*c(3,n)*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
            g = -g/c(4,n-1)
            c(4,n) = c(3,n-1)
         end if
      end if
      c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)

   end if

else if (ibcbeg < 1) then

   if (n > 2) then

      c(4,1) = c(3,3)
      c(3,1) = c(3,2) + c(3, 3)
      c(2,1) = ((c(3,2) + 2.0*c(3,1))*c(4,2)*c(3,3) + c(3,2)*c(3,2)*c(4,3))/c(3,1)
      do m = 2, l, 1
         g = -c(3,m+1)/c(4, m-1)
         c(2,m) = g*c(2,m-1) + 3.0*(c(3,m)*c(4,m+1) + c(3,m+1)*c(4,m))
         c(4,m) = g*c(3,m-1) + 2.0*(c(3,m) + c(3,m+1))
      enddo

      if (ibcend > 1) then
         c(2,n) = 3.0*c(4,n) + 0.50*c(3,n)*c(2,n)
         c(4,n) = 2.0
         g = -1.0/c(4,n-1)
      else if (ibcend < 1) then
         if ((3 == n) .and. (0 == ibcbeg)) then
            ! cannot arrive
            c(2,n) = 2.0*c(4,n)
            c(4,n) = 1.0
            g = -1.0/c(4,n-1)
         else
            g = c(3,n-1) + c(3,n)
            c(2,n) = ((c(3,n) + 2.0*g)*c(4,n)*c(3,n-1) + c(3,n)*c(3,n)*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
            g = -g/c(4,n-1)
            c(4,n) = c(3,n-1)
         end if
      end if
      c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)

   else

      c(4,1) = 1.0
      c(3,1) = 1.0
      c(2,1) = 2.0*c(4,2)
      c(2,n) = c(4,n)

      if (ibcend > 1) then
         c(2,n) = 3.0*c(4,n) + 0.50*c(3,n)*c(2,n)
         c(4,n) = 2.0
         g = -1.0/c(4,n-1)
         c(4,n) = g*c(3,n-1) + c(4,n)
         c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
      else if (ibcend < 1) then
         if (ibcbeg > 0) then
            ! cannot arrive
            c(2,n) = 2.0*c(4,n)
            c(4,n) = 1.0
            g = -1.d0/c(4,n-1)
            c(4,n) = g*c(3,n-1) + c(4,n)
            c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
         else
            c(2,n) = c(4,n)
         end if
      end if

   end if

end if


j = l
do while(j > 0)
   c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
   j = j-1
end do
do i = 2, n, 1
   dtau = c(3,i)
   divdf1 = (c(1,i) - c(1,i-1))/dtau
   divdf3 = c(2,i-1) + c(2,i) - 2.0*divdf1
   c(3,i-1) = 2.0*(divdf1 - c(2,i-1) - divdf3)/dtau
   c(4,i-1) = (divdf3/dtau)*(6.0/dtau)
enddo




! ! construct first equation from the boundary condition, of the form
! !             c(4,1)*s(1) + c(3,1)*s(2) = c(2,1)
! !MB if (ibcbeg-1) 11,15,16

! if (0 == (ibcbeg-1)) goto 15
! if ((ibcbeg-1) > 0) goto 16

! if (n > 2) goto 12
! ! no condition at left end and n = 2.
! c(4,1) = 1.d0
! c(3,1) = 1.d0
! c(2,1) = 2.d0*c(4,2)
! goto 25
! ! not-a-knot condition at left end and n > 2.
! 12 c(4,1) = c(3,3)
! c(3,1) = c(3,2) + c(3, 3)
! c(2,1) = ((c(3,2) + 2.d0*c(3,1))*c(4,2)*c(3,3) + &
!           c(3,2)**2*c(4,3))/c(3,1)
! goto 19
! ! slope prescribed at left end.

! 15 c(4,1) = 1.d0
! c(3,1) = 0.d0
! goto 18
! ! second derivative prescribed at left end.
! 16 c(4,1) = 2.d0
! c(3,1) = 1.d0
! c(2,1) = 3.d0*c(4,2) - c(3,2)/2.d0*c(2,1)
! 18 if (2 == n) goto 25
! ! if there are interior knots, generate the corresp. equations and
! ! carry out the forward pass of gauss elimination, after which the m-th
! ! equation reads  c(4,m)*s(m) + c(3,m)*s(m+1) = c(2,m).
! 19 do m = 2, l, 1
!   g = -c(3,m+1)/c(4, m-1)
!   c(2,m) = g*c(2,m-1) + 3.d0*(c(3,m)*c(4,m+1) + c(3,m+1)*c(4,m))
!   c(4,m) = g*c(3,m-1) + 2.d0*(c(3,m) + c(3,m+1))
! enddo
! ! construct last equation from the second boundary condition, of the form
! !           (-g*c(4,n-1))*s(n-1) + c(4,n)*s(n) = c(2,n)
! ! if slope is prescribed at right end, one can go directly to back-
! ! substitution, since c array happens to be set up just right for it
! ! at this point.
! !MB   if (ibcend-1) 21,30,24
! if (0 == (ibcend-1)) goto 30
! if ((ibcend-1) > 0) goto 24

! if ((3 == n) .and. (0 == ibcbeg)) goto 22
! ! not-a-knot and n >= 3, and either n>3 or  also not-a-knot at
! ! left end point.
! g = c(3,n-1) + c(3,n)
! c(2, n) = ((c(3,n) + 2.d0*g)*c(4,n)*c(3,n-1) + &
!            c(3,n)*c(3,n)*(c(1,n-1) - c(1,n-2))/c(3,n-1))/g
! g = -g/c(4,n-1)
! c(4,n) = c(3,n-1)
! goto 29
! ! either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
! ! knot at left end point).
! 22 c(2,n) = 2.d0*c(4,n)
! c(4,n) = 1.d0
! goto 28
! ! second derivative prescribed at right endpoint.
! 24 c(2,n) = 3.d0*c(4,n) + 0.50*c(3,n)*c(2,n)
! c(4, n) = 2.d0
! goto 28
! !MB25 if (ibcend-1) 26,30,24
! 25  if (0 == (ibcend-1)) goto 30
! if ((ibcend-1) > 0) goto 24
! if (ibcbeg > 0) goto 22
! ! not-a-knot at right endpoint and at left endpoint and n = 2.
! c(2,n) = c(4,n)
! goto 30
! 28  g = -1.d0/c(4,n-1)
! ! complete forward pass of gauss elimination.
! 29 c(4,n) = g*c(3,n-1) + c(4,n)
! c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
! ! ccarry out back substitution
! 30 j = l
! 40 c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
! j = j - 1
! if (j > 0) goto 40
! ! generate cubic coefficients in each interval, i.e., the deriv.s
! ! at its left endpoint, from value and slope at its endpoints.

! do i = 2, n, 1
!   dtau = c(3,i)
!   divdf1 = (c(1,i) - c(1,i-1))/dtau
!   divdf3 = c(2,i-1) + c(2,i) - 2.d0*divdf1
!   c(3,i-1) = 2.d0*(divdf1 - c(2,i-1) - divdf3)/dtau
!   c(4,i-1) = (divdf3/dtau)*(6.d0/dtau)
! enddo


return


end subroutine cubspl



subroutine free_mspline()

implicit none

deallocate(xx, sc, scc)

end subroutine free_mspline



real(4) function fsig(x1, x2, r)

real(4), intent(in) :: x1, x2, r

fsig = (x1 - r)*(x2 - r)

return

end function fsig


end module spline_m