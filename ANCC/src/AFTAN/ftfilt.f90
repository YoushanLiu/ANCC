!-------------------------------------------------
! FTAN filter: y = x*exp(-alpha*((om-om0)/om0)**2
!-------------------------------------------------
subroutine ftfilt(alpha, om0, dom, n, a, fs)

implicit none

integer(4), intent(in) :: n

real(4), intent(in) :: alpha, om0, dom

complex, intent(in) :: a(n)


complex, intent(out) :: fs(n)



integer(4) k

complex czero


real(4) ome, om2, b



czero = cmplx(0.0, 0.0)

do k = 1, n, 1

   fs(k) = czero
   b = 0.0
   ome = (k-1)*dom
   om2 = -alpha*(ome - om0)/om0*(ome - om0)/om0

   if (abs(om2) <= 40.0) then
      b = exp(om2)
      fs(k) = a(k)*b
   endif

enddo


return


end subroutine ftfilt
