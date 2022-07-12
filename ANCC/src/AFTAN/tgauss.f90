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
!-------------------------------------------------------------
! taper phase matched signal
!-------------------------------------------------------------
subroutine tgauss(fsnr, gt0, dw, dt, n, fmatch, seis, ss)

implicit none

integer(4), intent(in) :: n

real(8), intent(in) :: fsnr, gt0, dw, dt, fmatch

complex(8), intent(in) :: seis(n)


complex(8), intent(out) :: ss(n)


integer(4) i, ii, ism, nn, nnn
integer(4) nnl, nnnl, nnr, nnnr
integer(4) nc, nleft, nright, ier

real(8) freq, zero, dl, dr
real(8) sm, smw, tresh, tre, coef

complex(8), parameter :: czero = cmplx(0.d0,0.d0)

integer(4), dimension(:), allocatable :: left, right

real(8), dimension(:), allocatable :: smax, vleft, vright

allocate(left(1:n), right(1:n), stat=ier)
allocate(smax(1:n), vleft(1:n), vright(1:n), stat=ier)



ism = 1

zero = 0.d0
nc = nint(gt0/dt) + 1


! find global max, sm, and index ism
sm = 0.d0
do i = 1, n, 1
   smw = abs(seis(i))
   if (smw >= sm) then
      sm = smw
      ism = i
   endif
   smax(i) = smw
   ss(i) = seis(i)
enddo
coef = fsnr*sm


! find some local minima,# < 100 from left and right side of central max ism

! left side
nleft = 0
do i = ism-1, 2, -1

   dl = smax(i) - smax(i-1)
   dr = smax(i+1) - smax(i)

   if (((dl < zero) .and. (dr >= zero)) .or. &
       ((dl <= zero) .and. (dr > zero))) then
      nleft = nleft + 1
      left(nleft) = i
      vleft(nleft) = smax(i)
   endif

enddo


! right side
nright = 0
do i = ism+1, n-1, 1

   dl = smax(i) - smax(i-1)
   dr = smax(i+1) - smax(i)

   if (((dl < zero) .and. (dr >= zero)) .or. &
       ((dl <= zero) .and. (dr > zero))) then
      nright = nright + 1
      right(nright) = i
      vright(nright) = smax(i)
   endif

enddo


! left side, apply cutting
ii = 0
nnl = 0
nnnl = 0
if (0 /= nleft) then
   do i = 1, nleft, 1
      if (abs(ism - left(i))*dt > 5.d0) then
         if (vleft(i) < coef) then
            nnnl = left(i)
            ii = i
            exit
         endif
      endif
   enddo
end if
if (0 /= nnnl) then
   if (ii /= nleft) then
      nnl = left(ii + 1)
   else
      nnl = 1
   endif
endif


! right side, apply cutting
ii = 0
nnr = 0
nnnr = 0
if (0 /= nright) then
   do i = 1, nright, 1
      if (abs(ism - right(i))*dt > 5.d0) then
         if (vright(i) < coef) then
            nnr = right(i)
            ii = i
            exit
         endif
      endif
   enddo
end if
if (0 /= nnr) then
   if (ii /= nright) then
      nnnr = right(ii+1)
   else
      nnnr = n
   endif
endif


! ---
if ((0 /= nnnr) .and. (0 /= nnnl)) then
   nn = max(abs(ism-nnnl), abs(ism-nnr))
   nnn = max(abs(nnnl-nnl), abs(nnnr-nnr))
   nnnl = ism - nn
   nnl = nnnl - nnn
   nnr = ism + nn
   nnnr = nnr + nnn
endif

! setup cutting point for gaussian
tresh = log(sm) - 24.d0
if (((0 /= nnl) .and. (0 /= nnnl))) then
   ! expand left end by factor fmatch
   nnl = nint((nnl-ism)*fmatch) + ism
   nnnl = nint((nnnl-ism)*fmatch) + ism
   nnl = max(1, nnl)
   nnnl = max(1, nnnl)
   freq = (nnnl - nnl) + 1
   do i = 1, nnnl, 1
      tre = -0.5d0*(i-nnnl)/freq*(i-nnnl)/freq
      if (tre > tresh) then
         ss(i) = ss(i)*exp(tre)
      else
         ss(i) = czero
      endif
   enddo
end if


if (((0 /= nnr) .and. (0 /= nnnr))) then
   ! expand right end by factor fmatch
   nnr = nint((nnr - ism)*fmatch) + ism
   nnnr = nint((nnnr - ism)*fmatch) + ism
   nnr = min(n, nnr)
   nnnr = min(n, nnnr)
   freq = (nnnr - nnr) + 1
   do i = nnr, n, 1
      tre = -0.5d0*(i-nnr)/freq*(i-nnr)/freq
      if (tre > tresh) then
         ss(i) = ss(i)*exp(tre)
      else
         ss(i) = czero
      endif
   enddo
end if


deallocate(left, right)
deallocate(smax, vleft, vright)


return


end subroutine tgauss
