! This file is part of ANCC.
!
! AND is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! AND is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
! -----------------------------------------------------------------------------------------
! butterworth filter
! Author: Youshan Liu
! Date: 5/25/2016
!
! design the order of the butterworth filters
!
! buttordlp --- the order of the lowpass butterworth filter
! buttordhp --- the order of the highpass butterworth filter
! buttordbp --- the order of the bandpass butterworth filter
! buttordbs --- the order of the bandstop butterworth filter
!
! compute the filter coefficients of the butterworth filters
! buttlp --- the coefficients of the lowpass butterworth filter
! butthp --- the coefficients of the highpass butterworth filter
! buttbp --- the coefficients of the bandpass butterworth filter
! buttbs --- the coefficients of the bandstop butterworth filter
!
! buttflt --- IIR filter
!
! 3  dB--->0.5
! 10 dB--->1.d-1
! 20 dB--->1.d-2
! 30 dB--->1.d-3
! 40 dB--->1.d-4
! 50 dB--->1.d-5
! 60 dB--->1.d-30
! -----------------------------------------------------------------------------------------
!
! buttordlp
!
subroutine buttordlp(fp, fs, rp, rs, order)

implicit none

real(8), intent(inout) :: rp, rs
real(8), intent(inout) :: fp, fs

integer(4), intent(out) :: order


real(8) wa, w0


if (fs < fp) then
   write(*,"(A)") "Error: fs must be great than fp in the lowpass butterworth filter !"
   stop
end if

if (abs(fp) < 1.d-30) then
   write(*,"(A)") "Error: fp must be great than 0 in the lowpass butterworth filter !"
   stop
end if

rp = abs(rp)
rs = abs(rs)
if (abs(rp) < 1.d-10) then
   rp = 3.d0
end if

if (abs(rs) < 1.d-10) then
   rs = 40.d0
end if


wa = fs/fp

order = ceiling(abs(dlog10((10**(0.1d0*rs)-1.d0)/(10**(0.1d0*rp)-1.d0)) / &
                                                      (2.d0*dlog10(wa))))

w0 = wa/((10**(0.1d0*abs(rs)) - 1.d0)**(1.d0/(2.d0*(abs(order)))))
fp = fp*w0


end subroutine buttordlp
!
! buttordhp
!
subroutine buttordhp(fp, fs, rp, rs, order)

implicit none

real(8), intent(inout) :: rp, rs
real(8), intent(inout) :: fp, fs

integer(4), intent(out) :: order


real(8) wa, w0


if (fp < fs) then
   write(*,"(A)") "Error: fp must be great than fs in the highpass butterworth filter !"
   stop
end if

if (abs(fs) < 1.d-30) then
   write(*,"(A)") "Error: fs must be great than 0 in the highpass butterworth filter !"
   stop
end if

rp = abs(rp)
rs = abs(rs)
if (abs(rp) < 1.d-10) then
   rp = 3.d0
end if

if (abs(rs) < 1.d-10) then
   rs = 40.d0
end if

wa = fp/fs

order = ceiling(abs(log10((10**(0.1d0*rs)-1.d0)/(10**(0.1d0*rp)-1.d0)) / &
                                                     (2.d0*dlog10(wa))))

w0 = wa/((10**(0.1d0*abs(rs)) - 1.d0)**(1.d0/(2.d0*(abs(order)))))
fp = fp/w0


end subroutine buttordhp
!
! buttordbp
!
subroutine buttordbp(fp1, fp2, fs1, fs2, rp, rs, order)

implicit none

real(8), intent(inout) :: rp, rs
real(8), intent(inout) :: fp1, fp2
real(8), intent(inout) :: fs1, fs2


integer(4), intent(out) :: order


real(8) fa
real(8) fn(1:2)
real(8) f0(1:2)
real(8) fp(1:2)
real(8) fs(1:2)
real(8) ff(1:2)


if (fp2 < fp1) then
   write(*,"(A)") "Error: fp2 must be great than fp1 in the bandpass butterworth filter !"
   stop
end if

if (fs2 < fs1) then
   write(*,"(A)") "Error: fs2 must be great than fs1 in the bandpass butterworth filter !"
   stop
end if

if (fp1 < fs1) then
   write(*,"(A)") "Error: fp1 must be great than fs1 in the bandpass butterworth filter !"
   stop
end if

if (fs2 < fp2) then
   write(*,"(A)") "Error: fs2 must be great than fp2 in the bandpass butterworth filter !"
   stop
end if

rp = abs(rp)
rs = abs(rs)
if (abs(rp) < 1.d-10) then
   rp = 3.d0
end if

if (abs(rs) < 1.d-10) then
   rs = 40.d0
end if

fp(1) = fp1
fp(2) = fp2
fs(1) = fs1
fs(2) = fs2
ff = (fs**2-fp1*fp2)/(fs*(fp1-fp2))
fa = minval(abs(ff))

order = ceiling(abs(log10((10**(0.1d0*rs)-1.d0)/(10**(0.1d0*rp)-1.d0)) / &
                                                     (2.d0*dlog10(fa))))

fa = fa/((10**(0.1d0*abs(rs)) - 1.d0)**(1.d0/(2.d0*(abs(order)))))
f0(1) = -fa
f0(2) =  fa
fn = -0.5d0*f0*(fp2-fp1) + dsqrt(0.25d0*f0**2*(fp2-fp1)**2+fp1*fp2)
fp1 = minval(abs(fn))
fp2 = maxval(abs(fn))

end subroutine buttordbp
!
! buttordbs
!
subroutine buttordbs(fp1, fp2, fs1, fs2, rp, rs, order)

implicit none

real(8), intent(inout) :: rp, rs
real(8), intent(inout) :: fp1, fp2
real(8), intent(inout) :: fs1, fs2

integer(4), intent(out) :: order


real(8) fa, f0
real(8) fn(1:2)
real(8) fp(1:2)
real(8) fs(1:2)
real(8) ff(1:2)


if (fp2 < fp1) then
   write(*,"(A)") "Error: fp2 must be great than fp1 in the bandstop butterworth filter !"
   stop
end if

if (fs2 < fs1) then
   write(*,"(A)") "Error: fs2 must be great than fs1 in the bandstop butterworth filter !"
   stop
end if

if (fs1 < fp1) then
   write(*,"(A)") "Error: fs1 must be great than fp1 in the bandstop butterworth filter !"
   stop
end if

if (fp2 < fs2) then
   write(*,"(A)") "Error: fp2 must be great than fs2 in the bandstop butterworth filter !"
   stop
end if

rp = abs(rp)
rs = abs(rs)
if (abs(rp) < 1.d-10) then
   rp = 3.d0
end if

if (abs(rs) < 1.d-10) then
   rs = 40.d0
end if

fp(1) = fp1
fp(2) = fp2
fs(1) = fs1
fs(2) = fs2
ff = (fs*(fp1-fp2))/(fs**2-(fp1-fp2))
fa = maxval(abs(ff))
order = ceiling(abs(dlog10((10**(0.1d0*rs)-1.d0)/(10**(0.1d0*rp)-1.d0)) / &
                                                      (2.d0*dlog10(fa))))

f0 = fa/((10**(0.1d0*abs(rs)) - 1.d0)**(1.d0/(2.d0*(abs(order)))))
fn(1) = ((fp2-fp1)+dsqrt((fp2-fp1)**2+4.d0*f0**2*fp1*fp2))/(2.d0*f0)
fn(2) = ((fp2-fp1)-dsqrt((fp2-fp1)**2+4.d0*f0**2*fp1*fp2))/(2.d0*f0)

fp1 = minval(abs(fn))
fp2 = maxval(abs(fn))

end subroutine buttordbs
!
! buttlp
!
subroutine buttlp(order, dt, f, a, b)

implicit none

integer(4), intent(in) :: order

real(8), intent(in) :: dt, f
real(8), intent(inout) :: a(0:order)
real(8), intent(inout) :: b(0:order)


integer(4) k, n, m

real(8) r, r2, p, rp
real(8) dorder, den, PI

real(8), dimension(:,:), allocatable :: c, d


if (abs(f) < 1.d-30) then
   write(*,"(A)") "Error: f must be great than 0 in the lowpass butterworth filter !"
   stop
end if

m = int(order/2)
n = int((order+1)/2)
dorder = dble(order)
PI = 4.d0*datan(1.d0)
r = 1.d0/dtan(PI*f*dt)
r2 = r*r

allocate(c(1:3,0:n-1))
allocate(d(1:3,0:n-1))

a = 0.d0
b = 0.d0
c = 0.d0
d = 0.d0
do k = 0, m-1, 1
   p = -dsin(0.5d0*PI*(2*k+1)/dorder)
   rp = r*p
   den = r2-2.d0*rp+1.d0
   c(1,k) = (r2+2.0*rp+1.d0)/den
   c(2,k) = 2.d0*(1.d0-r2)/den
   c(3,k) = 1.d0
   d(1,k) = 1.d0/den
   d(2,k) = 2.d0/den
   d(3,k) = 1.d0/den
end do
do k = m, n-1, 1
   den = 1.d0+r
   c(1,k) = (1.d0-r)/den
   c(2,k) = 1.d0
   d(1,k) = 1.d0/den
   d(2,k) = 1.d0/den
end do

call polyval(order, n, m, c, d, a, b)

deallocate(c, d)

end subroutine buttlp
!
! butthp
!
subroutine butthp(order, dt, f, a, b)

implicit none

integer(4), intent(in)::order

real(8), intent(in) :: dt, f
real(8), intent(inout) :: a(0:order)
real(8), intent(inout) :: b(0:order)


integer(4) k, n, m

real(8) r, r2, p, rp
real(8) dorder, den, PI


real(8), dimension(:,:), allocatable :: c, d


if (abs(f) > 0.5d0/dt) then
   write(*,"(A)") "Error: fp must be less than the Nyquist frequency in the highpass butterworth filter !"
   stop
end if

m = int(order/2)
n = int((order+1)/2)
dorder = dble(order)
PI = 4.d0*datan(1.d0)
r = dtan(PI*f*dt)
r2 = r*r

allocate(c(1:3,0:n-1))
allocate(d(1:3,0:n-1))

a = 0.d0
b = 0.d0
c = 0.d0
d = 0.d0
do k = 0, m-1, 1
   p = -dsin(0.5d0*PI*(2*k+1)/dorder)
   rp = r*p
   den = r2-2.d0*rp+1.d0
   c(1,k) = (r2+2.0*rp+1.d0)/den
   c(2,k) = 2.d0*(r2-1.d0)/den
   c(3,k) = 1.d0
   d(1,k) = 1.d0/den
   d(2,k) = -2.d0/den
   d(3,k) = 1.d0/den
end do
do k = m, n-1, 1
   den = r+1.d0
   c(1,k) = (r-1.d0)/den
   c(2,k) = 1.d0
   d(1,k) = -1.d0/den
   d(2,k) = 1.d0/den
end do

call polyval(order, n, m, c, d, a, b)

deallocate(c, d)

end subroutine butthp
!
! buttbp
!
subroutine buttbp(order, dt, f1, f2, a, b)

implicit none

integer(4), intent(in) :: order

real(8), intent(in) :: dt, f1, f2
real(8), intent(inout) :: a(0:2*order)
real(8), intent(inout) :: b(0:2*order)


integer(4) k, n, m

real(8) u2, v2, vp
real(8) fi1, fi2, PI
real(8) dorder, den
real(8) u, v, p, p2, q2

real(8), dimension(:,:), allocatable :: c, d


if (abs(f1-f2) < 1.d-30) then
   write(*,"(A)") "f1 must be not equal to f2 in the butterworth bandpass filter !"
   stop
end if

if (min(abs(f1), abs(f2)) < 1.d-30) then
   write(*,"(A)") "Error: f1 must be great than 0 in the butterworth bandpass filter !"
   stop
end if

if (max(abs(f1), abs(f2)) > 0.5d0/dt) then
   write(*,"(A)") "Error: f2 must be less than the Nyquist frequency in the butterworth bandpass filter !"
   stop
end if


m = int(order/2)
n = int((order+1)/2)
dorder = dble(order)
PI = 4.d0*datan(1.d0)
fi1 = dtan(PI*f1*dt)
fi2 = dtan(PI*f2*dt)
p = 1.d0/abs(fi2-fi1)
p2 = p*p
q2 = p2*fi1*fi2
den = p2+q2
u = (q2-p2)/den
v = p/den
u2 = u*u
v2 = v*v

allocate(c(1:5,0:n-1))
allocate(d(1:5,0:n-1))

a = 0.d0
b = 0.d0
c = 0.d0
d = 0.d0
do k=0,m-1,1
   p = -dsin(0.5d0*PI*(2*k+1)/dorder)
   vp = v*p
   den = v2-2.d0*vp+1.d0
   c(1,k) = (v2+2.d0*vp+1.d0)/den
   c(2,k) = 4.d0*u*(1.d0+vp)/den
   c(3,k) = 2.d0*(2.d0*u2-v2+1.d0)/den
   c(4,k) = 4.d0*u*(1.d0-vp)/den
   c(5,k) = 1.d0
   d(1,k) = v2/den
   d(2,k) = 0.d0
   d(3,k) = -2.d0*v2/den
   d(4,k) = 0.d0
   d(5,k) = v2/den
end do
do k = m, n-1, 1
   den = 1.d0 + v
   c(1,k) = (1.d0-v)/den
   c(2,k) = 2.d0*u/den
   c(3,k) = 1.d0
   d(1,k) = -v/den
   d(2,k) = 0.d0
   d(3,k) = v/den
end do

call polyval2(order, n, m, c, d, a, b)

deallocate(c, d)

end subroutine buttbp
!
! buttbp
!
subroutine buttbs(order,dt,f1,f2,a,b)

implicit none

integer(4), intent(in) :: order

real(8), intent(in) :: dt, f1, f2
real(8), intent(inout) :: a(0:2*order)
real(8), intent(inout) :: b(0:2*order)


integer(4) k, n, m

real(8) u2, v2, vp
real(8) fi1, fi2, PI
real(8) dorder, den
real(8) u, v, p, p2, q2

real(8), dimension(:,:), allocatable :: c, d


if (abs(f1-f2) < 1.d-30) then
   write(*,"(A)") "f1 must be not equal to f2 in the butterworth bandstop filter !"
   stop
end if

if (min(abs(f1),abs(f2)) < 1.d-30) then
   write(*,"(A)") "Error: f1 must be great than 0 in the butterworth bandstop filter !"
   stop
end if

if (max(abs(f1),abs(f2)) > 0.5d0/dt) then
   write(*,"(A)") "Error: f2 must be less than the Nyquist frequency in the butterworth bandstop filter !"
   stop
end if


m = int(order/2)
n = int((order+1)/2)
dorder = dble(order)
PI = 4.d0*datan(1.d0)
fi1 = dtan(PI*f1*dt)
fi2 = dtan(PI*f2*dt)
p = 1.d0/abs(fi2-fi1)
p2 = p*p
q2 = p2*(fi1*fi2)
den = (p2+q2)
u = (q2-p2)/den
v = p/den
u2 = u*u
v2 = v*v

allocate(c(1:5,0:n-1))
allocate(d(1:5,0:n-1))

a = 0.d0
b = 0.d0
c = 0.d0
d = 0.d0
do k = 0, m-1, 1
   p = -dsin(0.5d0*PI*(2*k+1)/dorder)
   vp = v*p
   den = v2-2.d0*vp+1.d0
   c(1,k) = (v2+2.d0*vp+1.d0)/den
   c(2,k) = 4.d0*u*(1.d0+vp)/den
   c(3,k) = 2.d0*(2.d0*u2-v2+1.d0)/den
   c(4,k) = 4.d0*u*(1.d0-vp)/den
   c(5,k) = 1.d0
   d(1,k) = 1.d0/den
   d(2,k) = 4.d0*u/den
   d(3,k) = 2.d0*(2.d0*u2+1.d0)/den
   d(4,k) = 4.d0*u/den
   d(5,k) = 1.d0/den
end do
do k = m, n-1, 1
   den = 1.d0+v
   c(1,k) = (1.d0-v)/den
   c(2,k) = 2.d0*u/den
   c(3,k) = 1.d0
   d(1,k) = 1.d0/den
   d(2,k) = 2.d0*u/den
   d(3,k) = 1.d0/den
end do

call polyval2(order, n, m, c, d, a, b)

deallocate(c, d)

end subroutine buttbs
!
! polyval
!
subroutine polyval(order, n, m, c, d, a, b)

implicit none

integer(4), intent(in) :: n, m, order

real(8), intent(in) :: c(1:3,0:n-1)
real(8), intent(in) :: d(1:3,0:n-1)

real(8), intent(out) :: a(0:order)
real(8), intent(out) :: b(0:order)

integer(4) i, j, k
integer(4) nterm

real(8) ac(0:order)
real(8) bc(0:order)


a = 0.d0
b = 0.d0
if (1 == order) then
   ac(0:1) = c(1:2,0)
   bc(0:1) = d(1:2,0)
else
   a(0:2) = c(1:3,0)
   b(0:2) = d(1:3,0)

   nterm = 2
   do i = 1, m-1, 1
      ac = 0.d0
      bc = 0.d0
      do j = 0, nterm, 1
         do k = 1, 3, 1
            ac(j+k-1) = ac(j+k-1) + a(j)*c(k,i)
            bc(j+k-1) = bc(j+k-1) + b(j)*d(k,i)
         end do
      end do
      a(:) = ac(:)
      b(:) = bc(:)
      nterm = nterm + 2
   end do
   do i = m, n-1, 1
      ac = 0.d0
      bc = 0.d0
      do j = 0, nterm, 1
         do k = 1, 2, 1
            ac(j+k-1) = ac(j+k-1) + a(j)*c(k,i)
            bc(j+k-1) = bc(j+k-1) + b(j)*d(k,i)
         end do
      end do
      a(:) = ac(:)
      b(:) = bc(:)
      nterm = nterm + 1
   end do
   ac(:) = a(:)
   bc(:) = b(:)
end if

do j = 0, order, 1
   a(j) = ac(order-j)
   b(j) = bc(order-j)
end do

end subroutine polyval
!
! polyval
!
subroutine polyval2(order, n, m, c, d, a, b)

implicit none

integer(4), intent(in) :: n, m, order

real(8), intent(in) :: c(1:5,0:n-1)
real(8), intent(in) :: d(1:5,0:n-1)

real(8), intent(out) :: a(0:2*order)
real(8), intent(out) :: b(0:2*order)

integer(4) i, j, k
integer(4) nterm, order2


real(8) ac(0:2*order)
real(8) bc(0:2*order)


a = 0.d0
b = 0.d0
order2 = 2*order
if (1 == order) then
   ac(0:2) = c(1:3,0)
   bc(0:2) = d(1:3,0)
else
   a(0:4) = c(1:5,0)
   b(0:4) = d(1:5,0)

   nterm = 4
   do i = 1, m-1, 1
      ac = 0.d0
      bc = 0.d0
      do j = 0, nterm, 1
         do k = 1, 5, 1
            ac(j+k-1) = ac(j+k-1) + a(j)*c(k,i)
            bc(j+k-1) = bc(j+k-1) + b(j)*d(k,i)
         end do
      end do
      a(:) = ac(:)
      b(:) = bc(:)
      nterm = nterm + 4
   end do
   do i = m, n-1, 1
      ac = 0.d0
      bc = 0.d0
      do j = 0, nterm, 1
         do k = 1, 3, 1
            ac(j+k-1) = ac(j+k-1) + a(j)*c(k,i)
            bc(j+k-1) = bc(j+k-1) + b(j)*d(k,i)
         end do
      end do
      a(:) = ac(:)
      b(:) = bc(:)
      nterm = nterm + 2
   end do
   ac(:) = a(:)
   bc(:) = b(:)
end if

do j = 0, order2, 1
   a(j) = ac(order2-j)
   b(j) = bc(order2-j)
end do

end subroutine polyval2
!
! apply filter using the coefficients of butterworth filter
!
subroutine buttfilt(n, nt, a, b, x)

implicit none

integer, intent(in) :: n, nt

real(8), intent(in) :: a(0:n)
real(8), intent(in) :: b(0:n)

real(8), intent(inout) :: x(0:nt)


integer(4) j, k

real(8) mysum

real(8) z(-n:nt)
real(8) y(-n:nt)


z = 0.d0
y = 0.d0
z(0:nt) = x(0:nt)
do j = 0, nt, 1
   mysum = b(0) * z(j)
   do k = 1, n, 1
      mysum = mysum + (b(k)*z(j-k) - a(k)*y(j-k))
   end do
   y(j) = mysum
end do
x(0:nt) = y(0:nt)


end subroutine buttfilt


subroutine filtfilt(n, nt, a, b, x)

!Zero-phase forward and reverse digitial IIR filtering

implicit none

integer(4), intent(in) :: n, nt

real(8), intent(in) :: a(0:n)
real(8), intent(in) :: b(0:n)

real(8), intent(inout) :: x(0:nt)


integer(4) j, k

real(8) mysum

real(8) z(-n:nt)
real(8) y(-n:nt)


!forward
!call buttfilt(n, nt, a, b, x)
z = 0.d0
y = 0.d0
z(0:nt) = x(0:nt)
do j = 0, nt, 1
   mysum = b(0) * z(j)
   do k = 1, n, 1
      mysum = mysum + (b(k)*z(j-k) - a(k)*y(j-k))
   end do
   y(j) = mysum / a(0)
end do
z(0:nt) = y(nt:0:-1)


!backward
y = 0.d0
do j = 0, nt, 1
   mysum = b(0) * z(j)
   do k = 1, n, 1
      mysum = mysum + (b(k)*z(j-k) - a(k)*y(j-k))
   end do
   y(j) = mysum / a(0)
end do
z(0:nt) = y(nt:0:-1)
x(0:nt) = z(0:nt)


end subroutine filtfilt

!
!subroutine filtfilt(n, nt, a, b, x)
!
!!Zero-phase forward and reverse digitial IIR filtering
!
!implicit none
!
!integer, intent(in) :: n, nt
!
!real(8), intent(in) :: a(0:n)
!real(8), intent(in) :: b(0:n)
!
!real(8), intent(inout) :: x(0:nt)
!
!real(8) z(0:nt)
!
!
!!forward
!call buttfilt(n, nt, a, b, x)
!!backward
!z = 0.d0
!z(0:nt) = x(nt:0:-1)
!call buttfilt(n, nt, a, b, z)
!x(0:nt) = z(nt:0:-1)
!
!
!end subroutine filtfilt
