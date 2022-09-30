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
!--------------------------------------------------------------
! spectra tapering procedure
!--------------------------------------------------------------
subroutine tapers(omb, ome, dom, alpha, ns, omstart, inds, inde, omdom, ampdom)

implicit none

integer(4), intent(in) :: ns

real(8), intent(in) :: omb, ome, dom, alpha


integer(4), intent(out) :: inds, inde

real(8), intent(out) :: omstart

real(8), intent(out) :: ampdom(ns), omdom(ns)


integer(4) i, om1, om2, om3, om4

real(8), parameter :: PI = 4.d0*datan(1.d0)

real(8) om2d, om3d, wd, tresh



tresh = 0.5d0
ampdom(1:ns) = 0.d0

om2d = omb/dom
wd = max(16.d0, om2d*sqrt(tresh/alpha))
om1 = nint(max(1.d0, om2d-0.5d0*wd))
om2 = nint(min(real(ns), om1+wd))
do i = om1, om2, 1
   ampdom(i) = 0.5d0*(1.d0 - cos(PI*(i-om1)/(om2-om1)))
enddo


om3d = ome/dom
wd = max(16.d0, om3d*sqrt(tresh/alpha))
om4 = nint(min(real(ns), om3d+0.5d0*wd))
om3 = nint(max(1.d0, om4-wd))
do i = om3, om4, 1
   ampdom(i) = 0.5d0*(1.d0 + cos(PI*(i-om3)/(om4-om3)))
enddo


ampdom(om2:om3) = 1.d0


do i = 1, ns, 1
   omdom(i) = (i-1)*dom
enddo


omstart = omb

inds = om1
inde = om4


return


end subroutine tapers
