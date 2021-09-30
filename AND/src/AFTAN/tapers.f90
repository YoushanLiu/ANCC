!--------------------------------------------------------------
! spectra tapering procedure
!--------------------------------------------------------------
subroutine tapers(omb, ome, dom, alpha, ns, omstart, inds, inde, omdom, ampdom)

implicit none

integer(4), intent(in) :: ns

real(4), intent(in) :: omb, ome, dom, alpha


integer(4), intent(out) :: inds, inde

real(4), intent(out) :: omstart

real(4), intent(out) :: ampdom(ns), omdom(ns)


integer(4) i, om1, om2, om3, om4

real(8), parameter :: PI = 4.d0*datan(1.d0)

real(4) om2d, om3d, wd, tresh



tresh = 0.50
ampdom(1:ns) = 0.0

om2d = omb/dom
wd = max(16.0, om2d*sqrt(tresh/alpha))
om1 = nint(max(1.0, om2d-0.50*wd))
om2 = nint(min(real(ns), om1+wd))
do i = om1, om2, 1
   ampdom(i) = 0.50*(1.0 - cos(PI/(om2-om1)*(i-om1)))
enddo


om3d = ome/dom
wd = max(16.0, om3d*sqrt(tresh/alpha))
om4 = nint(min(real(ns), om3d+0.50*wd))
om3 = nint(max(1.0, om4-wd))
do i = om3, om4, 1
   ampdom(i) = 0.50*(1.0 + cos(PI/(om4-om3)*(i-om3)))
enddo


ampdom(om2:om3) = 1.0


do i = 1, ns, 1
   omdom(i) = (i-1)*dom
enddo


omstart = omb

inds = om1
inde = om4


return


end subroutine tapers
