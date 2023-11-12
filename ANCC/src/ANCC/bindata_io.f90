module bindata_io_m

use db_m
use sac_io_m


implicit none



contains



! ====================================================================
! write complex data into a local binary file
! filename: file name of the bianary file [input]
! n: number of points of the data [input]
! comdata: complex data array [input]
! ier: status indicator [output]
! ====================================================================
subroutine write_bindata(filename, n, comdata, ier)

implicit none

integer, intent(in) :: n

character(len=*), intent(in) :: filename

complex(SGL), dimension(:), intent(in) :: comdata


integer, intent(out) :: ier


integer iunit



iunit = 11 + nprocs + myrank

open(unit=iunit, file=filename, status='replace', action='write', &
                   iostat=ier, access='stream', form='unformatted')

   if (0 /= ier) then
      write(*,"(A,A)") "Error: Cannot open: ", trim(adjustl(filename))
      call flush(6)
      return
   end if

   write(iunit) n
   write(iunit) comdata(1:n)

close(unit=iunit)


end subroutine write_bindata


! ====================================================================
! read complex data from a local binary file
! filename: file name of the bianary file [input]
! n: number of points to be read [output]
! comdata: complex data array [output]
! ier: status indicator [output]
! ====================================================================
subroutine read_bindata(filename, n, comdata, ier)

implicit none

character(len=*), intent(in) :: filename


integer, intent(out) :: n, ier

complex(SGL), allocatable, dimension(:), intent(out) :: comdata


integer iunit



iunit = 11 + myrank


open(unit=iunit, file=filename, status='old', action='read', &
              iostat=ier, access='stream', form='unformatted')

   if (0 /= ier) then
      !write(*,*) trim(adjustl(filename))
      write(*,"(A,A)") "Error: Cannot read: ", trim(adjustl(filename))
      call flush(6)
      return
   end if

   read(iunit) n

   allocate(comdata(n), stat=ier)

   read(iunit) comdata

close(unit=iunit)


end subroutine read_bindata


! ====================================================================
! Convert cross-correlation from sac format to ascii format
! in order to generate input files for Huajian Yao' code
! filename: file name of the cross-correlation file in sac format [input]
! ====================================================================
subroutine sac2asc(filename)

implicit none

character(len=*), intent(in) :: filename


integer iunit
integer k, n, n0, ier

logical is_existing

real dt, stlo, stla, evlo, evla

type(sachead) shd

character(len=512) infile, outfile

real, allocatable, dimension(:) :: sacdata



iunit = 11 + nprocs + myrank



inquire(file=filename, exist=is_existing)
if (.not.(is_existing)) return


call sacio_readsac(filename, shd, sacdata, ier)


n = shd%npts
dt = shd%delta
stlo = shd%stlo
stla = shd%stla
evlo = shd%evlo
evla = shd%evla

n0 = n/2 + 1


infile = trim(adjustl(filename))
outfile = trim(adjustl(infile(1:len_trim(infile)-4)))//'.dat'


open(unit=iunit, file=trim(adjustl(outfile)), action='write', &
                                  status='replace', iostat=ier)

   if (0 /= ier) then
      write(*,"(A,A)") "Error: Cannot open: ", trim(adjustl(outfile))
      call flush(6)
      return
   end if

   write(iunit,"(F15.5,5x,F15.5,11x,F15.5)") evlo, evla, dt
   write(iunit,"(F15.5,5x,F15.5,11x,F15.5)") stlo, stla, dt

   do k = 1, n0, 1
      write(iunit,"(F15.5,4x,2ES25.10)") (k-1)*dt, sacdata(n0+k-1), sacdata(n0-k+1)
   end do

close(unit=iunit)


end subroutine sac2asc



end module bindata_io_m
