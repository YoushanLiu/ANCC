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
!* The sample of test driver for FTAN with phase match filter for
!* subroutines aftanpg and aftanipg
!
program AFTAN

use dispio_m
use aftanpg_m
use aftanipg_m

implicit none

integer(4) i, k, sac, nargc, iargc
integer(4) n, npoints, nfin, nfout1, nfout2
integer(4) nrow, ncol, npred, nprpv, ier, ioer

real(8) t0, dt, tresh, ffact1, ffact2
real(8) perc, taperl, fmatch, delta, tamp
real(8) vmin, vmax, tmin, tmax, snr, lambda
real(8) tmin_read, tmax_read

real(8) PIover4

real(8) x(1), y(1)

real(8), dimension(:), allocatable :: seis

real(8), dimension(:), allocatable :: prpvper, prpvvel
real(8), dimension(:,:), allocatable :: ampo, arr1, arr2, pred

character(1) isPerCut, isOutput, isVerbose

character(256) filename, parfile, outfile



sac = 1

!---  input command line arguments treatment
nargc = iargc()
if (1 /= nargc) then
   write(*,*) "Usage: AFTAN parameter_file \n"
   stop
endif


! read input.dat
open(unit=99, file='input.dat', status='old')

   do k = 1, 25, 1
      read(99,"(A)")
   end do
   read(99,*) PIover4
   read(99,*) vmin, vmax
   read(99,*) tmin_read, tmax_read
   read(99,*) isPerCut, lambda
   read(99,*) isOutput
   if (tmin_read < 1) then
      isOutput = 'Y'
   end if
   read(99,*) tresh
   read(99,*) ffact1, ffact2
   read(99,*) taperl
   read(99,*) snr
   read(99,*) fmatch
   do k = 1, 6, 1
      read(99,"(A)")
   end do
   read(99,*) isVerbose

close(unit=99)



! read ref1Dmod.dat
open(11, file='ref1Dmod.dat', status='old', iostat=ier)

   if (0 /= ier) then
      write(*,*)'Error: No phase velocity reference file!'
      stop
   end if

   nprpv = 0
   do
      read(11, *, iostat=ioer) x(1), y(1)
      if (0 /= ioer) exit
      nprpv = nprpv + 1
   end do

   rewind(11)

   allocate(prpvper(1:nprpv), prpvvel(1:nprpv), stat=ier)
   do k = 1, nprpv, 1
      read(11, *) prpvper(k), prpvvel(k)
   end do

close(11)



call getarg(1, parfile)
open(10, file=trim(adjustl(parfile)), status='old')


   do

      read(10, '(A)', iostat=ioer) filename
      if (0 /= ioer) exit


      filename = trim(adjustl(filename))
      outfile = trim(adjustl(filename(1:len_trim(filename)-4)))


      tmin = tmin_read
      tmax = tmax_read
      if ((isVerbose == 'Y') .or. (isVerbose == 'y')) then
         write(*,'(A, A)') 'AFTAN: ', trim(adjustl(filename))
      end if

      !
      ! read SAC or ascii data
      !

      call readhead(sac, filename, n, ier)
      allocate(seis(1:n), stat=ier)
      call readdata(sac, filename, dt, delta, t0, seis, ier)

      if ((isPerCut == 'Y') .or. (isPerCut == 'y')) then

         do k = 1, 100, 1

            x(1) = k
            call linear_interpo(prpvper(1:nprpv), prpvvel(1:nprpv), x, y, ier)

            if ((lambda*x(1)*y(1) > delta) .and. ((k+3) < tmax)) then
               tmax = k + 3
               exit
            end if

         end do

      end if


      nfin    = 64
      npoints = 5
      perc    = 50.0


      !---  FTAN with phase match filter. First Iteration.
      call aftanpg(PIover4, n, seis, t0, dt, delta, vmin, vmax, tmin, tmax, &
                tresh, ffact1, perc, npoints, taperl, nfin, nprpv, prpvper, &
                prpvvel, nfout1, arr1, nfout2, arr2, tamp, nrow, ncol, ampo, ier)

      if ((isOutput == 'Y') .or. (isOutput == 'y')) then
         call printres(dt, delta, nfout1, arr1, nfout2, arr2, &
                    tamp, nrow, ncol, ampo, ier, outfile, "_1")
      end if

      if (tmin >= 1) then
         call write_data(arr1, nfout1, arr2, nfout2, trim(adjustl(outfile))//'_1')
      endif
      if (allocated(ampo)) then
	 deallocate(ampo)
      end if
      if (allocated(arr1)) then
	 deallocate(arr1)
      end if


      if (0 == nfout2) then
         cycle
      end if


      !--- Make prediction based on the first iteration
      npred = nfout2
      tmin = arr2(2,1)
      tmax = arr2(2,nfout2)
      allocate(pred(1:nfout2,1:2))
      do i = 1, nfout2, 1
         pred(i,1) = arr2(2,i)
         pred(i,2) = arr2(3,i)
      enddo
      if (allocated(arr2)) then
	 deallocate(arr2)
      end if



      !--- FTAN with phase matching filter. Second Iteration.
      call aftanipg(PIover4, n, seis, t0, dt, delta, vmin, vmax, tmin, tmax, tresh, &
              ffact2, perc, npoints, taperl, nfin, snr, fmatch, npred, pred, nprpv, &
              prpvper, prpvvel, nfout1, arr1, nfout2, arr2, tamp, nrow, ncol, ampo, ier)

      if ((isOutput == 'Y') .or. (isOutput == 'y')) then
         call printres(dt, delta, nfout1, arr1, nfout2, arr2, &
                    tamp, nrow, ncol, ampo, ier, outfile, "_2")
      end if


      if (tmin >= 1) then
         call write_data(arr1, nfout1, arr2, nfout2, trim(adjustl(outfile))//'_2')
      end if


      deallocate(seis)
      if (allocated(arr1)) then
         deallocate(arr1)
      end if
      if (allocated(arr2)) then
         deallocate(arr2)
      end if
      if (allocated(ampo)) then
         deallocate(ampo)
      end if
      if (allocated(pred)) then
         deallocate(pred)
      end if

   end do


close(10)


deallocate(prpvper, prpvvel)


end program AFTAN
