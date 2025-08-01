#!/usr/bin/env python3
# -*- coding: UTF-8 -*-


'''
cutdata

This program cuts daily SAC files into segments using the ObsPy

Date: 25/10/2020
Author: Youshan Liu
Affiliation: Institute of Geology and Geophysics, Chinese Academy of Sciences


input data directory structure:
./top folder/stage folder/station folder/day folder/SAC files


After cutting, the files will be save in the following structures:
./data folder/yyyy/mm/yyyymmdd_hh00000/

for example:
./DATA/2007/10/20071003_000000/NCISP.NE00.BHZ.SAC


'''




import os
import re
import glob
from math import *
import numpy as np
from obspy import read
from obspy import UTCDateTime
from obspy.io.sac import SACTrace
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from multiprocessing.pool import ThreadPool



##############################################################
##############################################################
##############################################################
# some options to be set


# dryrun just for debug
# dryrun = True  => just print some directory information not to write sac files
# dryrun = False => do the actual files conversion
dryrun = False

# directory of the data
input_path = './DATA_Raw'
output_path = './DATA_cut'

# whether save input data
save_input_data = True


# segment length
segment_length = 3600*6


# component list to be cut
# three-components
#channels = ['HZ', 'HN', 'HE']
# only Z-component
channels = ['HZ']
# only N-component
#channels = ['HN']
# only E-component
#channels = ['HE']


# interpolation method: 'lanczos' is highest quality interpolation but expensive, the 'cubic' method may be a good choice
interpolation_method = 'cubic'
# You can have the following options
# 
# "lanczos": This offers the highest quality interpolation and should be chosen whenever possible. 
#            It is only due to legacy reasons that this is not the default method. The one downside it has is that it can be fairly expensive.
#            Essentially a fnite support version of sinc resampling (the ideal reconstruction filter).
#            For large values of a, it converges towards sinc resampling.
# if interpolation_method is lanczos, a value of 'a' that is the width of window in samples on either side shold be given.
# Values of a >= 20 show good results even for data that has energy close to the Nyquist frequency.
# Please see the https://docs.obspy.org/packages/autogen/obspy.signal.interpolation.lanczos_interpolation.html#obspy.signal.interpolation.lanczos_interpolation for details.
lanczos_winlen = 20 # this parameter is only valid for the 'lanczos' interpolation method
#
# "weighted_average_slopes": This is the interpolation method used by SAC. Refer to weighted_average_slopes() for more details.
#
# "slinear", "quadratic" and "cubic": spline interpolation of first, second or third order.
#
# "linear": Linear interpolation.
#
# "nearest": Nearest neighbour interpolation.


# nskip: a parameter is used to skip head and tail samples to check whether a segment data is zeros
if interpolation_method == 'lanczos':
	nskip = lanczos_winlen
else:
	nskip = 0



day_in_seconds = 24*3600

##############################################################



def create_sac_filename(stats):

	network = stats.sac.knetwk
	station = stats.sac.kstnm
	channel = stats.sac.kcmpnm

	filename = network + '.' + station + \
		'.' + channel + sac_suffix

	return filename



def create_sac_path(starttime):

	year  = '%4.4d' % starttime.year
	month = '%2.2d' % starttime.month
	day   = '%2.2d' % starttime.day
	hh    = '%2.2d' % starttime.hour
	mm    = '%2.2d' % starttime.minute
	ss    = '%2.2d' % starttime.second

	mydate = year + month + day

	month_folder = mydate + '_' + hh + mm + ss

	sac_path = year + '/' + month + \
				'/' + month_folder + '/'

	return sac_path



def date2str(starttime):

	year  = '%4.4d' % starttime.year
	month = '%2.2d' % starttime.month
	day   = '%2.2d' % starttime.day

	mydate = year + month + day

	return mydate



def merge_data(hour_files_list):

	st = Stream()
	data_len = 0.0
	for hour_file in hour_files_list:
		try:
			# faster for a specific file format
			st_tmp = read(hour_file, format='SAC', check_compression=False)
			#st_tmp = read(hour_file)
			if (not save_input_data):
				os.remove(hour_file)
		except:
			continue
		st += st_tmp
		data_len += (st_tmp[0].stats.npts-1)*st_tmp[0].stats.delta


	if (data_len < 0.50*segment_length):
		return []


	st.sort(['starttime'])

	# merge data
	if (len(st) > 1):
		st.merge(method=1, fill_value=0.0)
		#st.merge(method=1, fill_value='interpolate')


	#for i in range(0,len(st[0])):
	#	arr = st[0].data[i]
	#	if (np.ma.is_masked(arr)):
	#		st[0].data[i] = 0.0

	#if (np.ma.is_masked(st[0].data)):
	#	continue


	try:
		tr = st[0]
		# remove round error
		tr.stats.delta = round(tr.stats.delta*1e6)*1.e-6
	except:
		return


	stats = tr.stats
	dt = stats.delta
	df = stats.sampling_rate
	npts_daily = int(day_in_seconds*df)


	starttime = stats.starttime
	endtime = stats.endtime
	midtime = starttime + 0.50*stats.npts*dt
	starttime_daily = UTCDateTime(midtime.year, midtime.month, midtime.day, 0, 0, 0, 0)
	endtime_daily = starttime_daily + day_in_seconds - dt


	# indexes of the starttime and endtime in local temporal axis
	ibeg = ceil((max(starttime, starttime_daily) - starttime_daily)*df)
	iend = int((min(endtime, endtime_daily) - starttime_daily)*df)


	if (0 != starttime.microsecond % round(dt*1e6)):
		# interpolate to local temporal axis
		#tr.interpolate(df, 'cubic', starttime_daily+ibeg*dt, iend-ibeg+1, 0.0)
		#tr.interpolate(df, 'lanczos', starttime_daily+ibeg*dt, iend-ibeg+1, 0.0, a=21)
		tr.interpolate(df, interpolation_method, starttime_daily+ibeg*dt, iend-ibeg+1, 0.0, a=lanczos_winlen)
	else:
		tr = tr.slice(starttime_daily+ibeg*dt, starttime_daily+iend*dt)


	# pad array
	data = np.zeros(npts_daily)
	data[ibeg:iend+1] = tr.data


	# replace data
	tr.data = data
	tr.stats.starttime = starttime_daily


	del data, st, st_tmp

	return tr



def cutdata_daily(day_folder):

	day_path = station_path + day_folder + '/'
	print('Entering directory ' + day_path[len_topdir:-1])
	#print('\n')

	if (not os.path.isdir(day_path)):
		return

	for C in channels:

		hour_files_list = glob.glob(day_path + '*' + C + '*' + sac_suffix)

		if ([] == hour_files_list):
			continue


		try:
			tr = merge_data(hour_files_list)
			if ([] == tr):
				continue
			# remove round error
			tr.stats.delta = round(tr.stats.delta*1e6)*1.e-6
			dt = tr.stats.delta
			df = tr.stats.sampling_rate
			starttime_daily = tr.stats.starttime
		except:
			print('merge data failure')
			del hour_files_list
			continue


		# remove invalid data
		if ((max(tr.data[nwin:-nwin-1]) - min(tr.data[nskip:-1-nskip])) < 1.e-12):
			del tr, hour_files_list
			continue


		# cut data into segments
		nsegments = round(day_in_seconds / segment_length)
		#nsegments = int(day_in_seconds / segment_length)


		ibeg = 0
		tbeg = 0.0
		for i in range(nsegments):

			tend = min(tbeg + segment_length - dt, day_in_seconds)
			iend = int(tend * df)


			tr_out = tr.copy()
			tr_out.data = tr.data[ibeg:iend+1]
			starttime = starttime_daily + tbeg
			tr_out.stats.starttime = starttime


			# remove invalid data
			if ((max(tr_out.data[nwin:-nwin-1]) - min(tr_out.data[nskip:-1-nskip])) < 1.e-12):
				del tr_out
				ibeg = iend + 1
				tbeg = tend + dt
				if ((day_in_seconds - tend) < 0.4*segment_length):
					break
				else:
					continue


			sac_path = create_sac_path(starttime)
			sac_filename = create_sac_filename(tr_out.stats)

			mypath = output_path + '/' + sac_path

			# create the path if it does not exist
			if (not os.path.exists(mypath)):
				os.makedirs(mypath, exist_ok=True)


			outfile = mypath + sac_filename


			# write sac
			if (dryrun):

				if (os.path.exists(outfile) and os.path.isfile(outfile)):
					os.remove(outfile)

			else:

				sac = SACTrace.from_obspy_trace(tr_out)


				sac.nzyear = starttime.year
				sac.nzjday = starttime.julday
				sac.nzhour = starttime.hour
				sac.nzmin = starttime.minute
				#sac.nzsec = starttime.second
				##sac.nzmsec = int(starttime.microsecond*1.e-3)
				#sac.nzmsec = round(starttime.microsecond*1.e-3)
				sec = round((starttime.second + starttime.microsecond*1.e-6)*1000)
				sac.nzsec = int(sec*0.001)
				sac.nzmsec = int(sec - sac.nzsec*1000)
				sac.b = 0
				#sac.reftime += sac.b
				#sac.reftime = starttime


				sac.write(outfile)

				del sac


			del tr_out
			ibeg = iend + 1
			tbeg = tend + dt


			if ((day_in_seconds - tend) < 0.4*segment_length):
				break


		del tr, hour_files_list
		print('%s is done' % date2str(starttime_daily))

	print('Leaving directory ' + day_path[len_topdir:-1])
	#print('\n')

	return



def cutdata():

	global len_topdir, station_path, sac_suffix

	len_topdir = len(input_path) + 1

	stage_folders_list = os.listdir(input_path)


	sac_suffix = '.SAC'

	# convert reftek to sac
	for stage_folder in stage_folders_list:

		stage_path = input_path + '/' + stage_folder + '/'
		print('Entering directory ' + stage_path[len_topdir:-1])
		#print('\n')

		if (not os.path.isdir(stage_path)):
			continue

		station_folder_list = os.listdir(stage_path)

		for station_folder in station_folder_list:

			station_path = stage_path + station_folder + '/'
			print('Entering directory ' + station_path[len_topdir:-1])
			#print('\n')

			if (not os.path.isdir(station_path)):
				continue

			day_folders_list = os.listdir(station_path)

			pool = ThreadPool()
			pool.map(cutdata_daily, day_folders_list)
			pool.close()
			pool.join()

			del day_folders_list
			print('Leaving directory ' + station_path[len_topdir:-1])
			#print('\n')

		del station_folder_list
		print('Leaving directory ' + stage_path[len_topdir:-1])
		#print('\n')

	del stage_folders_list

	return



if __name__ == '__main__':

	#print('\n')
	print('cutdata: ')
	print('This program cuts daily data into segments using the ObsPy (parallel version)')
	print('Youshan Liu at Institute of Geology and Geophysics, Chinese Academy of Sciences')
	print('Welcome to send any bugs and suggestions to ysliu@mail.iggcas.ac.cn')
	print('\n')

	starttime = UTCDateTime()

	# absolution path
	#current_path = os.getcwd()

	# cutdata daily
	cutdata()

	endtime = UTCDateTime()

	elapsed_time = (endtime - starttime)

	print('\n')
	print('Start   time : %s' % starttime)
	print('End     time : %s' % endtime)
	print('Elapsed time : %f hours' % (elapsed_time / 3600.0))



