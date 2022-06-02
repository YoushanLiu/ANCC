#!/usr/bin/env python3
# -*- coding: UTF-8 -*-


'''
cutdata

This program cuts daily sac files into segments using the ObsPy

Date: 25/10/2020
Author: Youshan Liu
Affiliation: Institute of Geology and Geophysics, Chinese Academy of Sciences


After cutting, the files will be save in the following structures:
./input_folder/yyyy/mm/yyyymmdd_hh00000/

for example:
./DATA/2007/10/20071003_000000/NCISP.NE00.BHZ.SAC


'''




import os
import re
import glob
from math import *
import numpy as np
from obspy import read
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from obspy.core import UTCDateTime
from obspy.io.sac import SACTrace
from multiprocessing.dummy import Pool as ThreadPool



##############################################################
##############################################################
##############################################################
# some options to be set


# dryrun just for debug
# dryrun = True  => just print some direction information not to write sac files
# dryrun = False => do the actual files conversion
dryrun = False

# direction of the reftek data
input_folder = 'DATA_Raw'
output_folder = 'DATA_cut'

# whether save input data
save_input_data = True


# segment length
segment_length = 24*3600

# component list to be cut
# three-components
#component_list = ['Z', 'N', 'E']
# only Z-component
component_list = ['Z']
# only N-component
#component_list = ['N']
# only E-component
#component_list = ['E']


seconds_daily = 24*3600

##############################################################



def create_sac_filename(stats):

	network_name = stats.sac.knetwk
	station_name = stats.sac.kstnm
	channel_name = stats.sac.kcmpnm

	filename = network_name + '.' + station_name + \
		'.' + channel_name + sac_suffix

	return filename



def create_sac_path(stats):

	time  = UTCDateTime(stats.starttime)
	year  = '%4.4d' % time.year
	month = '%2.2d' % time.month
	day   = '%2.2d' % time.day
	hh    = '%2.2d' % time.hour
	mm    = '%2.2d' % time.minute
	ss    = '%2.2d' % time.second

	mydate = year + month + day

	month_folder = mydate + '_' + hh + mm + ss

	sac_path = year + '/' + month + \
				'/' + month_folder + '/'

	return sac_path



def date2str(starttime):

	time  = UTCDateTime(starttime)
	year  = '%4.4d' % time.year
	month = '%2.2d' % time.month
	day   = '%2.2d' % time.day

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
	except:
		return


	dt = tr.stats.delta
	df = tr.stats.sampling_rate
	npts_daily = int(seconds_daily*df)


	starttime = tr.stats.starttime
	endtime = tr.stats.endtime
	midtime = starttime + 0.50*tr.stats.npts*dt
	starttime_daily = UTCDateTime(midtime.year, midtime.month, midtime.day, 0, 0, 0, 0)
	endtime_daily = starttime_daily + seconds_daily - dt


	# indexes of the starttime and endtime in local temporal axis
	ibeg = ceil((max(starttime, starttime_daily) - starttime_daily)*df)
	iend = int((min(endtime, endtime_daily) - starttime_daily)*df)


	# interpolate to local temporal axis
	#tr.interpolate(df, 'cubic', starttime_daily+ibeg*dt, iend-ibeg+1, 0.0)
	tr.interpolate(df, 'lanczos', starttime_daily+ibeg*dt, iend-ibeg+1, 0.0, a=21)


	# pad array
	data = np.zeros(npts_daily)
	data[ibeg:iend+1] = tr.data


	# replace data
	tr.data = data
	tr.stats.starttime = starttime_daily


	del data, st, st_tmp

	return tr



def cutdata_daily(station_path):

	day_folders_list = os.listdir(station_path)

	for day_folder in day_folders_list:

		day_path = station_path + day_folder + '/'
		print('Entering directory ' + day_path[nrootdir:-1])
		print('\n')

		if (not os.path.isdir(day_path)):
			return

		for C in component_list:

			hour_files_list = glob.glob(day_path + '*' + C + sac_suffix)

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
			if ((max(tr.data) - min(tr.data)) < 1.e-12):
				del tr, hour_files_list
				continue


			# cut data into segments
			nsegments = round(seconds_daily / segment_length)
			#nsegments = int(seconds_daily / segment_length)


			ibeg = 0
			tbeg = 0.0
			for i in range(0,nsegments):

				tend = min(tbeg + segment_length - dt, seconds_daily)
				iend = int(tend * df)


				tr_out = tr.copy()
				tr_out.data = tr.data[ibeg:iend+1]
				tr_out.stats.starttime = starttime_daily + tbeg


				# remove invalid data
				if ((max(tr_out.data) - min(tr_out.data)) < 1.e-12):
					del tr_out
					ibeg = iend + 1
					tbeg = tend + dt
					continue


				sac_filename = create_sac_filename(tr_out.stats)
				sac_path = create_sac_path(tr_out.stats)

				mypath = output_path + sac_path

				# create the path if it does not exist
				if (not os.path.exists(mypath)):
					os.makedirs(mypath)


				outfile = mypath + sac_filename


				# write sac
				if (dryrun):
					if (os.path.exists(outfile) and os.path.isfile(outfile)):
						os.remove(outfile)
				else:
					sac = SACTrace.from_obspy_trace(tr_out)


					#sac.nzyear = tr_out.stats.starttime.year
					#sac.nzjday = tr_out.stats.starttime.julday
					#sac.nzhour = tr_out.stats.starttime.hour
					#sac.nzmin = tr_out.stats.starttime.minute
					#sac.nzsec = tr_out.stats.starttime.second
					#sac.nzmsec = int(tr_out.stats.starttime.microsecond*1.e-3)
					#sac.b = 0
					#sac.reftime += sac.b
					sac.reftime = tr_out.stats.starttime


					sac.write(outfile)

					del sac


				del tr_out
				ibeg = iend + 1
				tbeg = tend + dt


				if ((seconds_daily - tend) < 0.1*segment_length):
					break


			del tr, hour_files_list
			print('%s is done ... \n' % date2str(starttime_daily))

		print('Leaving directory ' + day_path[nrootdir:-1])
		print('\n')

	del day_folders_list

	return



def cutdata(current_path):

	global nrootdir, output_path, sac_suffix

	rootdir = current_path + '/' + input_folder + '/'
	nrootdir = len(rootdir)

	period_folders_list = os.listdir(rootdir)

	output_path = current_path + '/' + output_folder + '/'


	sac_suffix = '.SAC'

	# cut daily sac files into segments
	for period_folder in period_folders_list:

		period_path = rootdir + period_folder + '/'
		print('Entering directory ' + period_path[nrootdir:-1])
		print('\n')

		if (not os.path.exists(period_path)):
			continue

		station_folder_list = os.listdir(period_path)

		for station_folder in station_folder_list:

			station_path = period_path + station_folder + '/'
			print('Entering directory ' + station_path[nrootdir:-1])
			print('\n')

			if (not os.path.exists(station_path)):
				continue

			cutdata_daily(station_path)

			print('Leaving directory ' + station_path[nrootdir:-1])
			print('\n')

		del station_folder_list
		print('Leaving directory ' + period_path[nrootdir:-1])
		print('\n')

	del period_folders_list

	return



if __name__ == '__main__':

	print('\n')
	print('cutdata: ')
	print('This program cuts daily data into segments using the ObsPy (serial version)')
	print('Youshan Liu at Institute of Geology and Geophysics, Chinese Academy of Sciences')
	print('Welcome to send any bugs and suggestions to ysliu@mail.iggcas.ac.cn')
	print('\n\n')

	starttime = UTCDateTime()

	# get current path
	# absolution path
	#current_path = os.getcwd()
	# relative path
	current_path = '.'

	# cutdata daily
	cutdata(current_path)

	endtime = UTCDateTime()

	elapsed_time = (endtime - starttime)

	print('Start   time : %s' % starttime)
	print('End     time : %s' % endtime)
	print('Elapsed time : %f hours \n' % (elapsed_time / 3600.0))



