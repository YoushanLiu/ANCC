#!/usr/bin/env python3
# -*- coding: UTF-8 -*-


'''
reftek2sac

This program converts files from REFTEK to SAC format using the ObsPy

Date: 22/10/2020
Author: Youshan Liu
Affiliation: Institute of Geology and Geophysics, Chinese Academy of Sciences


directory structure:
./top folder/period folder/station folder/day folder/Reftek UnitID number/stream/reftek files

for example:
./DATA_Raw/2007_276_2008_005/NE00/2007276/9F78/1



Reftek 130 Disk Directory Structure

  Year
  |   Day of year
  |   |
  v   v
   __  _
  |  || |
 \2003032
 \2003033

		Unit ID number
		|
		v
		 __
   		|  |
  		\90F0

   		Datastream
   		|
   		v
  		\0
  		\1
  		\2

'''




import os
import re
from math import *
from obspy import read
from obspy import UTCDateTime
from obspy.io.sac import SACTrace
from multiprocessing.dummy import Pool as ThreadPool


##############################################################
##############################################################
##############################################################
# some options to be set


# dryrun just for debug
# dryrun = True  => just print some directory information not to write sac files
# dryrun = False => do the actual files conversion
dryrun = False

# directory of the reftek data
input_path = './DATA_Raw'

# station information file
station_list = './NEsta.lst'


# component list to be converted
# only Z-component
#components = ['0']
# only N-component
#components = ['1']
# only E-component
#components = ['2']
# three-components
components = ['0', '1', '2']

##############################################################
# preprocess options

# whether remove mean value or not
is_demean = True

# whethet remove linear and nonlinear trends or not
is_detrend = True
is_denonlinear = False

# whether apply bandpass filter or not
is_bandpass  = False
is_zerophase = True
flow  = 1.0/50.0
fhigh = 5.0

# whether downsampling seismograms
is_decimate = True
# the downsampling frequency
downsampling_rate = 10.0

##############################################################

# channel name, it consists of band code, instrument code, orientation code
channels = ['BHZ', 'BHN', 'BHE']


##############################################################
##############################################################
##############################################################
# read station list
def read_station_list(filename):

	class Station:
		def __init__(self):
			self.stnm = []
			self.stla = []
			self.stlo = []
			self.stel = []
			self.netwk = []


	try:
		with open(filename, 'r') as f:
			lines = f.readlines()
	except:
		raise Exception('Cannot open file %s !' % filename)
	sta = Station()
	netwk, stnm, stla, stlo, stel = None, None, None, None, None
	for line in lines[1:]:
		try:
			line_splited = line.split()
			if ([] == line_splited):
				continue
			netwk, stnm, stla, stlo, stel = line_splited
		except:
			raise Exception('Format error in %s !' % filename)
		sta.stnm.append(stnm)
		sta.netwk.append(netwk)
		sta.stla.append(float(stla))
		sta.stlo.append(float(stlo))
		sta.stel.append(float(stel))

	del lines

	return sta



def create_sac_filename(stats):

	time = UTCDateTime(stats.starttime)
	yyyy = '%4.4d' % time.year
	ddd  = '%3.3d' % time.julday
	hh   = '%2.2d' % time.hour
	mm   = '%2.2d' % time.minute
	ss   = '%2.2d' % time.second
	fff  = '%3.3d' % (0.001*time.microsecond)

	sac_filename = yyyy + '.' + ddd + '.' + hh + '.' + \
		   mm + '.' + ss + '.' + fff + '.' + \
		   stats.network + '.' + stats.station + '.' + \
		   stats.channel + sac_suffix

	return sac_filename



def findstr(str_src, str_target):

	'''
	findstr(str_src, str_target)
	find the index of a substring in a string
	str_src-> source string
	str_target -> substring
	'''

	n_src = len(str_src)
	n_target = len(str_target)
	i = 0
	ipos = []
	while str_target in str_src[i:]:
		idx = str_src.index(str_target, i, n_src)
		ipos.append(idx)
		i = (idx + n_target)

	return ipos



def convert_hourly(hour_files_path, day_path):

	hour_files_list = os.listdir(hour_files_path)

	idx = findstr(hour_files_path, '/')

	for hour_file in hour_files_list:

		#if ('EE' != hour_file[-4:-2]):
		#	continue

		# if this file is not a reftek format
		is_sacfile = False
		starts = [each.start() for each in re.finditer(sac_suffix, hour_file.upper())]
		ends = [start+len(sac_suffix) - 1 for start in starts]
		span = [(start, end) for start,end in zip(starts, ends)]
		is_sacfile = is_sacfile and (len(span) >= 1)

		if (is_sacfile):
			continue


		try:
			has_dot_in_filename = hour_file.index('.')
		except:
			has_dot_in_filename = False
		if (has_dot_in_filename):
			continue


		infile = hour_files_path + hour_file
		try:
			# faster for specific file format
			st = read(infile, format='REFTEK130', check_compression=False, component_codes=['0', '1', '2'])
			#st = read(infile, component_codes=['0', '1', '2'])
			#st = read(infile, component_codes=['Z', 'N', 'E'])
			# cancel out the sort in obspy
			#st.sort(reverse=True)
		except:
			continue


		if (len(st) > 3):
			st.sort(['starttime'])
			if ((st[nstream-1].stats.endtime - st[0].stats.starttime) > 129600):
				continue
			st.merge(method=1, fill_value=0)
			st.sort()



		for i in range(len(st)):

			try:
				tr = st[i]
				# remove round error
				tr.stats.delta = round(tr.stats.delta*1e6)*1.e-6
			except:
				continue

			if (tr.stats.channel[-1] not in components):
				continue
			channel = channels[i]



			#starttime = tr.stats.starttime
			#endtime = tr.stats.endtime
			#min2hour = int((starttime.hour + (starttime.second + starttime.microsecond*1.e-6)/60.0 + 40)/60.0)
			#if (0 == min2hour):
			#	dtinus = 1e6 / downsampling_rate
			#	#microsecond = ceil(microsecond / dtinus) * dtinus
			#	#starttime_first = starttime
			#	#starttime_first = UTCDateTime(starttime.year, starttime.month, starttime.day, starttime.hour, starttime.minute, starttime.second, microsecond, strict=False)
			#	sec = round(starttime.second*1e6 + starttime.microsecond)
			#	sec = ceil(sec / dtinus) * dtinus
			#	second = int(sec * 1.e-6)
			#	microsecond = int(sec - second*1e6)
			#	starttime_first = UTCDateTime(starttime.year, starttime.month, starttime.day, starttime.hour, starttime.minute, second, microsecond, strict=False)
			#else:
			#	starttime_first = UTCDateTime(starttime.year, starttime.month, starttime.day, starttime.hour + min2hour, 0, 0, 0)
			#	starttime_first = UTCDateTime(starttime.year, starttime.month, starttime.day, starttime.hour, 0, 0, 0) + 3600
			#if ((starttime.minute + (starttime.second + starttime.microsecond*1.e-6)/60.0) > 40.0):
			#	endtime_last = endtime
			#else:
			#	endtime_last = UTCDateTime(endtime.year, endtime.month, endtime.day, endtime.hour, 0, 0, 0)
			#if (starttime_first < endtime_last):
			#	tr = tr.slice(starttime=starttime_first, endtime=endtime_last, nearest_sample=False)
			#else:
			#	del tr
			#	continue




			# some preprocess
			if (is_demean):
				tr.detrend(type='demean')
			if (is_detrend):
				tr.detrend(type='linear')
			#if (is_denonlinear):
			#	#tr.detrend('polynomial', order=3)
			#	tr.detrend('spline', order=3, dspline=5)


			# bandpass
			if (is_bandpass):
				tr.filter('bandpass', freqmin=flow, freqmax=fhigh, corners=2, zerophase=is_zerophase)


			# downsampling
			if (is_decimate):
				df = tr.stats.sampling_rate
				if (downsampling_rate > df):
					print("Error: downsampling sampling rate cannot large than original sampling rate !")
					continue
				decimate_factor = int(df / downsampling_rate)
				if (abs(df - (decimate_factor*downsampling_rate)) > 0.0):
					print("Error: decimate factor can only be integer !")
					continue
				if (decimate_factor > 1):
					# Nyquist frequency of downsampling rate
					freq_lowpass = 0.49 * downsampling_rate
					#freq_lowpass = 0.49 * tr.stats.sampling_rate / decimate_factor
					if (not(is_bandpass and (fhigh <= freq_lowpass))):
						tr.filter('lowpass', freq=freq_lowpass, corners=2, zerophase=True)
					tr.decimate(factor=decimate_factor, strict_length=False, no_filter=True)


			# get the index of station name in station list
			#station_name = tr.stats.station
			#res = [sta.stnm.index(x) for x in sta.stnm if x.upper() == station_name.upper()]
			#if ([] != res):
			#	# first find station name in reftek head.
			#	# if the field "station" in reftek header is NULL, then try extract station name from folder name
			#	try:
			#		ipos = res[0]
			#	except:
			#		raise Exception('Error: station %s is not in the station list' % station_name)
			#else:
			#	#print("Warning: station field in reftek header is NULL")
			#	# try to extract station name based on folder name
			#	ipos = -1
			#	for j in range(len(sta.stnm)):
			#		res = findstr(hour_files_path[idx[-5]+1:idx[-4]], sta.stnm[j])
			#		if ([] != res):
			#			ipos = j
			#			break
			#	if (-1 == ipos):
			#		print("Error: station %s is not in the station list or field 'station' in reftek header is NULL" % station_name)
			#		return


			ipos = -1
			station_name = hour_files_path[idx[-5]+1:idx[-4]]
			for j in range(len(sta.stnm)):
				res = findstr(station_name, sta.stnm[j])
				#res = findstr(hour_files_path[len_topdir:-1], sta.stnm[j])
				if ([] != res):
					ipos = j
					break
			if (-1 == ipos):
				print("Error: station folder %s does not include the name of this station or this station is missing in the stainfo.lst" % station_name)
				return



			station = sta.stnm[ipos]
			network = sta.netwk[ipos]
			tr.stats.network = network
			tr.stats.station = station
			tr.stats.channel = channel
			starttime = tr.stats.starttime

			sac_filename = create_sac_filename(tr.stats)

			outfile = day_path + sac_filename


			if (dryrun):

				if (os.path.exists(outfile) and os.path.isfile(outfile)):
					os.remove(outfile)

			else:

				sac = SACTrace.from_obspy_trace(tr)


				# write stla, stlo, stel into sac header
				sac.stla = sta.stla[ipos]
				sac.stlo = sta.stlo[ipos]
				sac.stel = sta.stel[ipos]
				# write station name
				sac.kstnm = station
				# write network name
				sac.knetwk = network


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
				sac.reftime = tr.stats.starttime


				# write sac
				sac.write(outfile)
				#tr.write(outfile, format='sac')

				del sac

			del tr

		del st
		print('%s is done' % hour_file)

	del hour_files_list

	return



def convert_daily(day_folder):

	day_path = station_path + day_folder + '/'
	print('Entering directory ' + day_path[len_topdir:-1])
	#print('\n')

	if (not os.path.isdir(day_path)):
		return

	UnitID_folders_list = os.listdir(day_path)

	for UnitID in UnitID_folders_list:

		UnitID_path = day_path + UnitID + '/'

		if (not os.path.isdir(UnitID_path)):
			continue

		print('Entering directory ' + UnitID_path[len_topdir:-1])
		#print('\n')

		hour_files_path = UnitID_path + '1/'

		if (not os.path.exists(hour_files_path)):
			continue

		print('Entering directory ' + hour_files_path[len_topdir:-1])
		#print('\n')

		convert_hourly(hour_files_path, day_path)

		print('Leaving directory ' + hour_files_path[len_topdir:-1])
		#print('\n')

		print('Leaving directory ' + UnitID_path[len_topdir:-1])
		#print('\n')

	del UnitID_folders_list
	print('Leaving directory ' + day_path[len_topdir:-1])
	#print('\n')

	return



def reftek2sac():

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

		station_folders_list = os.listdir(stage_path)

		for station_folder in station_folders_list:

			station_path = stage_path + station_folder + '/'
			print('Entering directory ' + station_path[len_topdir:-1])
			#print('\n')

			if (not os.path.isdir(station_path)):
				continue

			day_folders_list = os.listdir(station_path)

			pool = ThreadPool()
			pool.map(convert_daily, day_folders_list)
			pool.close()
			pool.join()

			del day_folders_list
			print('Leaving directory ' + station_path[len_topdir:-1])
			#print('\n')

		del station_folders_list
		print('Leaving directory ' + stage_path[len_topdir:-1])
		#print('\n')

	del stage_folders_list

	return



if __name__ == '__main__':

	#print('\n')
	print('reftek2sac: ')
	print('This program converts files from reftek to sac format using the ObsPy (parallel version)')
	print('Youshan Liu at Institute of Geology and Geophysics, Chinese Academy of Sciences')
	print('Welcome to send any bugs and suggestions to ysliu@mail.iggcas.ac.cn')
	print('\n')

	starttime = UTCDateTime()

	# absolution path
	#current_path = os.getcwd()

	# read station information
	sta = read_station_list(current_path + '/' + station_list)

	# convert reftek file to sac format
	reftek2sac()

	endtime = UTCDateTime()

	elapsed_time = (endtime - starttime)

	print('\n')
	print('Start   time : %s' % starttime)
	print('End     time : %s' % endtime)
	print('Elapsed time : %f hours' % (elapsed_time / 3600.0))



