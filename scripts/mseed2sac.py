#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

'''
mseed2sac

This program converts files from MSEED to SAC format using the ObsPy

Date: 22/10/2020
Author: Youshan Liu
Affiliation: Institute of Geology and Geophysics, Chinese Academy of Sciences


directory structure:
./top folder/stage folder/station folder/DAS number/MSEED files

for example:
./DATA_Raw/DATA_1/5322/509



Directory Structure

 DAS number
 |
 v
  __
 |  |
 \509
     MSEED files
     |
     v
     \E???????.???

'''




import os
import re
import glob
from math import *
from obspy import read
from obspy import UTCDateTime
from obspy.io.sac import SACTrace



##############################################################
##############################################################
##############################################################
# some options to be set


# dryrun just for debug
# dryrun = True  => just print some directory information not to write sac files
# dryrun = False => do the actual files conversion
dryrun = False

# directory of the mseed data
input_path = './JD.Group5.EPS.Raw'
output_path = './SAC.Group5.EPS'

# station information file
station_list = './stainfo.lst'


# component list to be converted
# only Z-component
#components = ['Z']
# only N-component
#components = ['N']
# only E-component
#components = ['E']
# three-components
components = ['Z', 'N', 'E']

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
# the downsampling_rate, downsampling frequency
downsampling_rate = 20.0


segment_in_seconds = 3600*24

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
	netwk, name, stla, stlo, stel = None, None, None, None, None
	for line in lines[1:]:
		try:
			line_splited = line.split()
			if ([] == line_splited):
				continue
			netwk, name, stla, stlo, stel = line_splited
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

	day_path = yyyy + ddd

	return sac_filename, day_path



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



def convert_files(segment_files_path):

	segment_files_list = os.listdir(segment_files_path)
	#segment_files_list = glob.glob(segment_files_path + 'E*.*')

	idx = findstr(segment_files_path, '/')

	for segment_file in segment_files_list:

		if ('.LST' == segment_file[-4:] or '.LOG' == segment_file[-4:]):
			continue

		#if ('EE' != segment_file[-4:-2]):
		#	continue

		# if this file is not a reftek format
		is_sacfile = False
		starts = [each.start() for each in re.finditer(sac_suffix, segment_file.upper())]
		ends = [start+len(sac_suffix) - 1 for start in starts]
		span = [(start, end) for start,end in zip(starts, ends)]
		is_sac = is_sacfile and (len(span) >= 1)

		if (is_sacfile):
			continue


		#try:
		#	has_dot_in_filename = segment_file.index('.')
		#except:
		#	has_dot_in_filename = False
		#if (has_dot_in_filename):
		#	continue


		infile = segment_files_path + segment_file
		#st = read(infile, format='MSEED', check_compression=False, component_codes=['0', '1', '2'])
		try:
			# faster for specific file format
			st = read(infile, format='MSEED', check_compression=False, component_codes=['0', '1', '2'])
			#st = read(infile, component_codes=['0', '1', '2'])
			#st = read(infile, component_codes=['Z', 'N', 'E'])
			# cancel out the sort in obspy
			#st.sort(reverse=True)
		except:
			#print('Fail to read data')
			continue


		if (len(st) > 3):
			st.sort(['starttime'])
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


			#hours = tr.stats.starttime.hour*60 + tr.stats.starttime.minute
			#hour = tr.stats.starttime.hour + ceil(tr.stats.starttime.minute/60.0)
			#hour = tr.stats.starttime.hour + int((tr.stats.starttime.minute + 29.9)/60.0))
			starttime = tr.stats.starttime
			endtime = tr.stats.endtime
			#hour = starttime.hour + ceil((starttime.minute + (starttime.second + starttime.microsecond*1.e-6)/60.0)/60.0)
			#hour = starttime.hour + int((starttime.minute + (starttime.second + starttime.microsecond*1.e-6)/60.0 + 45)/60.0)
			#hour = starttime.hour + int((starttime.minute + (starttime.second + starttime.microsecond*1.e-6)/60.0 + 40)/60.0)
			min2hour = int((starttime.minute + (starttime.second + starttime.microsecond*1.e-6)/60.0 + 40)/60.0)
			if (0 == min2hour):
				#if (is_decimate):
				##df = tr.stats.sampling_rate
				##if (downsampling_rate > df):
				##print("Error: downsampling rate cannot large than original sampling rate !")
				##decimate_factor = int(df / downsampling_rate)
				##if (abs(df - (decimate_factor*downsampling_rate)) > 0.0):
				##print("Error: decimate factor can only be integer !")
				#dtinus = decimate_factor * dt * 1e6
				#else:
				#dtinus = dt * 1e6
				dtinus = 1e6 / downsampling_rate
				#microsecond = ceil(starttime.microsecond / dtinus) * dtinus
				##starttime_first = starttime
				#starttime_first = UTCDateTime(starttime.year, starttime.month, starttime.day, starttime.hour, starttime.minute, starttime.second, microsecond, strict=False)
				sec = round(starttime.second*1e6 + starttime.microsecond)
				sec = ceil(sec / dtinus) * dtinus
				second = int(sec * 1.e-6)
				microsecond = int(sec - second*1e6)
				#starttime_first = starttime
				starttime_first = UTCDateTime(starttime.year, starttime.month, starttime.day, starttime.hour, starttime.minute, second, microsecond, strict=False)
			else:
				starttime_first = UTCDateTime(starttime.year, starttime.month, starttime.day, starttime.hour + min2hour, 0, 0, 0)
			if ((endtime.minute + (endtime.second + endtime.microsecond*1.e-6)/60.0) > 40):
				endtime_last = endtime
			else:
				endtime_last = UTCDateTime(endtime.year, endtime.month, endtime.day, endtime.hour, 0, 0, 0)
			if (starttime_first < endtime_last):
				tr = tr.slice(starttime=starttime_first, endtime=endtime_last, nearest_sample=False)
			else:
				del tr
				continue



			# some preprocesses
			if (is_demean):
				tr.detrend(type='demean')
			if (is_detrend):
				tr.detrend(type='linear')
			#if (is_denonlinear):
			#	#tr.detrend('polynomial', order=50)
			#	tr.detrend('spline', order=3, dspline=5)


			# bandpass
			if (is_bandpass):
				tr.filter('bandpass', freqmin=flow, freqmax=fhigh, corners=2, zerophase=is_zerophase)


			# downsampling
			if (is_decimate):
				df = tr.stats.sampling_rate
				if (downsampling_rate > df):
					print("Error: downsampling rate cannot large than original sampling rate !")
					continue
				decimate_factor = int(df / downsampling_rate)
				if (abs(df - (decimate_factor*downsampling_rate)) > 0.0):
					print("Error: decimate factor can only be integer !")
					continue
				if (decimate_factor > 1):
					# Nyquist frequency of the downsampling rate
					freq_lowpass = 0.49 * tr.stats.sampling_rate / decimate_factor
					if (not(is_bandpass and (fhigh <= freq_lowpass))):
						tr.filter('lowpass', freq=freq_lowpass, corners=2, zerophase=True)
					tr.decimate(factor=decimate_factor, strict_length=False, no_filter=True)


			# get the index of station name in station list
			#station_path = tr.stats.station
			#res = [sta.stnm.index(x) for x in sta.stnm if x.upper() == station_path.upper()]
			#if ([] != res):
			#	# first find station name in mseed head.
			#	# if the field "station" in mseed header is NULL, then try extract station name from folder name
			#	try:
			#		ipos = res[0]
			#	except:
			#		raise Exception('Error: station %s is not in the station list' % station_path)
			#else:
			#	#print("Warning: station field in mseed header is NULL")
			#	# try to extract station name based on folder name
			#	ipos = -1
			#	for j in range(len(sta.stnm)):
			#		res = findstr(segment_files_path[idx[-3]+1:idx[-2]], sta.stnm[j])
			#		if ([] != res):
			#			ipos = j
			#			break
			#	if (-1 == ipos):
			#		print("Error: station %s is not in the station list or field 'station' in mseed header is NULL" % station_path)
			#		return


			ipos = -1
			station_path = segment_files_path[idx[-3]+1:idx[-2]]
			for j in range(len(sta.stnm)):
				res = findstr(station_path, sta.stnm[j])
				#res = findstr(segment_files_path[len_topdir:-1], sta.stnm[j])
				if ([] != res):
					ipos = j
					break
			if (-1 == ipos):
				print("Error: station folder %s does not include the name of this station or this station is missing in the stainfo.lst" % station_path)
				return


			station = sta.stnm[ipos]
			network = sta.netwk[ipos]
			tr.stats.station = sta.stnm[ipos]
			# set channel name
			#tr.stats.channel = channels[i]


			npts_org = tr.stats.npts
			dt = tr.stats.delta
			df = tr.stats.sampling_rate
			starttime_org = tr.stats.starttime
			endtime_org = tr.stats.endtime


			ibeg = 0
			starttime = starttime_org
			endtime = UTCDateTime(starttime.year, starttime.month, starttime.day, 23, 59, 59, 999999)


			iday1 = starttime_org.julday
			iday2 = endtime_org.julday
			#print('day1 = %d' % iday1)
			#print('day2 = %d' % iday2)
			for iday in range(iday1,iday2+1):

				#print('\niday = %d' % iday)
				#midtime = starttime + 0.5*segment_in_seconds
				#endtime = UTCDateTime(midtime.year, midtime.month, midtime.day, 23, 59, 59, 999999)
				#endtime = UTCDateTime(year=starttime.year, julday=starttime.julday, \
				# hour=23, minute=59, second=59, microsecond=999999)
				#time_duration = endtime - starttime
				#print(starttime)
				#print(endtime)
				#print(starttime.day)

				iend = min(ibeg + int((endtime-starttime)*df), npts_org)


				tr_out = tr.copy()
				tr_out.data = tr.data[ibeg:iend+1]
				tr_out.stats.starttime = starttime


				sac_filename, day_path = create_sac_filename(tr_out.stats)


				sac_path = output_path + '/' + sta.stnm[ipos] + '/' + day_path + '/'
				if (not os.path.exists(sac_path)):
					os.makedirs(sac_path, exist_ok=True)


				#print(sac_path)
				outfile = sac_path + sac_filename
				#print(outfile)



				if (dryrun):
					if (os.path.exists(outfile) and os.path.isfile(outfile)):
						os.remove(outfile)
				else:

					# set channel name
					#tr_out.stats.channel = channels[i]
					#print('channel = %s, %s, %s' %(tr.stats.channel, tr_out.stats.channel, channels[i]))


					sac = SACTrace.from_obspy_trace(tr_out)


					# write stla, stlo, stel into sac header
					sac.stla = sta.stla[ipos]
					sac.stlo = sta.stlo[ipos]
					sac.stel = sta.stel[ipos]
					# write station name
					sac.kstnm = station
					# write network name
					sac.knetwk = network


					sac.nzyear = tr_out.stats.starttime.year
					sac.nzjday = tr_out.stats.starttime.julday
					sac.nzhour = tr_out.stats.starttime.hour
					sac.nzmin = tr_out.stats.starttime.minute
					sac.nzsec = tr_out.stats.starttime.second
					#sac.nzmsec = int(tr_out.stats.starttime.microsecond*1.e-3)
					sac.nzmsec = round(tr_out.stats.starttime.microsecond*1.e-3)
					sac.b = 0
					#sac.reftime += sac.b
					#sac.reftime = tr_out.stats.starttime


					# write sac
					#print(infile)
					#print(outfile)
					#print('data_length = %d' % len(tr_out.data))
					sac.write(outfile)
					#tr.write(outfile, format='sac')

					del sac


				del tr_out
				starttime = starttime + (iend-ibeg+1)*dt
				endtime = starttime + segment_in_seconds - dt
				ibeg = iend + 1
				#print(tr_out.stats.npts)
				#print(tr_out.stats.starttime)
				#print(tr_out.stats.endtime)
				#print('\n')


				if (ibeg > npts_org):
					break

			del tr

		del st
		print('%s is done' % segment_file)

	del segment_files_list
	print('%s is done' % segment_files_path[idx[0]+1:idx[-1]])

	return



def convert_station(station_stage_path):

	station_folders_list = os.listdir(station_stage_path)

	for station_folder in station_folders_list:

		station_path = station_stage_path + station_folder + '/'
		print('Entering directory ' + station_path[len_topdir:-1])
		#print('\n')

		if (not os.path.exists(station_path)):
			continue

		UnitID_folders_list = os.listdir(station_path)

		for UnitID in UnitID_folders_list:

			UnitID_path = station_path + UnitID + '/'

			if (not os.path.isdir(UnitID_path)):
				continue

			print('Entering directory ' + UnitID_path[len_topdir:-1])
			#print('\n')

			#segment_files_path = UnitID_path + '1/'

			#if (not os.path.exists(segment_files_path)):
			#	continue

			#print('Entering directory ' + segment_files_path[len_topdir:-1])
			#print('\n')

			convert_files(UnitID_path)
			#convert_hourly(segment_files_path, day_path)

			#print('Leaving directory ' + segment_files_path[len_topdir:-1])
			#print('\n')

			print('Leaving directory ' + UnitID_path[len_topdir:-1])
			#print('\n')

		del UnitID_folders_list
		print('Leaving directory ' + station_path[len_topdir:-1])
		#print('\n')

	del station_folders_list

	return



def mseed2sac():

	global len_topdir, sac_suffix

	len_topdir = len(input_path) + 1

	if (not os.path.exists(output_path)):
		os.makedirs(output_path)


	stage_folders_list = os.listdir(input_path)


	sac_suffix = '.SAC'

	# convert mseed to SAC
	for station_stage_folder in stage_folders_list:

		station_stage_path = input_path + '/' + station_stage_folder + '/'
		print('Entering directory ' + station_stage_path[len_topdir:-1])
		#print('\n')

		if (not os.path.exists(station_stage_path)):
			continue

		convert_daily(station_stage_path)

		print('Leaving directory ' + station_stage_path[len_topdir:-1])
		#print('\n')

	del stage_folders_list

	return



if __name__ == '__main__':

	print('\n')
	print('mseed2sac: ')
	print('This program converts files from mseed to sac format using the ObsPy (serial version)')
	print('Youshan Liu at Institute of Geology and Geophysics, Chinese Academy of Sciences')
	print('Welcome to send any bugs and suggestions to ysliu@mail.iggcas.ac.cn')
	print('\n')

	starttime = UTCDateTime()

	# absolution path
	#current_path = os.getcwd()

	# read station information
	sta = read_station_list(station_list)

	# convert mseed file to SAC format
	mseed2sac()

	endtime = UTCDateTime()

	elapsed_time = (endtime - starttime)

	print('\n')
	print('Start   time : %s' % starttime)
	print('End     time : %s' % endtime)
	print('Elapsed time : %f hours' % (elapsed_time / 3600.0))



