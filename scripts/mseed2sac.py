#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

'''
mseed2sac

This program converts files from mseed to sac format using the ObsPy

Date: 22/10/2020
Author: Youshan Liu
Affiliation: Institute of Geology and Geophysics, Chinese Academy of Sciences


folders structure:
./your_data_folder/period folder/station folder/mseed UnitID number/day folder/stream/mseed files

for example:
./Raw/NE00_2007_276_2008_005/NE00/2007276/9F78/1



mseed 130 Disk Directory Structure

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
import glob
from obspy import read
from obspy.core import UTCDateTime
from obspy.io.sac import SACTrace



##############################################################
##############################################################
##############################################################
# some options to be set


# dryrun just for debug
# dryrun = True  => just print some direction information not to write sac files
# dryrun = False => do the actual files conversion
dryrun = False

# direction of the mseed data
#data_folder = 'test2'
#data_folder = 'JD5'
#output_folder = 'Group5_EPS'
data_folder = 'JD.Group5.EPS.Raw'
output_folder = 'SAC.Group5.EPS'

# station information file
station_list = 'stainfo.lst'


# component list to be converted
# only Z-component
#component_list = ['Z']
# only N-component
#component_list = ['N']
# only E-component
#component_list = ['E']
# three-components
component_list = ['Z', 'N', 'E']

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

seconds_segment = 24*3600
##############################################################
##############################################################
##############################################################


# channel name, it consists of band code, instrument code, orientation code
channel_name = ['BHZ', 'BHN', 'BHE']


##############################################################
##############################################################
##############################################################
# read station list
def read_station_list(filename):

	class station:
		def __init__(self):
			self.name = []
			self.stla = []
			self.stlo = []
			self.stel = []
			self.netwk = []


	name = ''
	stla = ''
	stlo = ''
	stel = ''
	netwk = ''
	sta = station()
	with open(filename, 'r') as f:
		lines = f.readlines()
	for line in lines[1:]:
		try:
			line_splited = line.split()
			if ([] == line_splited):
				continue
			netwk, name, stla, stlo, stel = line_splited
		except:
			raise Exception('Format error in %s ! ' % filename)
		sta.name.append(name)
		sta.netwk.append(netwk)
		sta.stla.append(float(stla))
		sta.stlo.append(float(stlo))
		sta.stel.append(float(stel))

	del lines

	return sta



def create_sac_filename(stats, network_name, channel_name, sac_suffix):

	time = UTCDateTime(stats.starttime)
	yyyy = '%4.4d' % time.year
	ddd  = '%3.3d' % time.julday
	hh   = '%2.2d' % time.hour
	mm   = '%2.2d' % time.minute
	ss   = '%2.2d' % time.second
	fff  = '%3.3d' % (0.001*time.microsecond)

	sac_filename = yyyy + '.' + ddd + '.' + hh + '.' + \
		           mm + '.' + ss + '.' + fff + '.' + \
		           network_name + '.' + stats.station + '.' + \
		           channel_name + sac_suffix

	day_path = yyyy + ddd

	return sac_filename, day_path



def findstr(str_src, str_target):

	'''
	findstr(str_src, str_target)
	find the index of a substring in a string
	str_src    -> source string
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



def convert_hourly(hour_files_path):

	hour_files_list = os.listdir(hour_files_path)
	#hour_files_list = glob.glob(hour_files_path + 'E*.*')

	for hour_file in hour_files_list:

		if ('.LST' == hour_file[-4:] or '.LOG' == hour_file[-4:]):
			continue

		#if ('EE' != hour_file[-4:-2]):
		#	continue

		# if this file is not a reftek format
		is_sacfile = False
		starts = [each.start() for each in re.finditer(sac_suffix, hour_file.upper())]
		ends = [start+len(sac_suffix) - 1 for start in starts]
		span = [(start, end) for start,end in zip(starts, ends)]
		is_sac = is_sacfile and (len(span) >= 1)

		if (is_sacfile):
			continue


		#try:
		#	has_dot_in_filename = hour_file.index('.')
		#except:
		#	has_dot_in_filename = False
		#if (has_dot_in_filename):
		#	continue


		infile = hour_files_path + hour_file
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



		for i in range(0,min(len(channel_name),len(st))):

			if (channel_name[i][-1] not in component_list):
				continue

			try:
				tr = st[i]
				# remove round error
				tr.stats.delta = round(tr.stats.delta*1e6)*1.e-6
			except:
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
				decimate_factor = int(df / downsampling_rate)
				if (abs(df - (decimate_factor*downsampling_rate)) > 0.0):
					print("Error: decimate factor can only be integer !")
				# Nyquist frequency of the downsampling rate
				freq_lowpass = 0.499 * tr.stats.sampling_rate / decimate_factor
				if (not(is_bandpass and (fhigh <= freq_lowpass))):
					tr.filter('lowpass', freq=freq_lowpass, corners=2, zerophase=True)
				tr.decimate(factor=decimate_factor, strict_length=False, no_filter=True)


			# get the index of station name in station list
			#station_name = tr.stats.station
			#res = [sta.name.index(x) for x in sta.name if x.upper() == station_name.upper()]
			#if ([] != res):
			#	# first find station name in mseed head.
			#	# if the field "station" in mseed header is NULL, then try extract station name from folder name
			#	try:
			#		ipos = res[0]
			#	except:
			#		raise Exception('Error: station %s is not in the station list' % station_name)
			#else:
			#	#print("Warning: station field in mseed header is NULL")
			#	# try to extract station name based on folder name
			#	ipos = -1
			#	for j in range(len(sta.name)):
			#		res = findstr(hour_files_path[nrootdir:-1], sta.name[j])
			#		if ([] != res):
			#			ipos = j
			#			break
			#	if (-1 == ipos):
			#		print("Error: station %s is not in the station list or filed 'staion' in mseed header is NULL" % station_name)
			#		return


			ipos = -1
			for j in range(len(sta.name)):
				res = findstr(hour_files_path[nrootdir:-1], sta.name[j])
				if ([] != res):
					ipos = j
					break
			if (-1 == ipos):
				print("Error: station %s is not in the station list or filed 'staion' in mseed header is NULL" % station_name)
				return


			network_name = sta.netwk[ipos]
			tr.stats.station = sta.name[ipos]
			# set channel name
			#tr.stats.channel = channel_name[i]


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
				#midtime = starttime + 0.5*seconds_segment
				#endtime = UTCDateTime(midtime.year, midtime.month, midtime.day, 23, 59, 59, 999999)
				#endtime = UTCDateTime(year=starttime.year, julday=starttime.julday, \
				#                 hour=23, minute=59, second=59, microsecond=999999)
				#time_duration = endtime - starttime
				#print(starttime)
				#print(endtime)
				#print(starttime.day)

				iend = min(ibeg + int((endtime-starttime)*df), npts_org)


				tr_out = tr.copy()
				tr_out.data = tr.data[ibeg:iend+1]
				tr_out.stats.starttime = starttime


				#sac_filename, day_path = create_sac_filename(tr_out.stats, network_name, channel_name[i], sac_suffix)
				sac_filename, day_path = create_sac_filename(tr_out.stats, network_name, tr.stats.channel, sac_suffix)


				sac_path = output_path + '/' + sta.name[ipos] + '/' + day_path + '/'
				if (not os.path.exists(sac_path)):
					os.makedirs(sac_path)


				#print(sac_path)
				outfile = sac_path + sac_filename
				#print(outfile)



				if (dryrun):
					if (os.path.exists(outfile) and os.path.isfile(outfile)):
						os.remove(outfile)
				else:

					# set channel name
					#tr_out.stats.channel = channel_name[i]
					#print('channel = %s, %s, %s' %(tr.stats.channel, tr_out.stats.channel, channel_name[i]))


					sac = SACTrace.from_obspy_trace(tr_out)


					# write stla, stlo, stel into sac header
					sac.stla = sta.stla[ipos]
					sac.stlo = sta.stlo[ipos]
					sac.stel = sta.stel[ipos]
					# write station name
					sac.kstnm = sta.name[ipos]
					# write network name
					sac.knetwk = network_name


					sac.nzyear = tr_out.stats.starttime.year
					sac.nzjday = tr_out.stats.starttime.julday
					sac.nzhour = tr_out.stats.starttime.hour
					sac.nzmin = tr_out.stats.starttime.minute
					sac.nzsec = tr_out.stats.starttime.second
					sac.nzmsec = int(tr_out.stats.starttime.microsecond*1.e-3)
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
				endtime = starttime + seconds_segment - dt
				ibeg = iend + 1
				#print(tr_out.stats.npts)
				#print(tr_out.stats.starttime)
				#print(tr_out.stats.endtime)
				##print('\n')


				if (ibeg > npts_org):
					break

			del tr

		del st
		print(hour_file + ' is done ... \n')

	del hour_files_list

	return



def convert_daily(station_period_path):

	day_folders_list = os.listdir(station_period_path)

	for day_folder in day_folders_list:

		day_path = station_period_path + day_folder + '/'
		print('Entering directory ' + day_path[nrootdir:-1])
		print('\n')

		if (not os.path.exists(day_path)):
			continue

		UnitID_folders_list = os.listdir(day_path)

		for UnitID in UnitID_folders_list:

			UnitID_path = day_path + UnitID + '/'

			if (not os.path.isdir(UnitID_path)):
				continue

			print('Entering directory ' + UnitID_path[nrootdir:-1])
			print('\n')

			#hour_files_path = UnitID_path + '1/'

			#if (not os.path.exists(hour_files_path)):
			#	continue

			#print('Entering directory ' + hour_files_path[nrootdir:-1])
			##print('\n')

			convert_hourly(UnitID_path)
			#convert_hourly(hour_files_path, day_path)

			#print('Leaving directory ' + hour_files_path[nrootdir:-1])
			##print('\n')

			print('Leaving directory ' + UnitID_path[nrootdir:-1])
			print('\n')

		del UnitID_folders_list
		print('Leaving directory ' + day_path[nrootdir:-1])
		print('\n')

	del day_folders_list

	return



def mseed2sac(current_path):

	global nrootdir, sta, sac_suffix, output_path

	rootdir = current_path + '/' + data_folder + '/'
	nrootdir = len(rootdir)

	output_path = current_path + '/' + output_folder
	if (not os.path.exists(output_path)):
		os.makedirs(output_path)


	period_folders_list = os.listdir(rootdir)


	sac_suffix = '.SAC'

	# convert mseed to sac
	for station_period_folder in period_folders_list:

		station_period_path = rootdir + station_period_folder + '/'
		print('Entering directory ' + station_period_path[nrootdir:-1])
		#print('\n')

		if (not os.path.exists(station_period_path)):
			continue

		convert_daily(station_period_path)

		print('Leaving directory ' + station_period_path[nrootdir:-1])
		#print('\n')

	del period_folders_list

	return



if __name__ == '__main__':

	#print('\n')
	print('mseed2sac: ')
	print('This program converts files from mseed to sac format using the ObsPy (serial version)')
	print('Youshan Liu at Institute of Geology and Geophysics, Chinese Academy of Sciences')
	print('Welcome to send any bugs and suggestions to ysliu@mail.iggcas.ac.cn')
	print('\n\n')

	starttime = UTCDateTime()

	# get current path
	# absolution path
	#current_path = os.getcwd()
	# relative path
	current_path = '.'

	# read station information
	sta = read_station_list(current_path + '/' + station_list)

	# convert mseed file to sac format
	mseed2sac(current_path)

	endtime = UTCDateTime()

	elapsed_time = (endtime - starttime)

	print('Start   time : %s' % starttime)
	print('End     time : %s' % endtime)
	print('Elapsed time : %f hours \n' % (elapsed_time / 3600.0))



