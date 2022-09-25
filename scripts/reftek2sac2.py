#!/usr/bin/env python3
# -*- coding: UTF-8 -*-


'''
reftek2sac

This program converts files from reftek to sac format using the ObsPy

Date: 22/10/2020
Author: Youshan Liu
Affiliation: Institute of Geology and Geophysics, Chinese Academy of Sciences


folders structure:
./your_data_folder/stage folder/station folder/Reftek UnitID number/day folder/stream/reftek files

for example:
./Raw2/2007_276_2008_005/NE00/2007276/9F78/1



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

# direction of the reftek data
data_folder = './Raw2'

# station information file
station_list = 'NEsta.lst'


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
is_demean = False

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

	return sac_filename



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



def convert_hourly(hour_files_path, day_path):

	hour_files_list = os.listdir(hour_files_path)

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
					print("Error: downsampling sampling rate cannot large than original sampling rate !")
				decimate_factor = int(df / downsampling_rate)
				if (abs(df - (decimate_factor*downsampling_rate)) > 0.0):
					print("Error: decimate factor can only be integer !")
				# Nyquist frequency of downsampling rate
				freq_lowpass = 0.499 * tr.stats.sampling_rate / decimate_factor
				if (not(is_bandpass and (fhigh <= freq_lowpass))):
					tr.filter('lowpass', freq=freq_lowpass, corners=2, zerophase=True)
				tr.decimate(factor=decimate_factor, strict_length=False, no_filter=True)


			# get the index of station name in station list
			#station_name = tr.stats.station
			#res = [sta.name.index(x) for x in sta.name if x.upper() == station_name.upper()]
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
			#	for j in range(len(sta.name)):
			#		res = findstr(hour_files_path[len_rootdir:-1], sta.name[j])
			#		if ([] != res):
			#			ipos = j
			#			break
			#	if (-1 == ipos):
			#		print("Error: station %s is not in the station list or filed 'staion' in reftek header is NULL" % station_name)
			#		return


			ipos = -1
			for j in range(len(sta.name)):
				res = findstr(hour_files_path[len_rootdir:-1], sta.name[j])
				if ([] != res):
					ipos = j
					break
			if (-1 == ipos):
				print("Error: station %s is not in the station list or filed 'staion' in reftek header is NULL" % station_name)
				return



			network_name = sta.netwk[ipos]
			tr.stats.station = sta.name[ipos]

			sac_filename = create_sac_filename(tr.stats, network_name, channel_name[i], sac_suffix)

			outfile = day_path + sac_filename


			if (dryrun):
				if (os.path.exists(outfile) and os.path.isfile(outfile)):
					os.remove(outfile)
			else:

				# set channel name
				tr.stats.channel = channel_name[i]


				sac = SACTrace.from_obspy_trace(tr)


				# write stla, stlo, stel into sac header
				sac.stla = sta.stla[ipos]
				sac.stlo = sta.stlo[ipos]
				sac.stel = sta.stel[ipos]
				# write station name
				sac.kstnm = sta.name[ipos]
				# write network name
				sac.knetwk = network_name


				sac.nzyear = tr.stats.starttime.year
				sac.nzjday = tr.stats.starttime.julday
				sac.nzhour = tr.stats.starttime.hour
				sac.nzmin = tr.stats.starttime.minute
				sac.nzsec = tr.stats.starttime.second
				sac.nzmsec = int(tr.stats.starttime.microsecond*1.e-3)
				sac.b = 0
				#sac.reftime += sac.b
				sac.reftime = tr.stats.starttime



				# write sac
				sac.write(outfile)
				#tr.write(outfile, format='sac')

				del sac

			del tr

		del st
		print(hour_file + ' is done ... \n')

	del hour_files_list

	return



def convert_daily(station_path):

	day_folders_list = os.listdir(station_path)

	for day_folder in day_folders_list:

		day_path = station_path + day_folder + '/'
		print('\t\tEntering directory ' + day_path[len_rootdir:-1])
		#print('\n')

		UnitID_folders_list = os.listdir(day_path)

		for UnitID in UnitID_folders_list:

			UnitID_path = day_path + UnitID + '/'

			if (not os.path.isdir(UnitID_path)):
				continue

			print('\t\t\tEntering directory ' + UnitID_path[len_rootdir:-1])
			#print('\n')

			hour_files_path = UnitID_path + '1/'

			if (not os.path.exists(hour_files_path)):
				continue

			print('\t\t\t\tEntering directory ' + hour_files_path[len_rootdir:-1])
			#print('\n')

			convert_hourly(hour_files_path, day_path)

			print('\t\t\t\tLeaving directory ' + hour_files_path[len_rootdir:-1])
			#print('\n')

			print('\t\t\tLeaving directory ' + UnitID_path[len_rootdir:-1])
			#print('\n')

		del UnitID_folders_list
		print('\t\tLeaving directory ' + day_path[len_rootdir:-1])
		#print('\n')

	del day_folders_list

	return



def reftek2sac(current_path):

	global len_rootdir, sta, sac_suffix

	rootdir = data_folder + '/'
	#rootdir = current_path + '/' + data_folder + '/'
	len_rootdir = len(rootdir)

	stage_folders_list = os.listdir(rootdir)


	sac_suffix = '.SAC'

	# convert reftek to sac
	for stage_folder in stage_folders_list:

		stage_path = rootdir + stage_folder + '/'
		print('Entering directory ' + stage_path[len_rootdir:-1])
		#print('\n')

		if (not os.path.exists(stage_path)):
			continue

		station_folders_list = os.listdir(stage_path)

		for station_folder in station_folders_list:

			station_path = stage_path + station_folder + '/'
			print('\tEntering directory ' + station_path[len_rootdir:-1])
			#print('\n')

			if (not os.path.exists(station_path)):
				continue

			convert_daily(station_path)

			print('\tLeaving directory ' + station_path[len_rootdir:-1])
			#print('\n')

		del station_folders_list
		print('Leaving directory ' + stage_path[len_rootdir:-1])
		#print('\n')

	del stage_folders_list

	return



if __name__ == '__main__':

	print('\n')
	print('reftek2sac: ')
	print('This program converts files from reftek to sac format using the ObsPy (serial version)')
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

	# convert reftek file to sac format
	reftek2sac(current_path)

	endtime = UTCDateTime()

	elapsed_time = (endtime - starttime)

	print('Start   time : %s' % starttime)
	print('End     time : %s' % endtime)
	print('Elapsed time : %f hours \n' % (elapsed_time / 3600.0))



