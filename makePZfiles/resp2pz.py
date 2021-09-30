#!/usr/bin/env python3
# -*- coding: UTF-8 -*-


'''

This program automatically generate PZ files based on the RESP files.


create PZ files for removing instrument response

RESP file format, please see the following website for details:
http://ds.iris.edu/ds/nodes/dmc/data/formats/resp/

RESP files can be downloaded from the following website:
http://ds.iris.edu/NRL/sensors/guralp/guralp_sensors.html


'''



import os
import re
import xlrd
import xdrlib, sys
from obspy.core import UTCDateTime



# path of the resp pool for sensor
sensor_resp_pool_path = './Sensor_RESP_pool'

# path of the resp pool for DAS
das_resp_pool_path = './DAS_RESP_pool'

# the used resp files list for sensor
sensor_resp_listname = './sensor_used.lst'

# the used resp files list for DAS
das_resp_listname = './das_used.lst'

# path of the generated PZ files
output_path = '../PZs_all'

# station log file
station_logfile = './NCISP6.xls'


# channel name, it consists of band code, instrument code, orientation code
channel_name = ['BHZ']


# the decimate rate must be interger so that the sampling_rate is divisible by it
downsampling_rate = 10.0



def read_resp(filename):

	'''
	read_resp(filename)
	read standard RESP file
	filename      -> RESP file name
	convert2displ -> tansfer to displacement
	'''

	# define regular expression
	A0_regexp = re.compile(r"^[ \t]*[^#]+[\w \t]+(A0)+[ \t]*(normalization)+[ \t]*((factor)+:[ \t]*([+-]*[0-9]*\.[0-9]+[EeDd]+[+-]*[0-9]+)||(^(\+)?\d(\.\d+)?]))+")
	gain_regexp = re.compile(r"^[ \t]*[^#]+[\w \t]+(Gain)+:[ \t]*([+-]*[0-9]*\.[0-9]+[EeDd]+[+-]*[0-9]+)+")
	zeros_regexp = re.compile(r"^[ \t]*[# \t]+(Complex)+[ \t]*(zeroes):+[ \w \t]*")
	poles_regexp = re.compile(r"^[ \t]*[# \t]+(Complex)+[ \t]*(poles):+[ \w \t]*")
	nzeros_regexp = re.compile(r"^[ \t]*[^#]+[\w \t]+(Number)+[ \t]*(of)+[ \t]*(zeroes)+:[0-9]*")
	npoles_regexp = re.compile(r"^[ \t]*[^#]+[\w \t]+(Number)+[ \t]*(of)+[ \t]*(poles)+:[0-9]*")
	input_unit_regexp = re.compile(r"^[ \t]*[^#]+[\w \t]+(Response)+[ \t]*(in)+[ \t]*(units)+[ \t]*(lookup)+:[ \w \t]*")
	sensitivity_regexp = re.compile(r"^[ \t]*[^#]+[\w \t]+(Sensitivity)+:[ \t]*([+-]*[0-9]*\.[0-9]+[EeDd]+[+-]*[0-9]+)+")


	with open(filename, 'r+', encoding='UTF-8') as f:
		lines = f.readlines()


	poles = []
	zeros = []
	npoles = 0
	nzeros = 0
	gain_all = 1.0
	A0_factor = 1.0
	iline_zeros = -1
	iline_poles = -1
	veloc2displ = False
	accel2displ = False
	for i in range(len(lines)):
		line = lines[i]
		# read input units
		match_obj = input_unit_regexp.search(line)
		if (None != match_obj):
			line_splited = line.split(':')
			if (findstr(line_splited[1], 'M/S')):
				veloc2displ = True
			if (findstr(line_splited[1], 'M/S**2')):
				accel2displ = True
		# read number of zeroes
		match_obj = nzeros_regexp.search(line)
		if (None != match_obj):
			line_splited = line.split(':')
			nzeros = int(line_splited[1])
		# read number of poles
		match_obj = npoles_regexp.search(line)
		if (None != match_obj):
			line_splited = line.split(':')
			npoles = int(line_splited[1])
		# find line of zeros and read zeros
		match_obj = zeros_regexp.search(line)
		if (None != match_obj):
			iline_zeros = i + 2
		if ((i >= iline_zeros) and (i <= iline_zeros + nzeros -1)):
			line_splited = line.split()
			_, _, realpart, imagpart = line_splited[0:4]
			zeros.append(float(realpart) + 1j*float(imagpart))
		# find line of poles and read poles
		match_obj = poles_regexp.search(line)
		if (None != match_obj):
			iline_poles = i+2
		if ((i >= iline_poles) and (i <= iline_poles + npoles - 1)):
			line_splited = line.split()
			_, _, realpart, imagpart = line_splited[0:4]
			poles.append(float(realpart) + 1j*float(imagpart))
		# read A0 normalization factor
		match_obj = A0_regexp.search(line)
		if (None != match_obj):
			line_splited = line.split(':')
			_, A0_factor = line_splited[0:2]
			A0_factor = float(A0_factor)
		# read Gain
		match_obj = gain_regexp.search(line)
		if (None != match_obj):
			line_splited = line.split(':')
			_, gain = line_splited[0:2]
			gain_all = gain_all * float(gain)
		# read sensitivity
		match_obj = sensitivity_regexp.search(line)
		if (None != match_obj):
			line_splited = line.split(':')
			_, sensitivity = line_splited[0:2]
			sensitivity = float(sensitivity)


	# convert to displacement
	if (veloc2displ and (nzeros >= 1)):
		zeros.append(0.0+0.0j)
	if (accel2displ and (nzeros >= 1)):
		zeros.append(0.0+0.0j)
		zeros.append(0.0+0.0j)


	CONSTANT = float(A0_factor) * gain_all
	pz = { \
		'poles': poles, \
		'zeros': zeros, \
		'A0': A0_factor, \
		'sensitivity': sensitivity, \
		'gain': gain_all, \
		'CONSTANT': CONSTANT}


	return pz



def read_sensor_resp_list(infile):

	'''
	read_sensor_resp_list(filename)
	read sensor resp list file
	filename -> filename of a resp file list
	'''

	class Sensor_resp_list:
		def __init__(self):
			self.filename = []
			self.sensor_type = []
			self.sensor_period = []


	with open(infile, 'r') as f:
		lines = f.readlines()

	filename = ''
	sensor_type = ''
	sensor_period = ''
	sensor_resp_list = Sensor_resp_list()
	for line in lines[1:]:
		try:
			line_splited = line.split()
			if ([] == line_splited):
				continue
			else:
				sensor_period = ''
				if (len(line_splited) >= 3):
					filename, sensor_type, sensor_period = line_splited
				else:
					filename, sensor_type = line_splited
		except:
			raise Exception('Format error in %s ' % filename)

		sensor_resp_list.filename.append(filename)
		sensor_resp_list.sensor_type.append(sensor_type)
		sensor_resp_list.sensor_period.append(sensor_period)

	return sensor_resp_list



def read_das_resp_list(infile):

	'''
	read_das_resp_list(filename)
	read das resp list file
	filename -> filename of a resp file list
	'''

	class DAS_resp_list:
		def __init__(self):
			self.filename = []
			self.das_type = []
			self.das_gain = []
			self.das_sampling_rate = []


	with open(infile, 'r') as f:
		lines = f.readlines()

	filename = ''
	das_type = ''
	das_gain = ''
	das_sampling_rate = ''
	das_resp_list = DAS_resp_list()
	for line in lines[1:]:
		try:
			line_splited = line.split()
			if ([] == line_splited):
				continue
			else:
				if (len(line_splited) >= 4):
					filename, das_type, das_gain, das_sampling_rate = line_splited
				else:
					raise Exception('Format error in %s ' % filename)
		except:
			raise Exception('Format error in %s ' % filename)

		das_resp_list.filename.append(filename)
		das_resp_list.das_type.append(das_type)
		das_resp_list.das_gain.append(das_gain)
		das_resp_list.das_sampling_rate.append(das_sampling_rate)

	return das_resp_list



def read_resp_pool(resp_pool_path, resp_list):

	'''
	read_resp_pool(resp_pool)
	read resp files pool
	resp_pool -> RESP pool filename
	'''

	resp_pool = []
	for file in resp_list.filename:
		filename = resp_pool_path + '/' + file
		resp = read_resp(filename)
		resp_pool.append(resp)

	return resp_pool



def open_xls(filename):

	'''
	open_xls(filename)
	open a excel workbook
	resp_pool -> RESP pool filename
	'''

	try:
		fid = xlrd.open_workbook(filename)
		return fid
	except:
		raise Exception('Cannot open %s \n' % filename)



def get_xls_title(filename, colname_idx=0, idx=0):

	'''
	get_xls_title(filename, colidx=0, idx=0)
	read a excel title
	filename    -> file name of excel
	colname_idx -> index of column name
	idx         -> index of excel sheet
	'''

	fid = open_xls(filename)
	table = fid.sheets()[idx]
	nrows = table.nrows
	ncols = table.ncols
	colnames = table.row_values(colname_idx)

	return colnames



def read_xls_byindex(filename, colname_idx=0, idx=0):

	'''
	read_xls_byindex(filename, colname_idx=0, idx=0)
	read a excel sheet by sheet index
	filename    -> file name of excel
	colname_idx -> index of column name
	idx         -> index of excel sheet
	'''

	data = open_xls(filename)
	table = data.sheets()[idx]
	nrows = table.nrows
	ncols = table.ncols
	colnames = table.row_values(colname_idx)
	table_data = []
	for irow in range(1, nrows):
		row = table.row_values(irow)
		if row:
			app = {}
			for i in range(len(colnames)):
				app[colnames[i]] = row[i]
			table_data.append(app)

	return table_data, nrows, ncols



def read_xls_byname(filename, colname_idx=0, sheet_name=u'Sheet1'):

	'''
	read_xls_byname
	read a excel sheet by sheet name
	filename    -> file name of excel
	colname_idx -> index of column name
	idx         -> index of excel sheet
	'''

	data = open_excel(filename)
	table = data.sheet_by_name(sheet_name)
	nrows = table.nrows
	colnames = table.row_values(colname_idx)
	table_data = []
	for irow in range(1, nrows):
		row = table.row_values(irow)
		if row:
			app = {}
			for i in range(len(colnames)):
				app[colnames[i]] = row[i]
			table_data.append(app)

	return table_data, nrows, ncols



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



def findstr_in_list(str_list, str_target):

	'''
	findstr(str_src, str_target)
	find the index of a substring in a string list
	str_list   -> source string list
	str_target -> substring
	'''

	ipos = []
	for i in range(len(str_list)):
		str_src = str_list[i]
		idxs = findstr(str_src, str_target)
		if [] != idxs:
			ipos.append(idxs[0] + i)

	return ipos



def time2str(time):

	'''
	time2str(time)
	convert a UTCDateTime to a string
	time -> UTCDateTime object
	'''

	year   = '%4.4d' % time.year
	month  = '%2.2d' % time.month
	day    = '%2.2d' % time.day
	hour   = '%2.2d' % time.hour
	minute = '%2.2d' % time.minute
	second = '%2.2d' % time.second

	return year + '-' + month + '-' + day + 'T' +\
		hour + ':' + minute + ':' + second



def reformat_time(timestr):

	'''
	reformat_time(timestr)
	convert a time string from one format to another format
	time -> UTCDateTime object
	'''

	time = UTCDateTime(timestr)

	return time2str(time)



def weehour_daily(starttime):

	time_tmp = UTCDateTime(starttime)
	time = UTCDateTime(time_tmp.year, time_tmp.month, time_tmp.day, 0, 0, 0, 000000)
	starttime = time2str(time)

	return starttime



# read excel
def find_colname_from_xls(station_logfile, network, station_name, starttime, endtime, \
			              latitude, longitude, elevation, sensor_type, sensor_depth, \
			              das_type, das_gain, das_sampling_rate):

	# read logfile title
	xls_title = get_xls_title(station_logfile)


	# get position of network
	try:
		network_idx = findstr_in_list(xls_title, network)[0]
	except:
		raise Exception('Cannot found Network in excel title')

	# get position of station_name
	try:
		station_name_idx = findstr_in_list(xls_title, station_name)[0]
	except:
		raise Exception('Cannot found Station code in excel title')

	# get position of starttime
	try:
		starttime_idx = findstr_in_list(xls_title, starttime)[0]
	except:
		raise Exception('Cannot found Starttime in excel title')

	# get position of endtime
	try:
		endtime_idx = findstr_in_list(xls_title, endtime)[0]
	except:
		raise Exception('Cannot found Endtime in excel title')

	# get position of latitude
	try:
		latitude_idx = findstr_in_list(xls_title, latitude)[0]
	except:
		raise Exception('Cannot found Latitude in excel title')

	# get position of longitude
	try:
		longitude_idx = findstr_in_list(xls_title, longitude)[0]
	except:
		raise Exception('Cannot found Longitude in excel title')

	# get position of elevation
	try:
		elevation_idx = findstr_in_list(xls_title, elevation)[0]
	except:
		raise Exception('Cannot found Elevation in excel title')

	# get position of sensor type
	try:
		sensor_type_idx = findstr_in_list(xls_title, sensor_type)[0]
	except:
		raise Exception('Cannot found Sensor type in excel title')

	# get position of sensor depth
	try:
		sensor_depth_idx = findstr_in_list(xls_title, sensor_depth)[0]
	except:
		raise Exception('Cannot found Sensor depth in excel title')

	# get position of DAS type
	try:
		das_type_idx = findstr_in_list(xls_title, das_type)[0]
	except:
		raise Exception('Cannot found DAS type in excel title')

	# get position of DAS gain
	try:
		das_gain_idx = findstr_in_list(xls_title, das_gain)[0]
	except:
		raise Exception('Cannot found Gain in excel title')

	# get position of DAS sample rate
	try:
		das_sampling_rate_idx = findstr_in_list(xls_title, das_sampling_rate)[0]
	except:
		raise Exception('Cannot found Sample rate in excel title')


	colnames = {network: xls_title[network_idx], \
		        station_name: xls_title[station_name_idx], \
		        starttime: xls_title[starttime_idx], \
		        endtime: xls_title[endtime_idx], \
		        latitude: xls_title[latitude_idx], \
		        longitude: xls_title[longitude_idx], \
		        elevation: xls_title[elevation_idx], \
		        sensor_type: xls_title[sensor_type_idx], \
		        sensor_depth: xls_title[sensor_depth_idx], \
		        das_type: xls_title[das_type_idx], \
		        das_gain: xls_title[das_gain_idx], \
		        das_sampling_rate: xls_title[das_sampling_rate_idx]}


	return colnames



def index_sensor_resp_pool(resp_list, sensor_type_in):

	iresp = []

	for i in range(len(resp_list.filename)):

		sensor_type = resp_list.sensor_type[i]
		sensor_period = resp_list.sensor_period[i]

		res1 = findstr(sensor_type_in, sensor_type)
		if (res1 != []):

			for ipos1 in res1:

				if ((sensor_period == '') or (sensor_period == [])):
					iresp = i
					break
				else:
					res2 = findstr(sensor_type_in, sensor_period)
					if (res2 != []):
						iresp = i

	return iresp



def index_das_resp_pool(resp_list, das_type_in, das_gain_in, das_sampling_rate_in):

	iresp = []

	ier = 0
	for i in range(len(resp_list.filename)):

		das_type = resp_list.das_type[i]
		das_gain = resp_list.das_gain[i]
		das_sampling_rate = resp_list.das_sampling_rate[i]

		res1 = findstr(str(das_type_in), str(das_type))
		if (res1 != []):

			for i in range(len(res1)):

				if (abs(float(das_gain_in) - float(das_gain)) < 1.e-9):

					if (abs(float(das_sampling_rate_in) - float(das_sampling_rate)) < 1.e-9):
						iresp = i
						break
					else:
						ier = -1

		else:
			ier = -1

		if (-1 == ier):
			print("Warning: Cannot find DAS of %s with gain %s and sampling rate %s " \
						     % das_type_in, das_gain_in, das_sampling_rate_in)
			print("Warning: DAS's instrument response will be ignored !")

	return iresp



def merge_paz(sensor_resp_pool, iresp_sensor, das_resp_pool, iresp_das):

	try:
		sensor_PZ = sensor_resp_pool[iresp_sensor]
	except:
		print("Error: Sensor's instrument response is mandatory ")
		return


	if ([] == iresp_das):

		return sensor_PZ

	else:

		das_PZ = das_resp_pool[iresp_das]
		poles = sensor_PZ.get('poles')
		das_poles = das_PZ.get('poles')
		if ([] != das_poles):
			poles += das_poles

		zeros = sensor_PZ.get('zeros')
		das_zeros = das_PZ.get('zeros')
		if ([] != das_zeros):
			zeros += das_zeros

		A0_factor = sensor_PZ.get('A0') * das_PZ.get('A0')
		sensitivity = sensor_PZ.get('sensitivity') * das_PZ.get('sensitivity')
		gain_all = sensor_PZ.get('gain') * das_PZ.get('gain')
		CONSTANT = sensor_PZ.get('CONSTANT') * das_PZ.get('CONSTANT')

		PZ = { \
			'poles': poles, \
			'zeros': zeros, \
			'A0': A0_factor, \
			'sensitivity': sensitivity, \
			'gain': gain_all, \
			'CONSTANT': CONSTANT}

	return PZ



def write_paz(f, PZ, network, station_name, C, starttime, endtime, \
	                        stla, stlo, stel, stdep, sampling_rate):

	f.write('* ***************************************************************\n')
	f.write('* NETWORK       : %s \n' % network)
	f.write('* STATION       : %s \n' % station_name)
	f.write('* LOCATION      : %s \n' % '')
	f.write('* CHANNEL       : %s \n' % C)
	f.write('* CREATED       : %s \n' % time2str(UTCDateTime()))
	f.write('* START         : %s \n' % reformat_time(starttime))
	f.write('* END           : %s \n' % reformat_time(endtime))
	f.write('* DESCRIPTION   : %s \n' % 'YSLIU')
	f.write('* LATITUDE      : %s \n' % stla)
	f.write('* LONGITUDE     : %s \n' % stlo)
	f.write('* ELEVATION     : %s \n' % stel)
	f.write('* DEPTH         : %f \n' % abs(float(stdep)))
	if ([] != findstr(C, 'Z')):
		f.write('* DIP (SEED)    : %f \n' % -90.0)
		f.write('* AZIMUTH       : %s \n' % 0.0)
	if ([] != findstr(C, 'N')):
		f.write('* DIP (SEED)    : %f \n' % 0.0)
		f.write('* AZIMUTH       : %s \n' % 0.0)
	if ([] != findstr(C, 'E')):
		f.write('* DIP (SEED)    : %f \n' % 0.0)
		f.write('* AZIMUTH       : %s \n' % 90.0)
	f.write('* SAMPLE RATE   : %f \n' % sampling_rate)
	f.write('* INPUT UNIT    : %s \n' % 'M')
	f.write('* OUTPUT UNIT   : %s \n' % 'COUNTS')
	f.write('* INSTTYPE      : %s \n' % 'None')
	f.write('* INSTGAIN      : %e (M/S)\n' % PZ.get('gain'))
	f.write('* SENSITIVITY   : %e (M/S)\n' % PZ.get('sensitivity'))
	f.write('* A0            : %f \n' % PZ.get('A0'))
	f.write('* ***************************************************************\n')
	zeros = PZ.get('zeros')
	nzeros = len(zeros)
	f.write('ZEROS %d \n' % nzeros)
	for i in range(0, nzeros):
		f.write('%f %f \n' % (zeros[i].real, zeros[i].imag))
	poles = PZ.get('poles')
	npoles = len(poles)
	f.write('POLES %d \n' % npoles)
	for i in range(0, npoles):
		f.write('%f %f \n' % (poles[i].real, poles[i].imag))
	f.write('CONSTANT %e \n' % PZ.get('CONSTANT'))

	return



def resp2pz(sensor_resp_pool_path, sensor_resp_listname, \
            das_resp_pool_path, das_resp_listname):

	# read resp pool of sensor
	sensor_resp_list = read_sensor_resp_list(sensor_resp_listname)
	sensor_resp_pool = read_resp_pool(sensor_resp_pool_path, sensor_resp_list)

	# read resp pool of DAS
	das_resp_list = read_das_resp_list(das_resp_listname)
	das_resp_pool = read_resp_pool(das_resp_pool_path, das_resp_list)



	# define marks
	network_key = 'Network'
	station_name_key = 'Station code'
	starttime_key = 'Starttime'
	endtime_key = 'Endtime'
	latitude_key = 'Latitude'
	longitude_key = 'Longitude'
	elevation_key = 'Elevation'
	sensor_type_key = 'Sensor type'
	sensor_depth_key = 'sensor depth'
	das_type_key = 'DAS type'
	das_sampling_rate_key = 'Sample rate'
	das_gain_key = 'Gain'


	colnames = find_colname_from_xls(station_logfile, network_key, station_name_key, \
				starttime_key, endtime_key, latitude_key, longitude_key, \
				elevation_key, sensor_type_key, sensor_depth_key, das_type_key, \
				das_gain_key, das_sampling_rate_key)


	# get column names for marks
	colname_of_network = colnames.get(network_key)
	colname_of_station_name = colnames.get(station_name_key)
	colname_of_starttime = colnames.get(starttime_key)
	colname_of_endtime = colnames.get(endtime_key)
	colname_of_latitude = colnames.get(latitude_key)
	colname_of_longitude = colnames.get(longitude_key)
	colname_of_elevation = colnames.get(elevation_key)
	colname_of_sensor_type = colnames.get(sensor_type_key)
	colname_of_sensor_depth = colnames.get(sensor_depth_key)
	colname_of_das_type = colnames.get(das_type_key)
	colname_of_das_gain = colnames.get(das_gain_key)
	colname_of_das_sampling_rate = colnames.get(das_sampling_rate_key)


	# read excel title
	tables, nrows, _ = read_xls_byindex(station_logfile)


	# create output path if it does not exist
	if (not os.path.exists(output_path)):
		os.makedirs(output_path)
	else:
		for file in os.listdir(output_path):
			filename = output_path + '/' + file
			if (os.path.exists(filename)):
				os.remove(filename)


	station_name_prev = ''
	for i in range(nrows):


		# get current station name
		try:
			station_name = tables[i].get(colname_of_station_name)
		except:
			continue


		iresp_sensor = []
		# find the index of the RESP file in the PZ pool
		if ((station_name != '') and (station_name != [])):

			sensor_type = tables[i].get(colname_of_sensor_type)

			# here, iresp_sensor is the index of the sensor resp file in the pool
			iresp_sensor = index_sensor_resp_pool(sensor_resp_list, sensor_type)


			das_type = tables[i].get(colname_of_das_type)
			das_gain = tables[i].get(colname_of_das_gain)
			das_sampling_rate = tables[i].get(colname_of_das_sampling_rate)

			# here, iresp_das is the index of the das resp file in the pool
			iresp_das = index_das_resp_pool(das_resp_list, das_type, das_gain, das_sampling_rate)


			PZ = merge_paz(sensor_resp_pool, iresp_sensor, das_resp_pool, iresp_das)


		# whether open a new file
		open_new_file = False
		if (station_name != station_name_prev):
			if (station_name_prev != '' or station_name == ''):
				f.close()
			open_new_file = True


		if ((station_name != '') and (station_name != [])):

			network = tables[i].get(colname_of_network)
			starttime = weehour_daily(tables[i].get(colname_of_starttime))
			endtime = tables[i].get(colname_of_endtime)
			stla = tables[i].get(colname_of_latitude)
			stlo = tables[i].get(colname_of_longitude)
			stel = tables[i].get(colname_of_elevation)
			stdep = tables[i].get(colname_of_sensor_depth)

			for C in channel_name:

				filename = output_path + '/' + \
					network + '.' + station_name + '.' + \
					'' + '.' + C + '.PZ'

				if (open_new_file):
					f = open(filename, 'w+')

				write_paz(f, PZ, network, station_name, C, starttime, endtime, \
					                 stla, stlo, stel, stdep, downsampling_rate)

		if (station_name != ''):
			print('station %s is done ... \n' % station_name)

		station_name_prev = station_name


	# now, close the last opened file object
	f.close()

	return



if __name__ == '__main__':

	print('\n')
	print('makePZs: ')
	print('This program automatically generate PZ files based on the RESP files.')
	print('Youshan Liu at Institute of Geology and Geophysics, Chinese Academy of Sciences')
	print('Wecome to send any bugs and suggestions to ysliu@mail.iggcas.ac.cn')
	print('Please downloaded the unavailable RESP files in RESP_pool from the following website.')
	print('http://ds.iris.edu/NRL/sensors/guralp/guralp_sensors.html')
	print('\n')

	resp2pz(sensor_resp_pool_path, sensor_resp_listname, das_resp_pool_path, das_resp_listname)

