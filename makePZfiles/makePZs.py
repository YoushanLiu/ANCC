#!/usr/bin/env python3


'''

This program automatically generate PZ files based on the RESP files.

install xlrd:
pip3 install xlrd

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



# path of the PZ files pool
RESP_pool_path = './Sensor_RESP_pool'

# the used PZ files list
RESP_listname = 'PZs_used.lst'

# path of the generated PZ files
output_path = '../PZs_all'

# station log file
station_logfile = './NCISP6.xls'


# channel name, it consists of band code, instrument code, orientation code
channels = ['BHZ']

# the decimate rate must be interger so that the sampling_rate is divisible by it
decimate_rate = 4



def read_resp(filename):

	'''
	read_resp(filename)
	read standard RESP file
	filename      -> RESP file name
	convert2displ -> tansfer to displacement
	'''

	# define regular expression
	A0_regexp = re.compile(r"^[ \t]*[^#]+[\w \t]+(A0)+[ \t]*(normalization)+[ \t]*(factor)+:[ \t]*([+-]*[0-9]*\.[0-9]+[EeDd]+[+-]*[0-9]+)+")
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
	if (veloc2displ):
		zeros.append(0.0+0.0j)
	if (accel2displ):
		zeros.append(0.0+0.0j)
		zeros.append(0.0+0.0j)


	CONSTANT = float(A0_factor) * gain_all
	resp = { \
		'poles': poles, \
		'zeros': zeros, \
		'A0': A0_factor, \
		'sensitivity': sensitivity, \
		'gain': gain_all, \
		'CONSTANT': CONSTANT}


	return resp



def read_resp_list(filename):

	'''
	read_resp_list(filename)
	read RESP list file
	filename -> filename of a RESP file list
	'''

	class RESP_list:
		def __init__(self):
			self.name = []
			self.keyword1 = []
			self.keyword2 = []


	name = ''
	keyword1 = ''
	keyword2 = ''
	resp_list = RESP_list()
	with open(filename, 'r') as f:
		lines = f.readlines()
	for line in lines[1:]:
		try:
			line_splited = line.split()
			if ([] == line_splited):
				continue
			else:
				keyword2 = ''
				if (len(line_splited) >= 3):
					name, keyword1, keyword2 = line_splited
				else:
					name, keyword1 = line_splited
		except:
			raise Exception('Format error in %s ' % filename)

		resp_list.name.append(name)
		resp_list.keyword1.append(keyword1)
		resp_list.keyword2.append(keyword2)

	return resp_list



def read_resp_pool(resp_list):

	'''
	read_resp_pool(resp_pool)
	read RESP files pool
	resp_pool -> RESP pool filename
	'''

	RESP_pool = []
	for file in resp_list.name:
		filename = RESP_pool_path + '/' + file
		resp = read_resp(filename)
		RESP_pool.append(resp)

	return RESP_pool



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



# read excel
def find_colname_from_xls(station_logfile, network, station, starttime, endtime, \
			latitude, longitude, elevation, sensor_type, sensor_depth, sample_rate):

	# read logfile title
	xls_title = get_xls_title(station_logfile)


	# get position of network
	try:
		network_idx = findstr_in_list(xls_title, network)[0]
	except:
		raise Exception('Cannot found Network in excel title')

	# get position of station
	try:
		station_idx = findstr_in_list(xls_title, station)[0]
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

	# get position of sample rate
	try:
		sample_rate_idx = findstr_in_list(xls_title, sample_rate)[0]
	except:
		raise Exception('Cannot found Sample rate in excel title')


	colnames = {network: xls_title[network_idx], \
		station: xls_title[station_idx], \
		starttime: xls_title[starttime_idx], \
		endtime: xls_title[endtime_idx], \
		latitude: xls_title[latitude_idx], \
		longitude: xls_title[longitude_idx], \
		elevation: xls_title[elevation_idx], \
		sensor_type: xls_title[sensor_type_idx], \
		sensor_depth: xls_title[sensor_depth_idx], \
		sample_rate: xls_title[sample_rate_idx]}


	return colnames



def index_resp_pool(RESP_list, sensor_type):

	iresp = []

	for i in range(len(RESP_list.name)):

		keyword1 = RESP_list.keyword1[i]
		keyword2 = RESP_list.keyword2[i]
	
		res1 = findstr(sensor_type, keyword1)
		if (res1 != []):

			for ipos1 in res1:

				if ((keyword2 == '') or (keyword2 == [])):
					iresp = i
					break
				else:
					res2 = findstr(sensor_type, keyword2)
					if (res2 != []):
						iresp = i

	return iresp



def weehour_daily(starttime):
	time_tmp = UTCDateTime(starttime)
	time = UTCDateTime(time_tmp.year, time_tmp.month, time_tmp.day, 0, 0, 0, 000000)
	starttime = time2str(time)

	return starttime



def write_paz(f, RESP_pool, iresp, network, station, C, \
		starttime, endtime, stla, stlo, stel, stdep, sample_rate):

	f.write('* **************************************************\n')
	f.write('* NETWORK       : %s\n' % network)
	f.write('* STATION       : %s\n' % station)
	f.write('* LOCATION      : %s\n' % '')
	f.write('* CHANNEL       : %s\n' % C)
	f.write('* CREATED       : %s\n' % time2str(UTCDateTime()))
	f.write('* START         : %s\n' % reformat_time(starttime))
	f.write('* END           : %s\n' % reformat_time(endtime))
	f.write('* DESCRIPTION   : %s.%s \n' % (network, station))
	f.write('* LATITUDE      : %s\n' % stla)
	f.write('* LONGITUDE     : %s\n' % stlo)
	f.write('* ELEVATION     : %s\n' % stel)
	f.write('* DEPTH         : %f\n' % abs(float(stdep)))
	if ('Z' == C[-1].upper()):
		f.write('* DIP (SEED)    : %f\n' % -90.0)
		f.write('* AZIMUTH       : %s\n' % 0.0)
	if ('N' == C[-1].upper()):
		f.write('* DIP (SEED)    : %f\n' % 0.0)
		f.write('* AZIMUTH       : %s\n' % 0.0)
	if ('E' == C[-1].upper()):
		f.write('* DIP (SEED)    : %f\n' % 0.0)
		f.write('* AZIMUTH       : %s\n' % 90.0)
	f.write('* SAMPLE RATE   : %f\n' % sample_rate)
	f.write('* INPUT UNIT    : %s\n' % 'M')
	f.write('* OUTPUT UNIT   : %s\n' % 'COUNTS')
	f.write('* INSTTYPE      : %s\n' % 'None')
	f.write('* INSTGAIN      : %e (M/S)\n' % RESP_pool[iresp].get('gain'))
	f.write('* SENSITIVITY   : %e (M/S)\n' % RESP_pool[iresp].get('sensitivity'))
	f.write('* A0            : %f\n' % RESP_pool[iresp].get('A0'))
	f.write('* **************************************************\n')
	zeros = RESP_pool[iresp].get('zeros')
	nzeros = len(zeros)
	f.write('ZEROS %d \n' % nzeros)
	for i in range(0, nzeros):
		f.write('%+f %+f \n' % (zeros[i].real, zeros[i].imag))
	poles = RESP_pool[iresp].get('poles')
	npoles = len(poles)
	f.write('POLES %d \n' % npoles)
	for i in range(0, npoles):
		f.write('%+f %+f \n' % (poles[i].real, poles[i].imag))
	f.write('CONSTANT %e \n' % RESP_pool[iresp].get('CONSTANT'))

	return



def resp2pz(RESP_pool_path, RESP_listname):

	RESP_list = read_resp_list(RESP_pool_path + '/' + RESP_listname)


	RESP_pool = read_resp_pool(RESP_list)


	# define marks
	network_mark = 'Network'
	station_mark = 'Station code'
	starttime_mark = 'Starttime'
	endtime_mark = 'Endtime'
	latitude_mark = 'Latitude'
	longitude_mark = 'Longitude'
	elevation_mark = 'Elevation'
	sensor_type_mark = 'Sensor type'
	sensor_depth_mark = 'sensor depth'
	sample_rate_mark = 'Sample rate'


	colnames = find_colname_from_xls(station_logfile, network_mark, station_mark, \
				starttime_mark, endtime_mark, latitude_mark, longitude_mark, \
				elevation_mark, sensor_type_mark, sensor_depth_mark, sample_rate_mark)


	# get column names for marks
	colname_of_network = colnames.get(network_mark)
	colname_of_station = colnames.get(station_mark)
	colname_of_starttime = colnames.get(starttime_mark)
	colname_of_endtime = colnames.get(endtime_mark)
	colname_of_latitude = colnames.get(latitude_mark)
	colname_of_longitude = colnames.get(longitude_mark)
	colname_of_elevation = colnames.get(elevation_mark)
	colname_of_sensor_type = colnames.get(sensor_type_mark)
	colname_of_sensor_depth = colnames.get(sensor_depth_mark)
	colname_of_sample_rate = colnames.get(sample_rate_mark)


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


	station_prev = ''
	for i in range(nrows):

		# get current station name
		try:
			station = tables[i].get(colname_of_station)
		except:
			continue

		iresp = []
		# find the index of the RESP file in the PZ pool
		if ((station != '') and (station != [])):

			sensor_type = tables[i].get(colname_of_sensor_type)

			# here, iresp is the index of the RESP file in the RESP_pool
			iresp = index_resp_pool(RESP_list, sensor_type)

		# whether open a new file
		open_new_file = False
		if (station != station_prev):
			if (station_prev != '' or station == ''):
				f.close()
			open_new_file = True


		if ((station != '') and (station != [])):

			network = tables[i].get(colname_of_network)
			starttime = weehour_daily(tables[i].get(colname_of_starttime))
			endtime = tables[i].get(colname_of_endtime)
			stla = tables[i].get(colname_of_latitude)
			stlo = tables[i].get(colname_of_longitude)
			stel = tables[i].get(colname_of_elevation)
			stdep = tables[i].get(colname_of_sensor_depth)
			sample_rate = float(tables[i].get(colname_of_sample_rate))/decimate_rate

			for C in channels:

				filename = output_path + '/' + \
					network + '.' + station + '.' + \
					'' + '.' + C + '.PZ'

				if (open_new_file):
					f = open(filename, 'w+')

				write_paz(f, RESP_pool, iresp, network, station, C, \
					starttime, endtime, stla, stlo, stel, stdep, sample_rate)

		if (station != ''):
			print('station %s is done ...' % station)

		station_prev = station


	# now, close the last opened file object
	f.close()

	return



if __name__ == '__main__':

	#print('\n')
	print('makePZs: ')
	print('This program automatically generates PZ files based on the RESP files.')
	print('Youshan Liu at Institute of Geology and Geophysics, Chinese Academy of Sciences')
	print('Wecome to send any bugs and suggestions to ysliu@mail.iggcas.ac.cn')
	print('Please downloaded the unavailable RESP files in RESP_pool from the following website.')
	print('http://ds.iris.edu/NRL/sensors/guralp/guralp_sensors.html')
	#print('\n')

	resp2pz(RESP_pool_path, RESP_listname)

