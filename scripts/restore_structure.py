#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

'''
reftek2sac

This program converts the files from reftek to sac format using the ObsPy

Date: 22/10/2020
Author: Youshan Liu
Affiliation：Institute of Geology and Geophysics, Chinese Academy of Sciences


folders structure:
./your_data_folder/period folder/station folder/Reftek UnitID number/day folder/stream/reftek files

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
data_folder = 'DATA_Raw'


##############################################################
##############################################################
##############################################################
def restore_daily(day_folder):

	day_path = station_path + day_folder + '/'
	print('\t\tEntering directory ' + day_path[nrootdir:-1])
	#print('\n')

	if (not os.path.isdir(day_path)):
		return

	if (not os.path.isdir(day_path + 'raw')):
		return

	UnitID_folders_list = os.listdir(day_path + 'raw')
	if ([] == UnitID_folders_list):
		return

	#print("rm -rf " + day_path + 'raw')
	#os.system("rm -rf " + day_path + 'raw')
	for UnitID in UnitID_folders_list:
		print("mv " + day_path + 'raw' + '/' + UnitID + ' ' + day_path)
		if (os.path.exists(day_path + 'raw' + '/' + UnitID)):
			os.system("mv " + day_path + 'raw' + '/' + UnitID + ' ' + day_path)
		print("rm -rf " + day_path + 'raw')
		os.system("rm -rf " + day_path + 'raw')
		os.system("rm -rf " + day_path + "R???.??")
		os.system("rm -rf " + day_path + "*.ref")

	del UnitID_folders_list
	print('\t\tLeaving directory ' + day_path[nrootdir:-1])
	#print('\n')

	return



def miniseed2reftek(current_path):

	global nrootdir, station_path, sta, sac_suffix

	rootdir = current_path + '/' + data_folder + '/'
	nrootdir = len(rootdir)

	period_folders_list = os.listdir(rootdir)


	sac_suffix = '.SAC'

	# convert reftek to sac
	for period_folder in period_folders_list:

		period_path = rootdir + period_folder + '/'
		print('Entering directory ' + period_path[nrootdir:-1])
		print('\n')

		if (not os.path.isdir(period_path)):
			continue

		station_folder_list = os.listdir(period_path)

		for station_folder in station_folder_list:

			station_path = period_path + station_folder + '/'
			print('Entering directory ' + station_path[nrootdir:-1])
			print('\n')

			if (not os.path.isdir(station_path)):
				continue

			day_folders_list = os.listdir(station_path)

			pool = ThreadPool()
			pool.map(restore_daily, day_folders_list)
			pool.close()
			pool.join()

			del day_folders_list
			print('Leaving directory ' + station_path[nrootdir:-1])
			print('\n')

		del station_folder_list
		print('Leaving directory ' + period_path[nrootdir:-1])
		print('\n')

	del period_folders_list

	return



if __name__ == '__main__':

	print('\n')
	print('reftek2sac: ')
	print('This program restore the files from miniseed to reftek format using the ObsPy (parallel version)')
	print('Youshan Liu at Institute of Geology and Geophysics, Chinese Academy of Sciences')
	print('Welcome to send any bugs and suggestions to ysliu@mail.iggcas.ac.cn')
	print('\n\n')

	starttime = UTCDateTime()

	# get current path
	# absolution path
	#current_path = os.getcwd()
	# relative path
	current_path = '.'

	# convert reftek file to sac format
	miniseed2reftek(current_path)

	endtime = UTCDateTime()

	elapsed_time = (endtime - starttime)

	print('Start   time : %s' % starttime)
	print('End     time : %s' % endtime)
	print('Elapsed time : %f hours \n' % (elapsed_time / 3600.0))



