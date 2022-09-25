#!/usr/bin/env python3
# -*- coding: UTF-8 -*-


'''
reftek2sac

This program convert files from reftek to sac format using the ObsPy

Date: 22/10/2020
Author: Youshan Liu
Affiliation: Institute of Geology and Geophysics, Chinese Academy of Sciences


folders structure:
./your_data_path/stage folder/station folder/Reftek UnitID number/day folder/stream/reftek files

for example:
./Raw/NE00_2007_276_2008_005/NE00/2007276/9F78/1



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
dryrun = True

# direction of the reftek data
data_path = './DATA_1/'


##############################################################
##############################################################
##############################################################
def reftek2sac(current_path):

	sac_suffix = '.SAC'

	rootdir = data_path + '/'
	#rootdir = current_path + '/' + data_path + '/'
	len_rootdir = len(rootdir)

	stage_folders_list = os.listdir(rootdir)

	#print(stage_folders_list)
	#return

	# convert reftek to sac
	for station_stage_folder in stage_folders_list:

		station_stage_path = rootdir + station_stage_folder + '/'
		print('Entering directory ' + station_stage_path[len_rootdir:-1])
		print('\n')

		if (not os.path.exists(station_stage_path)):
			continue

		day_folders_list = os.listdir(station_stage_path)

		#print(day_folders_list)
		#return

		for day_folder in day_folders_list:

			day_path = station_stage_path + day_folder + '/'
			print('\tEntering directory ' + day_path[len_rootdir:-1])
			print('\n')

			if (not os.path.exists(day_path)):
				continue

			# remove all sac files
			#print(day_path)
			os.system("rm -rf " + day_path + "./*" + sac_suffix)
			continue

			print('\tLeaving directory ' + day_path[len_rootdir:-1])
			print('\n')

		del day_folders_list
		print('Leaving directory ' + station_stage_path[len_rootdir:-1])
		print('\n')

	del stage_folders_list

	return



if __name__ == '__main__':

	print('\n')
	print('reftek2sac: ')
	print('This program convert files from reftek to sac format using the ObsPy (serial version)')
	print('Youshan Liu at Institute of Geology and Geophysics, Chinese Academy of Sciences')
	print('Welcome to send any bugs and suggestions to ysliu@mail.iggcas.ac.cn')
	print('\n\n')

	starttime = UTCDateTime()

	# get current path
	current_path = os.getcwd()

	# convert reftek file to sac format
	reftek2sac(current_path)

	endtime = UTCDateTime()

	elapsed_time = (endtime - starttime)

	print('Start   time : %s' % starttime)
	print('End     time : %s' % endtime)
	print('Elapsed time : %f hours \n' % (elapsed_time / 3600.0))



