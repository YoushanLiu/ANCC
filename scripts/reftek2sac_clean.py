#!/usr/bin/env python3
# -*- coding: UTF-8 -*-


'''
reftek2sac

This program converts files from REFTEK to SAC format using the ObsPy

Date: 22/10/2020
Author: Youshan Liu
Affiliation: Institute of Geology and Geophysics, Chinese Academy of Sciences


directory structure:
./top folder/station and stage folder/day folder/Reftek UnitID number/stream/reftek files

for example:
./DATA_Raw/NE00_2007_276_2008_005/2007276/9F78/1



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
from obspy import UTCDateTime
from obspy.io.sac import SACTrace



##############################################################
##############################################################
##############################################################
# some options to be set


# dryrun just for debug
# dryrun = True  => just print some directory information not to write sac files
# dryrun = False => do the actual files conversion
dryrun = True

# directory of the reftek data
data_path = './DATA_1'


##############################################################
##############################################################
##############################################################
def reftek2sac():

	sac_suffix = '.SAC'

	len_topdir = len(input_path) + 1

	stage_folders_list = os.listdir(input_path)

	# convert reftek to sac
	for station_stage_folder in stage_folders_list:

		station_stage_path = input_path + '/' + station_stage_folder + '/'
		print('Entering directory ' + station_stage_path[len_topdir:-1])
		#print('\n')

		if (not os.path.exists(station_stage_path)):
			continue

		day_folders_list = os.listdir(station_stage_path)

		#print(day_folders_list)
		#return

		for day_folder in day_folders_list:

			day_path = station_stage_path + day_folder + '/'
			print('Entering directory ' + day_path[len_topdir:-1])
			#print('\n')

			if (not os.path.exists(day_path)):
				continue

			# remove all sac files
			#print(day_path)
			os.system("rm -rf " + day_path + "*" + sac_suffix)
			continue

			print('Leaving directory ' + day_path[len_topdir:-1])
			#print('\n')

		del day_folders_list
		print('Leaving directory ' + station_stage_path[len_topdir:-1])
		#print('\n')

	del stage_folders_list

	return



if __name__ == '__main__':

	print('\n')
	print('reftek2sac_clean: ')
	print('This program convert files from reftek to sac format using the ObsPy (serial version)')
	print('Youshan Liu at Institute of Geology and Geophysics, Chinese Academy of Sciences')
	print('Welcome to send any bugs and suggestions to ysliu@mail.iggcas.ac.cn')
	print('\n')

	starttime = UTCDateTime()

	# get current path
	#current_path = os.getcwd()

	# convert reftek file to sac format
	reftek2sac()

	endtime = UTCDateTime()

	elapsed_time = (endtime - starttime)

	print('\n')
	print('Start   time : %s' % starttime)
	print('End     time : %s' % endtime)
	print('Elapsed time : %f hours' % (elapsed_time / 3600.0))



