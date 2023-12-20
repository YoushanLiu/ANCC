#!/usr/bin/env python3

#=======================================================================================
#=======================================================================================
# This script is the only one that you need to execute to process the seismic ambient
# noise data with parallel computing. It will invoke the <ANCC>, <TF_PWS> and
# <AFTAN> three executables to deal with the input seismic waveforms based on input
# parameters as defined in <input.dat>. Detailed description of each individual parameter
# can be found in <input.dat>. The output of this program is auto- or cross-correlation 
# functions or corresponding group and phase velocity dispersion data.
#=======================================================================================
# The PATH of input data folder is given in the <SACfolder> variable in the following
# script. The data structure should be */year/month/event/net.sta.cha.SAC
# e.g. DATA/2018/01/20180101_120000/NE.NE01.BHZ.SAC
#=======================================================================================
# Individual or combined full RESP or PZ file(s) should be put in <PZfolder>.
#=======================================================================================
# 1-D phase velocity reference model used in AFTAN is given by <ref1Dmod.dat>. It has
# the format of (period, reference phase velocity)
#=======================================================================================
# The PATH of output data folder is given in the <tarfolder> variable in the following
# script. The ouptput 'CC_AFTAN' folder contains all the cross-correlation waveforms and
# corresponding dispersion measurements. The output 'FINAL' folder contains all the final
# dispersion data. The output 'BOOTSTRAP' folder contains all the bootstrap data (optional).
#=======================================================================================
# Required auxiliary files and programs:
#   SAC
#   input.dat
#   ref1Dmod.dat
#=======================================================================================
# Notes: The algorithm of this program is modified from Barmin's code (http://ciei.color
#        ado.edu/Products/), including the cross-correlation part and the frequency-time
#        analysis part. The phase-weighted stacking program is modified from Guoliang Li'
#        code (guoleonlee@gmail.com)
#        Any bug reports, suggestions or comments related to this program can be directed
#        to Xingli Fan via <fanxldengcm@gmail.com> and <ysliu@mail.iggcas.ac.cn>.
#=======================================================================================
#=======================================================================================

import os
import re
import glob
import linecache
import multiprocessing
from obspy import UTCDateTime
from multiprocessing.dummy import Pool as ThreadPool


###########################################################
# Please modify the following three PATH varialbes with
# your own case and set the right core numbers for parallel
# processing according to your own computing resource.
###########################################################

SACfolder = "../DATA_cut"
PZfolder  = "../PZs_all"
tarfolder = "./Results_10Hz"

nprocs = multiprocessing.cpu_count()


###########################################################
def count_station_num(files_list):
	strs = files_list
	for i in range(0,len(files_list)):
		try:
			strs[i] = files_list[i].split('/')[-1].split('.')[-3]
		except:
			strs[i] = ''
	return len({}.fromkeys(strs).keys())


def check_autocorrelation(filename):
	try:
		line = linecache.getline(filename, 46)
	except:
		return False
	regexp = re.compile("r[^#]+[ \t\w{,3}]+")
	ans = regexp.search(line[0:8])
	if (ans is None):
		return 0
	if (ans.strip().upper() in "YES"):
		is_auto_correlation = True
	else:
		is_auto_correlation = False
	return is_auto_correlation


is_auto_correlation = check_autocorrelation("input.dat")
###########################################################
# Reomve empty folder(s) and file(s)
###########################################################
#os.system('find %s -depth -type "d" -empty -exec rmdir {} \;' %(SACfolder))
#os.system('find %s -name "*" -type f -size 0c | xargs -n 1 rm -f' % (SACfolder))

###########################################################
# Retrieve the station and event information.
###########################################################
os.system('rm -rf stations.tmp events.lst')


#awkstr = " | awk '{printf \"%8s  %8s  %8.4f  %8.4f  %5.3f\",$2,$3,$4,$5,$6}' >> stations.tmp"
awkstr = " | awk '{printf " + r'"%-8s  %-8s  %7.4f  %8.4f  %5.3f\n",$2,$3,$4,$5,$6}' + "' >> stations.tmp"
global month_folder
def scan_segment(segment):
    segment_folder = month_folder + '/' + segment
    sacfiles = segment_folder + '/*.SAC'
    file_list = glob.glob(sacfiles)
    if (not(is_auto_correlation) and (count_station_num(files_list) <= 1)):
        print("skip %s because of single station folder\n"%(segment_folder))
        return
    #os.system("saclst knetwk kstnm stla stlo delta f %s | awk '{print $2,$3,$4,$5,$6}' >> stations.tmp"%(sacfiles))
    os.system("saclst knetwk kstnm stla stlo delta f %s"%(sacfiles + awkstr))
    os.system("ls %s -d >> events.lst"%(segment_folder))


for year in os.listdir(SACfolder):
    year_folder = SACfolder + '/' + year
    for month in os.listdir(year_folder):
        month_folder = year_folder + '/' + month
        pool = ThreadPool()
        pool.map(scan_segment, os.listdir(month_folder))
        pool.close()
        pool.join()


##awkstr = " | awk '{printf \"%8s  %8s  %8.4f  %8.4f  %5.3f\",$2,$3,$4,$5,$6}' >> stations.tmp"
#awkstr = " | awk '{printf " + r'"%-8s  %-8s  %7.4f  %8.4f  %5.3f\n",$2,$3,$4,$5,$6}' + "' >> stations.tmp"
#for year in os.listdir(SACfolder):
#    year_folder = SACfolder + '/' + year
#    for month in os.listdir(year_folder):
#        month_folder = year_folder + '/' + month
#        for segment in os.listdir(month_folder):
#            segment_folder = month_folder + '/' + day
#            sacfiles = segment_folder + '/*.SAC'
#            files_list = glob.glob(sacfiles)
#            if (not(is_auto_correlation) and (count_station_num(files_list) <= 1)):
#                print("skip %s because of single station folder\n"%(segment_folder))
#                continue
#            #os.system("saclst knetwk kstnm stla stlo delta f %s | awk '{print $2,$3,$4,$5,$6}' >> stations.tmp"%(sacfiles))
#            os.system("saclst knetwk kstnm stla stlo delta f %s" % (sacfiles + awkstr))
#            os.system("ls %s -d >> events.lst"%(segment_folder))


os.system("sort stations.tmp | uniq > stations.lst")
os.system("rm -rf stations.tmp")

###########################################################
# Return if stations.lst contains duplicate stations
###########################################################
stations = []
duplicates = []
for line in open('stations.lst', 'r'):
    net, sta = line.split()[0:2]
    name = net+'.'+sta
    if name in stations:
        duplicates.append(name)
    else:
        stations.append(name)

if len(duplicates) > 0:
    print('Error: The following stations have multiple longitude or latitude or sampling rate.\n\
           Please check stations.lst for details and correct the corresponding sac \n\
           headers before running this program.')
    print(duplicates)
    exit()

###########################################################
# Execute the main program to process the ambient noise data.
###########################################################
starttime = UTCDateTime()
os.system("mpirun -np %d ANCC %s %s %s"%(nprocs, SACfolder, PZfolder, tarfolder))
endtime = UTCDateTime()

elapsed_time = endtime - starttime


print("\n")
print('Start   time : %s' % (starttime + 28800))
print('End     time : %s' % (endtime + 28800))
print('Elapsed time : %f hours \n' % (elapsed_time / 3600.0))
