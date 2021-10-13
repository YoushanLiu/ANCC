#!/usr/bin/env python3

#=======================================================================================
#=======================================================================================
# This script is the only one that you need to execute to process the seismic ambient
# noise data with parallel computing. It will invoke the <AND_Driver>, <TF_PWS> and
# <AFTAN_PROG> three executables to deal with the input seismic waveforms based on input
# parameters as defined in <input.dat>. Detailed description of each individual parameter
# can be found in <input.dat>. The output of this program is cross-correlation functions
# and corresponding group and phase velocity dispersion data.
#=======================================================================================
# The path of input data folder is given in the <SACfolder> variable in the following
# script. The data structure should be */year/month/event/net.sta.cha.SAC
# e.g. DATA/2018/01/20180101_120000/DB.EW01.LHZ.SAC
#=======================================================================================
# Individual or combined full RESP or PZ file(s) should be put in <PZfolder>.
#=======================================================================================
# 1-D phase velocity reference model used in AFTAN is given by <ref1Dmod.dat>. It has
# the format of (period, reference phase velocity)
#=======================================================================================
# The path of output data folder is given in the <tarfolder> variable in the following
# script. The ouptput 'CC_AFTAN' folder contains all the cross-correlation waveforms and
# corresponding dispersion measurements. The output 'FINAL' folder contains all the final
# dispersion data. The output 'BOOTSTRAP' folder contains all the bootstrap data (optional).
#=======================================================================================
# Required auxiliary files and programs:
#   sac
#   myTools.py
#   input.dat
#   ref1Dmod.dat
#=======================================================================================
# Notes: The algorithm of this program is modified from Barmin's code (http://ciei.color
#        ado.edu/Products/), including the cross-correlation part and the frequency-time
#        analysis part. The phase-weighted stacking program is modified from Guoliang Li'
#        code (guoleonlee@gmail.com)
#        Any bug reports, suggestions or comments related to this program can be directed
#        to Xingli Fan via <fanxldengcm@gmail.com>.
#=======================================================================================
#=======================================================================================

import os
import glob
import multiprocessing
from obspy.core import UTCDateTime
from multiprocessing.dummy import Pool as ThreadPool


###########################################################
# Please modify the following three path varialbes with
# your own case and set the right core numbers for parallel
# processing according to your own computing resource.
###########################################################

SACfolder = "../DATA_cut"
PZfolder  = "../PZs_all"
tarfolder = "./Result_10Hz"

nprocs = multiprocessing.cpu_count()


###########################################################

BIN_DIR = os.environ['HOME'] + '/bin'
# create links
if (not os.path.exists(os.getcwd() + '/AFTAN')):
	os.symlink(BIN_DIR + '/AFTAN', 'AFTAN')
else:
	os.remove('AFTAN')
	os.symlink(BIN_DIR + '/AFTAN', 'AFTAN')

if (not os.path.exists(os.getcwd() + '/AND')):
	os.symlink(BIN_DIR + '/AND', 'AND')
else:
	os.remove('AND')
	os.symlink(BIN_DIR + '/AND', 'AND')

if (not os.path.exists(os.getcwd() + '/TF_PWS')):
	os.symlink(BIN_DIR + '/TF_PWS', 'TF_PWS')
else:
	os.remove('TF_PWS')
	os.symlink(BIN_DIR + '/TF_PWS', 'TF_PWS')

###########################################################

# Reomve empty folder(s) and file(s)
###########################################################
os.system('find %s -depth -type "d" -empty -exec rmdir {} \;' %(SACfolder))
os.system('find %s -name "*" -type f -size 0c | xargs -n 1 rm -f' % (SACfolder))

# Retrive the station and event infomation.
###########################################################
os.system('rm -rf stations.junk events.lst')



#def scan_daily(day):
#    day_folder = month_folder + '/' + day
#    sacfiles = day_folder + '/*.SAC'
#    file_list = glob.glob(sacfiles)
#    if (len(file_list) <= 1):
#        print("skip %s because of single station folder\n"%(day_folder))
#        return
#    os.system("saclst knetwk kstnm stlo stla delta f %s | awk '{print $2,$3,$4,$5,$6}' >> stations.junk"%(sacfiles))
#    os.system("ls %s -d >> events.lst"%(day_folder))


#print("\n")
#for year in os.listdir(SACfolder):
#    year_folder = SACfolder + '/' + year
#    for month in os.listdir(year_folder):
#        month_folder = year_folder + '/' + month

#        day_folders_list = os.listdir(month_folder)

#        pool = ThreadPool()
#        pool.map(scan_daily, day_folders_list)
#        pool.close()
#        pool.join()

#os.system("sort stations.junk | uniq > stations.lst")
#os.system("rm -rf stations.junk")



print("\n")
for year in os.listdir(SACfolder):
    year_folder = SACfolder + '/' + year
    for month in os.listdir(year_folder):
        month_folder = year_folder + '/' + month
        for day in os.listdir(month_folder):
            day_folder = month_folder + '/' + day
            sacfiles = day_folder + '/*.SAC'
            file_list = glob.glob(sacfiles)
            #if (len(file_list) <= 1):
            #    print("skip %s because of single station folder\n"%(day_folder))
            #    continue
            os.system("saclst knetwk kstnm stlo stla delta f %s | awk '{print $2,$3,$4,$5,$6}' >> stations.junk"%(sacfiles))
            os.system("ls %s -d >> events.lst"%(day_folder))
os.system("sort stations.junk | uniq > stations.lst")
os.system("rm -rf stations.junk")


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


# Execute the main program to process the ambient noise data.
###########################################################
starttime = UTCDateTime()
os.system("mpirun -np %d ./AND %s %s %s"%(nprocs, SACfolder, PZfolder, tarfolder))
endtime = UTCDateTime()

elapsed_time = endtime - starttime


print("\n")
print('Start   time : %s' % (starttime + 28800))
print('End     time : %s' % (endtime + 28800))
print('Elapsed time : %f hours \n' % (elapsed_time / 3600.0))
