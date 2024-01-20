#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import os
import math
import numpy as np
from obspy import read
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from obspy import UTCDateTime
from obspy.io.sac import SACTrace
from obspy.geodetics import locations2degrees, degrees2kilometers
from multiprocessing.dummy import Pool as ThreadPool


# group velocity range [unit in m/s]
#vmin = 2700.0
vmin = 2900.0
vmax = 3450.0
# period range (unit in Hz.)
fmin = 1.0 / 25.0
fmax = 1.0
# window factor
win_factor_minus = 2.0
#win_factor_plus  = 2.0
win_factor_plus  = 2.0
# noise window factor
win_noise_factor = 4.0



def get_average(records):
	return sum(records) / len(records)


def get_variance(records):
	average = get_average(records)
	return sum([(x - average) ** 2 for x in records]) / len(records)


def get_standard_deviation(records):
	variance = get_variance(records)
	return math.sqrt(variance)


def get_rms(records):
	return math.sqrt(sum([x ** 2 for x in records]) / len(records))


def measure_snr(tr, tmin, tmax, tend, dt):

	n = len(tr)
	n0 = int(n/2)
	# index of signal and noise windows
	it1 = max(round(tmin/dt), 0)
	it2 = min(round(tmax/dt), n)
	it3 = min(round(tend/dt), n)
	#energy_signal = 0.50*(get_rms(tr[n0-it2:n0-it1]) + get_rms(tr[n0+it1:n0+it2]))
	#energy_noise = 0.50*(get_rms(tr[n0-it1+1:n0]) + get_rms(tr[n0:n0+it1-1]))
	#if (it3 < n):
	#	energy_noise += 0.50*(get_rms[n0-it3+1:n0] + get_rms(tr[n0+it2+1:n0+it3]))
	#tr_signal = []
	#tr_signal = 0.50(tr[n0+it1+1:n0+it2] + tr[n0-it1-1:n0-it2:-1])
	#tr_signal.append(tr[n0-it2:n0-it1+1])
	#tr_noise = []
	# precursory noise window
	#tr_noise = 0.50*(tr[n0:n0+it1] + tr[n0:n0-it1:-1])
	# trailing noise window
	#if ((n0 + it3) < n):
	#	tr_noise.append(0.50*(tr[n0+it2+1:n0+it3] + tr[n0-it2-1:n0-it3:-1]))

	# symmetric signal
	#print(n)
	#print(n0)
	#print(len(tr[n0:]))
	#print(len(tr[n0::-1]))
	tr_sym = 0.50*(tr[n0:] + tr[n0::-1])

	# signal window
	tr_signal = tr_sym[it1:it2]

	tr_noise = []
#	print(tr_sym[it2:it3])
#	tr_noise2 = tr_sym[it2:it3]
#	tr_noise.extend(tr_noise2)
#	print(tr_noise)
	# precursory noise window
	tr_noise.extend(tr_sym[0:it1])
#	tr_noise += tr_sym[0:it1]
#	print(tr_sym[0:it1])
#	print(tr_noise)
#	if (it1 >= 1):
#		tr_noise = tr_sym[0:it1]
	# trailing noise window
	tr_noise.extend(tr_sym[it2:it3])
#	print(tr_sym[it2:it3])
#	tr_noise += tr_sym[it2:it3]
#	print(tr_noise)
#	if (it2 < n0):
#		#tr_noise += tr_sym[it2:it3]
#		if ([] == tr_noise):
#			tr_noise = tr_sym[it2:it3]
#		else:
#			tr_noise += tr_sym[it2:it3]

	print("it1 = %d" % it1)
	print("it2 = %d" % it2)
	print("it3 = %d" % it3)
	print(len(tr_signal))
	print(len(tr_noise))
	print(get_rms(tr_signal))
	print(get_rms(tr_noise))

	rms_signal = get_rms(tr_signal)
	rms_noise = get_rms(tr_noise)

	del tr_signal, tr_noise, tr_sym

	return 20.0*math.log10(rms_signal / rms_noise)
	#return get_rms(tr_signal) / get_rms(tr_noise)



#file_path = './selective_stack_test/prestack'
file_path = './selective_stack_test/prestack_NE00_NE05'
file_list = os.listdir(file_path)


#print(file_list)
st = Stream()
for file in file_list:
	try:
		st_tmp = read(file_path + '/' + file)
	except:
		continue
	#print(file)
	st += st_tmp

#print(len(st))
#st.plot()
#st.show()
npts = st[0].stats.npts
dt = st[0].stats.delta
Tmax = (0.50*(npts+1) - 1)*dt
stla = st[0].stats.sac.stla
stlo = st[0].stats.sac.stlo
evla = st[0].stats.sac.evla
evlo = st[0].stats.sac.evlo

print("evla = %f, stla = %f" % (st[0].stats.sac.evla, st[0].stats.sac.stla))
print("evlo = %f, stlo = %f" % (st[0].stats.sac.evlo, st[0].stats.sac.stlo))


# linear stack
nstack = len(st)
tr_ls = np.zeros(npts)
for i in range(nstack):
	tr_ls += st[i].data


tr = st[0].copy()
sac = SACTrace.from_obspy_trace(tr)
tr.data = tr_ls
tr.write('./selective_stack_test/tr_ls.SAC', format='sac')
# compute distance
#dist = degrees2kilometers(locations2degrees(st[0].stats.sac.evla, st[0].stats.sac.evlo, \
#                                            st[0].stats.sac.stla, st[0].stats.sac.stlo)) * 1.e3
dist = degrees2kilometers(locations2degrees(evla, evlo, stla, stlo)) * 1.e3
tmin = max(dist / vmax - win_factor_minus/fmin, 0.0)
tmax = min(dist / vmin + win_factor_plus/fmin, Tmax)
tend = min(tmax + win_noise_factor/fmin, Tmax)

print("evla = %f, stla = %f" % (st[0].stats.sac.evla, st[0].stats.sac.stla))
print("evlo = %f, stlo = %f" % (st[0].stats.sac.evlo, st[0].stats.sac.stlo))
print("dist / vmax = %f" % (dist/vmax))
print("dist / vmin = %f" % (dist/vmin))
print("dist / 3000 = %f" % (dist/3000))
print("Tmax = %f" % Tmax)
print("dist = %f" % dist)
print("tmin = %f" % tmin)
print("tmax = %f" % tmax)
print("tend = %f" % tend)


# measure SNR
G = 1.0 + 1.0/nstack
G = 1.0
SNR_ls = measure_snr(tr_ls, tmin, tmax, tend, dt)


# selective stack
tr_ls2 = np.zeros(npts)
#stack_file_list = []
stack_file_list = ''
nstack2 = 0
flag = np.ones(nstack)
for i in range(nstack):
	tr = tr_ls - st[i].data
	SNR_tmp = measure_snr(tr, tmin, tmax, tend, dt)
	Q = SNR_tmp / SNR_ls
	if (Q > G):
		flag[i] = 0
		#tr_ls2 += st[i].data
	else:
		nstack2 += 1
		tr_ls2 += st[i].data
		#filename = file_path + '/' + file_list[i] + '\n'
		filename = ' < ' + file_path + '/' + file_list[i]
		stack_file_list += filename
		#stack_file_list.append(file_list[i])



tr_ls2 /= nstack2
SNR_ls2 = measure_snr(tr_ls2, tmin, tmax, tend, dt)

tr = st[0].copy()
sac = SACTrace.from_obspy_trace(tr)
tr.data = tr_ls2
tr.write('./selective_stack_test/tr_ls_selective.SAC', format='sac')



print(nstack)
print(nstack2)
print(len(stack_file_list))
print("Q = %f, G = %f \n" % (Q, G))
print(flag)
print("SNR of ls  %f" % SNR_ls)
print("SNR of ls2 %f" % SNR_ls2)
f2 = 100.0
f3 = 3.0
weight = 1
#cmd_str = os.listdir(file_path + '/') + ' | TF_PWS -B %f' % (1.0/f3)
#cmd_str = 'ls ' + file_path + '/*.SAC' + ' | TF_PWS -B %f' % (1.0/f3)
#cmd_str = 'ls ' + file_path + '/*.SAC | TF_PWS -B ' + str(f3) + \
#          ' -E ' + str(f2) + ' -W ' + str(weight) + \
#          ' -O ' + './selective_stack_test' + '/stack/NCISP6.NE00_NCISP6.NE05' + ' -P 1'

#print(stack_file_list)

cmd_str = stack_file_list + ' | TF_PWS -B ' + str(f3) + \
          ' -E ' + str(f2) + ' -W ' + str(weight) + \
          ' -O ' + './selective_stack_test' + '/stack/NCISP6.NE00_NCISP6.NE05' + ' -P 1'

print(cmd_str)
#os.system(cmd_str)


#files = os.system('ls ' + file_path + '/*.SAC')
#print(files)
#cmd_str = file_path + '/' + stack_file_list + ' | TF_PWS -B ' + str(f3) + \
#          ' -E ' + str(f2) + ' -W ' + str(weight) + \
#          ' -O ' + './selective_stack_test' + '/stack/NCISP6.NE00_NCISP6.NE05' + ' -P 1'
#print(cmd_str)
#os.system(cmd_str)
#print(os.listdir(file_path + '/') + ' | TF_PWS -B ' + str(1.0/f3) + \
#          ' -E ' + str(1.0/f2) + ' -W ' + str(weight) + \
#          ' -O ' + 'NCISP6.NE00_NCISP6.NE05' + ' -P 1')
#os.system(os.listdir(file_path + '/prestack/*.SAC') + ' | TF_PWS -B ' + str(1.0/f3) + \
#          ' -E ' + str(1.0/f2) + ' -W ' + str(weight) + \
#          ' -O ' + 'NCISP6.NE00_NCISP6.NE05' + ' -P ' + '1')
### create stack file list
##stack_file_list = []
##for
