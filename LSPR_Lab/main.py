#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from pprint import pprint

def find_fwhm(xdata, ydata):
    peak_idx = np.argmax(ydata)
    hm = ydata[peak_idx]/2
    min_fw_idx = np.argmin(np.abs(ydata[:peak_idx] - hm))
    max_fw_idx = np.argmin(np.abs(ydata[peak_idx:] - hm)) + peak_idx  

    fw = np.abs(xdata[min_fw_idx]-xdata[max_fw_idx])

    min_fw = (xdata[min_fw_idx], ydata[min_fw_idx])
    max_fw = (xdata[max_fw_idx], ydata[max_fw_idx])
    
    return min_fw, max_fw

def extract_from_xml(filename):
    X_set = []
    Y_set = []
    X = []
    Y = []
    in_data = False
    with open(filename, 'r') as f:
        for line in f:
            if '<data>' in line:
                in_data = True
                data_str = line.split('>')[1]
                values = data_str.split('\t')
                X.append(float(values[0]))
                Y.append(float(values[1]))
            elif '</data>' in line:
                in_data = False
                data_str = line.split('<')[0]
                values = data_str.split('\t')
                values = data_str.split('\t')
                X.append(float(values[0]))
                Y.append(float(values[1]))
                X_set.append(np.array(X))
                Y_set.append(np.array(Y))
                X = []
                Y = []
            elif in_data:
                values = line.split('\t')
                X.append(float(values[0]))
                Y.append(float(values[1]))
                
    return X_set, Y_set


print('Starting calculations...')

# Doing calculations for the unknown samples
# First measurement is a calibration run with only clean glass

# Load all measurement data, including the baseline run
trans_data = np.genfromtxt('data/trans_data.csv', delimiter=',')
trans_data = trans_data[:,:-1].T

# Remove the baseline and the two bonus measurements
trans_data = trans_data[2:-4,:]

trans_x = trans_data[::2,:]
trans_y = 100-trans_data[1::2,:]

# Find peaks and FWHM
peak_idx = [np.argmax(data) for data in trans_y]
norm_fwhm = []
fwhm_list = []
peak_list = []
amp_list = []

for x, y, i in zip(trans_x, trans_y, peak_idx):
    width = find_fwhm(x, y)
    fwhm = width[0][0] - width[1][0]
    fwhm_list.append(fwhm)
    norm_fwhm.append(fwhm/y[i])
    peak_list.append(x[i])
    amp_list.append(y[i])

print('Identified peaks:')
pprint(peak_list)
print('Identified peak amplitudes:')
pprint(peak_list)
print('Identified FWHM:')
pprint(fwhm_list)
print('FWHM/A:')
pprint(norm_fwhm)

# Plot stuff
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=14)

plt.figure()
plt.title('Extinction Spectrum')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Extiction [\%]')
ax = plt.axes()
n = ord('A')
for x, y, i in zip(trans_x, trans_y, peak_idx):
    plt.plot(x,y, label='Measurement {0}'.format(chr(n)), linewidth=2)
    fw = find_fwhm(x, y)
    ax.annotate('', xy=fw[0], xytext=fw[1], arrowprops=dict(arrowstyle='<->'))
    #plt.plot(fw[0], fw[1], 'b-o')
    n += 1


    
plt.plot([x[i] for x, i in zip(trans_x,peak_idx)],
         [y[i] for y, i in zip(trans_y,peak_idx)],
         'rd', label='Peaks')
    
plt.legend(loc='upper right')
plt.savefig('Images/measurement_data.png')


# Data processing for the liquid interface measurement
conc =  [0,
         10,
         20,
         30,
         40,
         50,
         60,
         70,
         80,
         90,
         100]

refr_idx = [1.33300,
            1.34242,
            1.35238,
            1.36253,
            1.37275,
            1.38313,
            1.39336,
            1.40340,
            1.41315,
            1.42262,
            1.43188]

abs_dens = [0.99707,
            1.0096,
            1.0229,
            1.0364,
            1.0495,
            1.0620,
            1.0740,
            1.0850,
            1.0945,
            1.1029,
            1.1099]

slope = (refr_idx[-1]-refr_idx[0])/(conc[-1]-conc[0])
r_idx = lambda c: slope*c + refr_idx[0]

levels = [755.365,
          756.302,
          756.914,
          757.989,
          758.662,
          759.859]

conc_data = [ 0,
             10,
             17,
             29,
             36,
             47]

ridx_data = [r_idx(c) for c in conc_data]

l_diff = [x2 - x1 for x1, x2 in zip(levels[:-1], levels[1:])]
n_diff = [x2 - x1 for x1, x2 in zip(ridx_data[:-1], ridx_data[1:])]
ln_ratio = [l/n for l, n in zip(l_diff, n_diff)]

print('Refractive index:')
pprint(ridx_data)
print('Delta lamda [nm]:')
pprint(l_diff)
print('Delta n:')
pprint(n_diff)
print('dl/dn [nm]:')
pprint(ln_ratio)
print('Mean resolution [nm]:')
pprint(np.mean(ln_ratio))

plt.figure()
plt.title('$n$ vs Concentration')
plt.xlabel('Concentration [wt\%]')
plt.ylabel('$n$')
plt.grid(True)
plt.xlim(0, 100)
plt.ylim(refr_idx[0], refr_idx[-1])
ax = plt.axes()
plt.plot(conc, refr_idx, linewidth=2, label='Interpolated data')
plt.plot(conc_data, ridx_data, 'r^', markersize=10, label='Experimental data')
plt.legend(loc='upper left')
for x, y in zip(conc_data, ridx_data):
    x0 = 0
    y0 = refr_idx[0]
    ax.annotate(s='', xy=(x,y0), xytext=(x,y), 
                arrowprops=dict(arrowstyle='->'))
    ax.annotate(s='', xy=(x0,y), xytext=(x,y), 
                arrowprops=dict(arrowstyle='->'))

plt.savefig('Images/refractive_vs_conc.png')

plt.figure()
plt.title('Sensitivity')
plt.xlabel('Step')
plt.ylabel(r'$\frac{\Delta n}{\Delta \lambda }$')
plt.grid(True)
plt.plot(range(1,len(ln_ratio)+1), ln_ratio, '-o', markersize=10)
plt.plot([1,5], [np.mean(ln_ratio), np.mean(ln_ratio)])

plt.savefig('Images/sensitivity.png')

philip_data = np.genfromtxt('data/peak_pos.txt')
philip_x = philip_data[:,0]
philip_y = philip_data[:,1]

plt.figure()
plt.title('Extiction Spectrum')
plt.xlabel('Time [s]')
plt.ylabel('Extiction [\%]')
plt.grid(True)
plt.plot(philip_x,philip_y, linewidth=2)

plt.savefig('Images/liquid_interface.png')

X_set, Y_set = extract_from_xml('data/skalman_lspr.inr')

# plt.figure()
# plt.title('Extiction')
# plt.xlabel('Wavelength [nm]')
# plt.ylabel('Transmission [\%]')
# plt.grid(True)
# for x, y in zip(X_set[0], Y_set[0]):
#     plt.plot(x,y,'b-.')

plt.show()
