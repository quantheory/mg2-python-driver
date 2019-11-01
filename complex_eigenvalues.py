#!/usr/bin/env python

import sys

import numpy as np
import scipy.linalg as la
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numdifftools as ndt

from mg2 import wv_sat_methods as wsm
from mg2 import micro_mg2_0 as mg

from mg2_constants import *

blocksize = 4096
num_files = 12
end_column = 48601

splits = [(blocksize*i, blocksize*(i+1)-1) for i in range(num_files)]
splits[-1] = (splits[-1][0], end_column)

EVALS_FILE_NAMES = ["/home/santos/Data/Jacobian_cutoff_{}-{}.nc".format(split[0], split[1])
                    for split in splits]

efiles = []
for name in EVALS_FILE_NAMES:
    efiles.append(nc4.Dataset(name, 'r'))

plt.autoscale(tight=True)
cmap = plt.get_cmap('Reds')

total_evalues = 0
for efile in efiles:
    total_evalues += 10*len(efile.dimensions['num_cell'])

all_evalues = np.zeros((total_evalues,), dtype='c16')

i = 0
for efile in efiles:
    num_cell = len(efile.dimensions['num_cell'])
    all_evalues[i:i+10*num_cell] = efile["eigenvalues"][:,:]["real"].flatten() + \
                                   efile["eigenvalues"][:,:]["imag"].flatten()*1.j
    i += 10*num_cell


# Plot values in the complex plane.
max_exp = 1.
cutoff_exp = -5.
cutoff = 10.**cutoff_exp
log_polar_evalues = np.zeros((total_evalues, 2))
for i in range(total_evalues):
    abs_eval = np.abs(all_evalues[i])
    if abs_eval > cutoff:
        log_polar_evalues[i,0] = np.log10(abs_eval) - cutoff_exp
log_polar_evalues[:,1] = np.angle(all_evalues)

cart_log_evalues = np.zeros((total_evalues, 2))
cart_log_evalues[:,0] = log_polar_evalues[:,0] * np.cos(log_polar_evalues[:,1])
cart_log_evalues[:,1] = log_polar_evalues[:,0] * np.sin(log_polar_evalues[:,1])

max_plotval = max_exp - cutoff_exp

nbins = 51
assert nbins % 2 == 1
bins = np.linspace(-max_plotval, max_plotval, nbins+1)

plot_data, _, _ = np.histogram2d(cart_log_evalues[:,1], cart_log_evalues[:,0],
                                 bins=[bins, bins])

plt.pcolor(bins, bins, plot_data, edgecolors='k', cmap=cmap)
circle_rad = np.log10(1./300.)-cutoff_exp
circle_x = circle_rad * np.cos(np.linspace(0., 2.*np.pi, 1001))
circle_y = circle_rad * np.sin(np.linspace(0., 2.*np.pi, 1001))
plt.plot(circle_x, circle_y)
converge_rad = 1./300.
converge_x = converge_rad * (np.cos(np.linspace(0., 2.*np.pi, 1001)) - 1.)
converge_y = converge_rad * np.sin(np.linspace(0., 2.*np.pi, 1001))
converge_abs = np.sqrt(converge_x**2 + converge_y**2)
converge_angle = np.arctan2(converge_y, converge_x)
plot_abs = np.maximum(np.log10(converge_abs) - cutoff_exp, 0.)
plot_x = plot_abs * np.cos(converge_angle)
plot_y = plot_abs * np.sin(converge_angle)
plt.plot(plot_x, plot_y)
plt.axvline(x=0., color='k', linewidth=2.)
plt.axhline(y=0., color='k', linewidth=2.)
plt.axis([-max_plotval, max_plotval, -max_plotval, max_plotval])
ticks = np.arange(-max_plotval, max_plotval+1., 1)
# The "0.01" hack below is there to ensure that "0" outputs a 10 with a positive sign.
tickvals = ["${}^{{{}}}$".format(int(np.sign(i+0.01)*10),int(abs(i)+cutoff_exp)) for i in ticks]
ax = plt.gca()
ax.set_xticks(ticks)
ax.set_xticklabels(tickvals)
ax.set_yticks(ticks)
ax.set_yticklabels(tickvals)
plt.clim(vmin=0., vmax=1.e4)
plt.colorbar()
plt.savefig('./complex_eigenvalues.png')
plt.close()

midpoint = nbins // 2
real_zeros = np.zeros((nbins,))
near_real_zeros = np.zeros((nbins,))
imag_zeros = np.zeros((nbins,))

for i in range(total_evalues):
    real_part = cart_log_evalues[i,0]
    if real_part < bins[0] or real_part >= bins[-1]:
        continue
    bin_index = nbins*2
    for j in range(nbins):
        if real_part < bins[j+1]:
            bin_index = j
            break
    if all_evalues[i].imag == 0:
        real_zeros[bin_index] += 1
    elif np.abs(all_evalues[i].imag) <= 1.e-5:
        near_real_zeros[bin_index] += 1
    else:
        imag_zeros[bin_index] += 1

all_zeros = real_zeros + near_real_zeros + imag_zeros
bin_centers = (bins[:-1] + bins[1:]) / 2

plt.bar(bin_centers, all_zeros,
        width=(bins[1]-bins[0]), color='r', label='All zeros')
plt.bar(bin_centers, real_zeros,
        width=(bins[1]-bins[0]), color='b', label='Real zeros')
plt.legend(loc='best')
ticks = np.arange(-max_plotval, max_plotval+1., 1)
tickvals = ["${}^{{{}}}$".format(int(np.sign(i+0.01)*10),int(abs(i)+cutoff_exp)) for i in ticks]
ax = plt.gca()
ax.set_xlim(left=bins[0], right=bins[-1])
ax.set_xticks(ticks)
ax.set_xticklabels(tickvals)
plt.savefig('./complex_vertint.png')
plt.close()

plt.plot(bin_centers, real_zeros / all_zeros, color='k')
ticks = np.arange(-max_plotval, max_plotval+1., 1)
tickvals = ["${}^{{{}}}$".format(int(np.sign(i+0.01)*10),int(abs(i)+cutoff_exp)) for i in ticks]
ax = plt.gca()
ax.set_xticks(ticks)
ax.set_xticklabels(tickvals)
plt.savefig('./complex_ratio.png')
plt.close()

plt.bar(bin_centers, all_zeros,
        width=(bins[1]-bins[0]), color='r', label='All zeros')
plt.bar(bin_centers, real_zeros + near_real_zeros,
        width=(bins[1]-bins[0]), color='b', label='Real zeros')
plt.legend(loc='best')
ticks = np.arange(-max_plotval, max_plotval+1., 1)
tickvals = ["${}^{{{}}}$".format(int(np.sign(i+0.01)*10),int(abs(i)+cutoff_exp)) for i in ticks]
ax = plt.gca()
ax.set_xlim(left=bins[0], right=bins[-1])
ax.set_xticks(ticks)
ax.set_xticklabels(tickvals)
plt.savefig('./complex_vertint_cutoff.png')
plt.close()

plt.plot(bin_centers, (real_zeros + near_real_zeros) / all_zeros, color='k')
ticks = np.arange(-max_plotval, max_plotval+1., 1)
tickvals = ["${}^{{{}}}$".format(int(np.sign(i+0.01)*10),int(abs(i)+cutoff_exp)) for i in ticks]
ax = plt.gca()
ax.set_xticks(ticks)
ax.set_xticklabels(tickvals)
plt.savefig('./complex_ratio_cutoff.png')
plt.close()
