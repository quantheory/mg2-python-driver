#!/usr/bin/env python

import numpy as np
import scipy.linalg as la
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

CLUSTER_FILE_NAME = "/home/santos/Data/MG2_data_collection.10_cluster_labels.0001-01-06-00000.nc"

cfile = nc4.Dataset(CLUSTER_FILE_NAME, 'r')

ncluster = len(cfile.dimensions['ncluster'])
nproc = len(cfile.dimensions['nproc'])

cluster_centers = np.abs(cfile.variables["cluster_centers"][:,:])
process_names = [
    "Rain Evap.",
    "Snow Subl.",
    "Vapor/Ice DMS",
    "Nucleation Dep.",
    "Berg. (Snow)",
    "Liq. Accr. Snow",
    "Immersion Frz.",
    "Contact Frz.",
    "Sec. Ice Prod.",
    "Het. Rain Frz.",
    "Rain Accr. Snow",
    "Berg. (Cloud)",
    "Autoconversion",
    "Accretion",
    "Ice Auto.",
    "Ice Accretion",
    "Rain Self-col.",
    "Snow Self-col.",
    "Drop. Activ.",
]
#process_names = cfile.variables["process_names"]

pind = np.arange(nproc+1)
cind = np.arange(ncluster+1)

for i in range(nproc):
    proc_max = np.amax(cluster_centers[:,i])
    if proc_max > 0.:
        cluster_centers[:,i] /= proc_max

cmap = plt.get_cmap('Reds')
max_val = np.abs(cluster_centers).max()
plt.autoscale(tight=True)
plt.pcolor(pind, cind, np.abs(cluster_centers), edgecolors='k', cmap=cmap)
#for i in range(nproc):
#    for j in range(ncluster):
#        plt.text(i,j+0.5,'{:0.1f}'.format(cluster_centers[j,i]),color='k',fontweight='demi',
#                 horizontalalignment='left',verticalalignment='center',fontsize=9)
plt.title("Average process rate for each cluster (scaled)")
plt.xlabel("Process")
plt.ylabel("Cluster index")
ax = plt.gca()
ax.set_xticks(pind)
ax.set_xticklabels(process_names,
                   size='small', rotation='vertical', wrap=True)
ax.tick_params('x', direction='out', pad=40)
plt.subplots_adjust(bottom=0.20)
ax.set_yticks(cind)
ax.set_yticklabels([str(i) for i in cind[:-1]],
                   fontdict={'verticalalignment': 'bottom'})
ax.tick_params('y', direction='out')
plt.clim(0, max_val)
plt.colorbar()
plt.savefig('./cluster_centers_2D_scaled.eps')
plt.close()
