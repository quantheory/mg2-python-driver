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
#blocksize = 4096
#num_files = 1
#end_column = 4095

timescale_type = 'eigenvalue'
#timescale_type = 'depletion_time'

splits = [(blocksize*i, blocksize*(i+1)-1) for i in range(num_files)]
splits[-1] = (splits[-1][0], end_column)

EVALS_FILE_NAMES = ["/home/santos/Data/Jacobian_cutoff_{}-{}.nc".format(split[0], split[1])
                    for split in splits]
CLUSTER_FILE_NAME = "/home/santos/Data/MG2_data_collection.10_cluster_labels.0001-01-06-00000.nc"

efiles = []
for name in EVALS_FILE_NAMES:
    efiles.append(nc4.Dataset(name, 'r'))
cfile = nc4.Dataset(CLUSTER_FILE_NAME, 'r')

if timescale_type == 'depletion_time':
    HIST_FILE_NAME = "/home/santos/Data/MG2_data_collection.cam.h1.0001-01-06-00000.nc"
    hfile = nc4.Dataset(HIST_FILE_NAME, 'r')

nproc = len(efiles[0].dimensions['nproc'])
ncluster = len(cfile.dimensions['ncluster'])

label = cfile.variables["label"]

timestep = 1.e-8
evap_col_steps = 1
evap_steps = 1
col_steps = 1
auto_accr_steps = 1
auto_steps = 1
accr_steps = 1

# Indices for each variable in the state array.
it = 0
iq = 1
iqc = 2
inc = 3
iqi = 4
ini = 5
iqr = 6
inr = 7
iqs = 8
ins = 9

if timescale_type == 'depletion_time':
    # Just a dumb way to get a 10-item list so we can change it.
    hvars = list(range(10))
    hvars[it] = hfile.variables["MG2IN_T"]
    hvars[iq] = hfile.variables["MG2IN_Q"]
    hvars[iqc] = hfile.variables["MG2IN_QC"]
    hvars[inc] = hfile.variables["MG2IN_NC"]
    hvars[iqi] = hfile.variables["MG2IN_QI"]
    hvars[ini] = hfile.variables["MG2IN_NI"]
    hvars[iqr] = hfile.variables["MG2IN_QR"]
    hvars[inr] = hfile.variables["MG2IN_NR"]
    hvars[iqs] = hfile.variables["MG2IN_QS"]
    hvars[ins] = hfile.variables["MG2IN_NS"]

def calc_twmd(vec):
    av = np.abs(vec)
    return (av[iq] + av[iqc] + av[iqi] + av[iqr] + av[iqs]) * 0.5

# Dictionary associating "short" to "pretty" names for MG2 processes.
short_names = [
    "cauto",
    "caccr",
    "revap",
    "rnagg",
    "racrs",
    "rfrez",
    "ssubl",
    "snagg",
    "vdepo",
    "cacwi",
    "cbrgi",
    "iauto",
    "iaccr",
    "cbrgs",
    "cacws",
    "cnact",
    "vnudp",
    "cnccc",
    "cncct",
    "anadj",
]
process_names = {
    "revap": "Rain Evap.",
    "ssubl": "Snow Subl.",
    "vdepo": "Vapor/Ice Transfer",
    "vnudp": "Nucleation Dep.",
    "cbrgs": "Berg. (Snow)",
    "cacws": "Liq. Accr. Snow",
    "cnccc": "Immersion Frz.",
    "cncct": "Contact Frz.",
    "cacwi": "Sec. Ice Prod.",
    "rfrez": "Het. Rain Frz.",
    "racrs": "Rain Accr. Snow",
    "cbrgi": "Berg. (Cloud)",
    "cauto": "Autoconversion",
    "caccr": "Accretion",
    "iauto": "Ice Auto.",
    "iaccr": "Ice Accretion",
    "rnagg": "Rain Self-col.",
    "snagg": "Snow Self-col.",
    "cnact": "Drop. Activ.",
    "anadj": "Size Limiters",
}

# Need a way to translate processes in a different order.
name_map = np.zeros((nproc,), dtype='u4')
for i in range(nproc):
    name_string = efiles[0]["process_names"][i]
    for j in range(nproc):
        if name_string == short_names[j]:
            name_map[j] = i

ind = np.arange(nproc)

min_t = 1.e-1
max_t = 1.e5
cutoff2 = 0.5
# process must be at least this active to ride this ride.
q_cutoff = 1.e-10
# equivalent to number of large rain particles made of q_cutoff mass.
n_cutoff = q_cutoff / (np.pi/6. * 1000. * (400.e-6)**3)
print("Number cutoff (1/kg/s): ", n_cutoff)

def process_is_relevant(tends):
    if calc_twmd(tends) >= q_cutoff:
        return True
    at = np.abs(tends)
    return (at[inc] > n_cutoff or at[ini] > n_cutoff or
            at[inr] > n_cutoff or at[ins] > n_cutoff)

plt.autoscale(tight=True)

evalues_all = [[] for i in range(ncluster)]
evalues_rel = [[] for i in range(ncluster)]
evalues2 = [dict() for i in range(ncluster)]
evalue_correlation = dict()
cluster_cases = [0 for i in range(ncluster)]
for name in short_names:
    for j in range(ncluster):
        evalues2[j][name] = []
    evalue_correlation[name] = dict()
    for name2 in short_names:
        evalue_correlation[name][name2] = 0.

for efile in efiles:
    num_cell = len(efile.dimensions['num_cell'])
    for ci in range(num_cell):
        column = efile["cell_coords"][ci,0]
        level = efile["cell_coords"][ci,1]
        tends = efile["process_rates"][ci,:,:]
        c = label[0,level,column]
        cluster_cases[c] += 1
        # Actually use the eigenvalues
        if timescale_type == 'eigenvalue':
            evals = efile["eigenvalues"][ci,:]["real"]
            associations = efile["associations"][ci,:,:]
            for i in range(10):
                evalues_all[c].append(evals[i])
                maxproc = np.argmax(associations[i,:])
                if process_is_relevant(tends[:,maxproc]):
                    evalues_rel[c].append(evals[i])
                for j in range(nproc):
                    if (associations[i,name_map[j]] >= cutoff2 and process_is_relevant(tends[:,name_map[j]])):
                        evalues2[c][short_names[j]].append(np.real(evals[i]))
                        for j2 in range(nproc):
                            evalue_correlation[short_names[j]][short_names[j2]] += associations[i,name_map[j2]]
        elif timescale_type == 'depletion_time':
            for i in range(nproc):
                if not process_is_relevant(tends[:,name_map[i]]):
                    continue
                timescale = 1.e100
                for j in range(10):
                    if tends[j,name_map[i]] < 0.:
                        this_timescale = hvars[j][0, level, column] / tends[j,name_map[i]]
                        if np.abs(this_timescale) < np.abs(timescale):
                            timescale = this_timescale
                evalues_all[c].append(1./timescale)
                evalues_rel[c].append(1./timescale)
                evalues2[c][short_names[i]].append(1./timescale)

#Convert the correlation sum to an average.
evalue_corr_array = np.zeros((nproc,nproc))
for i in range(nproc):
    num_evalues = 0
    for c in range(ncluster):
        num_evalues += len(evalues2[c][short_names[i]])
    if num_evalues == 0:
        continue
    for j in range(nproc):
        evalue_corr_array[i,j] = evalue_correlation[short_names[i]][short_names[j]] / num_evalues
        if i == j:
            evalue_corr_array[i,j] -= 0.5

if timescale_type == 'eigenvalue':
    suffix = ''
elif timescale_type == 'depletion_time':
    suffix = '_depletion'

# Number of bins for plotting purposes.
nbins = 50

for c in range(ncluster):
    print("Number of cases in cluster {} is {}.".format(c, cluster_cases[c]))

cmap = plt.get_cmap('Greys')

# Histogram of all positive eigenvalues.
bins = np.logspace(np.log10(1./max_t), np.log10(1./min_t), nbins+1)
evalues_all_clusters = []
evalues_all_rel = []
evalues_all_assoc = []
for c in range(ncluster):
    evalues_all_clusters += [t for t in evalues_all[c] if t > 1./max_t and t < 1./min_t]
    evalues_all_rel += [t for t in evalues_rel[c] if t > 1./max_t and t < 1./min_t]
    for name in short_names:
        evalues_all_assoc += [t for t in evalues2[c][name] if t > 1./max_t and t < 1./min_t]
plt.hist(evalues_all_clusters, bins=bins)
plt.hist(evalues_all_rel, bins=bins)
plt.hist(evalues_all_assoc, bins=bins)
plt.title("Number of positive eigenvalues \n(based on {} eigenvalues)"
          .format(len(evalues_all_clusters)))
plt.gca().set_xscale("log")
plt.xlabel("Eigenvalue (1/s)")
plt.xlim(1./max_t, 1./min_t)
plt.ylabel("Number of eigenvalues")
plt.savefig('./time_hist_all_values_pos{}.eps'.format(suffix))
plt.close()

# Histogram of all negative eigenvalues.
bins = -np.logspace(np.log10(1./min_t), np.log10(1./max_t), nbins+1)
evalues_all_clusters = []
evalues_all_rel = []
evalues_all_assoc = []
for c in range(ncluster):
    evalues_all_clusters += [t for t in evalues_all[c] if t > -1./min_t and t < -1./max_t]
    evalues_all_rel += [t for t in evalues_rel[c] if t > -1./min_t and t < -1./max_t]
    for name in short_names:
        evalues_all_assoc += [t for t in evalues2[c][name] if t > -1./min_t and t < -1./max_t]
plt.hist(evalues_all_clusters, bins=bins, color='b', edgecolor='k')
plt.hist(evalues_all_rel, bins=bins, color='g', edgecolor='k')
plt.hist(evalues_all_assoc, bins=bins, color='r', edgecolor='k')
plt.axvline(x=-1./300., color='k', linewidth=2.)
plt.title("Number of negative eigenvalues \n(based on {} eigenvalues)"
          .format(len(evalues_all_clusters)))
plt.gca().set_xscale("symlog", linthresh=1./max_t)
plt.xlabel("Eigenvalue (1/s)")
plt.xlim(-1./min_t, -1./max_t)
plt.ylabel("Number of eigenvalues")
plt.savefig('./time_hist_all_values_neg{}.eps'.format(suffix))
plt.close()

# Histogram of all near-zero eigenvalues.
bins = np.linspace(-1./max_t, 1./max_t, nbins+1)
evalues_all_clusters = []
evalues_all_rel = []
evalues_all_assoc = []
for c in range(ncluster):
    evalues_all_clusters += [t for t in evalues_all[c] if t <= 1./max_t and t >= -1./max_t]
    evalues_all_rel += [t for t in evalues_rel[c] if t <= 1./max_t and t >= -1./max_t]
    for name in short_names:
        evalues_all_assoc += [t for t in evalues2[c][name] if t <= 1./max_t and t >= -1./max_t]
plt.hist(evalues_all_clusters, bins=bins)
plt.hist(evalues_all_rel, bins=bins)
plt.hist(evalues_all_assoc, bins=bins)
plt.title("Number of near-zero eigenvalues \n(based on {} eigenvalues)"
          .format(len(evalues_all_clusters)))
plt.xlabel("Eigenvalue (1/s)")
plt.xlim(-1./max_t, 1./max_t)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel("Number of eigenvalues")
plt.savefig('./time_hist_all_values_n0{}.eps'.format(suffix))
plt.close()

# Now broken out by cluster.
cind = np.arange(ncluster+1)

# Histogram of all positive eigenvalues.
bins = np.logspace(np.log10(1./max_t), np.log10(1./min_t), nbins+1)
hist_values = np.zeros((ncluster, nbins))
eig_count = []
for c in range(ncluster):
    hist_values[c,:], _ = np.histogram([t for t in evalues_all[c] if t > 1./max_t and t < 1./min_t], bins=bins)
    row_norm = hist_values[c,:].sum()
    print("Cluster {} has {} positive eigenvalues.".format(c, row_norm))
    eig_count.append(int(row_norm))
    row_max = hist_values[c,:].max()
    if row_norm != 0:
        hist_values[c,:] /= row_max
plt.pcolor(bins, cind, hist_values, edgecolors='k', cmap=cmap)
plt.title("Distribution of positive eigenvalues by cluster")
plt.gca().set_xscale("log")
plt.xlabel("Eigenvalue (1/s)")
plt.xlim(1./max_t, 1./min_t)
plt.ylabel("Cluster Index")
ax = plt.gca()
ax.set_yticks(cind[:-1])
ax.set_yticklabels([str(i) for i in cind[:-1]],
                   fontdict={'verticalalignment': 'bottom'})
ax.tick_params('y', direction='out')
plt.clim(vmin=0., vmax=1.)
plt.colorbar(pad=0.1)
ax2 = ax.twinx()
ax2.set_yticks(cind)
ax2.set_yticklabels([str(i) for i in eig_count] + [""],
                    fontdict={'verticalalignment': 'bottom'})
plt.savefig('./time_hist_cluster_2D_pos{}.eps'.format(suffix))
plt.close()

# Histogram of all negative eigenvalues.
bins = -np.logspace(np.log10(1./min_t), np.log10(1./max_t), nbins+1)
hist_values = np.zeros((ncluster, nbins))
eig_count = []
for c in range(ncluster):
    hist_values[c,:], _ = np.histogram([t for t in evalues_all[c] if -t > 1./max_t and -t < 1./min_t], bins=bins)
    row_norm = hist_values[c,:].sum()
    print("Cluster {} has {} negative eigenvalues.".format(c, row_norm))
    eig_count.append(int(row_norm))
    row_max = hist_values[c,:].max()
    if row_norm != 0:
        hist_values[c,:] /= row_max
plt.pcolor(bins, cind, hist_values, edgecolors='k', cmap=cmap)
plt.axvline(x=-1./300., color='k', linewidth=2.)
plt.title("Distribution of negative eigenvalues by cluster")
plt.gca().set_xscale("symlog", linthresh=1./max_t)
plt.xlabel("Eigenvalue (1/s)")
plt.xlim(-1./min_t, -1./max_t)
plt.ylabel("Cluster Index")
ax = plt.gca()
ax.set_yticks(cind[:-1])
ax.set_yticklabels([str(i) for i in cind[:-1]],
                   fontdict={'verticalalignment': 'bottom'})
ax.tick_params('y', direction='out')
plt.clim(vmin=0., vmax=1.)
plt.colorbar(pad=0.15)
ax2 = ax.twinx()
ax2.set_yticks(cind)
ax2.set_yticklabels([str(i) for i in eig_count] + [""],
                    fontdict={'verticalalignment': 'bottom'})
plt.savefig('./time_hist_cluster_2D_neg{}.eps'.format(suffix))
plt.close()

# Histogram of all near-zero eigenvalues.
bins = np.linspace(-1./max_t, 1./max_t, nbins+1)
hist_values = np.zeros((ncluster, nbins))
eig_count = []
for c in range(ncluster):
    hist_values[c,:], _ = np.histogram([t for t in evalues_all[c] if t <= 1./max_t and t >= -1./max_t], bins=bins)
    row_norm = hist_values[c,:].sum()
    print("Cluster {} has {} near-zero eigenvalues.".format(c, row_norm))
    eig_count.append(int(row_norm))
    row_max = hist_values[c,:].max()
    if row_norm != 0:
        hist_values[c,:] /= row_max
plt.pcolor(bins, cind, hist_values, edgecolors='k', cmap=cmap)
plt.title("PDF of near-zero eigenvalues by cluster")
plt.xlabel("Eigenvalue (1/s)")
plt.xlim(-1./max_t, 1./max_t)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel("Cluster Index")
ax = plt.gca()
ax.set_yticks(cind[:-1])
ax.set_yticklabels([str(i) for i in cind[:-1]],
                   fontdict={'verticalalignment': 'bottom'})
ax.tick_params('y', direction='out')
plt.clim(vmin=0., vmax=1.)
plt.colorbar(pad=0.1)
ax2 = ax.twinx()
ax2.set_yticks(cind)
ax2.set_yticklabels([str(i) for i in eig_count] + [""],
                    fontdict={'verticalalignment': 'bottom'})
plt.savefig('./time_hist_cluster_2D_n0{}.eps'.format(suffix))
plt.close()

# Now broken up by process.
pind = np.arange(nproc, -1, -1)
bins = np.logspace(np.log10(1./max_t), np.log10(1./min_t), nbins+1)
hist_values = np.zeros((nproc, nbins))
eig_count = []
for i in range(nproc):
    pos_evalues = []
    for c in range(ncluster):
        pos_evalues += [t for t in evalues2[c][short_names[i]] if t > 1./max_t and t < 1./min_t]
    hist_values[i,:], _ = np.histogram(pos_evalues, bins=bins)
    row_norm = hist_values[i,:].sum()
    eig_count.append(int(row_norm))
    row_max = hist_values[i,:].max()
    print("Process {} has {} positive eigenvalues.".format(process_names[short_names[i]], row_norm))
    if row_norm != 0:
        hist_values[i,:] /= row_max
plt.pcolor(bins, pind, hist_values, edgecolors='k', cmap=cmap)
plt.title("Distribution of positive eigenvalues by process")
plt.gca().set_xscale("log")
plt.xlabel("Eigenvalue (1/s)")
plt.xlim(1./max_t, 1./min_t)
ax = plt.gca()
ax.set_yticks(pind[:-1])
ax.set_yticklabels((process_names[name] for name in short_names),
                   fontdict={'verticalalignment': 'top'},
                   size='small', wrap=True)
ax.tick_params('y', direction='out')
plt.subplots_adjust(left=0.25)
plt.ylabel("Process")
plt.clim(vmin=0., vmax=1.)
plt.colorbar(pad=0.16)
ax2 = ax.twinx()
ax2.set_yticks(pind[:-1])
ax2.set_yticklabels([str(i) for i in eig_count],
                    fontdict={'verticalalignment': 'top'})
plt.savefig('./time_hist_process_2D_pos{}.eps'.format(suffix))
plt.close()

# Histogram of all negative eigenvalues.
bins = -np.logspace(np.log10(1./min_t), np.log10(1./max_t), nbins+1)
hist_values = np.zeros((nproc, nbins))
eig_count = []
for i in range(nproc):
    neg_evalues = []
    for c in range(ncluster):
        neg_evalues += [t for t in evalues2[c][short_names[i]] if -t > 1./max_t and -t < 1./min_t]
    hist_values[i,:], _ = np.histogram(neg_evalues, bins=bins)
    row_sum = hist_values[i,:].sum()
    eig_count.append(int(row_sum))
    row_max = hist_values[i,:].max()
    print("Process {} has {} negative eigenvalues.".format(process_names[short_names[i]], row_norm))
    if row_norm != 0:
        hist_values[i,:] /= row_max
plt.pcolor(bins, pind, hist_values, edgecolors='k', cmap=cmap)
plt.axvline(x=-1./300., color='k', linewidth=2.)
plt.title("Distribution of negative eigenvalues by process")
plt.gca().set_xscale("symlog", linthresh=1./max_t)
plt.xlabel("Eigenvalue (1/s)")
plt.xlim(-1./min_t, -1./max_t)
ax = plt.gca()
ax.set_yticks(pind[:-1])
ax.set_yticklabels((process_names[name] for name in short_names),
                   fontdict={'verticalalignment': 'top'},
                   size='small', wrap=True)
ax.tick_params('y', direction='out')
plt.subplots_adjust(left=0.25)
plt.ylabel("Process")
plt.clim(vmin=0., vmax=1.)
plt.colorbar(pad=0.16)
ax2 = ax.twinx()
ax2.set_yticks(pind[:-1])
ax2.set_yticklabels([str(i) for i in eig_count],
                    fontdict={'verticalalignment': 'top'})
plt.savefig('./time_hist_process_2D_neg{}.eps'.format(suffix))
plt.close()

# Histogram of all near-zero eigenvalues.
bins = np.linspace(-1./max_t, 1./max_t, nbins+1)
hist_values = np.zeros((nproc, nbins))
eig_count = []
for i in range(nproc):
    n0_evalues = []
    for c in range(ncluster):
        n0_evalues += [t for t in evalues2[c][short_names[i]] if t <= 1./max_t and t >= -1./max_t]
    hist_values[i,:], _ = np.histogram(n0_evalues, bins=bins)
    row_norm = hist_values[i,:].sum()
    eig_count.append(int(row_norm))
    row_max = hist_values[i,:].max()
    print("Process {} has {} non-zero eigenvalues.".format(process_names[short_names[i]], row_norm))
    if row_norm != 0:
        hist_values[i,:] /= row_max
plt.pcolor(bins, pind, hist_values, edgecolors='k', cmap=cmap)
plt.title("Distribution of near-zero eigenvalues by process")
plt.xlabel("Eigenvalue (1/s)")
plt.xlim(-1./max_t, 1./max_t)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax = plt.gca()
ax.set_yticks(pind[:-1])
ax.set_yticklabels((process_names[name] for name in short_names),
                   fontdict={'verticalalignment': 'top'},
                   size='small', wrap=True)
ax.tick_params('y', direction='out')
plt.subplots_adjust(left=0.25)
plt.ylabel("Process")
plt.clim(vmin=0., vmax=1.)
plt.colorbar(pad=0.1)
ax2 = ax.twinx()
ax2.set_yticks(pind[:-1])
ax2.set_yticklabels([str(i) for i in eig_count],
                    fontdict={'verticalalignment': 'top'})
plt.savefig('./time_hist_process_2D_n0{}.eps'.format(suffix))
plt.close()

# Correlation array.
print("==========evalue_corr_array==========")
print(evalue_corr_array)
print("=====================================")
plt.pcolor(np.flipud(pind), pind, evalue_corr_array, edgecolors='k', cmap=cmap)
plt.title("Degree of association between given eigenvalues by primary process")
ax = plt.gca()
ax.set_xticks(np.flipud(pind[1:]))
ax.set_xticklabels((process_names[name] for name in short_names),
                   size='small', rotation='vertical', wrap=True)
ax.tick_params('x', direction='out', pad=40)
plt.subplots_adjust(bottom=0.25)
#plt.xlabel("Process with given cutoff")
ax.set_yticks(pind[:-1])
ax.set_yticklabels((process_names[name] for name in short_names),
                   fontdict={'verticalalignment': 'top'},
                   size='small', wrap=False)
ax.tick_params('y', direction='out')
plt.subplots_adjust(left=0.25)
plt.ylabel("Primary process")
plt.clim(vmin=0., vmax=0.4)
plt.colorbar()
plt.savefig('./process_association{}.eps'.format(suffix))
plt.close()


# Cluster and process.
for c in range(ncluster):
    pos_evalues2 = {}
    hist_values = np.zeros((nproc, nbins))
    bins = np.logspace(np.log10(1./max_t), np.log10(1./min_t), nbins+1)

    nmodes = 0

    for i in range(nproc):
        process = short_names[i]
        pos_evalues2[process] = np.array([t for t in evalues2[c][process] if t > 1./max_t and t < 1./min_t])
        nmodes += len(pos_evalues2[process])
        hist_values[i,:], _ = np.histogram(pos_evalues2[process], bins=bins)

    plt.pcolor(bins, pind, hist_values, edgecolors='k', cmap=cmap)
    plt.title("Number of positive eigenvalues for cluster {} \n(based on {} eigenvalues, cutoff={})"
              .format(c,
                      nmodes,
                      cutoff2))
    plt.gca().set_xscale("log")
    plt.xlabel("Eigenvalue (1/s)")
    plt.xlim(1./max_t, 1./min_t)
    ax = plt.gca()
    ax.set_yticks(pind[:-1])
    ax.set_yticklabels((process_names[name] for name in short_names),
                       fontdict={'verticalalignment': 'top'},
                       size='small', wrap=True)
    ax.tick_params('y', direction='out')
    plt.subplots_adjust(left=0.25)
    plt.colorbar()
    plt.savefig('./time_hist_2D_pos_c{}{}.eps'.format(c, suffix))
    plt.close()

    neg_evalues2 = {}
    hist_values = np.zeros((nproc, nbins))
    bins = np.logspace(np.log10(1./max_t), np.log10(1./min_t), nbins+1)

    nmodes = 0

    for i in range(nproc):
        process = short_names[i]
        neg_evalues2[process] = np.array([-t for t in evalues2[c][process] if -t > 1./max_t and -t < 1./min_t])
        nmodes += len(neg_evalues2[process])
        hist_values[i,:], _ = np.histogram(neg_evalues2[process], bins=bins)
        row_norm = hist_values[i,:].sum()
        if row_norm != 0:
            hist_values[i,:] /= row_norm

    plt.pcolor(bins, pind, hist_values, edgecolors='k', cmap=cmap)
    plt.axvline(x=1./300., color='k', linewidth=2.)
    plt.title("Number of negative eigenvalues for cluster {} \n(based on {} eigenvalues, cutoff={})"
              .format(c,
                      nmodes,
                      cutoff2))
    plt.gca().set_xscale("log")
    plt.xlabel("Eigenvalue (1/s)")
    plt.xlim(1./min_t, 1./max_t)
    ax = plt.gca()
    ax.set_yticks(pind[:-1])
    ax.set_yticklabels((process_names[name] for name in short_names),
                       fontdict={'verticalalignment': 'top'},
                       size='small', wrap=True)
    ax.tick_params('y', direction='out')
    plt.subplots_adjust(left=0.25)
    plt.ylabel("Process")
    plt.clim(vmin=0., vmax=0.2)
    plt.colorbar()
    plt.savefig('./time_hist_2D_neg_c{}{}.eps'.format(c, suffix))
    plt.close()

    n0_evalues2 = {}
    hist_values = np.zeros((nproc, nbins))
    bins = np.linspace(-1./max_t, 1./max_t, nbins+1)

    nmodes = 0

    for i in range(nproc):
        process = short_names[i]
        n0_evalues2[process] = np.array([t for t in evalues2[c][process] if t <= 1./max_t and t >= -1./max_t])
        nmodes += len(n0_evalues2[process])
        hist_values[i,:], _ = np.histogram(n0_evalues2[process], bins=bins)
        row_norm = hist_values[i,:].sum()
        if row_norm != 0:
            hist_values[i,:] /= row_norm

    plt.pcolor(bins, pind, hist_values, edgecolors='k', cmap=cmap)
    plt.title("Number of near-zero eigenvalues for cluster {} \n(based on {} eigenvalues, cutoff={})"
              .format(c,
                      nmodes,
                      cutoff2))
    plt.xlabel("Eigenvalue (1/s)")
    plt.xlim(-1./max_t, 1./max_t)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax = plt.gca()
    ax.set_yticks(pind[:-1])
    ax.set_yticklabels((process_names[name] for name in short_names),
                       fontdict={'verticalalignment': 'top'},
                       size='small', wrap=True)
    ax.tick_params('y', direction='out')
    plt.subplots_adjust(left=0.25)
    plt.ylabel("Process")
    plt.clim(vmin=0., vmax=0.2)
    plt.colorbar()
    plt.savefig('./time_hist_2D_n0_c{}{}.eps'.format(c, suffix))
    plt.close()
