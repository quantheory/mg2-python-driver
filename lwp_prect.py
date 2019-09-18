#!/usr/bin/env python

import numpy as np
import scipy.linalg as la
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numdifftools as ndt
import sklearn.cluster as clstr

from mg2_constants import *

HIST_FILE_NAME = "/p/lscratchh/santos36/ACME/short_timestep_ctrl_diags/run/short_timestep_ctrl_diags.cam.h1.0001-01-01-72000.nc"
SHORT_HIST_FILE_NAME = "/p/lscratchh/santos36/ACME/short_timestep_diags/run/short_timestep_diags.cam.h1.0001-01-01-72000.nc"

nbins = 20
qsmall = 1.e-18

cmap = plt.get_cmap('coolwarm')

file = nc4.Dataset(HIST_FILE_NAME, 'r')

ncol = len(file.dimensions['ncol'])
lev = len(file.dimensions['lev'])
ilev = len(file.dimensions['ilev'])

lwp = file.variables["MPLWPI"][0,:]
twp = file.variables["MPLWPI"][0,:] + file.variables["MPIWPI"][0,:]
precl = file.variables["PRECL"][0,:]
prect = file.variables["PRECT"][0,:]
qrsedten = file.variables["QRSEDTEN"][0,:,:]
evaprain = file.variables["EVAPPREC"][0,:,:] - file.variables["EVAPSNOW"][0,:,:]
prc = file.variables["PRCO"][0,:,:]
pra = file.variables["PRAO"][0,:,:]
pmids = file.variables['lev'][:] * 100.
pints = file.variables['ilev'][:] * 100.

precl_long = precl
prect_long = prect

PLT_FILE = "lwp_precl"

plt.plot(lwp, precl, '.')
plt.xlabel("Water path (kg/m^2)")
plt.ylabel("PRECL (m/s)")
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

plt.plot(lwp, precl, '.')
plt.xlabel("Water path (kg/m^2)")
plt.ylabel("PRECL (m/s)")
plt.axis([0., 6., 0., 8.e-7])
plt.savefig(PLT_FILE+"_scaled.eps")
plt.savefig(PLT_FILE+"_scaled.png")
plt.close()

PLT_FILE = "lwp_precl_hist"

plt.hist2d(lwp, precl, bins=nbins, range=((0.05, 1), (7.e-9, 1.e-7)), cmap=cmap)
plt.xlabel("Water path (kg/m^2)")
plt.ylabel("PRECL (m/s)")
plt.axis('tight')
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

qrsedtot = np.zeros((ncol,))
evapraintot = np.zeros((ncol,))
for i in range(lev):
    qrsedtot[:] = qrsedtot[:] + np.abs(qrsedten[i,:])*(pints[i] - pints[i+1])/(lev*gravit)
    evapraintot[:] = evapraintot[:] + evaprain[i,:]*(pints[i] - pints[i+1])/gravit

evaprainavg = np.zeros((lev,))
prcavg = np.zeros((lev,))
praavg = np.zeros((lev,))

for n in range(ncol):
    evaprainavg[:] = evaprainavg[:] + evaprain[:,n]*(pints[:-1] - pints[1:])/gravit
    prcavg[:] = prcavg[:] + prc[:,n]*(pints[:-1] - pints[1:])/gravit
    praavg[:] = praavg[:] + pra[:,n]*(pints[:-1] - pints[1:])/gravit

PLT_FILE = "evaprain_qrsedtot"

plt.plot(evapraintot, qrsedtot, '.')
plt.xlabel("Rain evaporation rate (kg/m^2)")
plt.ylabel("Sedimentation column average rate (kg/m^2)")
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

PLT_FILE = "precl_hist"

bins = np.linspace(qsmall, 2.e-8, nbins+1)
plt.hist(precl, bins=bins)
plt.xlabel("PRECL (m/s)")
plt.axis([0., 2.e-8, 0, 30000])
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

PLT_FILE = "evaprain_hist"

bins = np.linspace(-5.e-5, -qsmall, nbins+1)
plt.hist(evapraintot, bins=bins)
plt.xlabel("Rain evaporation rate (kg/m^2)")
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

file = nc4.Dataset(SHORT_HIST_FILE_NAME, 'r')

ncol = len(file.dimensions['ncol'])
lev = len(file.dimensions['lev'])
ilev = len(file.dimensions['ilev'])

lwp = file.variables["MPLWPI"][0,:]
twp = file.variables["MPLWPI"][0,:] + file.variables["MPIWPI"][0,:]
precl = file.variables["PRECL"][0,:]
prect = file.variables["PRECT"][0,:]
qrsedten = file.variables["QRSEDTEN"][0,:,:]
evaprain = file.variables["EVAPPREC"][0,:,:] - file.variables["EVAPSNOW"][0,:,:]
prc = file.variables["PRCO"][0,:,:]
pra = file.variables["PRAO"][0,:,:]
pmids = file.variables['lev'][:] * 100.
pints = file.variables['ilev'][:] * 100.

PLT_FILE = "lwp_precl_short"

plt.plot(lwp, precl, '.')
plt.xlabel("Water path (kg/m^2)")
plt.ylabel("PRECT (m/s)")
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

plt.plot(lwp, precl, '.')
plt.xlabel("Water path (kg/m^2)")
plt.ylabel("PRECT (m/s)")
plt.axis([0., 6., 0., 8.e-7])
plt.savefig(PLT_FILE+"_scaled.eps")
plt.savefig(PLT_FILE+"_scaled.png")
plt.close()

PLT_FILE = "lwp_precl_hist_short"

plt.hist2d(lwp, precl, bins=nbins, range=((0.05, 1.), (5.e-9, 1.e-7)), cmap=cmap)
plt.xlabel("Water path (kg/m^2)")
plt.ylabel("PRECL (m/s)")
plt.axis('tight')
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

qrsedtot = np.zeros((ncol,))
evapraintot = np.zeros((ncol,))
for i in range(lev):
    qrsedtot[:] = qrsedtot[:] + np.abs(qrsedten[i,:])*(pints[i] - pints[i+1])/(lev*gravit)
    evapraintot[:] = evapraintot[:] + evaprain[i,:]*(pints[i] - pints[i+1])/gravit

evaprainavg_short = np.zeros((lev,))
prcavg_short = np.zeros((lev,))
praavg_short = np.zeros((lev,))

for n in range(ncol):
    evaprainavg_short[:] = evaprainavg_short[:] + evaprain[:,n]*(pints[:-1] - pints[1:])/gravit
    prcavg_short[:] = prcavg_short[:] + prc[:,n]*(pints[:-1] - pints[1:])/gravit
    praavg_short[:] = praavg_short[:] + pra[:,n]*(pints[:-1] - pints[1:])/gravit

PLT_FILE = "evaprain_qrsedtot_short"

plt.plot(evapraintot, qrsedtot, '.')
plt.xlabel("Rain evaporation rate (kg/m^2)")
plt.ylabel("Sedimentation column average rate (kg/m^2)")
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

PLT_FILE = "evaprain_vertical"

plt.plot(evaprainavg, pmids, label='300s timestep')
plt.plot(evaprainavg_short, pmids, label='1s timestep')
plt.xlabel("Rain evaporation rate (kg/m^2)")
plt.ylabel("Pressure (Pa)")
plt.axis('tight')
plt.gca().invert_yaxis()
plt.legend(loc='best')
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

PLT_FILE = "prc_vertical"

plt.plot(prcavg, pmids, label='300s timestep')
plt.plot(prcavg_short, pmids, label='1s timestep')
plt.xlabel("Autoconversion rate (kg/m^2)")
plt.ylabel("Pressure (Pa)")
plt.axis('tight')
plt.gca().invert_yaxis()
plt.legend(loc='best')
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

PLT_FILE = "pra_vertical"

plt.plot(praavg, pmids, label='300s timestep')
plt.plot(praavg_short, pmids, label='1s timestep')
plt.xlabel("Accretion rate (kg/m^2)")
plt.ylabel("Pressure (Pa)")
plt.axis('tight')
plt.gca().invert_yaxis()
plt.legend(loc='best')
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

PLT_FILE = "precl_hist_short"

bins = np.linspace(qsmall, 2.e-8, nbins+1)
plt.hist(precl, bins=bins)
plt.xlabel("PRECL (m/s)")
plt.axis([0., 2.e-8, 0, 30000])
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

PLT_FILE = "precl_hist_overlap"

bins = np.linspace(qsmall, 2.e-8, nbins+1)
plt.hist(precl, bins=bins, label='1s timestep')
plt.hist(precl_long, bins=bins, label='300s timestep')
plt.xlabel("PRECL (m/s)")
plt.axis([0., 2.e-8, 0, 30000])
plt.legend(loc='best')
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

PLT_FILE = "evaprain_hist_short"

bins = np.linspace(-1.5e-5, -qsmall, nbins+1)
plt.hist(evapraintot, bins=bins)
plt.xlabel("Rain evaporation rate (kg/m^2)")
plt.savefig(PLT_FILE+".eps")
plt.savefig(PLT_FILE+".png")
plt.close()

print('Mean PRECL @ 300s = ', np.mean(precl_long))
print('Mean PRECL @ 1s = ', np.mean(precl))
print('Mean PRECT @ 300s = ', np.mean(prect_long))
print('Mean PRECT @ 1s = ', np.mean(prect))
