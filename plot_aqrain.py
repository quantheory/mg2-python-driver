#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits import basemap
import netCDF4 as nc4

cmap = plt.get_cmap('coolwarm')

REF_FILE_NAME = '/home/santos/Data/nsubr_control_ANN_climo.nc'
TEST_FILE_NAME = '/home/santos/Data/short_timestep_ANN_climo.nc'

rfile = nc4.Dataset(REF_FILE_NAME, 'r')
tfile = nc4.Dataset(TEST_FILE_NAME, 'r')

lat = rfile['lat']
lon = rfile['lon']
lev = rfile['lev']
ilev = rfile['ilev']
gw = rfile['gw']

print("Sum of gw:", gw[:].sum())
gw = gw[:]/gw[:].sum()

def global_mean(var):
    var_weighted = var.copy()
    for i in range(len(lat)):
        var_weighted[i,:] *= gw[i]
    return var_weighted.sum()/len(lon)

ref_rain = rfile['AQRAIN'][0,:,:,:]
test_rain = tfile['AQRAIN'][0,:,:,:]

ref_rain = np.mean(ref_rain, axis=2)
test_rain = np.mean(test_rain, axis=2)

i300 = 0
for level in ilev:
    if level > 300.:
        break
    i300+= 1

def forward(a):
    a = np.deg2rad(a)
    return np.sin(a)

def inverse(a):
    a = np.arcsin(a)
    return np.rad2deg(a)

plt.pcolor(lat[1:], ilev[i300:], ref_rain[i300:,1:-1]*1.e6)
plt.axis('tight')
ax = plt.gca()
ylim = ax.get_ylim()
ax.set_ylim([ylim[1], ylim[0]])
ax.set_yticks([1000., 800., 600., 400.])
ax.set_xlim([80., -80.])
ax.set_xticks([60., 30., 0., -30., -60.])
ax.set_xticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("Reference case average rain mass (mg/kg)")
plt.colorbar()
plt.clim([0., 8.])
plt.savefig('AQRAIN_ref.png')
plt.close()

plt.pcolor(lat[1:], ilev[i300:], test_rain[i300:,1:-1]*1.e6)
plt.axis('tight')
ax = plt.gca()
ylim = ax.get_ylim()
ax.set_ylim([ylim[1], ylim[0]])
ax.set_yticks([1000., 800., 600., 400.])
ax.set_xlim([80., -80.])
ax.set_xticks([60., 30., 0., -30., -60.])
ax.set_xticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("MG2 at 1s case average rain mass (mg/kg)")
plt.colorbar()
plt.clim([0., 8.])
plt.savefig('AQRAIN_short.png')
plt.close()

plt.pcolor(lat[1:], ilev[i300:], (test_rain[i300:,1:-1] - ref_rain[i300:,1:-1])*1.e6, cmap=cmap)
plt.axis('tight')
ax = plt.gca()
ylim = ax.get_ylim()
ax.set_ylim([ylim[1], ylim[0]])
ax.set_yticks([1000., 800., 600., 400.])
ax.set_xlim([80., -80.])
ax.set_xticks([60., 30., 0., -30., -60.])
ax.set_xticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("Difference between average rain masses (mg/kg)")
plt.colorbar()
plt.clim([-1.1, 1.1])
plt.savefig('AQRAIN_diff.png')
plt.close()


ref_rain = rfile['AQRAIN'][0,:,:,:]
test_rain = tfile['AQRAIN'][0,:,:,:]

ref_freqr = rfile['FREQR'][0,:,:,:]
test_freqr = tfile['FREQR'][0,:,:,:]

for i in range(len(lev)):
    for j in range(len(lat)):
        for k in range(len(lon)):
            if ref_freqr[i,j,k] == 0.:
                ref_rain[i,j,k] = 0.
            else:
                ref_rain[i,j,k] = ref_rain[i,j,k] / ref_freqr[i,j,k]
            if test_freqr[i,j,k] == 0.:
                test_rain[i,j,k] = 0.
            else:
                test_rain[i,j,k] = test_rain[i,j,k] / test_freqr[i,j,k]


bmap = basemap.Basemap(lon_0=180.)

plt.pcolormesh(lon[:], lat[:], ref_rain[50,:,:]*1.e6)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
plt.clim([0., 450.])
plt.savefig('AQRAIN_IA_681mb_ref.png')
plt.close()

plt.pcolormesh(lon[:], lat[:], test_rain[50,:,:]*1.e6)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
plt.clim([0., 450.])
plt.savefig('AQRAIN_IA_681mb_short.png')
plt.close()

plt.pcolormesh(lon[:], lat[:], (test_rain[50,:,:] - ref_rain[50,:,:])*1.e6, cmap=cmap)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
plt.clim([-400., 400.])
plt.savefig('AQRAIN_IA_681mb_diff.png')
plt.close()

ref_rain = np.mean(ref_rain, axis=2)
test_rain = np.mean(test_rain, axis=2)

plt.pcolor(lat[1:], ilev[i300:], ref_rain[i300:,1:-1]*1.e6)
plt.axis('tight')
ax = plt.gca()
ylim = ax.get_ylim()
ax.set_ylim([ylim[1], ylim[0]])
ax.set_yticks([1000., 800., 600., 400.])
ax.set_xscale('function', functions=(forward, inverse))
ax.set_xticks([80., 50., 30., 0., -30., -50., -80.])
ax.set_xticklabels(['80N', '50N', '30N', '0', '30S', '50S', '80S'])
plt.title("Reference case in-cloud rain mass (mg/kg)")
plt.ylabel("Pressure (mb)")
plt.colorbar()
plt.clim([0., 90.])
plt.savefig('AQRAIN_IA_ref.png')
plt.close()

plt.pcolor(lat[1:], ilev[i300:], test_rain[i300:,1:-1]*1.e6)
plt.axis('tight')
ax = plt.gca()
ylim = ax.get_ylim()
#ax.set_yscale('function', functions=(forward, inverse))
ax.set_ylim([ylim[1], ylim[0]])
ax.set_yticks([1000., 800., 600., 400.])
ax.set_xscale('function', functions=(forward, inverse))
ax.set_xticks([80., 50., 30., 0., -30., -50., -80.])
ax.set_xticklabels(['80N', '50N', '30N', '0', '30S', '50S', '80S'])
plt.title("MG2 at 1s case in-cloud rain mass (mg/kg)")
plt.ylabel("Pressure (mb)")
plt.colorbar()
plt.clim([0., 90.])
plt.savefig('AQRAIN_IA_short.png')
plt.close()

plt.pcolor(lat[1:], ilev[i300:], (test_rain[i300:,1:-1] - ref_rain[i300:,1:-1])*1.e6, cmap=cmap)
plt.axis('tight')
ax = plt.gca()
ylim = ax.get_ylim()
#ax.set_yscale('function', functions=(forward, inverse))
ax.set_ylim([ylim[1], ylim[0]])
ax.set_yticks([1000., 800., 600., 400.])
ax.set_xscale('function', functions=(forward, inverse))
ax.set_xticks([80., 50., 30., 0., -30., -50., -80.])
ax.set_xticklabels(['80N', '50N', '30N', '0', '30S', '50S', '80S'])
plt.title("Difference between in-cloud rain masses (mg/kg)")
plt.ylabel("Pressure (mb)")
plt.colorbar()
plt.clim([-25., 25.])
plt.savefig('AQRAIN_IA_diff.png')
plt.close()

print(lev[50], ilev[50:52])

ref_prect = rfile['PRECC'][0,:,:] + rfile['PRECL'][0,:,:]
test_prect = tfile['PRECC'][0,:,:] + tfile['PRECL'][0,:,:]

ref_prect *= 1000. * 86400.
test_prect *= 1000. * 86400.

print("Average PRECTs: ", global_mean(ref_prect), global_mean(test_prect),
      global_mean(test_prect - ref_prect))

plt.pcolormesh(lon[:], lat[:], ref_prect)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("Reference case average total precipitation (mm/day)")
plt.colorbar()
plt.clim([0.,35.])
plt.savefig('PRECT_ref.png')
plt.close()

plt.pcolormesh(lon[:], lat[:], test_prect)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("MG2 at 1s case average total precipitation (mm/day)")
plt.colorbar()
plt.clim([0.,35.])
plt.savefig('PRECT_short.png')
plt.close()

plt.pcolormesh(lon[:], lat[:], test_prect - ref_prect, cmap=cmap)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("Net difference in average total precipitation (mm/day)")
plt.colorbar()
plt.clim([-3., 3.])
plt.savefig('PRECT_diff.png')
plt.close()

ref_precl = rfile['PRECL'][0,:,:]
test_precl = tfile['PRECL'][0,:,:]

ref_precl *= 1000. * 86400.
test_precl *= 1000. * 86400.

print("Average PRECLs: ", global_mean(ref_precl), global_mean(test_precl),
      global_mean(test_precl - ref_precl))

plt.pcolormesh(lon[:], lat[:], ref_precl)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("Reference case average total precipitation (mm/day)")
plt.colorbar()
plt.clim([0.,35.])
plt.savefig('PRECL_ref.png')
plt.close()

plt.pcolormesh(lon[:], lat[:], test_precl)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("MG2 at 1s case average total precipitation (mm/day)")
plt.colorbar()
plt.clim([0.,35.])
plt.savefig('PRECL_short.png')
plt.close()

plt.pcolormesh(lon[:], lat[:], test_precl - ref_precl, cmap=cmap)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("Net difference in average total precipitation (mm/day)")
plt.colorbar()
plt.clim([-3., 3.])
plt.savefig('PRECL_diff.png')
plt.close()

ref_precc = rfile['PRECC'][0,:,:]
test_precc = tfile['PRECC'][0,:,:]

ref_precc *= 1000. * 86400.
test_precc *= 1000. * 86400.

print("Average PRECCs: ", global_mean(ref_precc), global_mean(test_precc),
      global_mean(test_precc - ref_precc))

plt.pcolormesh(lon[:], lat[:], ref_precc)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("Reference case average total precipitation (mm/day)")
plt.colorbar()
plt.clim([0.,35.])
plt.savefig('PRECC_ref.png')
plt.close()

plt.pcolormesh(lon[:], lat[:], test_precc)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("MG2 at 1s case average total precipitation (mm/day)")
plt.colorbar()
plt.clim([0.,35.])
plt.savefig('PRECC_short.png')
plt.close()

plt.pcolormesh(lon[:], lat[:], test_precc - ref_precc, cmap=cmap)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.title("Net difference in average total precipitation (mm/day)")
plt.colorbar()
plt.clim([-3., 3.])
plt.savefig('PRECC_diff.png')
plt.close()


ref_precc = rfile['SWCF'][0,:,:]
test_precc = tfile['SWCF'][0,:,:]

print("Average SWCFs: ", global_mean(ref_precc), global_mean(test_precc),
      global_mean(test_precc - ref_precc))

ref_precc = rfile['LWCF'][0,:,:]
test_precc = tfile['LWCF'][0,:,:]

print("Average LWCFs: ", global_mean(ref_precc), global_mean(test_precc),
      global_mean(test_precc - ref_precc))

ref_precc = rfile['TGCLDLWP'][0,:,:]
test_precc = tfile['TGCLDLWP'][0,:,:]

print("Average TGCLDLWPs: ", global_mean(ref_precc), global_mean(test_precc),
      global_mean(test_precc - ref_precc))

ref_precc = rfile['TGCLDIWP'][0,:,:]
test_precc = tfile['TGCLDIWP'][0,:,:]

print("Average TGCLDIWPs: ", global_mean(ref_precc), global_mean(test_precc),
      global_mean(test_precc - ref_precc))

ref_precc = rfile['CLDHGH'][0,:,:]
test_precc = tfile['CLDHGH'][0,:,:]

print("Average CLDHGHs: ", global_mean(ref_precc), global_mean(test_precc),
      global_mean(test_precc - ref_precc))

ref_precc = rfile['CLDLOW'][0,:,:]
test_precc = tfile['CLDLOW'][0,:,:]

print("Average CLDLOWs: ", global_mean(ref_precc), global_mean(test_precc),
      global_mean(test_precc - ref_precc))

ref_precc = rfile['CLDTOT'][0,:,:]
test_precc = tfile['CLDTOT'][0,:,:]

print("Average CLDTOTs: ", global_mean(ref_precc), global_mean(test_precc),
      global_mean(test_precc - ref_precc))

