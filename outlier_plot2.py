#!/usr/bin/env python

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import netCDF4 as nc4

from mg2 import wv_sat_methods as wsm
from mg2 import micro_mg2_0 as mg

from mg2_constants import *

HIST_FILE_NAME = "/g/g14/santos36/Data/MG2_data_collection.cam.h1.0001-01-06-00000.nc"

file = nc4.Dataset(HIST_FILE_NAME, 'r')

ncol = len(file.dimensions['ncol'])
lev = len(file.dimensions['lev'])
ilev = len(file.dimensions['ilev'])

errstring = wsm.wv_sat_methods_init(kind, tmelt, h2otrip, tboil, ttrice, epsilo)

assert errstring.decode().strip() == '', \
    "wv_sat_methods initialization error: "+errstring.decode()

errstring = mg.micro_mg_init(kind, gravit, rair, rh2o, cpair, tmelt, latvap,
                             latice, rhmini, dcs, dcs_tdep, uniform, do_cldice,
                             use_hetfrz_classnuc, precip_frac_method,
                             berg_eff_factor, allow_sed_supersat, ice_sed_ai,
                             prc_coef1, prc_exp, prc_exp1, cld_sed,
                             mg_prc_coeff_fix, alpha_grad, beta_grad)

assert errstring.decode().strip() == '', \
    "MG2 initialization error: "+errstring.decode()

mgncol = 1
t = file.variables["MG2IN_T"]
q = file.variables["MG2IN_Q"]
qc = file.variables["MG2IN_QC"]
qi = file.variables["MG2IN_QI"]
nc = file.variables["MG2IN_NC"]
ni = file.variables["MG2IN_NI"]
qr = file.variables["MG2IN_QR"]
qs = file.variables["MG2IN_QS"]
nr = file.variables["MG2IN_NR"]
ns = file.variables["MG2IN_NS"]
relvar = file.variables["MG2IN_RELVAR"]
accre_enhan = file.variables["MG2IN_ACCRE_ENHAN"]
p = file.variables["MG2IN_P"]
pdel = file.variables["MG2IN_PDEL"]
precipf = file.variables["MG2IN_PRECIP"]
liqcldf = file.variables["MG2IN_LIQCLDF"]
icecldf = file.variables["MG2IN_ICECLDF"]
naai = file.variables["MG2IN_NAAI"]
npccn = file.variables["MG2IN_NPCCN"]
rndst = np.empty((t.shape[0], t.shape[1], t.shape[2], 4))
rndst[:,:,:,0] = file.variables["MG2IN_RNDST1"][:]
rndst[:,:,:,1] = file.variables["MG2IN_RNDST2"][:]
rndst[:,:,:,2] = file.variables["MG2IN_RNDST3"][:]
rndst[:,:,:,3] = file.variables["MG2IN_RNDST4"][:]
nacon = np.empty((t.shape[0], t.shape[1], t.shape[2], 4))
nacon[:,:,:,0] = file.variables["MG2IN_NACON1"][:]
nacon[:,:,:,1] = file.variables["MG2IN_NACON2"][:]
nacon[:,:,:,2] = file.variables["MG2IN_NACON3"][:]
nacon[:,:,:,3] = file.variables["MG2IN_NACON4"][:]
frzimm = file.variables["MG2IN_FRZIMM"]
frzcnt = file.variables["MG2IN_FRZCNT"]
frzdep = file.variables["MG2IN_FRZDEP"]

t_loc = np.empty((mgncol, t.shape[1]), order='F')
q_loc = np.empty((mgncol, q.shape[1]), order='F')
qc_loc = np.empty((mgncol, qc.shape[1]), order='F')
qi_loc = np.empty((mgncol, qi.shape[1]), order='F')
nc_loc = np.empty((mgncol, nc.shape[1]), order='F')
ni_loc = np.empty((mgncol, ni.shape[1]), order='F')
qr_loc = np.empty((mgncol, qr.shape[1]), order='F')
qs_loc = np.empty((mgncol, qs.shape[1]), order='F')
nr_loc = np.empty((mgncol, nr.shape[1]), order='F')
ns_loc = np.empty((mgncol, ns.shape[1]), order='F')
relvar_loc = np.empty((mgncol, relvar.shape[1]), order='F')
accre_enhan_loc = np.empty((mgncol, accre_enhan.shape[1]), order='F')
p_loc = np.empty((mgncol, p.shape[1]), order='F')
pdel_loc = np.empty((mgncol, pdel.shape[1]), order='F')
precipf_loc = np.empty((mgncol, precipf.shape[1]), order='F')
liqcldf_loc = np.empty((mgncol, liqcldf.shape[1]), order='F')
icecldf_loc = np.empty((mgncol, icecldf.shape[1]), order='F')
naai_loc = np.empty((mgncol, naai.shape[1]), order='F')
npccn_loc = np.empty((mgncol, npccn.shape[1]), order='F')
rndst_loc = np.empty((mgncol, rndst.shape[1], 4), order='F')
nacon_loc = np.empty((mgncol, nacon.shape[1], 4), order='F')
frzimm_loc = np.empty((mgncol, frzimm.shape[1]), order='F')
frzcnt_loc = np.empty((mgncol, frzcnt.shape[1]), order='F')
frzdep_loc = np.empty((mgncol, frzdep.shape[1]), order='F')

final_time = 1800

timesteps = np.array([60, 120])
loc_arrays = {
    'T': t_loc,
    'Q': q_loc,
    'QC': qc_loc,
    'QI': qi_loc,
    'QR': qr_loc,
    'QS': qs_loc,
    'NR': nr_loc,
}
var_names = sorted(list(loc_arrays.keys()))
finals = {}
roughs = {}
for name in var_names:
    finals[name] = np.zeros((1 + final_time // timesteps[0], lev))
    roughs[name] = np.zeros((1 + final_time // timesteps[1], lev))

rain_evap_fine = np.zeros((final_time // timesteps[0], lev))
rain_evap_coarse = np.zeros((final_time // timesteps[1], lev))
rain_sed_fine = np.zeros((final_time // timesteps[0], lev))
rain_sed_coarse = np.zeros((final_time // timesteps[1], lev))
rain_auto_fine = np.zeros((final_time // timesteps[0], lev))
rain_auto_coarse = np.zeros((final_time // timesteps[1], lev))
rain_accr_fine = np.zeros((final_time // timesteps[0], lev))
rain_accr_coarse = np.zeros((final_time // timesteps[1], lev))

max_ind = 1532

for it in range(timesteps.size):
    print("Starting timestep=", timesteps[it])
    assert final_time % timesteps[it] == 0
    nsteps = final_time // timesteps[it]
    deltat = float(timesteps[it])

    t_loc[:,:] = t[0,:,max_ind].transpose()
    q_loc[:,:] = q[0,:,max_ind].transpose()
    qc_loc[:,:] = qc[0,:,max_ind].transpose()
    qi_loc[:,:] = qi[0,:,max_ind].transpose()
    nc_loc[:,:] = nc[0,:,max_ind].transpose()
    ni_loc[:,:] = ni[0,:,max_ind].transpose()
    qr_loc[:,:] = qr[0,:,max_ind].transpose()
    qs_loc[:,:] = qs[0,:,max_ind].transpose()
    nr_loc[:,:] = nr[0,:,max_ind].transpose()
    ns_loc[:,:] = ns[0,:,max_ind].transpose()
    relvar_loc[:,:] = relvar[0,:,max_ind].transpose()
    accre_enhan_loc[:,:] = accre_enhan[0,:,max_ind].transpose()
    p_loc[:,:] = p[0,:,max_ind].transpose()
    pdel_loc[:,:] = pdel[0,:,max_ind].transpose()
    precipf_loc[:,:] = precipf[0,:,max_ind].transpose()
    liqcldf_loc[:,:] = liqcldf[0,:,max_ind].transpose()
    icecldf_loc[:,:] = icecldf[0,:,max_ind].transpose()
    naai_loc[:,:] = naai[0,:,max_ind].transpose()
    npccn_loc[:,:] = npccn[0,:,max_ind].transpose()
    rndst_loc[:,:,:] = rndst[0,:,max_ind:max_ind+1,:].transpose([1, 0, 2])
    nacon_loc[:,:,:] = nacon[0,:,max_ind:max_ind+1,:].transpose([1, 0, 2])
    frzimm_loc[:,:] = frzimm[0,:,max_ind].transpose()
    frzcnt_loc[:,:] = frzcnt[0,:,max_ind].transpose()
    frzdep_loc[:,:] = frzdep[0,:,max_ind].transpose()

    if it == 0:
        for name in var_names:
            finals[name][0,:] = loc_arrays[name]
    else:
        for name in var_names:
            roughs[name][0,:] = loc_arrays[name]

    for n in range(nsteps):
        qcsinksum_rate1ord, tlat, qvlat, qctend, qitend, nctend, nitend, qrtend, \
            qstend, nrtend, nstend, effc, effc_fn, effi, prect, preci, nevapr, \
            evapsnow, prain, prodsnow, cmeout, deffi, pgamrad, lamcrad, qsout, dsout, \
            rflx, sflx, qrout, reff_rain, reff_snow, qcsevap, qisevap, qvres, cmeitot, \
            vtrmc, vtrmi, umr, ums, qcsedten, qisedten, qrsedten, qssedten, pratot, \
            prctot, mnuccctot, mnuccttot, msacwitot, psacwstot, bergstot, bergtot, \
            melttot, homotot, qcrestot, prcitot, praitot, qirestot, mnuccrtot, \
            pracstot, meltsdttot, frzrdttot, mnuccdtot, nrout, nsout, refl, arefl, \
            areflz, frefl, csrfl, acsrfl, fcsrfl, rercld, ncai, ncal, qrout2, qsout2, \
            nrout2, nsout2, drout2, dsout2, freqs, freqr, nfice, qcrat, errstring, \
            prer_evap \
            = mg.micro_mg_tend(deltat, t_loc, q_loc, qc_loc, qi_loc, nc_loc,
                               ni_loc, qr_loc, qs_loc, nr_loc, ns_loc,
                               relvar_loc, accre_enhan_loc, p_loc, pdel_loc,
                               precipf_loc, liqcldf_loc, icecldf_loc, naai_loc,
                               npccn_loc, rndst_loc, nacon_loc,
                               frzimm=frzimm_loc, frzcnt=frzcnt_loc,
                               frzdep=frzdep_loc, mgncol=mgncol, nlev=lev)

        # Should use geopotential!
        t_loc += tlat * deltat / cpair
        q_loc += qvlat * deltat
        q_loc[:,:] = np.where(q_loc < 1.e-12, 1.e-12, q_loc)
        qc_loc += qctend * deltat
        qc_loc[:,:] = np.where(qc_loc < 0., 0., qc_loc)
        qi_loc += qitend * deltat
        qi_loc[:,:] = np.where(qi_loc < 0., 0., qi_loc)
        qr_loc += qrtend * deltat
        qr_loc[:,:] = np.where(qr_loc < 0., 0., qr_loc)
        qs_loc += qstend * deltat
        qs_loc[:,:] = np.where(qs_loc < 0., 0., qs_loc)
        nc_loc += nctend * deltat
        nc_loc[:,:] = np.where(nc_loc > 1.e10, 1.e10, np.where(nc_loc < 1.e-12, 1.e-12, nc_loc))
        ni_loc += nitend * deltat
        ni_loc[:,:] = np.where(nc_loc > 1.e10, 1.e10, np.where(ni_loc < 1.e-12, 1.e-12, ni_loc))
        nr_loc += nrtend * deltat
        nr_loc[:,:] = np.where(nc_loc > 1.e10, 1.e10, np.where(nr_loc < 1.e-12, 1.e-12, nr_loc))
        ns_loc += nstend * deltat
        ns_loc[:,:] = np.where(nc_loc > 1.e10, 1.e10, np.where(ns_loc < 1.e-12, 1.e-12, ns_loc))

        if it == 0:
            for name in var_names:
                finals[name][n+1,:] = loc_arrays[name]
            rain_evap_fine[n,:] = prer_evap
            rain_sed_fine[n,:] = qrsedten
            rain_auto_fine[n,:] = prctot
            rain_accr_fine[n,:] = pratot
        else:
            for name in var_names:
                roughs[name][n+1,:] = loc_arrays[name]
            rain_evap_coarse[n,:] = prer_evap
            rain_sed_coarse[n,:] = qrsedten
            rain_auto_coarse[n,:] = prctot
            rain_accr_coarse[n,:] = pratot

        # Do something with final columns.

# Highest level to plot.
top = 0

levels = file.variables["ilev"][top:]
times_fine = np.linspace(0., final_time, finals['Q'].shape[0])
times_coarse = np.linspace(0., final_time, roughs['Q'].shape[0])

times_fine, lev_fine = np.meshgrid(times_fine, levels)
times_coarse, lev_coarse = np.meshgrid(times_coarse, levels)

assert timesteps[1] % timesteps[0] == 0
stride = timesteps[1] // timesteps[0]

cmap = plt.get_cmap('Greens')

for name in var_names:
    # This is to get the two plots to use the same color bar.
    # Note skipping the last values below, since they are not plotted.
    color_min = min(finals[name][:-1,top:].min(), roughs[name][:-1,top:].min())
    color_max = max(finals[name][:-1,top:].max(), roughs[name][:-1,top:].max())
    plt.pcolormesh(times_fine, lev_fine, finals[name].transpose()[top:,:], cmap=cmap)
    plt.xlabel("Time (s)")
    plt.ylabel("Reference pressure (mb)")
    plt.axis('tight')
    plt.gca().invert_yaxis()
    plt.clim(color_min, color_max)
    plt.colorbar()
    plt.savefig('./outlier_time{}s_{}.eps'.format(timesteps[0], name))
    plt.close()
    plt.pcolormesh(times_coarse, lev_coarse, roughs[name].transpose()[top:,:], cmap=cmap)
    plt.xlabel("Time (s)")
    plt.ylabel("Reference pressure (mb)")
    plt.axis('tight')
    plt.gca().invert_yaxis()
    plt.clim(color_min, color_max)
    plt.colorbar()
    plt.savefig('./outlier_time{}s_{}.eps'.format(timesteps[1], name))
    plt.close()
    plt.pcolormesh(times_coarse, lev_coarse, roughs[name].transpose()[top:,:] - finals[name][::stride,top:].transpose(), cmap=cmap)
    plt.xlabel("Time (s)")
    plt.ylabel("Reference pressure (mb)")
    plt.axis('tight')
    plt.gca().invert_yaxis()
    plt.colorbar()
    plt.savefig('./outlier_{}s-{}s_{}.eps'.format(timesteps[1], timesteps[0], name))
    plt.close()
