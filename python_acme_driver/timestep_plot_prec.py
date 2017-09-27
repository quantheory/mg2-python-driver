#!/usr/bin/env python

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import netCDF4 as nc4

from mg2 import wv_sat_methods as wsm
from mg2 import micro_mg2_0 as mg

HIST_FILE_NAME = "/global/homes/s/santos/project/MG2_data_collection/run/MG2_data_collection.cam.h1.0001-01-06-00000.nc"

file = nc4.Dataset(HIST_FILE_NAME, 'r')

ncol = len(file.dimensions['ncol'])
lev = len(file.dimensions['lev'])
ilev = len(file.dimensions['ilev'])

kind = 8
tmelt = 273.15
h2otrip = 273.16
tboil = 373.16
ttrice = 20.
mwwv = 18.016
mwdair = 28.966
epsilo = mwwv / mwdair
gravit = 9.80616
boltz = 1.38065e-23
avogad = 6.02214e26
rgas = boltz * avogad
rair = rgas / mwdair
rh2o = rgas / mwwv
cpair = 1.00464e3
latvap = 2.501e6
latice = 3.337e5
rhmini = 0.8
dcs = 195.e-6
dcs_tdep = True
uniform = False
do_cldice = True
use_hetfrz_classnuc = True
precip_frac_method = 'in_cloud'
berg_eff_factor = 0.1
allow_sed_supersat = False
ice_sed_ai = 500.
prc_coef1 = 30500.
prc_exp = 3.19
prc_exp1 = -1.2
cld_sed = 1.
mg_prc_coeff_fix = True


errstring = wsm.wv_sat_methods_init(kind, tmelt, h2otrip, tboil, ttrice, epsilo)
if str(errstring).strip() != '':
    print("wv_sat_methods initialization error: ", errstring)

errstring = mg.micro_mg_init(kind, gravit, rair, rh2o, cpair, tmelt, latvap,
                             latice, rhmini, dcs, dcs_tdep, uniform, do_cldice,
                             use_hetfrz_classnuc, precip_frac_method,
                             berg_eff_factor, allow_sed_supersat, ice_sed_ai,
                             prc_coef1, prc_exp, prc_exp1, cld_sed,
                             mg_prc_coeff_fix)
if str(errstring).strip() != '':
    print("MG2 initialization error: ", errstring)

mgncol = 128
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

prect_loc = np.empty((mgncol,), order='F')
preci_loc = np.empty((mgncol,), order='F')

total_columns = 48512
final_time = 1800

timesteps = np.array([1, 2, 5, 15, 30, 60, 90, 120, 180, 300, 360, 450, 600, 900, 1800])
loc_arrays = {
#    'T': t_loc,
#    'Q': q_loc,
#    'QC': qc_loc,
#    'QI': qi_loc,
#    'QR': qr_loc,
#    'QS': qs_loc,
    'PRECT': prect_loc,
    'PRECI': preci_loc,
}
var_names = sorted(list(loc_arrays.keys()))
norms = {}
diffs = {}
finals = {}
for name in var_names:
    norms[name] = np.zeros((total_columns, timesteps.size - 1))
    diffs[name] = np.zeros((total_columns, timesteps.size - 1))
    finals[name] = np.zeros((total_columns,))

for it in range(timesteps.size):
    assert final_time % timesteps[it] == 0
    nsteps = final_time / timesteps[it]
    deltat = float(timesteps[it])

    for offset in range(total_columns / mgncol):
        t_loc[:,:] = t[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        q_loc[:,:] = q[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        qc_loc[:,:] = qc[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        qi_loc[:,:] = qi[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        nc_loc[:,:] = nc[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        ni_loc[:,:] = ni[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        qr_loc[:,:] = qr[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        qs_loc[:,:] = qs[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        nr_loc[:,:] = nr[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        ns_loc[:,:] = ns[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        relvar_loc[:,:] = relvar[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        accre_enhan_loc[:,:] = accre_enhan[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        p_loc[:,:] = p[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        pdel_loc[:,:] = pdel[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        precipf_loc[:,:] = precipf[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        liqcldf_loc[:,:] = liqcldf[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        icecldf_loc[:,:] = icecldf[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        naai_loc[:,:] = naai[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        npccn_loc[:,:] = npccn[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        rndst_loc[:,:,:] = rndst[0,:,offset*mgncol:(offset+1)*mgncol,:].transpose([1, 0, 2])
        nacon_loc[:,:,:] = nacon[0,:,offset*mgncol:(offset+1)*mgncol,:].transpose([1, 0, 2])
        frzimm_loc[:,:] = frzimm[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        frzcnt_loc[:,:] = frzcnt[0,:,offset*mgncol:(offset+1)*mgncol].transpose()
        frzdep_loc[:,:] = frzdep[0,:,offset*mgncol:(offset+1)*mgncol].transpose()

        prect_loc[:] = 0.
        preci_loc[:] = 0.

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

            prect_loc += prect * deltat
            preci_loc += preci * deltat

        if it == 0:
            for name in var_names:
                finals[name][offset*mgncol:(offset+1)*mgncol] = loc_arrays[name]
        else:
            for name in var_names:
                diffs[name][offset*mgncol:(offset+1)*mgncol,it-1] \
                    = finals[name][offset*mgncol:(offset+1)*mgncol] - loc_arrays[name]

        # Do something with final columns.

for name in var_names:
    norms[name][:,:] = np.abs(diffs[name][:,:]) + 1.e-80

for name in var_names:
    medians = np.median(diffs[name], axis=0)
    means = diffs[name].mean(axis=0)
    plt.semilogx(timesteps[1:], medians, label='median')
    plt.semilogx(timesteps[1:], means, label='mean')
    plt.semilogx(timesteps[1:], diffs[name].max(axis=0), label='max')
    plt.semilogx(timesteps[1:], diffs[name].min(axis=0), label='min')
    plt.legend(loc='best')
    plt.xlabel('Seconds/timestep')
    plt.ylabel('Difference in {} from $\Delta$t={}'.format(name, timesteps[0]))
    plt.axis('tight')
    plt.savefig('./timesteps_{}_signed.eps'.format(name))
    plt.close()

    medians = np.median(norms[name], axis=0)
    means = norms[name].mean(axis=0)
    plt.loglog(timesteps[1:], medians, label='median')
    plt.loglog(timesteps[1:], means, label='mean')
    plt.loglog(timesteps[1:], norms[name].max(axis=0), label='max')
    plt.legend(loc='best')
    plt.xlabel('Seconds/timestep')
    plt.ylabel('Absolute value of difference in {} from $\Delta$t={}'.format(name, timesteps[0]))
    plt.axis('tight')
    plt.savefig('./timesteps_{}.eps'.format(name))
    plt.close()
    time_log_diff = np.log(timesteps[3]) - np.log(timesteps[1])
    print("Estimated median convergence rate for variable {}: {}".format(name, (np.log(medians[2]) - np.log(medians[0]))/time_log_diff))
    print("Estimated mean convergence rate for variable {}: {}".format(name, (np.log(means[2]) - np.log(means[0]))/time_log_diff))

    percent_5s = np.percentile(norms[name], 95., axis=0)
    norms_clipped = norms[name][:,:]
    for i in range(timesteps.size - 1):
        norms_clipped[:,i] = np.where(norms[name][:,i] > percent_5s[i], percent_5s[i], norms[name][:,i])
    medians_clipped = np.median(norms[name], axis=0)
    assert (medians_clipped == medians).all()
    means_clipped = norms[name].mean(axis=0)
    plt.loglog(timesteps[1:], medians, label='median')
    plt.loglog(timesteps[1:], means, label='mean')
    plt.loglog(timesteps[1:], norms_clipped.max(axis=0), label='max')
    plt.legend(loc='best')
    plt.xlabel('Seconds/timestep')
    plt.ylabel('Absolute value of difference in {} from $\Delta$t={}'.format(name, timesteps[0]))
    plt.axis('tight')
    plt.savefig('./timesteps_{}_clipped.eps'.format(name))
    plt.close()
    time_log_diff = np.log(timesteps[4]) - np.log(timesteps[2])
    print("Estimated clipped median convergence rate for variable {}: {}".format(name, (np.log(medians_clipped[3]) - np.log(medians_clipped[1]))/time_log_diff))
    print("Estimated clipped mean convergence rate for variable {}: {}".format(name, (np.log(means_clipped[3]) - np.log(means_clipped[1]))/time_log_diff))
