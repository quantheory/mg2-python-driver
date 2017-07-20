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

t_in_loc = np.empty((mgncol, t.shape[1]), order='F')
q_in_loc = np.empty((mgncol, q.shape[1]), order='F')
qc_in_loc = np.empty((mgncol, qc.shape[1]), order='F')
qi_in_loc = np.empty((mgncol, qi.shape[1]), order='F')
qr_in_loc = np.empty((mgncol, qr.shape[1]), order='F')
qs_in_loc = np.empty((mgncol, qs.shape[1]), order='F')

total_columns = 2048
final_time = 1800

use_col_num = 2

timesteps = np.array([5, 15, 30, 60, 120, 300, 900])
loc_arrays = {
    'T': t_loc,
    'Q': q_loc,
    'QC': qc_loc,
    'QI': qi_loc,
    'QR': qr_loc,
    'QS': qs_loc,
    'T_IN': t_in_loc,
    'Q_IN': q_in_loc,
    'QC_IN': qc_in_loc,
    'QI_IN': qi_in_loc,
    'QR_IN': qr_in_loc,
    'QS_IN': qs_in_loc,
}
var_names = sorted(list(loc_arrays.keys()))
out_vals = {}
out_fins = {}
finals = {}
final_means = {}
for name in var_names:
    out_vals[name] = []
    out_fins[name] = []
    for i in range(use_col_num):
        out_vals[name].append(np.zeros((lev, timesteps.size - 1)))
        out_fins[name].append(np.zeros((lev, timesteps.size - 1)))
    finals[name] = np.zeros((total_columns, lev))
    final_means[name] = np.zeros((lev,))

norm_diffs = [0.]*use_col_num

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

        t_in_loc[:,:] = t_loc
        q_in_loc[:,:] = q_loc
        qc_in_loc[:,:] = qc_loc
        qi_in_loc[:,:] = qi_loc
        qr_in_loc[:,:] = qr_loc
        qs_in_loc[:,:] = qs_loc

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
                finals[name][offset*mgncol:(offset+1)*mgncol,:] = loc_arrays[name]
        else:
            for i in range(mgncol):
                norm_diff = la.norm(pdel_loc[i,:] * (finals['Q'][offset*mgncol+i,:] - q_loc[i,:])) / gravit
                if norm_diff > norm_diffs[-1]:
                    norm_diffs.pop()
                    idx = -1
                    for j in range(use_col_num-1):
                        if norm_diff > norm_diffs[j]:
                            idx = j
                            break
                    norm_diffs.insert(idx, norm_diff)
                    for name in var_names:
                        for j in range(use_col_num-1, idx, -1):
                            out_vals[name][j][:,it-1] = out_vals[name][j-1][:,it-1]
                            out_fins[name][j][:,it-1] = out_fins[name][j-1][:,it-1]
                        out_vals[name][idx][:,it-1] = loc_arrays[name][i,:]
                        out_fins[name][idx][:,it-1] = finals[name][offset*mgncol+i,:]

        # Do something with final columns.

p_ref = file.variables['lev']

for name in var_names:
    final_means[name][:] = finals[name].mean(axis=0)

    for it in range(timesteps.size-1):
        mean_var = np.zeros((lev,))
        for j in range(use_col_num):
            mean_var += out_fins[name][j][:,it]
        mean_var /= use_col_num
        mean_var -= final_means[name]
        plt.plot(mean_var, p_ref, label='$\Delta$t={}'.format(timesteps[it+1]))
    plt.gca().invert_yaxis()
    plt.legend(loc='best')
    plt.xlabel('Anomaly in {}'.format(name))
    plt.ylabel('Reference pressure (hPa)')
    plt.axis('tight')
    plt.savefig('./outlier_t5_{}.eps'.format(name))
    plt.close()
