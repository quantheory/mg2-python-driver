#!/usr/bin/env python

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

HIST_FILE_NAME = "/g/g14/santos36/Data/MG2_data_collection.cam.h1.0001-01-06-00000.nc"
CLUSTER_FILE_NAME = "/g/g14/santos36/Data/MG2_data_collection.10_cluster_labels.0001-01-06-00000.nc"

hfile = nc4.Dataset(HIST_FILE_NAME, 'r')
cfile = nc4.Dataset(CLUSTER_FILE_NAME, 'r')

ncol = len(hfile.dimensions['ncol'])
lev = len(hfile.dimensions['lev'])
ilev = len(hfile.dimensions['ilev'])
ncluster = len(cfile.dimensions['ncluster'])

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

t = hfile.variables["MG2IN_T"]
q = hfile.variables["MG2IN_Q"]
qc = hfile.variables["MG2IN_QC"]
qi = hfile.variables["MG2IN_QI"]
nc = hfile.variables["MG2IN_NC"]
ni = hfile.variables["MG2IN_NI"]
qr = hfile.variables["MG2IN_QR"]
qs = hfile.variables["MG2IN_QS"]
nr = hfile.variables["MG2IN_NR"]
ns = hfile.variables["MG2IN_NS"]
relvar = hfile.variables["MG2IN_RELVAR"]
accre_enhan = hfile.variables["MG2IN_ACCRE_ENHAN"]
p = hfile.variables["MG2IN_P"]
pdel = hfile.variables["MG2IN_PDEL"]
precipf = hfile.variables["MG2IN_PRECIP"]
liqcldf = hfile.variables["MG2IN_LIQCLDF"]
icecldf = hfile.variables["MG2IN_ICECLDF"]
naai = hfile.variables["MG2IN_NAAI"]
npccn = hfile.variables["MG2IN_NPCCN"]
rndst = np.empty((t.shape[0], t.shape[1], t.shape[2], 4))
rndst[:,:,:,0] = hfile.variables["MG2IN_RNDST1"][:]
rndst[:,:,:,1] = hfile.variables["MG2IN_RNDST2"][:]
rndst[:,:,:,2] = hfile.variables["MG2IN_RNDST3"][:]
rndst[:,:,:,3] = hfile.variables["MG2IN_RNDST4"][:]
nacon = np.empty((t.shape[0], t.shape[1], t.shape[2], 4))
nacon[:,:,:,0] = hfile.variables["MG2IN_NACON1"][:]
nacon[:,:,:,1] = hfile.variables["MG2IN_NACON2"][:]
nacon[:,:,:,2] = hfile.variables["MG2IN_NACON3"][:]
nacon[:,:,:,3] = hfile.variables["MG2IN_NACON4"][:]
frzimm = hfile.variables["MG2IN_FRZIMM"]
frzcnt = hfile.variables["MG2IN_FRZCNT"]
frzdep = hfile.variables["MG2IN_FRZDEP"]

label = cfile.variables["label"]

t_loc = np.empty((1, t.shape[1]), order='F')
q_loc = np.empty((1, q.shape[1]), order='F')
qc_loc = np.empty((1, qc.shape[1]), order='F')
qi_loc = np.empty((1, qi.shape[1]), order='F')
nc_loc = np.empty((1, nc.shape[1]), order='F')
ni_loc = np.empty((1, ni.shape[1]), order='F')
qr_loc = np.empty((1, qr.shape[1]), order='F')
qs_loc = np.empty((1, qs.shape[1]), order='F')
nr_loc = np.empty((1, nr.shape[1]), order='F')
ns_loc = np.empty((1, ns.shape[1]), order='F')
relvar_loc = np.empty((1, relvar.shape[1]), order='F')
accre_enhan_loc = np.empty((1, accre_enhan.shape[1]), order='F')
p_loc = np.empty((1, p.shape[1]), order='F')
pdel_loc = np.empty((1, pdel.shape[1]), order='F')
precipf_loc = np.empty((1, precipf.shape[1]), order='F')
liqcldf_loc = np.empty((1, liqcldf.shape[1]), order='F')
icecldf_loc = np.empty((1, icecldf.shape[1]), order='F')
naai_loc = np.empty((1, naai.shape[1]), order='F')
npccn_loc = np.empty((1, npccn.shape[1]), order='F')
rndst_loc = np.empty((1, rndst.shape[1], 4), order='F')
nacon_loc = np.empty((1, nacon.shape[1], 4), order='F')
frzimm_loc = np.empty((1, frzimm.shape[1]), order='F')
frzcnt_loc = np.empty((1, frzcnt.shape[1]), order='F')
frzdep_loc = np.empty((1, frzdep.shape[1]), order='F')

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

# Dictionary associating "short" to "pretty" names for MG2 processes.
short_names = [
    "revap",
    "ssubl",
    "vdepo",
    "vnudp",
    "cbrgs",
    "cacws",
    "cnccc",
    "cncct",
    "cacwi",
    "rfrez",
    "racrs",
    "cbrgi",
    "cauto",
    "caccr",
    "iauto",
    "iaccr",
    "rnagg",
    "snagg",
    "anadj",
]
process_names = {
    "revap": "Rain Evap.",
    "ssubl": "Snow Subl.",
    "vdepo": "Vapor/Ice DMS",
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

def mg2_tendencies(level, t, q, qc, nc, qi, ni, qr, nr, qs, ns):
    """Extract the tendencies due to each process from an MG2 run.

    Returns a dictionary associating process names to an array of tendencies.
    """
    qcsinksum_rate1ord, tlat, qvlat, qctend, qitend, nctend, nitend, qrtend, \
        qstend, nrtend, nstend, effc, effc_fn, effi, prect, preci, nevapr, \
        evapsnow, prain, prodsnow, cmeout, deffi, pgamrad, lamcrad, qsout, dsout, \
        rflx, sflx, qrout, reff_rain, reff_snow, qcsevap, qisevap, qvres, cmeitot, \
        vtrmc, vtrmi, umr, ums, qcsedten, qisedten, qrsedten, qssedten, pratot, \
        prctot, mnuccctot, mnuccttot, msacwitot, psacwstot, bergstot, bergtot, \
        melttot, homotot, qcrestot, prcitot, praitot, qirestot, mnuccrtot, \
        pracstot, meltsdttot, frzrdttot, mnuccdtot, nrout, nsout, refl, arefl, \
        areflz, frefl, csrfl, acsrfl, fcsrfl, rercld, ncai, ncal, qrout2, qsout2, \
        nrout2, nsout2, drout2, dsout2, freqs, freqr, nfice, qcrat, \
        nadjtot, ncmeitot, mnudeptot, nnudeptot, npsacwstot, nnuccctot, nnuccttot, \
        nsacwitot, nnuccrtot, mnuccritot, nnuccritot, npracstot, npratot, nprctot, \
        nprc1tot, npraitot, nprcitot, nsubrtot, nraggtot, nsaggtot, subqr, subnr, \
        subdr, errstring, prer_evap \
        = mg.micro_mg_tend(evap_col_steps, evap_steps, col_steps, timestep, t[0,level],
                           q[0,level], qc[0,level], qi[0,level], nc[0,level],
                           ni[0,level], qr[0,level], qs[0,level], nr[0,level], ns[0,level],
                           relvar_loc[0,level], accre_enhan_loc[0,level], p_loc[0,level], pdel_loc[0,level],
                           precipf_loc[0,level], liqcldf_loc[0,level], icecldf_loc[0,level], naai_loc[0,level],
                           npccn_loc[0,level], rndst_loc[:,level:level+1,:], nacon_loc[:,level:level+1,:], precip_frac[0,level],
                           frzimm=frzimm_loc[0,level], frzcnt=frzcnt_loc[0,level],
                           frzdep=frzdep_loc[0,level], mgncol=1, nlev=1, do_sed=False,
                           do_inst=False)
    # Add total tendencies to a dictionary.
    tends = {}
    tends["Total"] = np.zeros((10,))
    tends["Total"][it] = tlat[0,0] / cpair
    tends["Total"][iq] = qvlat[0,0]
    tends["Total"][iqc] = qctend[0,0]
    tends["Total"][inc] = nctend[0,0]
    tends["Total"][iqi] = qitend[0,0]
    tends["Total"][ini] = nitend[0,0]
    tends["Total"][iqr] = qrtend[0,0]
    tends["Total"][inr] = nrtend[0,0]
    tends["Total"][iqs] = qstend[0,0]
    tends["Total"][ins] = nstend[0,0]
    # Divide up individual processes.
    for name in process_names:
        tends[name] = np.zeros((10,))
    # Rain evaporation.
    tends["revap"][it] = -prer_evap[0,0] * latvap / cpair
    tends["revap"][iq] = prer_evap[0,0]
    tends["revap"][iqr] = -prer_evap[0,0]
    tends["revap"][inr] = nsubrtot[0,0]
    # Snow sublimation.
    tends["ssubl"][it] = -evapsnow[0,0] * (latvap + latice) / cpair
    tends["ssubl"][iq] = evapsnow[0,0]
    tends["ssubl"][iqs] = -evapsnow[0,0]
    # [Snow number sublimation currently ignored.]
    # Vapor deposition.
    tends["vdepo"][it] = cmeitot[0,0] * (latvap + latice) / cpair
    tends["vdepo"][iq] = -cmeitot[0,0]
    tends["vdepo"][iqi] = cmeitot[0,0]
    tends["vdepo"][ini] = ncmeitot[0,0]
    # External scheme's nucleation deposition (as implemented in MG2).
    tends["vnudp"][it] = mnudeptot[0,0] * (latvap + latice) / cpair
    tends["vnudp"][iq] = -mnudeptot[0,0]
    tends["vnudp"][iqi] = mnudeptot[0,0]
    tends["vnudp"][ini] = nnudeptot[0,0]
    # Bergeron process over snow.
    tends["cbrgs"][it] = bergstot[0,0] * latice / cpair
    tends["cbrgs"][iqc] = -bergstot[0,0]
    tends["cbrgs"][iqs] = bergstot[0,0]
    # Collection of cloud water by snow.
    tends["cacws"][it] = psacwstot[0,0] * latice / cpair
    tends["cacws"][iqc] = -psacwstot[0,0]
    tends["cacws"][iqs] = psacwstot[0,0]
    tends["cacws"][inc] = -npsacwstot[0,0]
    # Immersion Freezing.
    tends["cnccc"][it] = mnuccctot[0,0] * latice / cpair
    tends["cnccc"][iqc] = -mnuccctot[0,0]
    tends["cnccc"][iqi] = mnuccctot[0,0]
    tends["cnccc"][inc] = -nnuccctot[0,0]
    if use_hetfrz_classnuc:
        tends["cnccc"][ini] = nnuccctot[0,0]
    # Contact Freezing.
    tends["cncct"][it] = mnuccttot[0,0] * latice / cpair
    tends["cncct"][iqc] = -mnuccttot[0,0]
    tends["cncct"][iqi] = mnuccttot[0,0]
    tends["cncct"][inc] = -nnuccttot[0,0]
    tends["cncct"][ini] = nnuccttot[0,0]
    # Secondary ice production.
    tends["cacwi"][it] = msacwitot[0,0] * latice / cpair
    tends["cacwi"][iqc] = -msacwitot[0,0]
    tends["cacwi"][iqi] = msacwitot[0,0]
    tends["cacwi"][ini] = nsacwitot[0,0]
    # Freezing of rain.
    tends["rfrez"][it] = (mnuccrtot[0,0] + mnuccritot[0,0]) * latice / cpair
    tends["rfrez"][iqi] = mnuccritot[0,0]
    tends["rfrez"][iqr] = - mnuccrtot[0,0] - mnuccritot[0,0]
    tends["rfrez"][iqs] = mnuccrtot[0,0]
    tends["rfrez"][ini] = nnuccritot[0,0]
    tends["rfrez"][inr] = - nnuccrtot[0,0] - nnuccritot[0,0]
    tends["rfrez"][ins] = nnuccrtot[0,0]
    # Accretion of rain onto snow.
    tends["racrs"][it] = pracstot[0,0] * latice / cpair
    tends["racrs"][iqr] = -pracstot[0,0]
    tends["racrs"][iqs] = pracstot[0,0]
    tends["racrs"][inr] = -npracstot[0,0]
    # Bergeron process over cloud ice.
    tends["cbrgi"][it] = bergtot[0,0] * latice / cpair
    tends["cbrgi"][iqc] = -bergtot[0,0]
    tends["cbrgi"][iqi] = bergtot[0,0]
    # Autoconversion.
    tends["cauto"][iqc] = -prctot[0,0]
    tends["cauto"][iqr] = prctot[0,0]
    tends["cauto"][inc] = -nprc1tot[0,0]
    tends["cauto"][inr] = nprctot[0,0]
    # Accretion.
    tends["caccr"][iqc] = -pratot[0,0]
    tends["caccr"][iqr] = pratot[0,0]
    tends["caccr"][inc] = -npratot[0,0]
    # Ice autoconversion.
    tends["iauto"][iqi] = -prcitot[0,0]
    tends["iauto"][iqs] = prcitot[0,0]
    tends["iauto"][ini] = -nprcitot[0,0]
    tends["iauto"][ins] = nprcitot[0,0]
    # Ice accretion.
    tends["iaccr"][iqi] = -praitot[0,0]
    tends["iaccr"][iqs] = praitot[0,0]
    tends["iaccr"][ini] = -npraitot[0,0]
    # Rain self-collection.
    tends["rnagg"][inr] = nraggtot[0,0]
    # Snow self-collection.
    tends["snagg"][ins] = nsaggtot[0,0]
    # Droplet activation.
    tends["cnact"][inc] = npccn_loc[0,level]
    # Number adjustment.
    tends["anadj"][inc] = nadjtot[0,0,0]
    tends["anadj"][ini] = nadjtot[0,0,1]
    tends["anadj"][inr] = nadjtot[0,0,2]
    tends["anadj"][ins] = nadjtot[0,0,3]
    # Check the above now.
    # TODO: fix substepped code process rates so that this works again.
    #budget_total = np.zeros((10,))
    #for name in process_names:
    #    budget_total += tends[name]
    #assert np.allclose(budget_total, tends["Total"])
    return (tends, drout2)

def update_state(level, tends, deltat):
    t_loc[0,level] += tends["Total"][it] * deltat
    q_loc[0,level] += tends["Total"][iq] * deltat
    q_loc[0,level] = max(1.e-12, q_loc[0,level])
    qc_loc[0,level] += tends["Total"][iqc] * deltat
    qc_loc[0,level] = max(0., qc_loc[0,level])
    qi_loc[0,level] += tends["Total"][iqi] * deltat
    qi_loc[0,level] = max(0., qi_loc[0,level])
    qr_loc[0,level] += tends["Total"][iqr] * deltat
    qr_loc[0,level] = max(0., qr_loc[0,level])
    qs_loc[0,level] += tends["Total"][iqs] * deltat
    qs_loc[0,level] = max(0., qs_loc[0,level])
    nc_loc[0,level] += tends["Total"][inc] * deltat
    nc_loc[0,level] = max(1.e-12, nc_loc[0,level])
    ni_loc[0,level] += tends["Total"][ini] * deltat
    ni_loc[0,level] = max(1.e-12, ni_loc[0,level])
    nr_loc[0,level] += tends["Total"][inr] * deltat
    nr_loc[0,level] = max(1.e-12, nr_loc[0,level])
    ns_loc[0,level] += tends["Total"][ins] * deltat
    ns_loc[0,level] = max(1.e-12, ns_loc[0,level])

def reset_state(level, column):
    t_loc[0,level] = t[0,level,column]
    q_loc[0,level] = q[0,level,column]
    qc_loc[0,level] = qc[0,level,column]
    qi_loc[0,level] = qi[0,level,column]
    nc_loc[0,level] = nc[0,level,column]
    ni_loc[0,level] = ni[0,level,column]
    qr_loc[0,level] = qr[0,level,column]
    qs_loc[0,level] = qs[0,level,column]
    nr_loc[0,level] = nr[0,level,column]
    ns_loc[0,level] = ns[0,level,column]

def diff_states(level, column):
    diffs = np.zeros((10,))
    diffs[it] = t_loc[0,level] - t[0,level,column]
    diffs[iq] = q_loc[0,level] - q[0,level,column]
    diffs[iqc] = qc_loc[0,level] - qc[0,level,column]
    diffs[inc] = nc_loc[0,level] - nc[0,level,column]
    diffs[iqi] = qi_loc[0,level] - qi[0,level,column]
    diffs[ini] = ni_loc[0,level] - ni[0,level,column]
    diffs[iqr] = qr_loc[0,level] - qr[0,level,column]
    diffs[inr] = nr_loc[0,level] - nr[0,level,column]
    diffs[iqs] = qs_loc[0,level] - qs[0,level,column]
    diffs[ins] = ns_loc[0,level] - ns[0,level,column]
    return diffs

def calc_twmt(diff1, diff2, deltat):
    diff_tends = (diff2 - diff1) / deltat
    return (abs(diff_tends[iq]) + abs(diff_tends[iqc]) + abs(diff_tends[iqi]) +
                 abs(diff_tends[iqr]) + abs(diff_tends[iqs])) / 2.

ind = np.arange(len(short_names))

num_columns = 2048

plt.autoscale(tight=True)

mean_1s_twmt = 0.
mean_300s_twmt = 0.
mean_together_twmt = 0.
mean_evap_twmt = 0.
mean_hybrid_twmt = 0.
mean_1s300s_twmt = 0.
mean_300stogether_twmt = 0.
mean_300shybrid_twmt = 0.
mean_300sevap_twmt = 0.
mean_300scol_twmt = 0.
mean_300sapart_twmt = 0.
mean_1stogether_twmt = 0.
mean_1shybrid_twmt = 0.
mean_1sevap_twmt = 0.
mean_1scol_twmt = 0.
mean_1sapart_twmt = 0.

mean_drout_1s = 0.
mean_drout_300s = 0.
mean_drout_together = 0.
mean_drout_hybrid = 0.
mean_drout_evap = 0.
mean_drout_col = 0.
mean_drout_apart = 0.

mean_nr_1s = 0.
mean_nr_300s = 0.
mean_nr_together = 0.
mean_nr_hybrid = 0.
mean_nr_evap = 0.
mean_nr_col = 0.
mean_nr_apart = 0.

tend_norms = []

number_grid_cells = 0

for column in range(num_columns):
    print("On column: ", column)
    t_loc[:,:] = t[0,:,column]
    q_loc[:,:] = q[0,:,column]
    qc_loc[:,:] = qc[0,:,column]
    qi_loc[:,:] = qi[0,:,column]
    nc_loc[:,:] = nc[0,:,column]
    ni_loc[:,:] = ni[0,:,column]
    qr_loc[:,:] = qr[0,:,column]
    qs_loc[:,:] = qs[0,:,column]
    nr_loc[:,:] = nr[0,:,column]
    ns_loc[:,:] = ns[0,:,column]
    relvar_loc[:,:] = relvar[0,:,column]
    accre_enhan_loc[:,:] = accre_enhan[0,:,column]
    p_loc[:,:] = p[0,:,column]
    pdel_loc[:,:] = pdel[0,:,column]
    precipf_loc[:,:] = precipf[0,:,column]
    liqcldf_loc[:,:] = liqcldf[0,:,column]
    icecldf_loc[:,:] = icecldf[0,:,column]
    naai_loc[:,:] = naai[0,:,column]
    npccn_loc[:,:] = npccn[0,:,column]
    rndst_loc[:,:,:] = rndst[0,:,column:column+1,:].transpose([1, 0, 2])
    nacon_loc[:,:,:] = nacon[0,:,column:column+1,:].transpose([1, 0, 2])
    frzimm_loc[:,:] = frzimm[0,:,column]
    frzcnt_loc[:,:] = frzcnt[0,:,column]
    frzdep_loc[:,:] = frzdep[0,:,column]
    precip_frac = mg.calc_precip_frac(qc_loc, qi_loc, qr_loc,
                                      qs_loc, precipf_loc, liqcldf_loc,
                                      icecldf_loc, mgncol=1, nlev=lev)
    for level in range(lev):
        c = label[0,level,column]
        if c != 9:
            continue
        number_grid_cells += 1
        timestep = 1.
        evap_col_steps = 1
        evap_steps = 1
        col_steps = 1
        tot_t_1s = np.zeros((10,))
        for m in range(300):
            tends, drout_1s = mg2_tendencies(level, t_loc, q_loc, qc_loc, nc_loc, qi_loc, ni_loc,
                                             qr_loc, nr_loc, qs_loc, ns_loc)
            tot_t_1s += tends["Total"]
            update_state(level, tends, timestep)
        mean_nr_1s += nr_loc[0,level]
        state_diff_1s = diff_states(level, column)
        reset_state(level, column)
        timestep = 300.
        evap_col_steps = 1
        evap_steps = 1
        col_steps = 1
        tends, drout_300s = mg2_tendencies(level, t_loc, q_loc, qc_loc, nc_loc, qi_loc, ni_loc,
                                           qr_loc, nr_loc, qs_loc, ns_loc)
        update_state(level, tends, timestep)
        mean_nr_300s += nr_loc[0,level]
        state_diff_300s = diff_states(level, column)
        reset_state(level, column)
        timestep = 300.
        evap_col_steps = 300
        evap_steps = 1
        col_steps = 1
        tends, drout_together = mg2_tendencies(level, t_loc, q_loc, qc_loc, nc_loc, qi_loc, ni_loc,
                                               qr_loc, nr_loc, qs_loc, ns_loc)
        update_state(level, tends, timestep)
        mean_nr_together += nr_loc[0,level]
        state_diff_substep_together = diff_states(level, column)
        reset_state(level, column)
        timestep = 300.
        evap_col_steps = 10
        evap_steps = 30
        col_steps = 1
        tends, drout_hybrid = mg2_tendencies(level, t_loc, q_loc, qc_loc, nc_loc, qi_loc, ni_loc,
                                             qr_loc, nr_loc, qs_loc, ns_loc)
        update_state(level, tends, timestep)
        mean_nr_hybrid += nr_loc[0,level]
        state_diff_substep_hybrid = diff_states(level, column)
        reset_state(level, column)
        timestep = 300.
        evap_col_steps = 1
        evap_steps = 300
        col_steps = 1
        tends, drout_evap = mg2_tendencies(level, t_loc, q_loc, qc_loc, nc_loc, qi_loc, ni_loc,
                                           qr_loc, nr_loc, qs_loc, ns_loc)
        update_state(level, tends, timestep)
        mean_nr_evap += nr_loc[0,level]
        state_diff_substep_evap = diff_states(level, column)
        reset_state(level, column)
        timestep = 300.
        evap_col_steps = 1
        evap_steps = 1
        col_steps = 300
        tends, drout_col = mg2_tendencies(level, t_loc, q_loc, qc_loc, nc_loc, qi_loc, ni_loc,
                                          qr_loc, nr_loc, qs_loc, ns_loc)
        update_state(level, tends, timestep)
        mean_nr_col += nr_loc[0,level]
        state_diff_substep_col = diff_states(level, column)
        reset_state(level, column)
        timestep = 300.
        evap_col_steps = 1
        evap_steps = 300
        col_steps = 300
        tends, drout_apart = mg2_tendencies(level, t_loc, q_loc, qc_loc, nc_loc, qi_loc, ni_loc,
                                            qr_loc, nr_loc, qs_loc, ns_loc)
        update_state(level, tends, timestep)
        mean_nr_apart += nr_loc[0,level]
        state_diff_substep_apart = diff_states(level, column)
        reset_state(level, column)
        
        mean_1s_twmt += calc_twmt(0., state_diff_1s, 300.)
        mean_300s_twmt += calc_twmt(0., state_diff_300s, 300.)
        mean_together_twmt += calc_twmt(0., state_diff_substep_together, 300.)
        mean_hybrid_twmt += calc_twmt(0., state_diff_substep_hybrid, 300.)
        mean_evap_twmt += calc_twmt(0., state_diff_substep_evap, 300.)
        mean_1s300s_twmt += calc_twmt(state_diff_300s, state_diff_1s, 300.)
        mean_300stogether_twmt += calc_twmt(state_diff_300s,
                                            state_diff_substep_together, 300.)
        mean_300shybrid_twmt += calc_twmt(state_diff_300s,
                                            state_diff_substep_hybrid, 300.)
        mean_300sevap_twmt += calc_twmt(state_diff_300s,
                                        state_diff_substep_evap, 300.)
        mean_300scol_twmt += calc_twmt(state_diff_300s,
                                       state_diff_substep_col, 300.)
        mean_300sapart_twmt += calc_twmt(state_diff_300s,
                                         state_diff_substep_apart, 300.)
        mean_1stogether_twmt += calc_twmt(state_diff_1s,
                                          state_diff_substep_together, 300.)
        mean_1shybrid_twmt += calc_twmt(state_diff_1s,
                                        state_diff_substep_hybrid, 300.)
        mean_1sevap_twmt += calc_twmt(state_diff_1s,
                                      state_diff_substep_evap, 300.)
        mean_1scol_twmt += calc_twmt(state_diff_1s,
                                     state_diff_substep_col, 300.)
        mean_1sapart_twmt += calc_twmt(state_diff_1s,
                                       state_diff_substep_apart, 300.)
        mean_drout_1s += drout_1s
        mean_drout_300s += drout_300s
        mean_drout_together += drout_together
        mean_drout_hybrid += drout_hybrid
        mean_drout_apart += drout_apart
        mean_drout_col += drout_col
        mean_drout_evap += drout_evap

mean_1s_twmt *= 1.e3 / number_grid_cells
mean_300s_twmt *= 1.e3 / number_grid_cells
mean_together_twmt *= 1.e3 / number_grid_cells
mean_hybrid_twmt *= 1.e3 / number_grid_cells
mean_evap_twmt *= 1.e3 / number_grid_cells
mean_1s300s_twmt *= 1.e3 / number_grid_cells
mean_300stogether_twmt *= 1.e3 / number_grid_cells
mean_300shybrid_twmt *= 1.e3 / number_grid_cells
mean_300sevap_twmt *= 1.e3 / number_grid_cells
mean_300scol_twmt *= 1.e3 / number_grid_cells
mean_300sapart_twmt *= 1.e3 / number_grid_cells
mean_1stogether_twmt *= 1.e3 / number_grid_cells
mean_1shybrid_twmt *= 1.e3 / number_grid_cells
mean_1sevap_twmt *= 1.e3 / number_grid_cells
mean_1scol_twmt *= 1.e3 / number_grid_cells
mean_1sapart_twmt *= 1.e3 / number_grid_cells

mean_drout_1s *= 1.e6 / number_grid_cells
mean_drout_300s *= 1.e6 / number_grid_cells
mean_drout_together *= 1.e6 / number_grid_cells
mean_drout_hybrid *= 1.e6 / number_grid_cells
mean_drout_apart *= 1.e6 / number_grid_cells
mean_drout_col *= 1.e6 / number_grid_cells
mean_drout_evap *= 1.e6 / number_grid_cells

mean_nr_1s *= 1. / number_grid_cells
mean_nr_300s *= 1. / number_grid_cells
mean_nr_together *= 1. / number_grid_cells
mean_nr_hybrid *= 1. / number_grid_cells
mean_nr_apart *= 1. / number_grid_cells
mean_nr_col *= 1. / number_grid_cells
mean_nr_evap *= 1. / number_grid_cells

print("mean 1s twmt:", mean_1s_twmt)
print("mean 300s twmt:", mean_300s_twmt)
print("mean together twmt:", mean_together_twmt)
print("mean hybrid twmt:", mean_hybrid_twmt)
print("mean evap twmt:", mean_evap_twmt)
print("mean 300s - 1s twmt diff:", mean_1s300s_twmt)
print("mean together - 300s twmt diff:", mean_300stogether_twmt)
print("mean hybrid - 300s twmt diff:", mean_300shybrid_twmt)
print("mean evap - 300s twmt diff:", mean_300sevap_twmt)
print("mean col - 300s twmt diff:", mean_300scol_twmt)
print("mean apart - 300s twmt diff:", mean_300sapart_twmt)
print("mean together - 1s twmt diff:", mean_1stogether_twmt)
print("mean hybrid - 1s twmt diff:", mean_1shybrid_twmt)
print("mean evap - 1s twmt diff:", mean_1sevap_twmt)
print("mean col - 1s twmt diff:", mean_1scol_twmt)
print("mean apart - 1s twmt diff:", mean_1sapart_twmt)

print("mean 1s drout:", mean_drout_1s)
print("mean 300s drout:", mean_drout_300s)
print("mean together drout:", mean_drout_together)
print("mean hybrid drout:", mean_drout_hybrid)
print("mean apart drout:", mean_drout_apart)
print("mean col drout:", mean_drout_col)
print("mean evap drout:", mean_drout_evap)

print("mean 1s nr:", mean_nr_1s)
print("mean 300s nr:", mean_nr_300s)
print("mean together nr:", mean_nr_together)
print("mean hybrid nr:", mean_nr_hybrid)
print("mean apart nr:", mean_nr_apart)
print("mean col nr:", mean_nr_col)
print("mean evap nr:", mean_nr_evap)
