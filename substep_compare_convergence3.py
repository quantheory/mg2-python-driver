#!/usr/bin/env python

import csv

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

import accumulator as acc

matplotlib.rcParams['lines.linewidth'] = 2.0

HIST_FILE_NAME = "/home/santos/Data/MG2_data_collection.cam.h1.0001-01-06-00000.nc"
CLUSTER_FILE_NAME = "/home/santos/Data/MG2_data_collection.10_cluster_labels.0001-01-06-00000.nc"

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
        = mg.micro_mg_tend(evap_col_steps, evap_steps, col_steps, auto_accr_steps, auto_steps, accr_steps, timestep, t[0,level],
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
    budget_total = np.zeros((10,))
    for name in process_names:
        budget_total += tends[name]
    assert np.allclose(budget_total, tends["Total"])
    water_mass_change = tends["Total"][iq] + tends["Total"][iqc] + \
                        tends["Total"][iqi] + tends["Total"][iqr] + \
                        tends["Total"][iqs]
    assert np.allclose(0., water_mass_change)
    return (tends, drout2, subqr, subnr, subdr)

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

flag = False
ice_flag = False

def run_mg2_substepped(num_steps):
    # Icky global usage.
    global timestep
    timestep = 300./num_steps
    global flag
    flag = False
    global ice_flag
    ice_flag = False
    for m in range(num_steps):
        tends, drout, subqr, subnr, subdr = mg2_tendencies(level, t_loc, q_loc, qc_loc, nc_loc, qi_loc, ni_loc,
                                                              qr_loc, nr_loc, qs_loc, ns_loc)
        other_tends = tends["Total"] - tends["anadj"] - tends["cnact"] - tends["cauto"] - tends["caccr"]# - tends["rnagg"] - tends["revap"]
        if abs(other_tends[iqi]) > 1.e-10 or abs(other_tends[iqs]) > 1.e-10:
            ice_flag = True
        if abs(tends["revap"][iqr]) > 1.e-10:
            flag = True
        update_state(level, tends, timestep)
    #drout = (qr_loc[0,level]/(np.pi * 1000. * nr_loc[0,level]))**(1./3.)
    state_diff = diff_states(level, column)
    reset_state(level, column)
    return state_diff

ind = np.arange(len(short_names))

#num_columns = 48602
num_columns = 1024
cluster = 6
max_power = 11

plt.autoscale(tight=True)

mean_twnt_ref = 0.
mean_delqr_ref = 0.
mean_twnt_coarse = 0.
mean_delqr_coarse = 0.

coarse_diff_acc = acc.Accumulator([acc.mean, acc.median, acc.max])
coarse_dqr_acc = acc.Accumulator([acc.mean, acc.median, acc.max])

make_accs = lambda n: [acc.Accumulator([acc.mean, acc.median, acc.max])
                       for i in range(n)]

all_diff_acc = make_accs(max_power-1)
together_diff_acc = make_accs(max_power-1)
apart_diff_acc = make_accs(max_power-1)
auto_diff_acc = make_accs(max_power-1)
accr_diff_acc = make_accs(max_power-1)

all_dqr_acc = make_accs(max_power-1)
together_dqr_acc = make_accs(max_power-1)
apart_dqr_acc = make_accs(max_power-1)
auto_dqr_acc = make_accs(max_power-1)
accr_dqr_acc = make_accs(max_power-1)

coarse_qrs = []

subsample_file = open('subsample_list_icefilter.csv', 'w', encoding='ascii', newline='')
subsample_writer = csv.writer(subsample_file)
subsample_writer.writerow(['column', 'level'])

number_grid_cells = 0

evap_col_steps = 1
evap_steps = 1
col_steps = 1
for column in range(num_columns):
    if column % 10 == 0:
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
        if c != cluster:
            continue
        #print('precip_frac is ', precip_frac[0,level])
        #print('t is ', t_loc[0,level])
        #print('p is ', p_loc[0,level])
        #print('qr is ', qr_loc[0,level])
        #print('nr is ', nr_loc[0,level])
        #print('q is ', q_loc[0,level])
        #print('liqcldf is ', liqcldf_loc[0,level])
        auto_accr_steps = 1
        auto_steps = 1
        accr_steps = 1
        state_diff_ref = run_mg2_substepped(2**max_power)
        if ice_flag:
            continue
        mean_twnt_ref += calc_twmt(state_diff_ref, 0., 300.)
        mean_delqr_ref += state_diff_ref[iqr] / 300.
        state_diff_coarse = run_mg2_substepped(1)
        if flag:
            print("Evap rate significant for coarse results.")
            print("column, level: ", column, ",", level)
        state_diff_coarse_ish = run_mg2_substepped(2)
        if flag:
            print("Evap rate significant for coarse-ish results.")
            print("column, level: ", column, ",", level)
        number_grid_cells += 1
        subsample_writer.writerow([str(column), str(level)])
        mean_twnt_coarse += calc_twmt(state_diff_coarse, 0., 300.)
        mean_delqr_coarse += state_diff_coarse[iqr] / 300.
        coarse_qrs.append(1.e3*(state_diff_coarse[iqr] + qr_loc[0,level]))
        state_diffs_all = np.zeros((max_power-1, 10))
        state_diffs_all[0,:] = state_diff_coarse_ish
        for i in range(2, max_power):
            state_diffs_all[i-1,:] = run_mg2_substepped(2**i)
        coarse_error = calc_twmt(state_diff_coarse, state_diff_ref, 300.)
        fine_error = calc_twmt(state_diffs_all[0,:], state_diff_ref, 300.)
        if coarse_error < fine_error:
            print("150 second error greater than 300 second error")
            print("column, level: ", column, ",", level)
            print("coarse, fine: ", coarse_error, ",", fine_error)
            print(state_diff_coarse)
            print(state_diffs_all[0,:])
        state_diffs_together = np.zeros((max_power-1, 10))
        for i in range(1, max_power):
            auto_accr_steps = 2**i
            state_diffs_together[i-1,:] = run_mg2_substepped(1)
        auto_accr_steps = 1
        state_diffs_apart = np.zeros((max_power-1, 10))
        for i in range(1, max_power):
            auto_steps = 2**i
            accr_steps = 2**i
            state_diffs_apart[i-1,:] = run_mg2_substepped(1)
        auto_steps = 1
        accr_steps = 1
        state_diffs_auto = np.zeros((max_power-1, 10))
        for i in range(1, max_power):
            auto_steps = 2**i
            state_diffs_auto[i-1,:] = run_mg2_substepped(1)
        auto_steps = 1
        state_diffs_accr = np.zeros((max_power-1, 10))
        for i in range(1, max_power):
            accr_steps = 2**i
            state_diffs_accr[i-1,:] = run_mg2_substepped(1)
        accr_steps = 1

        coarse_diff_acc.push(calc_twmt(state_diff_coarse, state_diff_ref, 300.))
        for i in range(max_power-1):
            all_diff_acc[i].push(calc_twmt(state_diffs_all[i,:], state_diff_ref, 300.))
            together_diff_acc[i].push(calc_twmt(state_diffs_together[i,:], state_diff_ref, 300.))
            apart_diff_acc[i].push(calc_twmt(state_diffs_apart[i,:], state_diff_ref, 300.))
            auto_diff_acc[i].push(calc_twmt(state_diffs_auto[i,:], state_diff_ref, 300.))
            accr_diff_acc[i].push(calc_twmt(state_diffs_accr[i,:], state_diff_ref, 300.))

        coarse_dqr_acc.push(np.abs(state_diff_coarse[iqr] - state_diff_ref[iqr]) / 300.)
        for i in range(max_power-1):
            all_dqr_acc[i].push(np.abs(state_diffs_all[i,iqr] - state_diff_ref[iqr]) / 300.)
            together_dqr_acc[i].push(np.abs(state_diffs_together[i,iqr] - state_diff_ref[iqr]) / 300.)
            apart_dqr_acc[i].push(np.abs(state_diffs_apart[i,iqr] - state_diff_ref[iqr]) / 300.)
            auto_dqr_acc[i].push(np.abs(state_diffs_auto[i,iqr] - state_diff_ref[iqr]) / 300.)
            accr_dqr_acc[i].push(np.abs(state_diffs_accr[i,iqr] - state_diff_ref[iqr]) / 300.)

subsample_file.close()

mean_all_diffs = np.zeros((max_power,))
mean_together_diffs = np.zeros((max_power,))
mean_apart_diffs = np.zeros((max_power,))
mean_auto_diffs = np.zeros((max_power,))
mean_accr_diffs = np.zeros((max_power,))

median_all_diffs = np.zeros((max_power,))
median_together_diffs = np.zeros((max_power,))
median_apart_diffs = np.zeros((max_power,))
median_auto_diffs = np.zeros((max_power,))
median_accr_diffs = np.zeros((max_power,))

max_all_diffs = np.zeros((max_power,))
max_together_diffs = np.zeros((max_power,))
max_apart_diffs = np.zeros((max_power,))
max_auto_diffs = np.zeros((max_power,))
max_accr_diffs = np.zeros((max_power,))

mean_all_dqrs = np.zeros((max_power,))
mean_together_dqrs = np.zeros((max_power,))
mean_apart_dqrs = np.zeros((max_power,))
mean_auto_dqrs = np.zeros((max_power,))
mean_accr_dqrs = np.zeros((max_power,))

mean_coarse_diff, median_coarse_diff, max_coarse_diff = coarse_diff_acc.output()
mean_coarse_diff *= 1.e3
median_coarse_diff *= 1.e3
max_coarse_diff *= 1.e3
mean_all_diffs[0] = mean_coarse_diff
mean_together_diffs[0] = mean_coarse_diff
mean_apart_diffs[0] = mean_coarse_diff
mean_auto_diffs[0] = mean_coarse_diff
mean_accr_diffs[0] = mean_coarse_diff
median_all_diffs[0] = median_coarse_diff
median_together_diffs[0] = median_coarse_diff
median_apart_diffs[0] = median_coarse_diff
median_auto_diffs[0] = median_coarse_diff
median_accr_diffs[0] = median_coarse_diff
max_all_diffs[0] = max_coarse_diff
max_together_diffs[0] = max_coarse_diff
max_apart_diffs[0] = max_coarse_diff
max_auto_diffs[0] = max_coarse_diff
max_accr_diffs[0] = max_coarse_diff

mean_coarse_dqr = coarse_dqr_acc.output()[0]
mean_coarse_dqr *= 1.e3
mean_all_dqrs[0] = mean_coarse_dqr
mean_together_dqrs[0] = mean_coarse_dqr
mean_apart_dqrs[0] = mean_coarse_dqr
mean_auto_dqrs[0] = mean_coarse_dqr
mean_accr_dqrs[0] = mean_coarse_dqr

mean_twnt_ref *= 1.e3 / number_grid_cells
mean_delqr_ref *= 1.e3 / number_grid_cells
mean_twnt_coarse *= 1.e3 / number_grid_cells
mean_delqr_coarse *= 1.e3 / number_grid_cells

for i in range(1, max_power):
    accums = all_diff_acc[i-1].output()
    mean_all_diffs[i] = accums[0] * 1.e3
    median_all_diffs[i] = accums[1] * 1.e3
    max_all_diffs[i] = accums[2] * 1.e3
    accums = together_diff_acc[i-1].output()
    mean_together_diffs[i] = accums[0] * 1.e3
    median_together_diffs[i] = accums[1] * 1.e3
    max_together_diffs[i] = accums[2] * 1.e3
    accums = apart_diff_acc[i-1].output()
    mean_apart_diffs[i] = accums[0] * 1.e3
    median_apart_diffs[i] = accums[1] * 1.e3
    max_apart_diffs[i] = accums[2] * 1.e3
    accums = accr_diff_acc[i-1].output()
    mean_accr_diffs[i] = accums[0] * 1.e3
    median_accr_diffs[i] = accums[1] * 1.e3
    max_accr_diffs[i] = accums[2] * 1.e3
    accums = auto_diff_acc[i-1].output()
    mean_auto_diffs[i] = accums[0] * 1.e3
    median_auto_diffs[i] = accums[1] * 1.e3
    max_auto_diffs[i] = accums[2] * 1.e3
    accums = all_dqr_acc[i-1].output()
    mean_all_dqrs[i] = accums[0] * 1.e3
    accums = together_dqr_acc[i-1].output()
    mean_together_dqrs[i] = accums[0] * 1.e3
    accums = apart_dqr_acc[i-1].output()
    mean_apart_dqrs[i] = accums[0] * 1.e3
    accums = auto_dqr_acc[i-1].output()
    mean_auto_dqrs[i] = accums[0] * 1.e3
    accums = accr_dqr_acc[i-1].output()
    mean_accr_dqrs[i] = accums[0] * 1.e3

print("Reference mean twnt and dqr/dt: ", mean_twnt_ref, mean_delqr_ref)

print("Coarse run mean twnt and dqr/dt: ", mean_twnt_coarse, mean_delqr_coarse)

suffix = "_icefilter"

timesteps = 300. / (2 ** np.arange(max_power))

plt.loglog(timesteps, mean_all_diffs, label='MG2', color='b')
plt.loglog(timesteps, mean_auto_diffs, label='Autoconversion only', color='y')
plt.loglog(timesteps, mean_accr_diffs, label='Accretion only', color='r')
plt.loglog(timesteps, mean_together_diffs, label='Coupled Auto/Accr', color='g')
#plt.loglog(timesteps, mean_apart_diffs, label='Auto/Accr Separately')
ref_line = timesteps*mean_all_diffs[-4]/timesteps[-4]
plt.loglog(timesteps, ref_line, 'k--', label='1st-order reference')
plt.title('Mean error in Cluster {} subsample'.format(cluster))
plt.xlabel('Timestep (s)')
plt.ylabel('Total water mass difference (g/kg/s)')
min_y = min(mean_all_diffs[-1], mean_auto_diffs[-1], mean_accr_diffs[-1],
            mean_together_diffs[-1], ref_line[-1])
max_y = max(ref_line[0], mean_all_diffs[0])
plt.axis([timesteps[-1], timesteps[0], min_y, 2.*max_y])
plt.legend(loc='lower right')
plt.savefig('substep_convergence_mean_c{}{}.eps'.format(cluster,suffix))
plt.close()

plt.loglog(timesteps, median_all_diffs, label='MG2', color='b')
plt.loglog(timesteps, median_auto_diffs, label='Autoconversion only', color='y')
plt.loglog(timesteps, median_accr_diffs, label='Accretion only', color='r')
plt.loglog(timesteps, median_together_diffs, label='Coupled Auto/Accr', color='g')
#plt.loglog(timesteps, median_apart_diffs, label='Auto/Accr Separately')
plt.loglog(timesteps, timesteps*median_all_diffs[-4]/timesteps[-4], 'k--',
           label='1st-order reference')
plt.title('Median error in Cluster {} subsample'.format(cluster))
plt.xlabel('Timestep (s)')
plt.ylabel('Total water mass difference (g/kg/s)')
plt.axis('tight')
plt.legend(loc='lower right')
plt.savefig('substep_convergence_median_c{}{}.eps'.format(cluster,suffix))
plt.close()

plt.loglog(timesteps, max_all_diffs, label='MG2', color='b')
plt.loglog(timesteps, max_auto_diffs, label='Autoconversion only', color='y')
plt.loglog(timesteps, max_accr_diffs, label='Accretion only', color='r')
plt.loglog(timesteps, max_together_diffs, label='Coupled Auto/Accr', color='g')
#plt.loglog(timesteps, max_apart_diffs, label='Auto/Accr Separately')
plt.loglog(timesteps, timesteps*max_all_diffs[-4]/timesteps[-4], 'k--',
           label='1st-order reference')
plt.title('Maximum error in Cluster {} subsample'.format(cluster))
plt.xlabel('Timestep (s)')
plt.ylabel('Total water mass difference (g/kg/s)')
plt.axis('tight')
plt.legend(loc='lower right')
plt.savefig('substep_convergence_max_c{}{}.eps'.format(cluster,suffix))
plt.close()

plt.loglog(timesteps, mean_all_dqrs, label='MG2', color='b')
#plt.loglog(timesteps, mean_auto_dqrs, label='Autoconversion only')
plt.loglog(timesteps, mean_accr_dqrs, label='Accretion only', color='r')
plt.loglog(timesteps, mean_together_dqrs, label='Coupled Auto/Accr', color='g')
#plt.loglog(timesteps, mean_apart_dqrs, label='Auto/Accr Separately')
plt.loglog(timesteps, timesteps*mean_all_dqrs[-4]/timesteps[-4], 'k--',
           label='1st-order reference')
plt.title('Mean rain mass error in Cluster {} subsample'.format(cluster))
plt.xlabel('Timestep (s)')
plt.ylabel('Rain tendency diff (g/kg/s)')
plt.axis('tight')
plt.legend(loc='lower right')
plt.savefig('substep_convergence_qr_c{}{}.eps'.format(cluster,suffix))
plt.close()

plt.hist(coarse_qrs, bins=200)
plt.xlabel('Final QR (g/kg)')
plt.ylabel('Number of cases ({} total)'.format(number_grid_cells))
plt.axis('tight')
plt.savefig('post_300s_qr_c{}{}.eps'.format(cluster,suffix))
plt.close()
