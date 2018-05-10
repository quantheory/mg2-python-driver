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

timestep = 1.e-8

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

# Approximate expected size of variables, used to scale to around O(1). Since we
# just need something in the ballpark, the below list somewhat arbitrarily picks
# the triple point of water, the mass of a modest amount of hydrometeor, and
# a number concentration of 1/kg.
scale = np.zeros((10,))
scale[it] = 273.16
scale[iq] = 6.11657e-3
scale[iqc] = 1.e-6
scale[inc] = 1.
scale[iqi] = 1.e-6
scale[ini] = 1.
scale[iqr] = 1.e-6
scale[inr] = 1.
scale[iqs] = 1.e-6
scale[ins] = 1.

# This matrix scales the (hopefully O(1)) dummy variables to get back variables
# in the original units.
s = np.diag(scale)

# We also want to nudge things slightly so that when some number is added,
# mass is added as well. This is mainly because the number really doesn't mean
# anything when the mass is zero, and we want to avoid having to account for the
# physically meaningless limiters that come into play when the mass is zero.

# To avoid the minimum number limiters, we similarly add small amounts of number
# when mass is added.

# For cloud liquid and ice, add particles at 20 (low) and 45 (high) microns.
# Density of water is 1000 kg/m^3.
s[iqc,inc] = scale[inc] * np.pi/6. * 1000. * (20.e-6)**3
s[inc,iqc] = scale[iqc] / (np.pi/6. * 1000. * (45.e-6)**3)
# Density of cloud ice is taken to be 500 kg/m^3.
s[iqi,ini] = scale[ini] * np.pi/6. * 500. * (20.e-6)**3
s[ini,iqi] = scale[iqi] / (np.pi/6. * 500. * (45.e-6)**3)
# For rain and snow, add particles at 100 (low) and 400 (high) microns.
s[iqr,inr] = scale[inr] * np.pi/6. * 1000. * (100.e-6)**3
s[inr,iqr] = scale[iqr] / (np.pi/6. * 1000. * (400.e-6)**3)
# Density of snow is 250 kg/m^3.
s[iqs,ins] = scale[ins] * np.pi/6. * 250. * (100.e-6)**3
s[ins,iqs] = scale[iqs] / (np.pi/6. * 250. * (400.e-6)**3)

# Note that the condition number of s typically seems horrendous, but this is
# only because typical values for the different variables (e.g. qc vs nc) tend
# to be of wildly different orders of magnitude. In practice the error in
# inverting s seems to be fairly negligible.
# QR factorization is nonetheless useful for reducing issues with conditioning.
sq, sr = la.qr(s)

# Left-multiply by the inverse of s.
def s_inv_mult(state):
    return la.solve_triangular(sr, sq.T.dot(state))

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
    "revap": "Rain Evaporation",
    "ssubl": "Snow Sublimation",
    "vdepo": "Vapor/Ice DMS",
    "vnudp": "Nucleation\nDeposition",
    "cbrgs": "Bergeron (Snow)",
    "cacws": "Liq. Collection\nby Snow",
    "cnccc": "Immersion Freezing",
    "cncct": "Contact Freezing",
    "cacwi": "Secondary Ice\nProduction",
    "rfrez": "Hetero. Rain\nFreezing",
    "racrs": "Rain Accretion\nby Snow",
    "cbrgi": "Bergeron\n(Cloud ice)",
    "cauto": "Autoconversion",
    "caccr": "Accretion",
    "iauto": "Ice Autoconversion",
    "iaccr": "Ice Accretion",
    "rnagg": "Rain\nSelf-collection",
    "snagg": "Snow\nSelf-collection",
    "cnact": "Droplet\nActivation",
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
        nrout2, nsout2, drout2, dsout2, freqs, freqr, nfice, qcrat, precip_frac, \
        nadjtot, ncmeitot, mnudeptot, nnudeptot, npsacwstot, nnuccctot, nnuccttot, \
        nsacwitot, nnuccrtot, mnuccritot, nnuccritot, npracstot, npratot, nprctot, \
        nprc1tot, npraitot, nprcitot, nsubrtot, nraggtot, nsaggtot, errstring, \
        prer_evap \
        = mg.micro_mg_tend(timestep, t, q, qc, qi, nc,
                           ni, qr, qs, nr, ns,
                           relvar_loc, accre_enhan_loc, p_loc, pdel_loc,
                           precipf_loc, liqcldf_loc, icecldf_loc, naai_loc,
                           npccn_loc, rndst_loc, nacon_loc,
                           frzimm=frzimm_loc, frzcnt=frzcnt_loc,
                           frzdep=frzdep_loc, mgncol=1, nlev=lev, do_sed=False,
                           do_inst=False)
    # Add total tendencies to a dictionary.
    tends = {}
    tends["Total"] = np.zeros((10,))
    tends["Total"][it] = tlat[0,level]
    tends["Total"][iq] = qvlat[0,level]
    tends["Total"][iqc] = qctend[0,level]
    tends["Total"][inc] = nctend[0,level]
    tends["Total"][iqi] = qitend[0,level]
    tends["Total"][ini] = nitend[0,level]
    tends["Total"][iqr] = qrtend[0,level]
    tends["Total"][inr] = nrtend[0,level]
    tends["Total"][iqs] = qstend[0,level]
    tends["Total"][ins] = nstend[0,level]
    # Divide up individual processes.
    for name in process_names:
        tends[name] = np.zeros((10,))
    # Rain evaporation.
    tends["revap"][it] = -prer_evap[0,level] * latvap
    tends["revap"][iq] = prer_evap[0,level]
    tends["revap"][iqr] = -prer_evap[0,level]
    tends["revap"][inr] = nsubrtot[0,level]
    # Snow sublimation.
    tends["ssubl"][it] = -evapsnow[0,level] * (latvap + latice)
    tends["ssubl"][iq] = evapsnow[0,level]
    tends["ssubl"][iqs] = -evapsnow[0,level]
    # [Snow number sublimation currently ignored.]
    # Vapor deposition.
    tends["vdepo"][it] = cmeitot[0,level] * (latvap + latice)
    tends["vdepo"][iq] = -cmeitot[0,level]
    tends["vdepo"][iqi] = cmeitot[0,level]
    tends["vdepo"][ini] = ncmeitot[0,level]
    # External scheme's nucleation deposition (as implemented in MG2).
    tends["vnudp"][it] = mnudeptot[0,level] * (latvap + latice)
    tends["vnudp"][iq] = -mnudeptot[0,level]
    tends["vnudp"][iqi] = mnudeptot[0,level]
    tends["vnudp"][ini] = nnudeptot[0,level]
    # Bergeron process over snow.
    tends["cbrgs"][it] = bergstot[0,level] * latice
    tends["cbrgs"][iqc] = -bergstot[0,level]
    tends["cbrgs"][iqs] = bergstot[0,level]
    # Collection of cloud water by snow.
    tends["cacws"][it] = psacwstot[0,level] * latice
    tends["cacws"][iqc] = -psacwstot[0,level]
    tends["cacws"][iqs] = psacwstot[0,level]
    tends["cacws"][inc] = -npsacwstot[0,level]
    # Immersion Freezing.
    tends["cnccc"][it] = mnuccctot[0,level] * latice
    tends["cnccc"][iqc] = -mnuccctot[0,level]
    tends["cnccc"][iqi] = mnuccctot[0,level]
    tends["cnccc"][inc] = -nnuccctot[0,level]
    if use_hetfrz_classnuc:
        tends["cnccc"][ini] = nnuccctot[0,level]
    # Contact Freezing.
    tends["cncct"][it] = mnuccttot[0,level] * latice
    tends["cncct"][iqc] = -mnuccttot[0,level]
    tends["cncct"][iqi] = mnuccttot[0,level]
    tends["cncct"][inc] = -nnuccttot[0,level]
    tends["cncct"][ini] = nnuccttot[0,level]
    # Secondary ice production.
    tends["cacwi"][it] = msacwitot[0,level] * latice
    tends["cacwi"][iqc] = -msacwitot[0,level]
    tends["cacwi"][iqi] = msacwitot[0,level]
    tends["cacwi"][ini] = nsacwitot[0,level]
    # Freezing of rain.
    tends["rfrez"][it] = (mnuccrtot[0,level] + mnuccritot[0,level]) * latice
    tends["rfrez"][iqi] = mnuccritot[0,level]
    tends["rfrez"][iqr] = - mnuccrtot[0,level] - mnuccritot[0,level]
    tends["rfrez"][iqs] = mnuccrtot[0,level]
    tends["rfrez"][ini] = nnuccritot[0,level]
    tends["rfrez"][inr] = - nnuccrtot[0,level] - nnuccritot[0,level]
    tends["rfrez"][ins] = nnuccrtot[0,level]
    # Accretion of rain onto snow.
    tends["racrs"][it] = pracstot[0,level] * latice
    tends["racrs"][iqr] = -pracstot[0,level]
    tends["racrs"][iqs] = pracstot[0,level]
    tends["racrs"][inr] = -npracstot[0,level]
    # Bergeron process over cloud ice.
    tends["cbrgi"][it] = bergtot[0,level] * latice
    tends["cbrgi"][iqc] = -bergtot[0,level]
    tends["cbrgi"][iqi] = bergtot[0,level]
    # Autoconversion.
    tends["cauto"][iqc] = -prctot[0,level]
    tends["cauto"][iqr] = prctot[0,level]
    tends["cauto"][inc] = -nprc1tot[0,level]
    tends["cauto"][inr] = nprctot[0,level]
    # Accretion.
    tends["caccr"][iqc] = -pratot[0,level]
    tends["caccr"][iqr] = pratot[0,level]
    tends["caccr"][inc] = -npratot[0,level]
    # Ice autoconversion.
    tends["iauto"][iqi] = -prcitot[0,level]
    tends["iauto"][iqs] = prcitot[0,level]
    tends["iauto"][ini] = -nprcitot[0,level]
    tends["iauto"][ins] = nprcitot[0,level]
    # Ice accretion.
    tends["iaccr"][iqi] = -praitot[0,level]
    tends["iaccr"][iqs] = praitot[0,level]
    tends["iaccr"][ini] = -npraitot[0,level]
    # Rain self-collection.
    tends["rnagg"][inr] = nraggtot[0,level]
    # Snow self-collection.
    tends["snagg"][ins] = nsaggtot[0,level]
    # Droplet activation.
    tends["cnact"][inc] = npccn_loc[0,level]
    # Number adjustment.
    tends["anadj"][inc] = nadjtot[0,level,0]
    tends["anadj"][ini] = nadjtot[0,level,1]
    tends["anadj"][inr] = nadjtot[0,level,2]
    tends["anadj"][ins] = nadjtot[0,level,3]
    # Check the above now.
    budget_total = np.zeros((10,))
    for name in process_names:
        budget_total += tends[name]
    assert np.allclose(budget_total, tends["Total"])
    return tends

def run_mg2(state, level):
    """Run MG2 on a single column, and return state tendencies.

    This function is specialized for use in finding the Jacobian at a single
    level.
    """
    state = s.dot(state)
    # Note that we have to deal with the possibility that floating-point error
    # changes the sign of the hydrometeor masses/numbers.
    my_t = t_loc.copy()
    my_q = q_loc.copy()
    my_qc = qc_loc.copy()
    my_nc = nc_loc.copy()
    my_qi = qi_loc.copy()
    my_ni = ni_loc.copy()
    my_qr = qr_loc.copy()
    my_nr = nr_loc.copy()
    my_qs = qs_loc.copy()
    my_ns = ns_loc.copy()
    my_t[0,level] = state[it]
    my_q[0,level] = state[iq]
    my_qc[0,level] = max(state[iqc], 0.)
    my_nc[0,level] = max(state[inc], 1.e-12)
    my_qi[0,level] = max(state[iqi], 0.)
    my_ni[0,level] = max(state[ini], 1.e-12)
    my_qr[0,level] = max(state[iqr], 0.)
    my_nr[0,level] = max(state[inr], 1.e-12)
    my_qs[0,level] = max(state[iqs], 0.)
    my_ns[0,level] = max(state[ins], 1.e-12)
    tends = mg2_tendencies(level, my_t, my_q, my_qc, my_nc, my_qi, my_ni, my_qr,
                           my_nr, my_qs, my_ns)
    return s_inv_mult(tends["Total"])

ind = np.arange(len(short_names))

column = 856

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

def get_modes_at_level(level):
    state = np.zeros((10,))
    state[it] = t_loc[0,level]
    state[iq] = q_loc[0,level]
    state[iqc] = qc_loc[0,level]
    state[inc] = nc_loc[0,level]
    state[iqi] = qi_loc[0,level]
    state[ini] = ni_loc[0,level]
    state[iqr] = qr_loc[0,level]
    state[inr] = nr_loc[0,level]
    state[iqs] = qs_loc[0,level]
    state[ins] = ns_loc[0,level]
    tends = mg2_tendencies(level, t_loc, q_loc, qc_loc, nc_loc, qi_loc, ni_loc,
                           qr_loc, nr_loc, qs_loc, ns_loc)

    j_mg2 = ndt.Jacobian(run_mg2, method='forward', order=2)
    j_loc = s.dot(sq.dot(la.solve_triangular(sr, j_mg2(s_inv_mult(state),
                                                       level).T, trans='T')).T)
    evals, evecs = la.eig(j_loc)
    # Sometimes la.eig spits out complex eigenvalues unnecessarily.
    if all(np.imag(evals) == 0.):
        evals = np.real(evals)
    if (np.imag(evecs) == 0.).all():
        evecs = np.real(evecs)
    return tends, evals, evecs

def calculate_indices(tends, evals, evecs):
    strengths = la.solve(evecs, tends["Total"])
    strengths /= la.norm(strengths, 1)
    participations = np.zeros((len(evals), len(ind)))
    for j in range(len(ind)):
        participations[:,j] = la.solve(evecs, tends[short_names[j]])
    for i in range(len(evals)):
        part_norm = la.norm(participations[i,:], 1)
        str_sign = np.sign(np.real(strengths[i]))
        # Skip normalization for "inactive" modes.
        if part_norm == 0.:
            continue
        if str_sign == 0.:
            str_sign = 1.
        participations[i,:] /= part_norm * str_sign
    return (strengths, participations)

def evec_bar_graph(evals, strengths, participations):
    width = 0.08
    bar_plots = []
    for i in range(len(evals)):
        bar_plots.append(plt.bar(ind + width*i, participations[i,:], width))

    ax = plt.gca()
    ax.set_xticks(ind)
    ax.set_xticklabels((process_names[name] for name in short_names),
                       size='xx-small', rotation='vertical', wrap=True)
    ax.tick_params('x', direction='out', pad=25)
    plt.subplots_adjust(bottom=0.18)
    plt.legend(bar_plots, ("Timescale={: .2e}s,Strength={: .2e}".format(1./evalue, np.abs(strength))
                           for evalue, strength in zip(evals, strengths)),
               loc='best', fontsize='xx-small')
    plt.savefig('./evecs.eps')
    plt.close()

# For examining a given point
#tends, evals, evecs = get_modes_at_level(54)
#strengths, participations = calculate_indices(tends, evals, evecs)

#evec_bar_graph(evals, strengths, participations)

min_t = 1.e-2
max_t = 1.e5
cutoff1 = 0.1
cutoff2 = 0.5
cutoff3 = 0.9
tend_cutoff = 1.e-10
num_columns = 2048

plt.autoscale(tight=True)

evalues1 = [dict() for i in range(ncluster)]
evalues2 = [dict() for i in range(ncluster)]
evalues3 = [dict() for i in range(ncluster)]
for i in range(len(short_names)):
    for j in range(ncluster):
        evalues1[j][short_names[i]] = []
        evalues2[j][short_names[i]] = []
        evalues3[j][short_names[i]] = []

tend_norms = []

singular_counter = 0

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
    for level in range(lev):
        c = label[0,level,column]
        tends, evals, evecs = get_modes_at_level(level)
        tot_t = tends["Total"]
        tend_norm = (abs(tot_t[iq]) + abs(tot_t[iqc]) + abs(tot_t[iqi]) +
                     abs(tot_t[iqr]) + abs(tot_t[iqs])) / 2.
        tend_norms.append(tend_norm)
        if tend_norm < tend_cutoff:
            continue
        try:
            strengths, participations = calculate_indices(tends, evals, evecs)
        except la.LinAlgError:
            singular_counter += 1
            print("Singular matrix or something. Ahhhh!")
            continue
        for i in range(len(evals)):
            for j in range(len(short_names)):
                if participations[i,j] >= cutoff1:
                    evalues1[c][short_names[j]].append(np.real(evals[i]))
                    if participations[i,j] >= cutoff2:
                        evalues2[c][short_names[j]].append(np.real(evals[i]))
                        if participations[i,j] >= cutoff3:
                            evalues3[c][short_names[j]].append(np.real(evals[i]))

# Number of bins for plotting purposes.
nbins = 50

cmap = plt.get_cmap('Greys')

for c in range(ncluster):
    pos_evalues1 = {}
    pos_evalues2 = {}
    pos_evalues3 = {}
    hist_values = np.zeros((len(short_names)-1, nbins))
    bins = np.logspace(np.log10(1./max_t), np.log10(1./min_t), nbins+1)

    nmodes = 0

    for i in range(len(short_names)):
        process = short_names[i]
        if process == 'anadj':
            continue
        pos_evalues1[process] = np.array([t for t in evalues1[c][process] if t > 1./max_t and t < 1./min_t])
        pos_evalues2[process] = np.array([t for t in evalues2[c][process] if t > 1./max_t and t < 1./min_t])
        pos_evalues3[process] = np.array([t for t in evalues3[c][process] if t > 1./max_t and t < 1./min_t])
        #if len(pos_evalues1[process]) == 0:
        #    continue
        nmodes += len(pos_evalues2[process])
        hist_values[i,:], _ = np.histogram(pos_evalues2[process], bins=bins)
        #plt.hist(pos_evalues1[process], bins=bins)
        #if len(pos_evalues2[process]) != 0:
        #    plt.hist(pos_evalues2[process], bins=bins)
        #    if len(pos_evalues3[process]) != 0:
        #        plt.hist(pos_evalues3[process], bins=bins)

    plt.pcolor(bins, ind, hist_values, edgecolors='k', cmap=cmap)
    plt.title("Number of positive eigenvalues for cluster {} \n(based on {} modes, cutoff={})"
              .format(c,
                      nmodes,
                      cutoff2))
    plt.gca().set_xscale("log")
    plt.xlabel("Eigenvalue (1/s)")
    plt.xlim(1./max_t, 1./min_t)
    ax = plt.gca()
    ax.set_yticks(ind)
    # ANADJ kludge: leave out last again.
    ax.set_yticklabels((process_names[name] for name in short_names[:-1]),
                       size='xx-small', wrap=True)
    ax.tick_params('y', direction='out', pad=25)
    plt.subplots_adjust(left=0.18)
    plt.colorbar()
    plt.savefig('./time_hist_2D_pos_c{}.eps'.format(c))
    plt.close()

    neg_evalues2 = {}
    hist_values = np.zeros((len(short_names)-1, nbins))
    bins = np.logspace(np.log10(1./max_t), np.log10(1./min_t), nbins+1)

    nmodes = 0

    for i in range(len(short_names)):
        process = short_names[i]
        if process == 'anadj':
            continue
        neg_evalues2[process] = np.array([-t for t in evalues2[c][process] if -t > 1./max_t and -t < 1./min_t])
        nmodes += len(neg_evalues2[process])
        hist_values[i,:], _ = np.histogram(neg_evalues2[process], bins=bins)

    meshind, meshbins = np.meshgrid(ind, bins)
    plt.pcolor(bins, ind, hist_values, edgecolors='k', cmap=cmap)
    plt.title("Number of negative eigenvalues for cluster {} \n(based on {} modes, cutoff={})"
              .format(c,
                      nmodes,
                      cutoff2))
    plt.gca().set_xscale("log")
    plt.xlabel("Eigenvalue (1/s)")
    plt.xlim(1./min_t, 1./max_t)
    ax = plt.gca()
    ax.set_yticks(ind)
    # ANADJ kludge: leave out last again.
    ax.set_yticklabels((process_names[name] for name in short_names[:-1]),
                       size='xx-small', wrap=True)
    ax.tick_params('y', direction='out', pad=25)
    plt.subplots_adjust(left=0.18)
    plt.colorbar()
    plt.savefig('./time_hist_2D_neg_c{}.eps'.format(c))
    plt.close()

    n0_evalues2 = {}
    hist_values = np.zeros((len(short_names)-1, nbins))
    bins = np.linspace(-1./max_t, 1./max_t, nbins+1)

    nmodes = 0

    for i in range(len(short_names)):
        process = short_names[i]
        if process == 'anadj':
            continue
        n0_evalues2[process] = np.array([t for t in evalues2[c][process] if t <= 1./max_t and t >= -1./max_t])
        nmodes += len(n0_evalues2[process])
        hist_values[i,:], _ = np.histogram(n0_evalues2[process], bins=bins)

    meshind, meshbins = np.meshgrid(ind, bins)
    plt.pcolor(bins, ind, hist_values, edgecolors='k', cmap=cmap)
    plt.title("Number of near-zero eigenvalues for cluster {} \n(based on {} modes, cutoff={})"
              .format(c,
                      nmodes,
                      cutoff2))
    plt.xlabel("Eigenvalue (1/s)")
    plt.xlim(-1./max_t, 1./max_t)
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax = plt.gca()
    ax.set_yticks(ind)
    # ANADJ kludge: leave out last again.
    ax.set_yticklabels((process_names[name] for name in short_names[:-1]),
                       size='xx-small', wrap=True)
    ax.tick_params('y', direction='out', pad=25)
    plt.subplots_adjust(left=0.18)
    plt.colorbar()
    plt.savefig('./time_hist_2D_n0_c{}.eps'.format(c))
    plt.close()


# Tend plot

num_points = len(tend_norms)
tend_norms = np.array([t for t in tend_norms if t > 0.])
num_non0_points = len(tend_norms)
bins = np.linspace(0., np.percentile(tend_norms, 80.), 50)
plt.hist(tend_norms, bins=bins)
plt.title("Number of grid points with given net water mass tendency")
plt.xlabel("Net water tendency")
plt.ylabel("Number of grid points")
plt.savefig('./tend_hist.eps')
plt.close()

# Log version of plot

bins = np.logspace(np.log10(tend_norms.min()), np.log10(tend_norms.max()), 50)
plt.hist(tend_norms, bins=bins)
plt.gca().set_xscale("log")
plt.title("Number of grid points with given net water mass tendency")
plt.xlabel("Net water tendency")
plt.ylabel("Number of grid points")
plt.savefig('./tend_hist_log.eps')
plt.close()

cutoff_score = stats.percentileofscore(tend_norms, tend_cutoff)
print("Cutoff was {} ({}% of non-zero tendencies were below this).".format(tend_cutoff, cutoff_score))
print("{}/{} grid points had non-zero water mass tendency.".format(num_non0_points, num_points))
print("Minimum/Maximum tendency: ", tend_norms.min(), tend_norms.max())
print("Number of grid points with singular matrix issues: ", singular_counter)
