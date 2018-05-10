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
nproc = len(short_names)
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

tend_cutoff = 1.e-10
num_columns = ncol

tendencies = np.zeros((lev*num_columns, nproc))

def collect_tendencies(tends):
    tends_out = np.zeros((nproc,))
    for i in range(nproc):
        name = short_names[i]
        if name in ("revap", "ssubl", "vdepo", "vnudp"):
            tends_out[i] = tends[name][iq]
        elif name in ("cbrgs", "cacws", "cnccc", "cncct", "cacwi", "cbrgi", "cauto", "caccr"):
            tends_out[i] = tends[name][iqc]
        elif name in ("rfrez", "racrs"):
            tends_out[i] = tends[name][iqr]
        elif name in ("iauto", "iaccr"):
            tends_out[i] = tends[name][iqi]
        elif name in ("rnagg",):
            tends_out[i] = tends[name][inr]
        elif name in ("snagg",):
            tends_out[i] = tends[name][ins]
        elif name in ("cnact",):
            tends_out[i] = tends[name][inc]
        elif name in ("anadj",):
            tends_out[i] = tends[name][inc] + tends[name][ini] + \
                           tends[name][ins] + tends[name][inr]
        else:
            assert False, "Unrecognized process short name "+name+"."
    return tends_out

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
        tends = mg2_tendencies(level, t_loc, q_loc, qc_loc, nc_loc, qi_loc, ni_loc,
                               qr_loc, nr_loc, qs_loc, ns_loc)
        #tot_t = tends["Total"]
        #tend_norm = (abs(tot_t[iq]) + abs(tot_t[iqc]) + abs(tot_t[iqi]) +
        #             abs(tot_t[iqr]) + abs(tot_t[iqs])) / 2.
        #if tend_norm < tend_cutoff:
        #    continue
        tendencies[level + column*lev,:] = collect_tendencies(tends)

# Normalize tendencies to treat processes equally.
process_cutoffs = np.zeros((nproc,))
mean_abs = np.zeros((nproc,))
count = np.zeros((nproc,))
process_cutoffs[short_names.index('revap')] = 2.e-7
process_cutoffs[short_names.index('ssubl')] = 1.e-8
process_cutoffs[short_names.index('vdepo')] = 1.e-8
process_cutoffs[short_names.index('cbrgs')] = 2.e-9
process_cutoffs[short_names.index('cacws')] = 1.e-8
process_cutoffs[short_names.index('cnccc')] = 2.e-9
process_cutoffs[short_names.index('cncct')] = 1.e-11
process_cutoffs[short_names.index('cacwi')] = 1.e-11
process_cutoffs[short_names.index('rfrez')] = 1.e-7
process_cutoffs[short_names.index('racrs')] = 1.e-8
process_cutoffs[short_names.index('cbrgi')] = 1.e-5
process_cutoffs[short_names.index('cauto')] = 1.e-8
process_cutoffs[short_names.index('caccr')] = 1.e-8
process_cutoffs[short_names.index('iauto')] = 2.e-8
process_cutoffs[short_names.index('iaccr')] = 2.e-9
process_cutoffs[short_names.index('rnagg')] = 200.
process_cutoffs[short_names.index('snagg')] = 0.5
for i in range(nproc):
    abs_tendencies = np.abs(tendencies[:,i])
    abs_tendencies = np.where(abs_tendencies > process_cutoffs[i], abs_tendencies, 0.)
    count[i] = np.count_nonzero(abs_tendencies)
    if count[i] == 0:
        string = "Process "+short_names[i]+" was not active in any grid cell.\n"
        print(string)
        with open('./process_rates.txt', 'a+') as proc_rate_file:
            proc_rate_file.write(string)
        continue
    mean_abs[i] = abs_tendencies.sum()/count[i]
    string = "Mean magnitude for process "+short_names[i]+" is "+repr(mean_abs[i])+".\n" + \
             "Count for process "+short_names[i]+" is "+repr(count[i])+"/"+repr(lev*num_columns)+".\n"
    print(string)
    with open('./process_rates.txt', 'a+') as proc_rate_file:
        proc_rate_file.write(string)
    tendencies[:,i] /= mean_abs[i]

# Furthermore, clip the scaled tendencies to remove the worst outliers.
#tendencies = np.where(tendencies > 100., 100., np.where(tendencies < -100., -100., tendencies))

plt.autoscale(tight=True)

n_clusters = 10
# ANADJ kludge: we leave anadj out because we know that it's last.
clusters = clstr.KMeans(n_clusters=n_clusters).fit(tendencies[:,:-1])

# Below is necessary for non-KMeans clustering algorithms.
# cluster_centers = np.zeros((n_clusters,nproc))
# points_in_cluster = np.zeros((n_clusters,))
# for j in range(lev*num_columns):
#     cluster_centers[clusters.labels_[j],:] += tendencies[j,:]
#     points_in_cluster[clusters.labels_[j]] += 1
# for i in range(n_clusters):
#     cluster_centers[i,:] /= points_in_cluster[i]

OUT_FILE_NAME = '/g/g14/santos36/Data/MG2_data_collection.10_cluster_labels.0001-01-06-00000.nc'

# Copy over dimensions from the history file.
out_file = nc4.Dataset(OUT_FILE_NAME, 'w')
out_file.createDimension("ncol", ncol)
out_file.createDimension("time", None)
out_file.createDimension("lev", lev)
out_file.createDimension("nproc", nproc)
out_file.createDimension("name_len", 64)
out_file.createDimension("ncluster", n_clusters)

out_file.createVariable("time", "f8", ("time",))
out_file.variables['time'][:] = file.variables['time'][:]
out_file.createVariable("lev", "f8", ("lev",))
out_file.variables['lev'][:] = file.variables['lev'][:]
out_file.createVariable("label", "u1", ("time", "lev", "ncol"))
for column in range(num_columns):
    for level in range(lev):
        out_file.variables['label'][0,level,column] = clusters.labels_[level + column*lev]
out_file.createVariable("process_names", "S1", ("nproc", "name_len"))
out_file["process_names"]._Encoding = "utf-8"
for i in range(len(short_names)):
    out_file["process_names"][i] = process_names[short_names[i]]
out_file.createVariable("inertia", "f8")
out_file["inertia"][...] = clusters.inertia_
out_file.createVariable("cluster_centers", "f8", ("ncluster", "nproc"))
out_file["cluster_centers"][:,:-1] = clusters.cluster_centers_
# anadj currently not output since it isn't part of the clustering.
out_file["cluster_centers"][:,-1] = 0.
out_file.createVariable("process_cutoffs", "f8", ("nproc",))
out_file["process_cutoffs"][:] = process_cutoffs
out_file.createVariable("mean_abs", "f8", ("nproc",))
out_file["mean_abs"][:] = mean_abs
out_file.createVariable("count", "f8", ("nproc",))
out_file["count"][:] = count
out_file.source_file = HIST_FILE_NAME

out_file.close()

width = 0.1
ind = np.arange(nproc-1)

# Plot with all clusters.
bar_plots = []
for i in range(n_clusters):
    bar_plots.append(plt.bar(ind + width*i, clusters.cluster_centers_[i,:], width))

# Filename prefix for plotting outputs.
PLOT_PREFIX = "clusters10_scaled"

ax = plt.gca()
ax.set_xticks(ind)
# ANADJ kludge: leave out last again.
ax.set_xticklabels((process_names[name] for name in short_names[:-1]),
                   size='xx-small', rotation='vertical', wrap=True)
ax.tick_params('x', direction='out', pad=25)
plt.subplots_adjust(bottom=0.18)
plt.legend(bar_plots, ("Cluster {}".format(i)
                       for i in range(n_clusters)),
           loc='best', fontsize='xx-small')
plt.savefig('./{}.eps'.format(PLOT_PREFIX))
plt.close()

# Plots for each cluster.
width = 0.8
for i in range(n_clusters):
    count = 0
    for label in clusters.labels_:
        if label == i:
            count += 1
    plt.bar(ind, clusters.cluster_centers_[i,:], width)
    with open('./cluster_centers10.txt', 'a+') as cluster_file:
        string = "Cluster {} has inertia {!r} and centroid {!r}\n"
        string = string.format(i, clusters.inertia_, clusters.cluster_centers_[i,:])
        cluster_file.write(string)
    ax = plt.gca()
    ax.set_xticks(ind)
    # ANADJ kludge: leave out last again.
    ax.set_xticklabels((process_names[name] for name in short_names[:-1]),
                       size='xx-small', rotation='vertical', wrap=True)
    ax.tick_params('x', direction='out', pad=25)
    plt.subplots_adjust(bottom=0.18)
    plt.title("Center of mass of cluster {}, containing {} points".format(i, count))
    plt.savefig('./{}_{}.eps'.format(PLOT_PREFIX, i))
    plt.close()
