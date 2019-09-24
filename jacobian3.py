#!/usr/bin/env python

import sys

import numpy as np
import scipy.linalg as la
import scipy.stats as stats
import netCDF4 as nc4
import numdifftools as ndt

from mg2 import wv_sat_methods as wsm
from mg2 import micro_mg2_0 as mg

from mg2_constants import *

start_column = int(sys.argv[1])
end_column = int(sys.argv[2])


HIST_FILE_NAME = "/home/santos/Data/MG2_data_collection.cam.h1.0001-01-06-00000.nc"
OUT_FILE_NAME = "/home/santos/Data/Jacobian_cutoff_{}-{}.nc".format(start_column, end_column)

hfile = nc4.Dataset(HIST_FILE_NAME, 'r')

ncol = len(hfile.dimensions['ncol'])
lev = len(hfile.dimensions['lev'])
ilev = len(hfile.dimensions['ilev'])

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
    "cnact",
    "anadj",
]
nproc = len(short_names)
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
    tends["Total"][it] = tlat[0,0]
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
    tends["revap"][it] = -prer_evap[0,0] * latvap
    tends["revap"][iq] = prer_evap[0,0]
    tends["revap"][iqr] = -prer_evap[0,0]
    tends["revap"][inr] = nsubrtot[0,0]
    # Snow sublimation.
    tends["ssubl"][it] = -evapsnow[0,0] * (latvap + latice)
    tends["ssubl"][iq] = evapsnow[0,0]
    tends["ssubl"][iqs] = -evapsnow[0,0]
    # [Snow number sublimation currently ignored.]
    # Vapor deposition.
    tends["vdepo"][it] = cmeitot[0,0] * (latvap + latice)
    tends["vdepo"][iq] = -cmeitot[0,0]
    tends["vdepo"][iqi] = cmeitot[0,0]
    tends["vdepo"][ini] = ncmeitot[0,0]
    # External scheme's nucleation deposition (as implemented in MG2).
    tends["vnudp"][it] = mnudeptot[0,0] * (latvap + latice)
    tends["vnudp"][iq] = -mnudeptot[0,0]
    tends["vnudp"][iqi] = mnudeptot[0,0]
    tends["vnudp"][ini] = nnudeptot[0,0]
    # Bergeron process over snow.
    tends["cbrgs"][it] = bergstot[0,0] * latice
    tends["cbrgs"][iqc] = -bergstot[0,0]
    tends["cbrgs"][iqs] = bergstot[0,0]
    # Collection of cloud water by snow.
    tends["cacws"][it] = psacwstot[0,0] * latice
    tends["cacws"][iqc] = -psacwstot[0,0]
    tends["cacws"][iqs] = psacwstot[0,0]
    tends["cacws"][inc] = -npsacwstot[0,0]
    # Immersion Freezing.
    tends["cnccc"][it] = mnuccctot[0,0] * latice
    tends["cnccc"][iqc] = -mnuccctot[0,0]
    tends["cnccc"][iqi] = mnuccctot[0,0]
    tends["cnccc"][inc] = -nnuccctot[0,0]
    if use_hetfrz_classnuc:
        tends["cnccc"][ini] = nnuccctot[0,0]
    # Contact Freezing.
    tends["cncct"][it] = mnuccttot[0,0] * latice
    tends["cncct"][iqc] = -mnuccttot[0,0]
    tends["cncct"][iqi] = mnuccttot[0,0]
    tends["cncct"][inc] = -nnuccttot[0,0]
    tends["cncct"][ini] = nnuccttot[0,0]
    # Secondary ice production.
    tends["cacwi"][it] = msacwitot[0,0] * latice
    tends["cacwi"][iqc] = -msacwitot[0,0]
    tends["cacwi"][iqi] = msacwitot[0,0]
    tends["cacwi"][ini] = nsacwitot[0,0]
    # Freezing of rain.
    tends["rfrez"][it] = (mnuccrtot[0,0] + mnuccritot[0,0]) * latice
    tends["rfrez"][iqi] = mnuccritot[0,0]
    tends["rfrez"][iqr] = - mnuccrtot[0,0] - mnuccritot[0,0]
    tends["rfrez"][iqs] = mnuccrtot[0,0]
    tends["rfrez"][ini] = nnuccritot[0,0]
    tends["rfrez"][inr] = - nnuccrtot[0,0] - nnuccritot[0,0]
    tends["rfrez"][ins] = nnuccrtot[0,0]
    # Accretion of rain onto snow.
    tends["racrs"][it] = pracstot[0,0] * latice
    tends["racrs"][iqr] = -pracstot[0,0]
    tends["racrs"][iqs] = pracstot[0,0]
    tends["racrs"][inr] = -npracstot[0,0]
    # Bergeron process over cloud ice.
    tends["cbrgi"][it] = bergtot[0,0] * latice
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
    output = np.zeros((10*(nproc+1),))
    output[0:10] = s_inv_mult(tends["Total"])
    for i in range(nproc):
        output[10*(i+1):10*(i+2)] = s_inv_mult(tends[short_names[i]])
    return output

ind = np.arange(nproc)

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
    tot_t = tends["Total"]
    tend_norm = (abs(tot_t[iq]) + abs(tot_t[iqc]) + abs(tot_t[iqi]) +
                 abs(tot_t[iqr]) + abs(tot_t[iqs])) / 2.
    if tend_norm < tend_cutoff:
        return tends, tend_norm, None, None, None
    j_mg2 = ndt.Jacobian(run_mg2, method='forward', order=2)
    j_locs_dimless = j_mg2(s_inv_mult(state), level)
    j_loc = s.dot(sq.dot(la.solve_triangular(sr, j_locs_dimless[0:10,:].T, trans='T')).T)
    evals, evecs = la.eig(j_loc)
    # Sometimes la.eig spits out complex eigenvalues unnecessarily.
    if all(np.imag(evals) == 0.):
        evals = np.real(evals)
    if (np.imag(evecs) == 0.).all():
        evecs = np.real(evecs)
    associations = np.zeros((10, nproc))
    for i in range(nproc):
        j_proc_loc = s.dot(sq.dot(la.solve_triangular(sr, j_locs_dimless[10*(i+1):10*(i+2),:].T, trans='T')).T)
        j_proc_j_basis = la.solve(evecs, j_proc_loc.dot(evecs))
        associations[:,i] = np.abs(np.diag(j_proc_j_basis))
    for i in range(10):
        assoc_norm = associations[i,:].sum()
        if assoc_norm != 0.:
            associations[i,:] /= assoc_norm
    return tends, tend_norm, evals, evecs, associations

def calculate_indices(tends, evals, evecs):
    strengths = la.solve(evecs, tends["Total"])
    strengths /= la.norm(strengths, 1)
    participations = np.zeros((len(evals), nproc))
    for j in range(nproc):
        participations[:,j] = np.real(la.solve(evecs, tends[short_names[j]]))
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

min_t = 1.e-2
max_t = 1.e5
cutoff2 = 0.5
tend_cutoff = 1.e-10
#num_columns = 48602

grid_cells = []

tend_norms = []

singular_counter = 0

num_above_cutoff = 0

for column in range(start_column, end_column+1):
    print("On column: ", column, file=sys.stderr)
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
        tends, tend_norm, evals, evecs, associations = get_modes_at_level(level)
        tend_norms.append(tend_norm)
        if tend_norm < tend_cutoff:
            continue
        num_above_cutoff += 1
        try:
            strengths, participations = calculate_indices(tends, evals, evecs)
        except la.LinAlgError:
            singular_counter += 1
            print("Singular matrix or something. Ahhhh!")
            continue
        this_cell = dict()
        this_cell["coords"] = np.array((column, level))
        this_cell["tends"] = tends
        this_cell["evals"] = evals
        this_cell["participations"] = participations
        this_cell["associations"] = associations
        grid_cells.append(this_cell)

num_cell = len(grid_cells)

out_file = nc4.Dataset(OUT_FILE_NAME, 'w')
out_file.createDimension("num_cell", num_cell)
out_file.createDimension("ncoords", 2)
out_file.createDimension("nproc", nproc)
out_file.createDimension("nvars", 10)
out_file.createDimension("name_len", 5)

out_file.createVariable("process_names", "S1", ("nproc", "name_len"))
out_file["process_names"]._Encoding = "ascii"
process_name_data = np.array(short_names, dtype="S5")
out_file["process_names"][:] = nc4.stringtochar(process_name_data)

out_file.createVariable("cell_coords", "u4", ("num_cell", "ncoords"))
out_file.createVariable("total_rates", "f8", ("num_cell", "nvars"))
out_file.createVariable("process_rates", "f8", ("num_cell", "nvars", "nproc"))

# Build complex type.
complex128 = np.dtype([("real", np.float64), ("imag", np.float64)])
c128 = out_file.createCompoundType(complex128, "c128")

out_file.createVariable("eigenvalues", c128, ("num_cell", "nvars"))
out_file.createVariable("associations", "f8", ("num_cell", "nvars", "nproc"))
out_file.createVariable("participations", "f8", ("num_cell", "nvars", "nproc"))

for i in range(num_cell):
    cell = grid_cells[i]
    out_file["cell_coords"][i,:] = cell["coords"]
    out_file["total_rates"][i,:] = cell["tends"]["Total"]
    for j in range(nproc):
        out_file["process_rates"][i,:,j] = cell["tends"][short_names[j]]
    evals = np.empty((10,), complex128)
    evals["real"] = cell["evals"].real
    evals["imag"] = cell["evals"].imag
    out_file["eigenvalues"][i,:] = evals
    out_file["associations"][i,:,:] = cell["associations"]
    out_file["participations"][i,:,:] = cell["participations"]

out_file.tendency_cutoff = tend_cutoff

out_file.close()


# Tend plot

num_points = len(tend_norms)
tend_norms = np.array([t for t in tend_norms if t > 0.])
num_non0_points = len(tend_norms)

cutoff_score = stats.percentileofscore(tend_norms, tend_cutoff)
print("Cutoff was {} ({}% of non-zero tendencies were below this, {} points were above).".format(tend_cutoff, cutoff_score, num_above_cutoff))
print("{}/{} grid points had non-zero water mass tendency.".format(num_non0_points, num_points))
print("Minimum/Maximum tendency: ", tend_norms.min(), tend_norms.max())
print("Number of grid points with singular matrix issues: ", singular_counter)
