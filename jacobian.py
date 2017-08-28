#!/usr/bin/env python

import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numdifftools as ndt

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

column = 856
# About 800 mb.
level = 54

state = np.array((
    t[0,level,column],
    q[0,level,column],
    qc[0,level,column],
    nc[0,level,column],
    qi[0,level,column],
    ni[0,level,column],
    qr[0,level,column],
    nr[0,level,column],
    qs[0,level,column],
    ns[0,level,column],
))

# Approximate expected size of variables, used to scale to around O(1). Since we
# just need something in the ballpark, the below list somewhat arbitrarily picks
# the triple point of water, the mass of a modest amount of hydrometeor, and
# whatever the actual number concentrations are as scaling factors.
scale = np.array((
    273.16,
    6.11657e-3,
    1.e-6,
    max(nc[0,level,column], 1.),
    1.e-6,
    max(ni[0,level,column], 1.),
    1.e-6,
    max(nr[0,level,column], 1.),
    1.e-6,
    max(ns[0,level,column], 1.),
))

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
s[2,3] = scale[3] * np.pi/6. * 1000. * (20.e-6)**3
s[3,2] = scale[2] / (np.pi/6. * 1000. * (45.e-6)**3)
# Density of cloud ice is taken to be 500 kg/m^3.
s[4,5] = scale[5] * np.pi/6. * 500. * (20.e-6)**3
s[5,4] = scale[4] / (np.pi/6. * 500. * (45.e-6)**3)
# For rain and snow, add particles at 100 (low) and 400 (high) microns.
s[6,7] = scale[7] * np.pi/6. * 1000. * (100.e-6)**3
s[7,6] = scale[6] / (np.pi/6. * 1000. * (400.e-6)**3)
# Density of snow is 250 kg/m^3.
s[8,9] = scale[9] * np.pi/6. * 250. * (100.e-6)**3
s[9,8] = scale[8] / (np.pi/6. * 250. * (400.e-6)**3)

# Note that the condition number of s typically seems horrendous, but this is
# only because typical values for the different variables (e.g. qc vs nc) tend
# to be of wildly different orders of magnitude. In practice the error in
# inverting s seems to be fairly negligible.
# QR factorization is nonetheless useful for reducing issues with conditioning.
sq, sr = la.qr(s)

# Left-multiply by the inverse of s.
def s_inv_mult(state):
    return la.solve_triangular(sr, sq.T.dot(state))

print(state)
print(s_inv_mult(state))

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

def run_mg2(state):
    """Run MG2 on a single column, and return state tendencies.

    This function is specialized for use in finding the Jacobian at a single
    level.
    Order of variables is t, q, qc, nc, qi, ni, qr, nr, qs, ns.

    Note: This modifies global state (t_loc, q_loc, etc.) because I am lazy.
    """
    state = s.dot(state)
    # Note that we have to deal with the possibility that floating-point error
    # changes the sign of the hydrometeor masses/numbers.
    t_loc[0,level] = state[0]
    q_loc[0,level] = state[1]
    qc_loc[0,level] = max(state[2], 0.)
    nc_loc[0,level] = max(state[3], 1.e-12)
    qi_loc[0,level] = max(state[4], 0.)
    ni_loc[0,level] = max(state[5], 1.e-12)
    qr_loc[0,level] = max(state[6], 0.)
    nr_loc[0,level] = max(state[7], 1.e-12)
    qs_loc[0,level] = max(state[8], 0.)
    ns_loc[0,level] = max(state[9], 1.e-12)
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
        = mg.micro_mg_tend(timestep, t_loc, q_loc, qc_loc, qi_loc, nc_loc,
                           ni_loc, qr_loc, qs_loc, nr_loc, ns_loc,
                           relvar_loc, accre_enhan_loc, p_loc, pdel_loc,
                           precipf_loc, liqcldf_loc, icecldf_loc, naai_loc,
                           npccn_loc, rndst_loc, nacon_loc,
                           frzimm=frzimm_loc, frzcnt=frzcnt_loc,
                           frzdep=frzdep_loc, mgncol=1, nlev=lev, do_sed=False, do_inst=False)
    tendencies = np.array((tlat[0,level], qvlat[0,level], qctend[0,level], \
                           nctend[0,level], qitend[0,level], nitend[0,level], \
                           qrtend[0,level], nrtend[0,level], qstend[0,level], \
                           nstend[0,level]))
    return s_inv_mult(tendencies)

j_mg2 = ndt.Jacobian(run_mg2, method='forward', order=2)

j_loc = s.dot(sq.dot(la.solve_triangular(sr, j_mg2(s_inv_mult(state)).T, trans='T')).T)

evals, evecs = la.eig(j_loc)

print(j_loc)
print(',\n'.join([str(evecs[:,i]) for i in range(len(state))]))
# Sometimes la.eig spits out complex eigenvalues unnecessarily.
if all(np.imag(evals) == 0.):
    evals = np.real(evals)
print(evals)

with np.errstate(divide='ignore'):
    print(1./evals)
