"""
mg2 is a "shadow module" for the _mg2 module, which cleans up the _mg2 interface
and adds functionality.
"""

import _mg2
import numpy             as np
import numpy.linalg      as la
import matplotlib.pyplot as plt
import netCDF4           as nc4
import os
import pycurl

from   mg2_constants import *

################################################################################

wsm      = _mg2.wv_sat_methods
mg       = _mg2.micro_mg2_0
mg_utils = _mg2.micro_mg_utils

################################################################################

def obtain_datafile(filename, url):
    """
    Check to see if filename exists. If not, download the file from the given
    URL. If it does exist, assume the file has already been downloaded and skip
    the download.
    """
    if not os.path.isfile(filename):
        with open(filename, 'wb') as f:
            c = pycurl.Curl()
            c.setopt(c.URL, url)
            c.setopt(c.WRITEDATA, f)
            print('Downloading', filename, '...', end='')
            c.perform()
            c.close()
            print('Done')

################################################################################

class tendency(object):

    def __init__(self, qcsinksum_rate1ord, tlat, qvlat, qctend, qitend, nctend,
                 nitend, qrtend, qstend, nrtend, nstend, effc, effc_fn, effi,
                 prect, preci, nevapr, evapsnow, prain, prodsnow, cmeout, deffi,
                 pgamrad, lamcrad, qsout, dsout, rflx, sflx, qrout, reff_rain,
                 reff_snow, qcsevap, qisevap, qvres, cmeitot, vtrmc, vtrmi, umr,
                 ums, qcsedten, qisedten, qrsedten, qssedten, pratot, prctot,
                 mnuccctot, mnuccttot, msacwitot, psacwstot, bergstot, bergtot,
                 melttot, homotot, qcrestot, prcitot, praitot, qirestot,
                 mnuccrtot, pracstot, meltsdttot, frzrdttot, mnuccdtot, nrout,
                 nsout, refl, arefl, areflz, frefl, csrfl, acsrfl, fcsrfl,
                 rercld, ncai, ncal, qrout2, qsout2, nrout2, nsout2, drout2,
                 dsout2, freqs, freqr, nfice, qcrat, errstring, prer_evap):
        self.qcsinksum_rate1ord = qcsinksum_rate1ord
        self.tlat               = tlat
        self.qvlat              = qvlat
        self.qctend             = qctend
        self.qitend             = qitend
        self.nctend             = nctend
        self.nitend             = nitend
        self.qrtend             = qrtend
        self.qstend             = qstend
        self.nrtend             = nrtend
        self.nstend             = nstend
        self.effc               = effc
        self.effc_fn            = effc_fn
        self.effi               = effi
        self.prect              = prect
        self.preci              = preci
        self.nevapr             = nevapr
        self.evapsnow           = evapsnow
        self.prain              = prain
        self.prodsnow           = prodsnow
        self.cmeout             = cmeout
        self.deffi              = deffi
        self.pgamrad            = pgamrad
        self.lamcrad            = lamcrad
        self.qsout              = qsout
        self.dsout              = dsout
        self.rflx               = rflx
        self.sflx               = sflx
        self.qrout              = qrout
        self.reff_rain          = reff_rain
        self.reff_snow          = reff_snow
        self.qcsevap            = qcsevap
        self.qisevap            = qisevap
        self.qvres              = qvres
        self.cmeitot            = cmeitot
        self.vtrmc              = vtrmc
        self.vtrmi              = vtrmi
        self.umr                = umr
        self.ums                = ums
        self.qcsedten           = qcsedten
        self.qisedten           = qisedten
        self.qrsedten           = qrsedten
        self.qssedten           = qssedten
        self.pratot             = pratot
        self.prctot             = prctot
        self.mnuccctot          = mnuccctot
        self.mnuccttot          = mnuccttot
        self.msacwitot          = msacwitot
        self.psacwstot          = psacwstot
        self.bergstot           = bergstot
        self.bergtot            = bergtot
        self.melttot            = melttot
        self.homotot            = homotot
        self.qcrestot           = qcrestot
        self.prcitot            = prcitot
        self.praitot            = praitot
        self.qirestot           = qirestot
        self.mnuccrtot          = mnuccrtot
        self.pracstot           = pracstot
        self.meltsdttot         = meltsdttot
        self.frzrdttot          = frzrdttot
        self.mnuccdtot          = mnuccdtot
        self.nrout              = nrout
        self.nsout              = nsout
        self.refl               = refl
        self.arefl              = arefl
        self.areflz             = areflz
        self.frefl              = frefl
        self.csrfl              = csrfl
        self.acsrfl             = acsrfl
        self.fcsrfl             = fcsrfl
        self.rercld             = rercld
        self.ncai               = ncai
        self.ncal               = ncal
        self.qrout2             = qrout2
        self.qsout2             = qsout2
        self.nrout2             = nrout2
        self.nsout2             = nsout2
        self.drout2             = drout2
        self.dsout2             = dsout2
        self.freqs              = freqs
        self.freqr              = freqr
        self.nfice              = nfice
        self.qcrat              = qcrat
        self.errstring          = errstring
        self.prer_evap          = prer_evap

################################################################################

def micro_mg_tendency(deltat, t_loc, q_loc, qc_loc, qi_loc, nc_loc, ni_loc,
                      qr_loc, qs_loc, nr_loc, ns_loc, relvar_loc,
                      accre_enhan_loc, p_loc, pdel_loc, precipf_loc,
                      liqcldf_loc, icecldf_loc, naai_loc, npccn_loc, rndst_loc,
                      nacon_loc, frzimm, frzcnt, frzdep, mgncol, nlev):
    out_args = mg.micro_mg_tend(deltat, t_loc, q_loc, qc_loc, qi_loc, nc_loc,
                                ni_loc, qr_loc, qs_loc, nr_loc, ns_loc,
                                relvar_loc, accre_enhan_loc, p_loc, pdel_loc,
                                precipf_loc, liqcldf_loc, icecldf_loc, naai_loc,
                                npccn_loc, rndst_loc, nacon_loc, frzimm=frzimm,
                                frzcnt=frzcnt, frzdep=frzdep, mgncol=mgncol,
                                nlev=nlev)
    return tendency(*out_args)

################################################################################

class convergence_test(object):

    ############################################################################

    def __init__(self, filename, mgncol=128):
        # Open the NetCDF file
        self.file   = nc4.Dataset(filename, 'r')
        self.mgncol = mgncol
        self.ncol   = len(self.file.dimensions['ncol'])
        self.lev    = len(self.file.dimensions['lev' ])
        self.ilev   = len(self.file.dimensions['ilev'])

        # Define the variables
        self.t              = self.file.variables["MG2IN_T"]
        self.q              = self.file.variables["MG2IN_Q"]
        self.qc             = self.file.variables["MG2IN_QC"]
        self.qi             = self.file.variables["MG2IN_QI"]
        self.nc             = self.file.variables["MG2IN_NC"]
        self.ni             = self.file.variables["MG2IN_NI"]
        self.qr             = self.file.variables["MG2IN_QR"]
        self.qs             = self.file.variables["MG2IN_QS"]
        self.nr             = self.file.variables["MG2IN_NR"]
        self.ns             = self.file.variables["MG2IN_NS"]
        self.relvar         = self.file.variables["MG2IN_RELVAR"]
        self.accre_enhan    = self.file.variables["MG2IN_ACCRE_ENHAN"]
        self.p              = self.file.variables["MG2IN_P"]
        self.pdel           = self.file.variables["MG2IN_PDEL"]
        self.precipf        = self.file.variables["MG2IN_PRECIP"]
        self.liqcldf        = self.file.variables["MG2IN_LIQCLDF"]
        self.icecldf        = self.file.variables["MG2IN_ICECLDF"]
        self.naai           = self.file.variables["MG2IN_NAAI"]
        self.npccn          = self.file.variables["MG2IN_NPCCN"]
        self.rndst          = np.empty((self.t.shape[0], self.t.shape[1], self.t.shape[2], 4))
        self.rndst[:,:,:,0] = self.file.variables["MG2IN_RNDST1"][:]
        self.rndst[:,:,:,1] = self.file.variables["MG2IN_RNDST2"][:]
        self.rndst[:,:,:,2] = self.file.variables["MG2IN_RNDST3"][:]
        self.rndst[:,:,:,3] = self.file.variables["MG2IN_RNDST4"][:]
        self.nacon          = np.empty((self.t.shape[0], self.t.shape[1], self.t.shape[2], 4))
        self.nacon[:,:,:,0] = self.file.variables["MG2IN_NACON1"][:]
        self.nacon[:,:,:,1] = self.file.variables["MG2IN_NACON2"][:]
        self.nacon[:,:,:,2] = self.file.variables["MG2IN_NACON3"][:]
        self.nacon[:,:,:,3] = self.file.variables["MG2IN_NACON4"][:]
        self.frzimm         = self.file.variables["MG2IN_FRZIMM"]
        self.frzcnt         = self.file.variables["MG2IN_FRZCNT"]
        self.frzdep         = self.file.variables["MG2IN_FRZDEP"]

        # Define the local variables
        self.t_loc           = np.empty((self.mgncol, self.t.shape[1]),           order='F')
        self.q_loc           = np.empty((self.mgncol, self.q.shape[1]),           order='F')
        self.qc_loc          = np.empty((self.mgncol, self.qc.shape[1]),          order='F')
        self.qi_loc          = np.empty((self.mgncol, self.qi.shape[1]),          order='F')
        self.nc_loc          = np.empty((self.mgncol, self.nc.shape[1]),          order='F')
        self.ni_loc          = np.empty((self.mgncol, self.ni.shape[1]),          order='F')
        self.qr_loc          = np.empty((self.mgncol, self.qr.shape[1]),          order='F')
        self.qs_loc          = np.empty((self.mgncol, self.qs.shape[1]),          order='F')
        self.nr_loc          = np.empty((self.mgncol, self.nr.shape[1]),          order='F')
        self.ns_loc          = np.empty((self.mgncol, self.ns.shape[1]),          order='F')
        self.relvar_loc      = np.empty((self.mgncol, self.relvar.shape[1]),      order='F')
        self.accre_enhan_loc = np.empty((self.mgncol, self.accre_enhan.shape[1]), order='F')
        self.p_loc           = np.empty((self.mgncol, self.p.shape[1]),           order='F')
        self.pdel_loc        = np.empty((self.mgncol, self.pdel.shape[1]),        order='F')
        self.precipf_loc     = np.empty((self.mgncol, self.precipf.shape[1]),     order='F')
        self.liqcldf_loc     = np.empty((self.mgncol, self.liqcldf.shape[1]),     order='F')
        self.icecldf_loc     = np.empty((self.mgncol, self.icecldf.shape[1]),     order='F')
        self.naai_loc        = np.empty((self.mgncol, self.naai.shape[1]),        order='F')
        self.npccn_loc       = np.empty((self.mgncol, self.npccn.shape[1]),       order='F')
        self.rndst_loc       = np.empty((self.mgncol, self.rndst.shape[1], 4),    order='F')
        self.nacon_loc       = np.empty((self.mgncol, self.nacon.shape[1], 4),    order='F')
        self.frzimm_loc      = np.empty((self.mgncol, self.frzimm.shape[1]),      order='F')
        self.frzcnt_loc      = np.empty((self.mgncol, self.frzcnt.shape[1]),      order='F')
        self.frzdep_loc      = np.empty((self.mgncol, self.frzdep.shape[1]),      order='F')

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

        self.loc_arrays = {'T' : self.t_loc,
                           'Q' : self.q_loc,
                           'QC': self.qc_loc,
                           'QI': self.qi_loc,
                           'QR': self.qr_loc,
                           'QS': self.qs_loc
                           }
        self.var_names = sorted(list(self.loc_arrays.keys()))

    ############################################################################

    def initialize(self, total_columns, timesteps):
        self.total_columns = total_columns
        self.timesteps     = np.array(timesteps)
        self.final_time    = self.timesteps[-1]
        self.norms = {}
        self.finals = {}
        for name in self.var_names:
            self.norms[name]  = np.zeros((self.total_columns, self.timesteps.size - 1))
            self.finals[name] = np.zeros((self.total_columns, self.lev))

        # Check (annoying) limitation on mgncol.
        assert self.total_columns % self.mgncol == 0, \
               "total columns ({}) does not divide MG2 batch size ({})".format(self.total_columns,
                                                                               self.mgncol)

        # Check to make sure timestep sizes are valid before entering loop.
        for timestep in self.timesteps:
            assert self.final_time % timestep == 0., \
                   "timestep ({}) does not divide final time ({})".format(timestep,
                                                                          self.final_time)

    ############################################################################

    def run_tests(self):
        for it in range(self.timesteps.size):
            print("Starting timestep =", self.timesteps[it])
            nsteps = int(round(self.final_time / self.timesteps[it]))
            deltat = float(self.timesteps[it])

            for offset in range(self.total_columns // self.mgncol):
                self.t_loc[:,:]           = self.t[          0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.q_loc[:,:]           = self.q[          0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.qc_loc[:,:]          = self.qc[         0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.qi_loc[:,:]          = self.qi[         0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.nc_loc[:,:]          = self.nc[         0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.ni_loc[:,:]          = self.ni[         0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.qr_loc[:,:]          = self.qr[         0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.qs_loc[:,:]          = self.qs[         0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.nr_loc[:,:]          = self.nr[         0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.ns_loc[:,:]          = self.ns[         0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.relvar_loc[:,:]      = self.relvar[     0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.accre_enhan_loc[:,:] = self.accre_enhan[0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.p_loc[:,:]           = self.p[          0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.pdel_loc[:,:]        = self.pdel[       0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.precipf_loc[:,:]     = self.precipf[    0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.liqcldf_loc[:,:]     = self.liqcldf[    0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.icecldf_loc[:,:]     = self.icecldf[    0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.naai_loc[:,:]        = self.naai[       0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.npccn_loc[:,:]       = self.npccn[      0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.rndst_loc[:,:,:]     = self.rndst[      0,:,offset*self.mgncol:(offset+1)*self.mgncol,:].transpose([1, 0, 2])
                self.nacon_loc[:,:,:]     = self.nacon[      0,:,offset*self.mgncol:(offset+1)*self.mgncol,:].transpose([1, 0, 2])
                self.frzimm_loc[:,:]      = self.frzimm[     0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.frzcnt_loc[:,:]      = self.frzcnt[     0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()
                self.frzdep_loc[:,:]      = self.frzdep[     0,:,offset*self.mgncol:(offset+1)*self.mgncol].transpose()

                for n in range(nsteps):
                    tend = micro_mg_tendency(deltat, self.t_loc, self.q_loc,
                                             self.qc_loc, self.qi_loc,
                                             self.nc_loc, self.ni_loc,
                                             self.qr_loc, self.qs_loc,
                                             self.nr_loc, self.ns_loc,
                                             self.relvar_loc,
                                             self.accre_enhan_loc, self.p_loc,
                                             self.pdel_loc, self.precipf_loc,
                                             self.liqcldf_loc, self.icecldf_loc,
                                             self.naai_loc, self.npccn_loc,
                                             self.rndst_loc, self.nacon_loc,
                                             frzimm=self.frzimm_loc,
                                             frzcnt=self.frzcnt_loc,
                                             frzdep=self.frzdep_loc,
                                             mgncol=self.mgncol, nlev=self.lev)

                    # Should use geopotential!
                    self.t_loc      += tend.tlat  * deltat / cpair
                    self.q_loc      += tend.qvlat * deltat
                    self.q_loc[:,:]  = np.where(self.q_loc < 1.e-12, 1.e-12, self.q_loc)
                    self.qc_loc     += tend.qctend * deltat
                    self.qc_loc[:,:] = np.where(self.qc_loc < 0., 0., self.qc_loc)
                    self.qi_loc     += tend.qitend * deltat
                    self.qi_loc[:,:] = np.where(self.qi_loc < 0., 0., self.qi_loc)
                    self.qr_loc     += tend.qrtend * deltat
                    self.qr_loc[:,:] = np.where(self.qr_loc < 0., 0., self.qr_loc)
                    self.qs_loc     += tend.qstend * deltat
                    self.qs_loc[:,:] = np.where(self.qs_loc < 0., 0., self.qs_loc)
                    self.nc_loc     += tend.nctend * deltat
                    self.nc_loc[:,:] = np.where(self.nc_loc > 1.e10, 1.e10, np.where(self.nc_loc < 1.e-12, 1.e-12, self.nc_loc))
                    self.ni_loc     += tend.nitend * deltat
                    self.ni_loc[:,:] = np.where(self.nc_loc > 1.e10, 1.e10, np.where(self.ni_loc < 1.e-12, 1.e-12, self.ni_loc))
                    self.nr_loc     += tend.nrtend * deltat
                    self.nr_loc[:,:] = np.where(self.nc_loc > 1.e10, 1.e10, np.where(self.nr_loc < 1.e-12, 1.e-12, self.nr_loc))
                    self.ns_loc     += tend.nstend * deltat
                    self.ns_loc[:,:] = np.where(self.nc_loc > 1.e10, 1.e10, np.where(self.ns_loc < 1.e-12, 1.e-12, self.ns_loc))

                if it == 0:
                    for name in self.var_names:
                        self.finals[name][offset*self.mgncol:(offset+1)*self.mgncol,:] = self.loc_arrays[name]
                else:
                    for name in self.var_names:
                        self.norms[name][offset*self.mgncol:(offset+1)*self.mgncol,it-1] \
                            = la.norm(self.finals[name][offset*self.mgncol:(offset+1)*self.mgncol,:] - self.loc_arrays[name],
                                      axis=1)

    ############################################################################

    def plot_variable_convergence(self, name, percentiles, estimate_slopes=True):
        """
        Create a convergence plot for a given variable.

        Inputs:
        name            - Variable name (must be one of the keys of loc_arrays).
        percentiles     - Percentiles of error to plot (e.g. 50. is the median).
        estimate_slopes - Print out an estimated slope for each percentile.
        """
        fig = plt.figure(figsize=(10,7.5))
        ax  = fig.add_subplot(1, 1, 1)
        for p in percentiles:
            norms_p = np.percentile(self.norms[name], p, axis=0)
            ax.loglog(self.timesteps[1:], norms_p, label='Percentile={:.1f}'.format(p))
            if estimate_slopes:
                coefs = np.polyfit(np.log(self.timesteps[1:]), np.log(norms_p), 1)
                print("Estimated slope for percentile {:.1f} is {}.".format(p, coefs[0]))
        ax.set_xlabel('Timestep (sec)', fontsize='x-large')
        ax.set_ylabel('Norm of difference in {} from $\Delta$t={}'.format(name, self.timesteps[0]),
                      fontsize='x-large')
        ax.legend(loc='best', fontsize='x-large')

    ############################################################################

    def plot_column_variable_convergence(self, name, percentiles, timestep, estimate_slopes=True):
        """
        Create a convergence plot for a given variable, tracking particular columns.

        Inputs:
        name            - Variable name (must be one of the keys of loc_arrays).
        percentiles     - Percentiles of error to plot (e.g. 50. is the median).
        timestep        - Which timestep to evaluate the percentiles at (must be in
                          the "timesteps" list).
        estimate_slopes - Print out an estimated slope for each column.
        """
        fig = plt.figure(figsize=(10,7.5))
        ax  = fig.add_subplot(1, 1, 1)
        indices = np.argsort(self.norms[name][:,list(self.timesteps).index(timestep)-1])
        for p in percentiles:
            norms_p = self.norms[name][indices[round(p/100.*(self.total_columns-1))]]
            ax.loglog(self.timesteps[1:], norms_p, label='Percentile={:.1f}'.format(p))
            if estimate_slopes:
                coefs = np.polyfit(np.log(self.timesteps[1:]), np.log(norms_p), 1)
                print("Estimated slope for percentile {:.1f} is {}.".format(p, coefs[0]))
        ax.set_xlabel('Timestep (sec)', fontsize='x-large')
        ax.set_ylabel('Norm of difference in {} from $\Delta$t={}'.format(name, self.timesteps[0]),
                      fontsize='x-large')
        ax.legend(loc='best', fontsize='x-large')
