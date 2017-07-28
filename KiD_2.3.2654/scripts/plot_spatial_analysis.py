from netCDF4 import Dataset
import matplotlib.pyplot as pyplot
import numpy as np
from analysis_tools import *

T = 3600
Z = 120
L = 120
cmap = 'bwr'

#limiterList = getDefaultLimiterList()
limiterList = [Limiter('qr_conserv','qr_conserv_mag','rain_mass','rescale'),
               Limiter('nr_conserv','nr_conserv_mag','rain_number','rescale'),
               Limiter('qc_conserv','qc_conserv_mag','cloud_mass','rescale'),
               Limiter('nc_conserv','nc_conserv_mag','cloud_number','rescale')]
#runList = getDefaultRunList(limiterList)
runList = [Run('warm1_mg2_acme_v1beta_dt1.0_mstep30.nc',30,limiterList),
           Run('warm1_mg2_acme_v1beta_dt1.0_mstep120.nc',120,limiterList),
           Run('warm1_mg2_acme_v1beta_dt1.0_mstep300.nc',300,limiterList)]

refrun = Run('warm1_mg2_acme_v1beta_dt1.0_mstep1.nc',1,limiterList)

for run in runList:
    f1, axarray1 = pyplot.subplots(3,len(limiterList),sharex=True,sharey=True,
                                    num=run.dt,figsize=(16,9))
    f2, axarray2 = pyplot.subplots(3,len(limiterList),sharex=True,sharey=True,
                                    num=run.dt+1,figsize=(16,9))
    for k in range(3):
      axarray1[k,0].set_ylabel('height')
      axarray2[k,0].set_ylabel('height')
    for k,limiter in enumerate(limiterList):
      axarray1[2,k].set_xlabel('time (s)')
      axarray2[2,k].set_xlabel('time (s)')
      [t,z] = np.meshgrid(np.arange(0,T,run.dt),np.linspace(0,Z,L))
      q = limiter.getQ(run)
      qRef = limiter.getQRef(run,refrun)
      limiterError = limiter.getLimiterError(run)
      cumulativeLimiterError = limiter.computeCumulativeLimiterError(run)
      qError = limiter.computeQError(run,refrun)
      # figure for reference solution, limited solution, limiter error
      cax = axarray1[0,k].pcolor(t,z,qRef,cmap=cmap,
            vmin=-np.amax(abs(qRef)),vmax=np.amax(abs(qRef)))
      axarray1[0,k].set_title(limiter.qName + ' (reference)')
      f1.colorbar(cax,ax=axarray1[0,k])
      cax = axarray1[1,k].pcolor(t,z,q,cmap=cmap,
            vmin=-np.amax(abs(qRef)),vmax=np.amax(abs(qRef)))
      axarray1[1,k].set_title(limiter.qName + ' (' + limiter.name + ')')
      f1.colorbar(cax,ax=axarray1[1,k])
      cax = axarray1[2,k].pcolor(t,z,limiterError,cmap=cmap,
            vmin=-np.amax(abs(limiterError)),vmax=np.amax(abs(limiterError)))
      axarray1[2,k].set_title('limiter error')
      f1.colorbar(cax,ax=axarray1[2,k])
      # figure for limited solution, q error, cumulative error ratio
      cax = axarray2[0,k].pcolor(t,z,q,cmap=cmap,
            vmin=-np.amax(abs(qRef)),vmax=np.amax(abs(qRef)))
      axarray2[0,k].set_title(limiter.qName + ' (' + limiter.name + ')')
      f2.colorbar(cax,ax=axarray2[0,k])
      cax = axarray2[1,k].pcolor(t,z,qError,cmap=cmap,
            vmin=-np.amax(qError),vmax=np.amax(qError))
      axarray2[1,k].set_title(limiter.qName + ' error')
      f2.colorbar(cax,ax=axarray2[1,k])
      value = cumulativeLimiterError/(qError+1e-12)
      maxvalue = np.amax(np.abs(value))
      cax = axarray2[2,k].pcolor(t,z,value,cmap=cmap,
            vmin=-min(1.0,maxvalue),vmax=min(1.0,maxvalue))
      axarray2[2,k].set_title('accumulated limiter error')
      f2.colorbar(cax,ax=axarray2[2,k])

pyplot.show()
