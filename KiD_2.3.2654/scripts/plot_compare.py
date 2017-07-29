from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pyplot
from analysis_tools import *

limiterDict = {'qr_conserv':Limiter('qr_conserv','qr_conserv_mag','rain_mass','rescale'),
               'nr_conserv':Limiter('nr_conserv','nr_conserv_mag','rain_number','rescale'),
               'qc_conserv':Limiter('qc_conserv','qc_conserv_mag','cloud_mass','rescale'),
               'nc_conserv':Limiter('nc_conserv','nc_conserv_mag','cloud_number','rescale')
               }

runDict = {'mstep1':Run('/home/vogl2/Desktop/temp/warm1_v0/warm1_mg2_acme_v1beta_dt1.0_mstep1.nc',1,limiterDict.values()),
           'mstep1_adaptive':Run('/home/vogl2/Desktop/temp/warm1_v1_R/warm1_mg2_acme_v1beta_dt1.0_mstep1.nc',1,limiterDict.values()),
           'mstep1_adaptive2':Run('/home/vogl2/Desktop/temp/warm1_v2_RC/warm1_mg2_acme_v1beta_dt1.0_mstep1.nc',1,limiterDict.values()),
           'mstep5':Run('/home/vogl2/Desktop/temp/warm1_v0/warm1_mg2_acme_v1beta_dt1.0_mstep5.nc',5,limiterDict.values()),
           'mstep5_adaptive':Run('/home/vogl2/Desktop/temp/warm1_v1_R/warm1_mg2_acme_v1beta_dt1.0_mstep5.nc',5,limiterDict.values()),
           'mstep5_adaptive2':Run('/home/vogl2/Desktop/temp/warm1_v2_RC/warm1_mg2_acme_v1beta_dt1.0_mstep5.nc',5,limiterDict.values()),
           'mstep10':Run('/home/vogl2/Desktop/temp/warm1_v0/warm1_mg2_acme_v1beta_dt1.0_mstep10.nc',10,limiterDict.values()),
           'mstep10_adaptive':Run('/home/vogl2/Desktop/temp/warm1_v1_R/warm1_mg2_acme_v1beta_dt1.0_mstep10.nc',10,limiterDict.values()),
           'mstep10_adaptive2':Run('/home/vogl2/Desktop/temp/warm1_v2_RC/warm1_mg2_acme_v1beta_dt1.0_mstep10.nc',10,limiterDict.values()),
           'mstep15':Run('/home/vogl2/Desktop/temp/warm1_v0/warm1_mg2_acme_v1beta_dt1.0_mstep15.nc',15,limiterDict.values()),
           'mstep15_adaptive':Run('/home/vogl2/Desktop/temp/warm1_v1_R/warm1_mg2_acme_v1beta_dt1.0_mstep15.nc',15,limiterDict.values()),
           'mstep15_adaptive2':Run('/home/vogl2/Desktop/temp/warm1_v2_RC/warm1_mg2_acme_v1beta_dt1.0_mstep15.nc',15,limiterDict.values()),
           'mstep30':Run('/home/vogl2/Desktop/temp/warm1_v0/warm1_mg2_acme_v1beta_dt1.0_mstep30.nc',30,limiterDict.values()),
           'mstep30_adaptive':Run('/home/vogl2/Desktop/temp/warm1_v1_R/warm1_mg2_acme_v1beta_dt1.0_mstep30.nc',30,limiterDict.values()),
           'mstep30_adaptive2':Run('/home/vogl2/Desktop/temp/warm1_v2_RC/warm1_mg2_acme_v1beta_dt1.0_mstep30.nc',30,limiterDict.values()),
           'mstep60':Run('/home/vogl2/Desktop/temp/warm1_v0/warm1_mg2_acme_v1beta_dt1.0_mstep60.nc',60,limiterDict.values()),
           'mstep60_adaptive':Run('/home/vogl2/Desktop/temp/warm1_v1_R/warm1_mg2_acme_v1beta_dt1.0_mstep60.nc',60,limiterDict.values()),
           'mstep60_adaptive2':Run('/home/vogl2/Desktop/temp/warm1_v2_RC/warm1_mg2_acme_v1beta_dt1.0_mstep60.nc',60,limiterDict.values()),
           'mstep120':Run('/home/vogl2/Desktop/temp/warm1_v0/warm1_mg2_acme_v1beta_dt1.0_mstep120.nc',120,limiterDict.values()),
           'mstep120_adaptive':Run('/home/vogl2/Desktop/temp/warm1_v1_R/warm1_mg2_acme_v1beta_dt1.0_mstep120.nc',120,limiterDict.values()),
           'mstep120_adaptive2':Run('/home/vogl2/Desktop/temp/warm1_v2_RC/warm1_mg2_acme_v1beta_dt1.0_mstep120.nc',120,limiterDict.values()),
           'mstep300':Run('/home/vogl2/Desktop/temp/warm1_v0/warm1_mg2_acme_v1beta_dt1.0_mstep300.nc',300,limiterDict.values()),
           'mstep300_adaptive':Run('/home/vogl2/Desktop/temp/warm1_v1_R/warm1_mg2_acme_v1beta_dt1.0_mstep300.nc',300,limiterDict.values()),
           'mstep300_adaptive2':Run('/home/vogl2/Desktop/temp/warm1_v2_RC/warm1_mg2_acme_v1beta_dt1.0_mstep300.nc',300,limiterDict.values()),
           'mstep600':Run('/home/vogl2/Desktop/temp/warm1_v0/warm1_mg2_acme_v1beta_dt1.0_mstep600.nc',600,limiterDict.values()),
           'mstep600_adaptive':Run('/home/vogl2/Desktop/temp/warm1_v1_R/warm1_mg2_acme_v1beta_dt1.0_mstep600.nc',600,limiterDict.values()),
           'mstep600_adaptive2':Run('/home/vogl2/Desktop/temp/warm1_v2_RC/warm1_mg2_acme_v1beta_dt1.0_mstep600.nc',600,limiterDict.values()),
           'mstep900':Run('/home/vogl2/Desktop/temp/warm1_v0/warm1_mg2_acme_v1beta_dt1.0_mstep900.nc',900,limiterDict.values()),
           'mstep900_adaptive':Run('/home/vogl2/Desktop/temp/warm1_v1_R/warm1_mg2_acme_v1beta_dt1.0_mstep900.nc',900,limiterDict.values()),
           'mstep900_adaptive2':Run('/home/vogl2/Desktop/temp/warm1_v2_RC/warm1_mg2_acme_v1beta_dt1.0_mstep900.nc',900,limiterDict.values()),
           'mstep1200':Run('/home/vogl2/Desktop/temp/warm1_v0/warm1_mg2_acme_v1beta_dt1.0_mstep1200.nc',1200,limiterDict.values()),
           'mstep1200_adaptive':Run('/home/vogl2/Desktop/temp/warm1_v1_R/warm1_mg2_acme_v1beta_dt1.0_mstep1200.nc',1200,limiterDict.values()),
           'mstep1200_adaptive2':Run('/home/vogl2/Desktop/temp/warm1_v2_RC/warm1_mg2_acme_v1beta_dt1.0_mstep1200.nc',1200,limiterDict.values())
           }

comparisonSizes = (1,30,120,300)
convergenceSizes = (5,10,15,30,60,120,300,600,900,1200)

qcABS = np.amax(abs(limiterDict['qc_conserv'].getQ(runDict['mstep1'])))
qrABS = np.amax(abs(limiterDict['qr_conserv'].getQ(runDict['mstep1'])))

if 0:
  for size in comparisonSizes:
    original = runDict['mstep'+str(size)]
    adaptive = runDict['mstep'+str(size)+'_adaptive']
    adaptive2 = runDict['mstep'+str(size)+'_adaptive2']
    qc = limiterDict['qc_conserv'].getQ(original)
    qr = limiterDict['qr_conserv'].getQ(original)
    qc_adaptive = limiterDict['qc_conserv'].getQ(adaptive)
    qr_adaptive = limiterDict['qr_conserv'].getQ(adaptive)
    qc_adaptive2 = limiterDict['qc_conserv'].getQ(adaptive2)
    qr_adaptive2 = limiterDict['qr_conserv'].getQ(adaptive2)
    [t,z] = np.meshgrid(np.arange(0,3600,original.dt),np.linspace(0,120,120))
    f, axarray = pyplot.subplots(2,3,figsize=(18,10))
    cax = axarray[0,0].pcolor(t,z,qc,cmap='bwr',vmin=-qcABS,vmax=qcABS)
    f.colorbar(cax,ax=axarray[0,0])
    cax = axarray[0,1].pcolor(t,z,qc_adaptive,cmap='bwr',vmin=-qcABS,vmax=qcABS)
    f.colorbar(cax,ax=axarray[0,1])
    cax = axarray[0,2].pcolor(t,z,qc_adaptive2,cmap='bwr',vmin=-qcABS,vmax=qcABS)
    f.colorbar(cax,ax=axarray[0,2])
    cax = axarray[1,0].pcolor(t,z,qr,cmap='bwr',vmin=-qrABS,vmax=qrABS)
    f.colorbar(cax,ax=axarray[1,0])
    cax = axarray[1,1].pcolor(t,z,qr_adaptive,cmap='bwr',vmin=-qrABS,vmax=qrABS)
    f.colorbar(cax,ax=axarray[1,1])
    cax = axarray[1,2].pcolor(t,z,qr_adaptive2,cmap='bwr',vmin=-qrABS,vmax=qrABS)
    f.colorbar(cax,ax=axarray[1,2])
    axarray[0,0].set_title('Cloud Mass (none)')
    axarray[0,1].set_title('Cloud Mass (rain)')
    axarray[0,2].set_title('Cloud Mass (rain+cloud)')
    axarray[1,0].set_title('Rain Mass (none)')
    axarray[1,1].set_title('Rain Mass (rain)')
    axarray[1,2].set_title('Rain Mass (rain+cloud)')
    for i in range(2):
      for j in range(3):
        axarray[i,j].set_xlabel('time (s)')
        axarray[i,j].set_ylabel('height')
    pyplot.savefig('mstep'+str(size)+'.png')

# Final Time Error Holders
qcL2ErrorFinalTime = np.zeros((3,len(convergenceSizes)))
qrL2ErrorFinalTime = np.zeros((3,len(convergenceSizes)))
qcLIErrorFinalTime = np.zeros((3,len(convergenceSizes)))
qrLIErrorFinalTime = np.zeros((3,len(convergenceSizes)))
# All Time Error Holders
qcL2ErrorAllTime = np.zeros((3,len(convergenceSizes)))
qrL2ErrorAllTime = np.zeros((3,len(convergenceSizes)))
qcLIErrorAllTime = np.zeros((3,len(convergenceSizes)))
qrLIErrorAllTime = np.zeros((3,len(convergenceSizes)))
dtlist = np.empty(len(convergenceSizes))
for j,size in enumerate(convergenceSizes):
  original = runDict['mstep'+str(size)]
  adaptive1 = runDict['mstep'+str(size)+'_adaptive']
  adaptive2 = runDict['mstep'+str(size)+'_adaptive2']
  qc0 = limiterDict['qc_conserv'].getQ(original)
  qr0 = limiterDict['qr_conserv'].getQ(original)
  qc1 = limiterDict['qc_conserv'].getQ(adaptive1)
  qr1 = limiterDict['qr_conserv'].getQ(adaptive1)
  qc2 = limiterDict['qc_conserv'].getQ(adaptive2)
  qr2 = limiterDict['qr_conserv'].getQ(adaptive2)
  qcRef = limiterDict['qc_conserv'].getQRef(original,runDict['mstep1'])
  qrRef = limiterDict['qr_conserv'].getQRef(original,runDict['mstep1'])
  dtlist[j] = original.dt
  # Final time errors
  qcL2ErrorFinalTime[0,j] = np.sqrt(np.sum(np.square(qc0[:,-1]-qcRef[:,-1]))/120)
  qcL2ErrorFinalTime[1,j] = np.sqrt(np.sum(np.square(qc1[:,-1]-qcRef[:,-1]))/120)
  qcL2ErrorFinalTime[2,j] = np.sqrt(np.sum(np.square(qc2[:,-1]-qcRef[:,-1]))/120)
  qrL2ErrorFinalTime[0,j] = np.sqrt(np.sum(np.square(qr0[:,-1]-qrRef[:,-1]))/120)
  qrL2ErrorFinalTime[1,j] = np.sqrt(np.sum(np.square(qr1[:,-1]-qrRef[:,-1]))/120)
  qrL2ErrorFinalTime[2,j] = np.sqrt(np.sum(np.square(qr2[:,-1]-qrRef[:,-1]))/120)
  qcLIErrorFinalTime[0,j] = np.amax(abs(qc0[:,-1]-qcRef[:,-1]))
  qcLIErrorFinalTime[1,j] = np.amax(abs(qc1[:,-1]-qcRef[:,-1]))
  qcLIErrorFinalTime[2,j] = np.amax(abs(qc2[:,-1]-qcRef[:,-1]))
  qrLIErrorFinalTime[0,j] = np.amax(abs(qr0[:,-1]-qrRef[:,-1]))
  qrLIErrorFinalTime[1,j] = np.amax(abs(qr1[:,-1]-qrRef[:,-1]))
  qrLIErrorFinalTime[2,j] = np.amax(abs(qr2[:,-1]-qrRef[:,-1]))
  # All time errors
  qcL2ErrorAllTime[0,j] = np.sqrt(np.sum(np.square(qc0-qcRef))/120/3601)
  qcL2ErrorAllTime[1,j] = np.sqrt(np.sum(np.square(qc1-qcRef))/120/3601)
  qcL2ErrorAllTime[2,j] = np.sqrt(np.sum(np.square(qc2-qcRef))/120/3601)
  qrL2ErrorAllTime[0,j] = np.sqrt(np.sum(np.square(qr0-qrRef))/120/3601)
  qrL2ErrorAllTime[1,j] = np.sqrt(np.sum(np.square(qr1-qrRef))/120/3601)
  qrL2ErrorAllTime[2,j] = np.sqrt(np.sum(np.square(qr2-qrRef))/120/3601)
  qcLIErrorAllTime[0,j] = abs(np.amax(qc0-qcRef))
  qcLIErrorAllTime[1,j] = abs(np.amax(qc1-qcRef))
  qcLIErrorAllTime[2,j] = abs(np.amax(qc2-qcRef))
  qrLIErrorAllTime[0,j] = abs(np.amax(qr0-qrRef))
  qrLIErrorAllTime[1,j] = abs(np.amax(qr1-qrRef))
  qrLIErrorAllTime[2,j] = abs(np.amax(qr2-qrRef))

f,axarray = pyplot.subplots(2,2,num=1,figsize=(10,10))
coeff = qcL2ErrorFinalTime[0,0]/dtlist[0]
axarray[0,0].loglog(dtlist,qcL2ErrorFinalTime[0,:],'-o',label='no adaptive stepping')
axarray[0,0].loglog(dtlist,qcL2ErrorFinalTime[1,:],'-o',label='adaptive rain')
axarray[0,0].loglog(dtlist,qcL2ErrorFinalTime[2,:],'-o',label='adaptive rain & cloud')
axarray[0,0].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
axarray[0,0].axis('equal')
axarray[0,0].set_xlabel('dt')
axarray[0,0].set_title('Cloud Mass Error (L2, Final Time)')
axarray[0,0].legend()
coeff = qcL2ErrorAllTime[0,0]/dtlist[0]
axarray[1,0].loglog(dtlist,qcL2ErrorAllTime[0,:],'-o',label='no adaptive stepping')
axarray[1,0].loglog(dtlist,qcL2ErrorAllTime[1,:],'-o',label='adaptive rain')
axarray[1,0].loglog(dtlist,qcL2ErrorAllTime[2,:],'-o',label='adaptive rain & cloud')
axarray[1,0].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
axarray[1,0].axis('equal')
axarray[1,0].set_xlabel('dt')
axarray[1,0].set_title('Cloud Mass Error (L2, All Time)')
axarray[1,0].legend()
coeff = qcLIErrorFinalTime[0,0]/dtlist[0]
axarray[0,1].loglog(dtlist,qcLIErrorFinalTime[0,:],'-o',label='no adaptive stepping')
axarray[0,1].loglog(dtlist,qcLIErrorFinalTime[1,:],'-o',label='adaptive rain')
axarray[0,1].loglog(dtlist,qcLIErrorFinalTime[2,:],'-o',label='adaptive rain & cloud')
axarray[0,1].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
axarray[0,1].axis('equal')
axarray[0,1].set_xlabel('dt')
axarray[0,1].set_title('Cloud Mass Error (Linf, Final Time)')
axarray[0,1].legend()
coeff = qcLIErrorAllTime[0,0]/dtlist[0]
axarray[1,1].loglog(dtlist,qcLIErrorAllTime[0,:],'-o',label='no adaptive stepping')
axarray[1,1].loglog(dtlist,qcLIErrorAllTime[1,:],'-o',label='adaptive rain')
axarray[1,1].loglog(dtlist,qcLIErrorAllTime[2,:],'-o',label='adaptive rain & cloud')
axarray[1,1].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
axarray[1,1].axis('equal')
axarray[1,1].set_xlabel('dt')
axarray[1,1].set_title('Cloud Mass Error (Linf, All Time)')
axarray[1,1].legend()
f.savefig('cloud_mass_convergence.png')


f,axarray = pyplot.subplots(2,2,num=2,figsize=(10,10))
coeff = qrL2ErrorFinalTime[0,0]/dtlist[0]
axarray[0,0].loglog(dtlist,qrL2ErrorFinalTime[0,:],'-o',label='no adaptive stepping')
axarray[0,0].loglog(dtlist,qrL2ErrorFinalTime[1,:],'-o',label='adaptive rain')
axarray[0,0].loglog(dtlist,qrL2ErrorFinalTime[2,:],'-o',label='adaptive rain & cloud')
axarray[0,0].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
axarray[0,0].axis('equal')
axarray[0,0].set_xlabel('dt')
axarray[0,0].set_title('Rain Mass Error (L2, Final Time)')
axarray[0,0].legend()
coeff = qrL2ErrorAllTime[0,0]/dtlist[0]
axarray[1,0].loglog(dtlist,qrL2ErrorAllTime[0,:],'-o',label='no adaptive stepping')
axarray[1,0].loglog(dtlist,qrL2ErrorAllTime[1,:],'-o',label='adaptive rain')
axarray[1,0].loglog(dtlist,qrL2ErrorAllTime[2,:],'-o',label='adaptive rain & cloud')
axarray[1,0].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
axarray[1,0].axis('equal')
axarray[1,0].set_xlabel('dt')
axarray[1,0].set_title('Rain Mass Error (L2, All Time)')
axarray[1,0].legend()
coeff = qrLIErrorFinalTime[0,0]/dtlist[0]
axarray[0,1].loglog(dtlist,qrLIErrorFinalTime[0,:],'-o',label='no adaptive stepping')
axarray[0,1].loglog(dtlist,qrLIErrorFinalTime[1,:],'-o',label='adaptive rain')
axarray[0,1].loglog(dtlist,qrLIErrorFinalTime[2,:],'-o',label='adaptive rain & cloud')
axarray[0,1].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
axarray[0,1].axis('equal')
axarray[0,1].set_xlabel('dt')
axarray[0,1].set_title('Rain Mass Error (Linf, Final Time)')
axarray[0,1].legend()
coeff = qrLIErrorAllTime[0,0]/dtlist[0]
axarray[1,1].loglog(dtlist,qrLIErrorAllTime[0,:],'-o',label='no adaptive stepping')
axarray[1,1].loglog(dtlist,qrLIErrorAllTime[1,:],'-o',label='adaptive rain')
axarray[1,1].loglog(dtlist,qrLIErrorAllTime[2,:],'-o',label='adaptive rain & cloud')
axarray[1,1].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
axarray[1,1].axis('equal')
axarray[1,1].set_xlabel('dt')
axarray[1,1].set_title('Rain Mass Error (Linf, All Time)')
axarray[1,1].legend()
f.savefig('rain_mass_convergence.png')

pyplot.show()
