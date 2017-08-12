from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pyplot
import socket
import os
from analysis_tools import *

hostname = socket.gethostname()
if ("tux" in hostname):
  wrkdir = os.getenv('HOME') + '/workspace/micro_physics/KiD_2.3.2654/output'

runDirectory = wrkdir + '/warm1_v1'
runDirectory_adaptive = wrkdir + '/warm1_v1'
runDirectory_adaptive2 = wrkdir + '/warm1_v1'

quantityDict = {'Rain Mass':Quantity('rain_mass'),
                'Cloud Mass':Quantity('cloud_mass')}

runDict = {'mstep1':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep1.nc',1,quantityList=quantityDict.values()),
           'mstep1_adaptive':Run(runDirectory_adaptive + '/warm1_mg2_acme_v1beta_dt1.0_mstep1.nc',1,quantityList=quantityDict.values()),
           'mstep1_adaptive2':Run(runDirectory_adaptive2 + '/warm1_mg2_acme_v1beta_dt1.0_mstep1.nc',1,quantityList=quantityDict.values()),
           'mstep5':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep5.nc',5,quantityList=quantityDict.values()),
           'mstep5_adaptive':Run(runDirectory_adaptive + '/warm1_mg2_acme_v1beta_dt1.0_mstep5.nc',5,quantityList=quantityDict.values()),
           'mstep5_adaptive2':Run(runDirectory_adaptive2 + '/warm1_mg2_acme_v1beta_dt1.0_mstep5.nc',5,quantityList=quantityDict.values()),
           'mstep10':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep10.nc',10,quantityList=quantityDict.values()),
           'mstep10_adaptive':Run(runDirectory_adaptive + '/warm1_mg2_acme_v1beta_dt1.0_mstep10.nc',10,quantityList=quantityDict.values()),
           'mstep10_adaptive2':Run(runDirectory_adaptive2 + '/warm1_mg2_acme_v1beta_dt1.0_mstep10.nc',10,quantityList=quantityDict.values()),
           'mstep15':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep15.nc',15,quantityList=quantityDict.values()),
           'mstep15_adaptive':Run(runDirectory_adaptive + '/warm1_mg2_acme_v1beta_dt1.0_mstep15.nc',15,quantityList=quantityDict.values()),
           'mstep15_adaptive2':Run(runDirectory_adaptive2 + '/warm1_mg2_acme_v1beta_dt1.0_mstep15.nc',15,quantityList=quantityDict.values()),
           'mstep30':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep30.nc',30,quantityList=quantityDict.values()),
           'mstep30_adaptive':Run(runDirectory_adaptive + '/warm1_mg2_acme_v1beta_dt1.0_mstep30.nc',30,quantityList=quantityDict.values()),
           'mstep30_adaptive2':Run(runDirectory_adaptive2 + '/warm1_mg2_acme_v1beta_dt1.0_mstep30.nc',30,quantityList=quantityDict.values()),
           'mstep60':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep60.nc',60,quantityList=quantityDict.values()),
           'mstep60_adaptive':Run(runDirectory_adaptive + '/warm1_mg2_acme_v1beta_dt1.0_mstep60.nc',60,quantityList=quantityDict.values()),
           'mstep60_adaptive2':Run(runDirectory_adaptive2 + '/warm1_mg2_acme_v1beta_dt1.0_mstep60.nc',60,quantityList=quantityDict.values()),
           'mstep120':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep120.nc',120,quantityList=quantityDict.values()),
           'mstep120_adaptive':Run(runDirectory_adaptive + '/warm1_mg2_acme_v1beta_dt1.0_mstep120.nc',120,quantityList=quantityDict.values()),
           'mstep120_adaptive2':Run(runDirectory_adaptive2 + '/warm1_mg2_acme_v1beta_dt1.0_mstep120.nc',120,quantityList=quantityDict.values()),
           'mstep300':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep300.nc',300,quantityList=quantityDict.values()),
           'mstep300_adaptive':Run(runDirectory_adaptive + '/warm1_mg2_acme_v1beta_dt1.0_mstep300.nc',300,quantityList=quantityDict.values()),
           'mstep300_adaptive2':Run(runDirectory_adaptive2 + '/warm1_mg2_acme_v1beta_dt1.0_mstep300.nc',300,quantityList=quantityDict.values()),
           'mstep600':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep600.nc',600,quantityList=quantityDict.values()),
           'mstep600_adaptive':Run(runDirectory_adaptive + '/warm1_mg2_acme_v1beta_dt1.0_mstep600.nc',600,quantityList=quantityDict.values()),
           'mstep600_adaptive2':Run(runDirectory_adaptive2 + '/warm1_mg2_acme_v1beta_dt1.0_mstep600.nc',600,quantityList=quantityDict.values()),
           'mstep900':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep900.nc',900,quantityList=quantityDict.values()),
           'mstep900_adaptive':Run(runDirectory_adaptive + '/warm1_mg2_acme_v1beta_dt1.0_mstep900.nc',900,quantityList=quantityDict.values()),
           'mstep900_adaptive2':Run(runDirectory_adaptive2 + '/warm1_mg2_acme_v1beta_dt1.0_mstep900.nc',900,quantityList=quantityDict.values()),
           'mstep1200':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep1200.nc',1200,quantityList=quantityDict.values()),
           'mstep1200_adaptive':Run(runDirectory_adaptive + '/warm1_mg2_acme_v1beta_dt1.0_mstep1200.nc',1200,quantityList=quantityDict.values()),
           'mstep1200_adaptive2':Run(runDirectory_adaptive2 + '/warm1_mg2_acme_v1beta_dt1.0_mstep1200.nc',1200,quantityList=quantityDict.values())
           }

comparisonSizes = (1,30,120,300)
convergenceSizes = (5,10,15,30,60,120,300,600,900,1200)



for l,size in enumerate(comparisonSizes):
  f, axarray = pyplot.subplots(len(quantityDict.keys()),3,figsize=(18,10),num=l)

  for k,qname in enumerate(quantityDict.keys()):
    qABS = np.amax(abs(quantityDict[qname].getQ(runDict['mstep1'])))
    original = runDict['mstep'+str(size)]
    adaptive = runDict['mstep'+str(size)+'_adaptive']
    adaptive2 = runDict['mstep'+str(size)+'_adaptive2']
    q = quantityDict[qname].getQ(original)
    q_adaptive = quantityDict[qname].getQ(adaptive)
    q_adaptive2 = quantityDict[qname].getQ(adaptive2)
    [t,z] = np.meshgrid(np.arange(0,3600,original.dt),np.linspace(0,120,120))
    cax = axarray[k,0].pcolor(t,z,q,cmap='bwr',vmin=-qABS,vmax=qABS)
    f.colorbar(cax,ax=axarray[k,0])
    cax = axarray[k,1].pcolor(t,z,q_adaptive,cmap='bwr',vmin=-qABS,vmax=qABS)
    f.colorbar(cax,ax=axarray[k,1])
    cax = axarray[k,2].pcolor(t,z,q_adaptive2,cmap='bwr',vmin=-qABS,vmax=qABS)
    f.colorbar(cax,ax=axarray[k,2])
    axarray[k,0].set_title(qname + ' (none)')
    axarray[k,1].set_title(qname + ' (rain)')
    axarray[k,2].set_title(qname + ' (rain+cloud)')

    for i in range(len(quantityDict.keys())):
      for j in range(3):
        axarray[i,j].set_xlabel('time (s)')
        axarray[i,j].set_ylabel('height')
    pyplot.savefig('mstep'+str(size)+'.png')


for k,qname in enumerate(quantityDict.keys()):
  # Final Time Error Holders
  l2ErrorFinalTime = np.zeros((3,len(convergenceSizes)))
  lIErrorFinalTime = np.zeros((3,len(convergenceSizes)))
  # All Time Error Holders
  l2ErrorAllTime = np.zeros((3,len(convergenceSizes)))
  lIErrorAllTime = np.zeros((3,len(convergenceSizes)))
  dtlist = np.empty(len(convergenceSizes))
  for j,size in enumerate(convergenceSizes):
    original = runDict['mstep'+str(size)]
    adaptive1 = runDict['mstep'+str(size)+'_adaptive']
    adaptive2 = runDict['mstep'+str(size)+'_adaptive2']
    q0 = quantityDict[qname].getQ(original)
    q1 = quantityDict[qname].getQ(adaptive1)
    q2 = quantityDict[qname].getQ(adaptive2)
    qRef = quantityDict[qname].getQRef(original,runDict['mstep1'])
    dtlist[j] = original.dt
    # Final time errors
    l2ErrorFinalTime[0,j] = np.sqrt(np.sum(np.square(q0[:,-1]-qRef[:,-1]))/120)
    l2ErrorFinalTime[1,j] = np.sqrt(np.sum(np.square(q1[:,-1]-qRef[:,-1]))/120)
    l2ErrorFinalTime[2,j] = np.sqrt(np.sum(np.square(q2[:,-1]-qRef[:,-1]))/120)
    lIErrorFinalTime[0,j] = np.amax(abs(q0[:,-1]-qRef[:,-1]))
    lIErrorFinalTime[1,j] = np.amax(abs(q1[:,-1]-qRef[:,-1]))
    lIErrorFinalTime[2,j] = np.amax(abs(q2[:,-1]-qRef[:,-1]))
    # All time errors
    l2ErrorAllTime[0,j] = np.sqrt(np.sum(np.square(q0-qRef))/120*original.dt)
    l2ErrorAllTime[1,j] = np.sqrt(np.sum(np.square(q1-qRef))/120*original.dt)
    l2ErrorAllTime[2,j] = np.sqrt(np.sum(np.square(q2-qRef))/120*original.dt)
    lIErrorAllTime[0,j] = np.amax(abs(q0-qRef))
    lIErrorAllTime[1,j] = np.amax(abs(q1-qRef))
    lIErrorAllTime[2,j] = np.amax(abs(q2-qRef))

  f,axarray = pyplot.subplots(2,2,figsize=(10,10),num=len(comparisonSizes)+k)
  coeff = l2ErrorFinalTime[0,0]/dtlist[0]
  axarray[0,0].loglog(dtlist,l2ErrorFinalTime[0,:],'-o',label='no adaptive stepping')
  axarray[0,0].loglog(dtlist,l2ErrorFinalTime[1,:],'-o',label='adaptive rain')
  axarray[0,0].loglog(dtlist,l2ErrorFinalTime[2,:],'-o',label='adaptive rain & cloud')
  axarray[0,0].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
  axarray[0,0].axis('equal')
  axarray[0,0].set_xlabel('dt')
  axarray[0,0].set_title(qname + ' Error (L2, Final Time)')
  axarray[0,0].legend()
  coeff = l2ErrorAllTime[0,0]/dtlist[0]
  axarray[1,0].loglog(dtlist,l2ErrorAllTime[0,:],'-o',label='no adaptive stepping')
  axarray[1,0].loglog(dtlist,l2ErrorAllTime[1,:],'-o',label='adaptive rain')
  axarray[1,0].loglog(dtlist,l2ErrorAllTime[2,:],'-o',label='adaptive rain & cloud')
  axarray[1,0].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
  axarray[1,0].axis('equal')
  axarray[1,0].set_xlabel('dt')
  axarray[1,0].set_title(qname + ' Error (L2, All Time)')
  axarray[1,0].legend()
  coeff = lIErrorFinalTime[0,0]/dtlist[0]
  axarray[0,1].loglog(dtlist,lIErrorFinalTime[0,:],'-o',label='no adaptive stepping')
  axarray[0,1].loglog(dtlist,lIErrorFinalTime[1,:],'-o',label='adaptive rain')
  axarray[0,1].loglog(dtlist,lIErrorFinalTime[2,:],'-o',label='adaptive rain & cloud')
  axarray[0,1].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
  axarray[0,1].axis('equal')
  axarray[0,1].set_xlabel('dt')
  axarray[0,1].set_title(qname + ' Error (Linf, Final Time)')
  axarray[0,1].legend()
  coeff = lIErrorAllTime[0,0]/dtlist[0]
  axarray[1,1].loglog(dtlist,lIErrorAllTime[0,:],'-o',label='no adaptive stepping')
  axarray[1,1].loglog(dtlist,lIErrorAllTime[1,:],'-o',label='adaptive rain')
  axarray[1,1].loglog(dtlist,lIErrorAllTime[2,:],'-o',label='adaptive rain & cloud')
  axarray[1,1].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
  axarray[1,1].axis('equal')
  axarray[1,1].set_xlabel('dt')
  axarray[1,1].set_title(qname + ' Error (Linf, All Time)')
  axarray[1,1].legend()
  f.savefig(quantityDict[qname].name + '_convergence.png')

pyplot.show()
