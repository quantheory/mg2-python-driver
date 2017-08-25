from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pyplot
import socket
import os
from analysis_tools import *

hostname = socket.gethostname()
if ("tux" in hostname):
  wrkdir = os.getenv('HOME') + '/workspace/micro_physics/KiD_2.3.2654/output'

# Dictionary for sedimentation methods (name: directory)
methodDict = {'v0: original': wrkdir + '/warm1_v0',
              'v1: time-varying speed': wrkdir + '/warm1_v1',
              'v2: v1 with algebra mod': wrkdir + '/warm1_v2'}

# Dictionary for quantities of interest (name: Quantity object)
quantityDict = {'Rain Mass':Quantity('rain_mass'),
                'Vapour':Quantity('vapour'),
                'Cloud Mass': Quantity('cloud_mass') }

# Dictionary for runs (method name: individual run dictionary)
runDict = {}
for method in methodDict.keys():
  runDirectory = methodDict[method]
  # dictionary for runs (run name: run directory)
  runs = {'mstep1':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep1.nc',1,quantityList=quantityDict.values()),
             'mstep5':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep5.nc',5,quantityList=quantityDict.values()),
             'mstep10':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep10.nc',10,quantityList=quantityDict.values()),
             'mstep15':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep15.nc',15,quantityList=quantityDict.values()),
             'mstep30':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep30.nc',30,quantityList=quantityDict.values()),
             'mstep60':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep60.nc',60,quantityList=quantityDict.values()),
             'mstep120':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep120.nc',120,quantityList=quantityDict.values()),
             'mstep300':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep300.nc',300,quantityList=quantityDict.values()),
             'mstep600':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep600.nc',600,quantityList=quantityDict.values()),
             'mstep900':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep900.nc',900,quantityList=quantityDict.values()),
             'mstep1200':Run(runDirectory + '/warm1_mg2_acme_v1beta_dt1.0_mstep1200.nc',1200,quantityList=quantityDict.values()),
             }
  runDict[method] = runs

# Dictionary for timings (quantity name: list of timings)
timingDict = {'rain': ['Sedimentation']}

# List of resolutions for plotting solutions across methods
comparisonSizes = (30,120,300) #1,30,120,300)

# List of resolutions for convergence plot
convergenceSizes = (5,10,15,30,60,120,300,600,900,1200)

# Plot timing information
if 1:
  for l, quantity in enumerate(timingDict.keys()):
    f, axarray = pyplot.subplots(2,2,figsize=(18,10),num=l)
    timingList = timingDict[quantity]

    for method in methodDict.keys():
      runDirectory = methodDict[method]
      currentRunDict = runDict[method]
      callNumbers = np.zeros((len(convergenceSizes),len(timingList)))
      wallClocks = np.zeros((len(convergenceSizes),len(timingList)))
      ratios = np.zeros((len(convergenceSizes),len(timingList)))
      dtlist = np.empty(len(convergenceSizes))

      for j,size in enumerate(convergenceSizes):
        dtlist[j] = currentRunDict['mstep'+str(size)].dt
        timingFile = open(runDirectory+'/warm1_mg2_acme_v1beta_dt1.0_mstep'+str(size)+'_timing.txt')
        lines = list(timingFile)

        for line in lines:
          for k,currentTiming in enumerate(timingList):
            if (line.strip().startswith(currentTiming+' ('+quantity+')')):
              words = line.split('('+quantity+')')
              tmp = words[1].split()
              callNumbers[j,k] = int(tmp[0])
              wallClocks[j,k] = float(tmp[2])
              ratios[j,k] = wallClocks[j,k]/callNumbers[j,k]
      axarray[0,0].plot(dtlist,wallClocks[:,0],'-o',label=method)
      axarray[0,1].plot(dtlist,callNumbers[:,0],'-o',label=method)
      axarray[1,0].plot(dtlist,ratios[:,0],'-o',label=method)

    axarray[0,0].set_ylabel('Wall Clock')
    axarray[0,0].legend()
    axarray[0,1].set_ylabel('Call #')
    axarray[0,1].legend()
    axarray[1,0].set_ylabel('Ratio')
    axarray[1,0].legend()

  pyplot.show()
  exit()


# Only plot solutions across methods if there aren't many methods
numMethods = len(methodDict.keys())
if (numMethods <= 3):

  for l,size in enumerate(comparisonSizes):
    f, axarray = pyplot.subplots(len(quantityDict.keys()),numMethods,figsize=(18,10),num=l)

    for j,qname in enumerate(quantityDict.keys()):
      currentRunDict = runDict[methodDict.keys()[0]]
      qABS = np.amax(abs(quantityDict[qname].getQ(currentRunDict['mstep1'])))

      for k,method in enumerate(methodDict.keys()):
        currentRunDict = runDict[method]
        currentRun = currentRunDict['mstep'+str(size)]
        q = quantityDict[qname].getQ(currentRun)
        [t,z] = np.meshgrid(np.arange(0,3600,currentRun.dt),np.linspace(0,120,120))
        cax = axarray[j,k].pcolor(t,z,q,cmap='bwr',vmin=-qABS,vmax=qABS)
        f.colorbar(cax,ax=axarray[j,k])
        axarray[j,k].set_title(method)

    for i in range(len(quantityDict.keys())):
      for j in range(numMethods):
        axarray[i,j].set_xlabel('time (s)')
        axarray[i,j].set_ylabel('height')
    pyplot.savefig('mstep'+str(size)+'.png')

# Plot convergence of Final Time / All Time errors in L2 / LInf norms
for k,qname in enumerate(quantityDict.keys()):

  f,axarray = pyplot.subplots(2,2,figsize=(10,10),num=len(comparisonSizes)+k)
  # Final Time Error Holders
  l2ErrorFinalTime = np.zeros((len(methodDict.keys()),len(convergenceSizes)))
  lIErrorFinalTime = np.zeros((len(methodDict.keys()),len(convergenceSizes)))
  # All Time Error Holders
  l2ErrorAllTime = np.zeros((len(methodDict.keys()),len(convergenceSizes)))
  lIErrorAllTime = np.zeros((len(methodDict.keys()),len(convergenceSizes)))
  dtlist = np.empty(len(convergenceSizes))

  for i,method in enumerate(methodDict.keys()):
    currentRunDict = runDict[method]
    for j,size in enumerate(convergenceSizes):
      currentRun = currentRunDict['mstep'+str(size)]
      q = quantityDict[qname].getQ(currentRun)
      qRef = quantityDict[qname].getQRef(currentRun,currentRunDict['mstep1'])
      dtlist[j] = currentRun.dt
      # Final time errors
      l2ErrorFinalTime[i,j] = np.sqrt(np.sum(np.square(q[:,-1]-qRef[:,-1]))/120)
      lIErrorFinalTime[i,j] = np.amax(abs(q[:,-1]-qRef[:,-1]))
      # All time errors
      l2ErrorAllTime[i,j] = np.sqrt(np.sum(np.square(q-qRef))/120*currentRun.dt)
      lIErrorAllTime[i,j] = np.amax(abs(q-qRef))
      # Plots
    axarray[0,0].loglog(dtlist,l2ErrorFinalTime[i,:],'-o',label=method)
    axarray[1,0].loglog(dtlist,l2ErrorAllTime[i,:],'-o',label=method)
    axarray[0,1].loglog(dtlist,lIErrorFinalTime[i,:],'-o',label=method)
    axarray[1,1].loglog(dtlist,lIErrorAllTime[i,:],'-o',label=method)

  coeff = l2ErrorFinalTime[0,0]/dtlist[0]
  #axarray[0,0].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
  axarray[0,0].axis('equal')
  axarray[0,0].set_xlabel('dt')
  axarray[0,0].set_title(qname + ' Error (L2, Final Time)')
  axarray[0,0].legend()
  coeff = l2ErrorAllTime[0,0]/dtlist[0]
  #axarray[1,0].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
  axarray[1,0].axis('equal')
  axarray[1,0].set_xlabel('dt')
  axarray[1,0].set_title(qname + ' Error (L2, All Time)')
  axarray[1,0].legend()
  coeff = lIErrorFinalTime[0,0]/dtlist[0]
  #axarray[0,1].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
  axarray[0,1].axis('equal')
  axarray[0,1].set_xlabel('dt')
  axarray[0,1].set_title(qname + ' Error (Linf, Final Time)')
  axarray[0,1].legend()
  coeff = lIErrorAllTime[0,0]/dtlist[0]
  #axarray[1,1].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
  axarray[1,1].axis('equal')
  axarray[1,1].set_xlabel('dt')
  axarray[1,1].set_title(qname + ' Error (Linf, All Time)')
  axarray[1,1].legend()
  #f.savefig(quantityDict[qname].name + '_convergence.png')
    # tmp = np.concatenate(([dtlist],l2ErrorFinalTime),axis=0)
    # np.savetxt(quantityDict[qname].name + '_l2ErrorFinalTime.txt', tmp.T)
    # tmp = np.concatenate(([dtlist],l2ErrorAllTime),axis=0)
    # np.savetxt(quantityDict[qname].name + '_l2ErrorAllTime.txt', tmp.T)
    # tmp = np.concatenate(([dtlist],lIErrorFinalTime),axis=0)
    # np.savetxt(quantityDict[qname].name + '_lIErrorFinalTime.txt', tmp.T)
    # tmp = np.concatenate(([dtlist],lIErrorAllTime),axis=0)
    # np.savetxt(quantityDict[qname].name + '_lIErrorAllTime.txt', tmp.T)

pyplot.show()
