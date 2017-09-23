from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pyplot
import socket
import os
import sys
from analysis_tools import *

wrkdir = os.getenv('PWD')
if (len(sys.argv) > 1):
  case = sys.argv[1]
else:
  case = 'warm1'


# Flags for what to plot
PLOT_PROPAGATION_SPEED = False
PLOT_TIMING_INFO = False
PLOT_SOLUTIONS = True
PLOT_CONVERGENCE = True
# Dictionary for sedimentation methods (name: directory)
methodDict = {#'v0: original': wrkdir + '/' + case + '_v0',
#             'v1: time-varying speed': wrkdir + '/' + case + '_v1',
             'v2: v1 with algebra mod': wrkdir + '/' + case + '_v2',
#             'v3: v2 with nonlin rain': wrkdir + '/' + case + '_v3',
#             'v4: v1 with comp flag': wrkdir + '/' + case +'_v4',
#             'v5: v4 with algebra mod': wrkdir + '/' + case + '_v5',
#             'v6: v5 with nonlin rain': wrkdir + '/' + case + '_v6'}
             'v7: v2 with WPA': wrkdir + '/' + case + '_v7',
             'v?': wrkdir + '/' + case}
             #'v8: v7 with nonlin rain': wrkdir + '/' + case + '_v8'}

# Dictionary for quantities of interest (name: Quantity object)
quantityDict = {'Rain Mass':Quantity('rain_mass'),
                'Vapour':Quantity('vapour'),
                'Cloud Mass': Quantity('cloud_mass') }

# Dictionary for runs (method name: individual run dictionary)
runDict = {}
for method in methodDict.keys():
  runDirectory = methodDict[method]
  # dictionary for runs (run name: run directory)
  runs = {'mstep1':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep1.nc',1,quantityList=quantityDict.values()),
             'mstep5':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep5.nc',5,quantityList=quantityDict.values()),
             'mstep10':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep10.nc',10,quantityList=quantityDict.values()),
             'mstep15':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep15.nc',15,quantityList=quantityDict.values()),
             'mstep30':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep30.nc',30,quantityList=quantityDict.values()),
             'mstep60':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep60.nc',60,quantityList=quantityDict.values()),
             'mstep120':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep120.nc',120,quantityList=quantityDict.values()),
             'mstep300':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep300.nc',300,quantityList=quantityDict.values()),
             'mstep400':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep600.nc',400,quantityList=quantityDict.values()),
             'mstep600':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep600.nc',600,quantityList=quantityDict.values()),
             'mstep900':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep900.nc',900,quantityList=quantityDict.values()),
             'mstep1200':Run(runDirectory + '/' + case + '_mg2_acme_v1beta_dt1.0_mstep1200.nc',1200,quantityList=quantityDict.values()),
             }
  runDict[method] = runs

# Dictionary for timings (quantity name: list of timings)
timingDict = {'rain': ['Lambda Calculation', 'Fall Speed Calculation', 'Flux Calculation']}
gptlLoopNumbers = [10000, 10000, 10000]

# List of resolutions for plotting solutions across methods
comparisonSizes = (5,120,300,600) #1,30,120,300)

# List of resolutions for convergence plot
convergenceSizes = (5,10,15,30,60,120,300,400,600,900,1200)

# Set simulation time based on case
if (case == 'warm1' or case == 'warm3'):
  T = 3600
elif (case == 'warm2'):
  T = 7200

# Plot propagation speed information
if (PLOT_PROPAGATION_SPEED):
  f1, axarray1 = pyplot.subplots(2,2,figsize=(10,10))
  f2, axarray2 = pyplot.subplots(2,len(comparisonSizes),figsize=(12,8))

  for method in methodDict.keys():
    runDirectory = methodDict[method]
    currentRunDict = runDict[method]
    dtlist = np.empty(len(convergenceSizes))
    fallSpeedMaxs = np.empty((4,len(convergenceSizes)))
    for j,size in enumerate(convergenceSizes):
      dtlist[j] = currentRunDict['mstep'+str(size)].dt
      sedInfo = np.loadtxt(runDirectory+'/'+case+'_mg2_acme_v1beta_dt1.0_mstep'+str(size)+'_sedinfo.txt',
                                skiprows=1)
      fallSpeedMaxs[0,j] = np.max(sedInfo[:,2])
      fallSpeedMaxs[1,j] = np.max(sedInfo[:,4])
      fallSpeedMaxs[2,j] = np.max(sedInfo[:,6])
      fallSpeedMaxs[3,j] = np.max(sedInfo[:,8])

    axarray1[0,0].semilogx(dtlist,fallSpeedMaxs[0,:],'-o',label=method)
    axarray1[0,1].semilogx(dtlist,fallSpeedMaxs[1,:],'-o',label=method)
    axarray1[1,0].semilogx(dtlist,fallSpeedMaxs[2,:],'-o',label=method)
    axarray1[1,1].semilogx(dtlist,fallSpeedMaxs[3,:],'-o',label=method)

    for j,size in enumerate(comparisonSizes):
      dt = currentRunDict['mstep'+str(size)].dt
      sedInfo = np.loadtxt(runDirectory+'/'+case+'_mg2_acme_v1beta_dt1.0_mstep'+str(size)+'_sedinfo.txt',
                                skiprows=1)
      nstep = T/size + 1
      t = np.linspace(0,T,nstep)
      axarray2[0,j].plot(t,sedInfo[:,4],'-o',label=method)
      axarray2[1,j].plot(t,sedInfo[:,6],'-o',label=method)


  axarray1[0,0].set_title('Max Propagation Speed (Ice)')
  axarray1[0,0].legend()
  axarray1[0,1].set_title('Max Propagation Speed (Cloud)')
  axarray1[0,1].legend()
  axarray1[1,0].set_title('Max Propagation Speed (Rain)')
  axarray1[1,0].legend()
  axarray1[1,1].set_title('Max Propagation Speed (Snow)')
  axarray1[1,1].legend()
  for j,size in enumerate(comparisonSizes):
    axarray2[0,j].set_title(str(size))
    axarray2[0,j].legend()
    axarray2[1,j].set_title(str(size))
    axarray2[1,j].legend()

# Plot timing information
if (PLOT_TIMING_INFO):
  for quantity in timingDict.keys():
    timingList = timingDict[quantity]
    f, axarray = pyplot.subplots(3,len(timingList),figsize=(18,10))

    for method in methodDict.keys():
      runDirectory = methodDict[method]
      currentRunDict = runDict[method]
      callNumbers = np.empty((len(convergenceSizes),len(timingList)))
      wallClocks = np.empty((len(convergenceSizes),len(timingList)))
      totalClock = np.empty(len(convergenceSizes))
      dtlist = np.empty(len(convergenceSizes))

      for j,size in enumerate(convergenceSizes):
        dtlist[j] = currentRunDict['mstep'+str(size)].dt
        timingFile = open(runDirectory+'/'+case+'_mg2_acme_v1beta_dt1.0_mstep'+str(size)+'_timing.txt')
        lines = list(timingFile)

        for line in lines:
          for k,currentTiming in enumerate(timingList):
            if (line.strip().startswith(currentTiming+' ('+quantity+')')):
              words = line.split('('+quantity+')')
              tmp = words[1].split()
              callNumbers[j,k] = int(tmp[0])
              wallClocks[j,k] = float(tmp[2])/gptlLoopNumbers[k]

      for k in range(len(timingList)):
        axarray[0,k].semilogx(dtlist,wallClocks[:,k],'-o',label=method)
        totalClock += wallClocks[:,k]
      axarray[1,0].semilogx(dtlist,callNumbers[:,0],'-o',label=method)
      axarray[1,1].semilogx(dtlist,callNumbers[:,2],'-o',label=method)
      axarray[1,2].semilogx(dtlist,totalClock,'-o',label=method)
      axarray[2,0].semilogx(dtlist,wallClocks[:,0]/callNumbers[:,0],'-o',label=method)
      axarray[2,1].semilogx(dtlist,wallClocks[:,1]/callNumbers[:,1],'-o',label=method)
      axarray[2,2].semilogx(dtlist,totalClock/callNumbers[:,2],'-o',label=method)

    for k, timing in enumerate(timingList):
      axarray[0,k].set_title(quantity + ' ' + timing + ' (Wall Clock)')
      axarray[0,k].legend()
    axarray[1,0].set_title('Number of CFL updates')
    axarray[1,0].legend()
    axarray[1,1].set_title('Sed. Subcycle Steps')
    axarray[1,1].legend()
    axarray[1,2].set_title('Total Sed. Wall Clock')
    axarray[1,2].legend()
    axarray[2,0].set_title('Lambda Wall Clock Per Update')
    axarray[2,0].legend()
    axarray[2,1].set_title('Fall Speed  Clock Per Update')
    axarray[2,1].legend()
    axarray[2,2].set_title('Total Wall Clock Per Subcycle Step')
    axarray[2,2].legend()

# Plot solutions and comparisons across methods
if (PLOT_SOLUTIONS):
  numMethods = len(methodDict.keys())
  for size in comparisonSizes:
    f1, axarray1 = pyplot.subplots(len(quantityDict.keys()),numMethods,figsize=(18,10))
    f2, axarray2 = pyplot.subplots(len(quantityDict.keys()),numMethods,figsize=(18,10))

    for j,qname in enumerate(quantityDict.keys()):
      currentRunDict = runDict[methodDict.keys()[0]]
      refRun = currentRunDict['mstep1']

      for k,method in enumerate(methodDict.keys()):
        currentRunDict = runDict[method]
        currentRun = currentRunDict['mstep'+str(size)]
        q = quantityDict[qname].getQ(currentRun)
        qRef = quantityDict[qname].getQRef(currentRun,refRun)
        qABS = np.amax(abs(qRef))
        [t,z] = np.meshgrid(np.arange(0,T,currentRun.dt),np.linspace(0,120,120))
        cax = axarray1[j,k].pcolor(t,z,q,cmap='bwr',vmin=-qABS,vmax=qABS)
        f1.colorbar(cax,ax=axarray1[j,k])
        axarray1[j,k].set_title(method)
        diff = (q-qRef)/qABS
        diffABS = np.amax(abs(diff))
        cax = axarray2[j,k].pcolor(t,z,diff,cmap='bwr',vmin=-diffABS,vmax=diffABS)
        f2.colorbar(cax,ax=axarray2[j,k])
        axarray2[j,k].set_title(method + '(difference)')

    #f1.savefig('mstep'+str(size)+'.png')
    #f2.savefig('mstep'+str(size)+'_error.png')

# Plot convergence of Final Time / All Time errors in L2 / LInf norms
if (PLOT_CONVERGENCE):
  for qname in quantityDict.keys():

    f,axarray = pyplot.subplots(2,2,figsize=(10,10))
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
    axarray[0,0].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
    axarray[0,0].axis('equal')
    axarray[0,0].set_xlabel('dt')
    axarray[0,0].set_title(qname + ' Error (L2, Final Time)')
    axarray[0,0].legend()
    coeff = l2ErrorAllTime[0,0]/dtlist[0]
    axarray[1,0].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
    axarray[1,0].axis('equal')
    axarray[1,0].set_xlabel('dt')
    axarray[1,0].set_title(qname + ' Error (L2, All Time)')
    axarray[1,0].legend()
    coeff = lIErrorFinalTime[0,0]/dtlist[0]
    axarray[0,1].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
    axarray[0,1].axis('equal')
    axarray[0,1].set_xlabel('dt')
    axarray[0,1].set_title(qname + ' Error (Linf, Final Time)')
    axarray[0,1].legend()
    coeff = lIErrorAllTime[0,0]/dtlist[0]
    axarray[1,1].loglog(dtlist,dtlist*coeff,'--',label='reference: first order')
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
