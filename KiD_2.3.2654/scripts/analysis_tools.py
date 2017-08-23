from netCDF4 import Dataset
import numpy as np

class Quantity():
  def __init__(self, name):
    self.name = name
    self.runDictionary = {}
    self.refrunDictionary = {}
    self.errorDictionary = {}

  def getQ(self,run):
    if (run.name not in self.runDictionary.keys()):
      q = run.variables[self.name][:]
      q = np.ma.array(q,mask=run.mask)
      self.runDictionary[run.name] = np.ma.compress_cols(q)
    return self.runDictionary[run.name]

  def getQRef(self,run,refrun):
    if (run.name not in self.refrunDictionary.keys()):
      qRef = refrun.variables[self.name][:]
      qRef = np.ma.array(qRef,mask=run.mask)
      self.refrunDictionary[run.name] = np.ma.compress_cols(qRef)
    return self.refrunDictionary[run.name]

  def computeQError(self,run,refrun=None):
    if (run.name not in self.errorDictionary.keys()):
      if (run.name not in self.runDictionary.keys()):
        q = getQ(run)
      else:
        q = self.runDictionary[run.name]
      if (run.name not in self.refrunDictionary.keys()):
        qRef = getQRef(run,refrun)
      else:
        qRef = self.refrunDictionary[run.name]
      self.errorDictionary[run.name] = np.abs(q-qRef)
    return self.errorDictionary[run.name]

class Limiter():

  def __init__(self, name, magnitudeName, qName, limiterType='other', description=None):
    self.name = name
    self.magnitudeName = magnitudeName
    self.q = Quantity(qName)
    self.limiterType = limiterType
    self.description = description
    self.limiterErrorDictionary = {}
    self.cumulativeLimiterErrorDictionary = {}
    self.qErrorDictionary = {}

  def getLimiterError(self,run):
    if (run.name not in self.limiterErrorDictionary.keys()):
      limiterError = run.variables[self.magnitudeName]
      self.limiterErrorDictionary[run.name] = np.ma.compress_cols(limiterError)
    return self.limiterErrorDictionary[run.name]

  def computeCumulativeLimiterError(self,run):
    if (run.name not in self.cumulativeLimiterErrorDictionary.keys()):
      limiterError = self.getLimiterError(run)
      cumulativeError = np.zeros(limiterError.shape)
      for j in range(1,limiterError.shape[1]):
        cumulativeError[:,j] = cumulativeError[:,j-1] + limiterError[:,j]
      self.cumulativeLimiterErrorDictionary[run.name] = cumulativeError
    return self.cumulativeLimiterErrorDictionary[run.name]

  def getQ(self,run): # for backwards compatibility
    return self.q.getQ(run)

  def getQRef(self,run,refrun): # for backwards compatibility
    return self.q.getQRef(run, refrun)

  def computeQError(self,run,refrun): # for backwards compatibility
    return self.q.computeQError(run,refrun)

class Run():

    def __init__(self,runName,dt,limiterList=None,quantityList=None):
      self.name = runName
      self.dt = dt
      self.variables = {}
      print('Loading ' + runName)
      dataset = Dataset(runName,mode='r')
      if (quantityList is None):
        for limiter in limiterList:
          if (limiter.q.name not in self.variables.keys()):
            self.variables[limiter.q.name] = dataset.variables[limiter.q.name][:]
          self.variables[limiter.magnitudeName] = dataset.variables[limiter.magnitudeName][:]
        self.mask = self.variables[limiterList[0].magnitudeName].mask
      elif (limiterList is None):
        for quantity in quantityList:
          if (quantity.name not in self.variables.keys()):
            self.variables[quantity.name] = dataset.variables[quantity.name][:]
        # TODO: Make this more dynamic
        self.mask = dataset.variables['qc_conserv'][:].mask

      dataset.close()

################################################################################

def getDefaultLimiterList():

  limiterList = []

  # Rescale-type limiters
  lim = Limiter('qc_conserv','qc_conserv_mag','cloud_mass','rescale')
  limiterList.append(lim)
  lim = Limiter('nc_conserv','nc_conserv_mag','cloud_number','rescale')
  limiterList.append(lim)
  lim = Limiter('ice_number','ni_conserv','ni_conserv_mag','rescale')
  limiterList.append(lim)
  lim = Limiter('qr_conserv','qr_conserv_mag','rain_mass','rescale')
  limiterList.append(lim)
  lim = Limiter('nr_conserv','nr_conserv_mag','rain_number','rescale')
  limiterList.append(lim)
  lim = Limiter('snow_number','ns_conserv','ns_conserv_mag','rescale')
  limiterList.append(lim)
  lim = Limiter('ice_mass','qi_conserv','qi_conserv_mag','rescale')
  limiterList.append(lim)
  lim = Limiter('snow_mass','qs_conserv','qs_conserv_mag','rescale')
  limiterList.append(lim)
  lim = Limiter('ice_number','ni_tend_lim','ni_tend_lim_mag','rescale')
  limiterList.append(lim)

  return limiterList

################################################################################

def getDefaultRunList(limiterList):

  runList = []

  run = Run('warm1_mg2_acme_v1beta_dt1.0_mstep5.nc',5,limiterList)
  runList.append(run)
  run = Run('warm1_mg2_acme_v1beta_dt1.0_mstep10.nc',10,limiterList)
  runList.append(run)
  run = Run('warm1_mg2_acme_v1beta_dt1.0_mstep15.nc',15,limiterList)
  runList.append(run)
  run = Run('warm1_mg2_acme_v1beta_dt1.0_mstep30.nc',30,limiterList)
  runList.append(run)
  run = Run('warm1_mg2_acme_v1beta_dt1.0_mstep60.nc',60,limiterList)
  runList.append(run)
  run = Run('warm1_mg2_acme_v1beta_dt1.0_mstep120.nc',120,limiterList)
  runList.append(run)
  run = Run('warm1_mg2_acme_v1beta_dt1.0_mstep300.nc',300,limiterList)
  runList.append(run)
  run = Run('warm1_mg2_acme_v1beta_dt1.0_mstep600.nc',600,limiterList)
  runList.append(run)
  run = Run('warm1_mg2_acme_v1beta_dt1.0_mstep900.nc',900,limiterList)
  runList.append(run)
  run = Run('warm1_mg2_acme_v1beta_dt1.0_mstep1200.nc',1200,limiterList)
  runList.append(run)

  return runList