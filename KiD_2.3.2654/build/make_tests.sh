#!/usr/bin/env python
# ==============================================================================
# This python script will create a set of namelist input files for KiD and save 
# them in the directory specified by rundir.
#
# The naming convention is: case_mphys_dt#.nml
# For exmaple: Using the warm1 case with mg2 microphsysics and a 1s time step 
# the name list is warm1_mg2_dt1.nml
#
# for dt<1, (given the assumption of powers of 2), powers of 2 are used instead
# of time in seconds:
#     for dt = 0.25 -> warm1_mg2_dt-2.nml
#     for dt = 4.0  -> warm1_mg2_dt4.nml
#
# D.J. Gardner @ LLNL
# Dec 2016
# ==============================================================================

import subprocess  # launch subprocesses
import numpy as np # math functions
import os          # operating system functions
import sys
from shutil import copyfile

wrkdir = '/p/lscratchd/dgardner/KiD'

# test case: warm1, warm7, mixed1
casename = 'warm1'

mphys='mg2_acme_v1beta'

# dynamics time step sizes
dtvals = [ 1.0 ]

# dtm = dt * mstep, physics time step sizes
mstepvals = [ 1, 5, 10, 15, 30, 60, 300, 600, 900, 1200]

# ------------------------------------------------------------------------------
# make KiD
command = "make_kid.sh 8 " + wrkdir + '/' + casename
ierr = subprocess.call(command, shell=True)

if (ierr != 0):
       print "ERROR: Build Failed"
       sys.exit()

# ------------------------------------------------------------------------------
# make name lists

for dt in dtvals:
       for mstep in mstepvals:
              
              infile  = open(casename+'.nml','r')  
              
              nmlname = casename+'_'+mphys+'_dt'+str(dt)+'_mstep'+str(mstep)+'.nml'
              outfile = open(wrkdir+'/'+casename+'/'+nmlname,'w')
              
              for line in infile:
                     if 'MPHYS' in line:
                            line = line.replace('MPHYS', '\''+str(mphys)+'\'')
                     if 'DT' in line:
                            line = line.replace('DT', str(dt))
                     if 'MSTEP' in line:
                            line = line.replace('MSTEP', str(mstep))
                     outfile.write(line)
              outfile.close()
              infile.close()
