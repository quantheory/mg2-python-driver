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

import numpy as np # math functions
import os          # operating system functions

wrkdir   = '/p/lscratchd/dgardner/KiD'
casename = 'warm1_test'

icase   = 101     # KiD case
mphys   = 'mg2'   # Microphysics option
dg_dt   = 1.0     # diagnostic output time
final_t = 14400.0 # final simulation time (tctrl(1) in KiD)

# time step sizes
dt_vals = [0.25, 0.5, 1.0, 5.0, 10.0,
           15.0, 30.0, 60.0, 90.0, 120.0,
           300.0, 600.0, 900.0, 1200.0, 1800.0]

# iterate over step sizes dt
for dt in dt_vals:

       if (dt < 1):
              dt_str = str(int(np.log2(dt)))
       else:
              dt_str = str(int(dt))
              
       nmlname  = 'KiD_'+mphys+'_dt'+dt_str+'.nml'
       filename = wrkdir+'/'+casename+'/'+nmlname
                            
       f = os.open(filename,os.O_CREAT|os.O_WRONLY)
       
       # Microphysics Namelist
       nml = ('&mphys \n'
              '! hydrometeor names \n'
              'h_names = \'cloud\',  \'rain\',  \'ice\',  \'snow\',  \'graupel\' \n'
              '! number of moments for each species \n'
              'num_h_moments = 2,2,2,2,0 \n'
              'num_h_bins    = 1,1,1,1,1 \n'
              '! Background values for each moment (assumed the same for all species) \n'
              'mom_init = 0,0,0 \n'
              '! Aerosol initialization \n'
              'num_aero_moments = 0,0,0 \n'
              'num_aero_bins    = 1 \n'
              'aero_N_init      = 0.0d0, 50.0d6, 0.0d0 \n'
              'aero_sig_init    = 0.0d0, 1.4d0, 0.0d0 \n'
              'aero_rd_init     = 0.0d0, 0.05d-6, 0.0d0 \n'
              '/ \n\n')
       
       # Case Namelist    
       nml = nml + ('&case \n'
                    'icase = '+str(icase)+' \n'
                    '/ \n\n')
       
       # Control Namelist 
       nml = nml + ('&control \n'
                    'mphys_scheme = \''+mphys+'\'\n'
                    'dt      = '+str(dt)+'d0 \n'
                    'dgstart = 0.0d0 \n'
                    'dg_dt   = '+ str(dg_dt)+'d0 \n'
                    '! scaling factor for updraft \n'
                    'wctrl(1) = 0.0d0 \n'
                    '! final time \n'
                    'tctrl(1) = '+str(final_t)+'\n'
                    '! half period of forcing \n'
                    'tctrl(2) = 0.0d0 \n'
                    '/ \n\n')
       
       # Switch Namelist 
       nml = nml + ('&switch \n'
                    'l_advect            =.false. \n'
                    'l_diverge           =.false. \n'
                    'l_fix_theta         =.true.  \n'
                    'l_diverge_advection =.false. \n'
                    'l_fix_aerosols      =.true.  \n'
                    'l_periodic_bound    =.true.  \n'
                    '/ \n\n')
       
       nml = nml + ('&addcontrol \n'
                    'iiwarm = .true. \n'
                    '/ \n\n')
       
       os.write(f,nml)
       os.close(f)
