from netCDF4 import Dataset
import numpy as np
import pylab as pyplot

# Specify limiter names as they appear in the NetCDF (.nc)file
limitercount_list = ['qric_lim','qric_qsmall','nric_qsmall','nric_neg', \
    'qsic_lim', 'qsic_qsmall', 'nsic_qsmall', 'nsic_neg', \
    'qc_conserv', 'ice_nuc_lim', 'ice_dep_lim', \
    'nc_conserv', 'qr_conserv', 'nr_conserv', 'qi_conserv', \
    'ni_conserv', 'qs_conserv', 'ni_tend_lim']

limitermag_list = ['rain_evap_lim_mag', 'snow_sub_lim_mag', 'ice_sub_lim_mag']

#limitermag_list = ['rain_evap_lim_mag', 'snow_sub_lim_mag', 'ice_sub_lim_mag']
limitermag_list = ['qr_conserv_mag', 'nr_conserv_mag']
q_list = ['rain_mass','rain_number']

# Specify the file names and corresponding dt values for each of the runs
run_dict = {'warm1_mg2_acme_v1beta_dt1.0_mstep1.nc':1.0, \
    'warm1_mg2_acme_v1beta_dt1.0_mstep5.nc':5.0, \
    'warm1_mg2_acme_v1beta_dt1.0_mstep10.nc':10.0, \
    'warm1_mg2_acme_v1beta_dt1.0_mstep15.nc':15.0, \
    'warm1_mg2_acme_v1beta_dt1.0_mstep30.nc':30.0, \
    'warm1_mg2_acme_v1beta_dt1.0_mstep60.nc':60.0, \
    'warm1_mg2_acme_v1beta_dt1.0_mstep120.nc':120.0, \
    'warm1_mg2_acme_v1beta_dt1.0_mstep300.nc':300.0, \
    'warm1_mg2_acme_v1beta_dt1.0_mstep600.nc':600.0, \
    'warm1_mg2_acme_v1beta_dt1.0_mstep900.nc':900.0, \
    'warm1_mg2_acme_v1beta_dt1.0_mstep1200.nc':1200.0}

# Reference run
reference_run = 'warm1_mg2_acme_v1beta_dt1.0_mstep1.nc'

################################################################################
# Plot the ratio of limited timsteps to total timesteps
################################################################################

# Initialize variables to hold ratios and grab dt values
num_runs = len(run_dict.keys())
num_limiters = len(limitercount_list)
spatial_ratio = np.zeros((num_limiters,num_runs))
global_ratio = np.zeros((num_limiters,num_runs))
dt_list = sorted(run_dict.values())

# Compute ratios
for j,limiter in enumerate(limitercount_list):
    for k,run in enumerate(sorted(run_dict,key=run_dict.get)):
        data = Dataset(run,mode='r')
        var = data.variables[limiter][:]
        num_levels = np.shape(var)[0];
        num_timesteps = np.shape(var)[1];

        # Spatial ratio is sum_n (sum_n var_i^n) / (# timesteps * # levels)
        spatial_ratio[j,k] = np.sum(var)/(num_levels*num_timesteps)

        # Global ratio is sum_n ( max_i y_i^n ) / (# timesteps)
        tmp = np.max(var,0)
        global_ratio[j,k] = np.sum(tmp)/num_timesteps

    # If limiters was ever active, add plot to corresponding figure
    if (np.sum(spatial_ratio[j,:]) > 0):
        if (np.max(spatial_ratio[j,:]) > 0.1):
            F = 0
        else:
            F = 1
        if (np.max(global_ratio[j,:]) > 0.1):
            G = 2
        else:
            G = 3

        pyplot.figure(F)
        pyplot.plot(dt_list,spatial_ratio[j,:],'-o',label=limiter)
        pyplot.figure(G)
        pyplot.plot(dt_list,global_ratio[j,:],'-o',label=limiter)

# Create plots and labels
for j in range(4):
    pyplot.figure(j)
    pyplot.xlabel('dt',fontsize=15)
    pyplot.ylabel('limited timesteps / timesteps',fontsize=15)
    pyplot.legend()

pyplot.show()


################################################################################
# Plot the ratio of limiter error to total error
################################################################################

# Initialize variables to hold ratios
num_limiters = len(limitermag_list)
ratio = np.zeros((num_limiters,num_runs))

# Compute ratios
for j,limiter in enumerate(limitermag_list):
    # Get reference solution
    data = Dataset(reference_run, mode='r')
    qref = data.variables[q_list[j]][:,-1]

    for k,run in enumerate(sorted(run_dict,key=run_dict.get)):
        data = Dataset(run,mode='r')
        var = data.variables[limiter][:]
        q = data.variables[q_list[j]][:,-1]

        ratio[j,k] = np.max(np.abs(np.sum(var,1)/(q-qref)))

    pyplot.figure(10)
    pyplot.plot(dt_list,ratio[j,:],'-o',label=limiter)

# Create plots and labels
pyplot.figure(10)
pyplot.xlabel('dt',fontsize=15)
pyplot.ylabel('limited timesteps / timesteps',fontsize=15)
pyplot.legend()

pyplot.show()
