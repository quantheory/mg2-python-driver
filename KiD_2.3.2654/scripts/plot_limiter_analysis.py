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

# Reference run information
reference_run = 'warm1_mg2_acme_v1beta_dt1.0_mstep1.nc'
T = 3600
num_levels = 120

# Figures to save (empty dictionary to show all plots)
save_dict = {}
#save_dict = {'number_ratio_global1':2, \
#    'number_ratio_global2':3, \
#    'number_ratio_spatial1':0, \
#    'number_ratio_spatial2':1, \
#    'error_ratio':10, \
#    'error_normalized1':30, \
#    'error_normalized2':31}

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
        num_timesteps = T/dt_list[k]

        # Spatial ratio is sum_n (sum_n var_i^n) / (# timesteps * # levels)
        spatial_ratio[j,k] = np.sum(var)/(num_levels*num_timesteps)

        # Global ratio is sum_n ( max_i y_i^n ) / (# timesteps)
        tmp = np.max(var,axis=0)
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
        pyplot.semilogx(dt_list,spatial_ratio[j,:],'-o',label=limiter)
        pyplot.figure(G)
        pyplot.semilogx(dt_list,global_ratio[j,:],'-o',label=limiter)

# Create plots and labels
for j in range(2):
    pyplot.figure(j)
    pyplot.xlabel('dt',fontsize=15)
    pyplot.ylabel('Sum_z Sum_t limiter_bool(z,t) / (levels*timesteps)',fontsize=15)
    pyplot.legend()
for j in range(2,4):
    pyplot.figure(j)
    pyplot.xlabel('dt',fontsize=15)
    pyplot.ylabel('Sum_t Max_z limiter_bool(z,t) / (timesteps)',fontsize=15)
    pyplot.legend()

################################################################################
# Plot the ratio of limiter error to total error
################################################################################

# Initialize variables to hold ratios
num_limiters = len(limitermag_list)
ratio = np.zeros((num_limiters,num_runs,num_levels))

# Compute ratios
for j,limiter in enumerate(limitermag_list):
    # Get reference solution
    data = Dataset(reference_run, mode='r')
    qref = data.variables[q_list[j]][:,-1]

    for k,run in enumerate(sorted(run_dict,key=run_dict.get)):
        if (run is not reference_run):
            data = Dataset(run,mode='r')
            var = data.variables[limiter][:]
            q = data.variables[q_list[j]][:,-1]

            ratio[j,k,:] = np.sum(var,axis=1)/(q-qref)

    pyplot.figure(10)
    maxratio = np.max(np.abs(ratio),axis=2)
    pyplot.semilogx(dt_list,maxratio[j,:],'-o',label=limiter)

# Create plots and labels
pyplot.figure(10)
pyplot.xlabel('dt',fontsize=15)
pyplot.ylabel('dt Max_z |Sum_t (a-1) f2(y,t) / (y-yref)|',fontsize=15)
pyplot.legend()

################################################################################
# Plot the ratio of per-step limiter error to per-step change
################################################################################

# Initialize variables to hold ratios
num_limiters = len(limitermag_list)

# Compute ratios
for j,limiter in enumerate(limitermag_list):
    for k,run in enumerate(sorted(run_dict,key=run_dict.get)):
        data = Dataset(run,mode='r')

        # get limitermag values, save mask, and compress
        tmp = data.variables[limiter][:]
        mask = np.ma.getmask(tmp)
        var = np.ma.compress_cols(tmp)

        # get y values, apply mask, and compress
        tmp = data.variables[q_list[j]][:]
        tmp = np.ma.array(tmp,mask=mask)
        q = np.ma.compress_cols(tmp)

        N = np.shape(q)[1]
        tmp = var[:,0:(N-1)] / (q[:,1:N] - q[:,0:(N-1)] + 1e-16)
        ratio = np.max(np.abs(tmp),axis=0)

        pyplot.figure(10*(j+2) + (k < 5))
        label = 'dt = %f' % dt_list[k]
        pyplot.plot(dt_list[k]*np.arange(N-1),ratio,'-o',label=label)

    for k in range(len(run_dict.keys())/5):
        pyplot.figure(10*(j+2) + k)
        pyplot.title(limiter)
        pyplot.xlabel('t',fontsize=15)
        pyplot.ylabel('dt (a_n-1) f2_n / (y_np1-y_n)',fontsize=15)
        pyplot.legend()

################################################################################
# Save the appropriate plots or show all plots
################################################################################
if (len(save_dict.keys()) > 0):
    for name in save_dict.keys():
        pyplot.figure(save_dict[name])
        pyplot.savefig(name + '.png')
else:
    pyplot.show()
