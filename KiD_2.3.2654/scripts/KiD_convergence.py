#!/usr/bin/env python
# -------------------------------------------------------------------------------
# Programmer(s):  David J. Gardner @ LLNL
# -------------------------------------------------------------------------------
# Script to plot convergence in KiD tests
# -------------------------------------------------------------------------------
# Last updated: 11 May 2017
# -------------------------------------------------------------------------------

def main():

    import argparse
    import sys
    import os
    import subprocess
    import operator
    
    from netCDF4 import Dataset # netcdf reader
    import FileHelper
    
    import numpy as np
    
    import matplotlib        as mpl
    import matplotlib.pyplot as plt
    import matplotlib.cm     as cm   # colormap for graphics
    
    parser = argparse.ArgumentParser(description = 'Plot KiD NetCDF files')
    
    # required arguments
    parser.add_argument('Label', type=str,
                        help='Name for saving figure')

    parser.add_argument('reffile', type=str,
                        help='Reference solution output')
        
    parser.add_argument('--outfiles', dest='outfiles', type=str, nargs='+', 
                        required=True,
                        help='Test solution outputs')

    parser.add_argument('--stepsizes', dest='stepsizes', type=float, nargs='+', 
                        required=True,
                        help='Test solution outputs')
    
    # plotting options
    parser.add_argument('--Warm', dest='Warm', action='store_true', 
                        default=False,
                        help='Signal Warm test case')

    parser.add_argument('--SetColors', type=int, nargs='+', dest='SetColors',
                        help='Set line colors')
    
    parser.add_argument('--show', dest='show', action='store_true', 
                        default=False,
                        help='Display plot on the screen rather than writing to file')
    
    parser.add_argument('--debug', dest='debug', action='store_true', default=False,
                        help='Turn on debugging ouput')

    parser.add_argument('--WriteArrays', dest='WriteArrays', action='store_true', 
                        default=False,
                        help='Flag to write convergence arrays to file')

        
    # parse command line args
    args = parser.parse_args()

    # number of output files and step sizes
    noutfiles  = len(args.outfiles)
    nstepsizes = len(args.stepsizes)

    if (noutfiles != nstepsizes):
        print "ERROR: number of steps does not equal number of flies"
        sys.exit()

    # line colors (qualitative colors from colorbrewer2.org)
    if (args.SetColors):    

        if (len(args.SetColors) != noutfiles):
            msg = "ERROR: len(SetColors) != len(outfiles)"
            raise Exception(msg)

        totalcolors = max(args.SetColors)+1

    else:
        totalcolors = noutfiles

        args.SetColors = range(noutfiles)

    if (totalcolors < 10):
        colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3',
                  '#ff7f00','#a65628','#f781bf','#999999','#000000']
    elif (totalcolors < 14):
        colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c',
                  '#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',
                  '#6a3d9a','#ffff99','#b15928','#000000']
    else:
        colors = cm.Set1(np.linspace(0, 1, totalcolors))

    Cvalue = []
    for i in args.SetColors:
        Cvalue.append(colors[i])

    # plot results for final time
    # ---------------------------------------------------------------------------

    # reference data file
    data = Dataset(args.reffile, mode='r')

    # solutions at final time
    qv_ref = data.variables['vapour'][:,-1]
    qc_ref = data.variables['cloud_mass'][:,-1]
    qr_ref = data.variables['rain_mass'][:,-1]
    qi_ref = data.variables['ice_mass'][:,-1]
    qs_ref = data.variables['snow_mass'][:,-1]

    nc_ref = data.variables['cloud_number'][:,-1]
    nr_ref = data.variables['rain_number'][:,-1]
    ni_ref = data.variables['ice_number'][:,-1]
    ns_ref = data.variables['snow_number'][:,-1]
               
    data.close()

    # if (args.Warm):
    #     state_ref = np.concatenate((qv_ref, qc_ref, qr_ref))
    # else:
    #     state_ref = np.concatenate((qv_ref, qc_ref, qr_ref, qi_ref, qs_ref))
        
    # print "Shape: ",np.shape(state_ref)

    err_qv = np.empty([noutfiles]) 

    err_qc = np.empty([noutfiles])  
    err_qr = np.empty([noutfiles])  
    err_qi = np.empty([noutfiles])  
    err_qs = np.empty([noutfiles])  

    err_nc = np.empty([noutfiles])  
    err_nr = np.empty([noutfiles])  
    err_ni = np.empty([noutfiles])  
    err_ns = np.empty([noutfiles])  

    # err_state = np.empty([noutfiles]) 

    for i in range(noutfiles):

        print "Reading:",args.outfiles[i]

        # data file
        data = Dataset(args.outfiles[i], mode='r')   

        qv = data.variables['vapour'][:,-1]

        qc = data.variables['cloud_mass'][:,-1]
        qr = data.variables['rain_mass'][:,-1]
        qi = data.variables['ice_mass'][:,-1]
        qs = data.variables['snow_mass'][:,-1]

        nc = data.variables['cloud_number'][:,-1]
        nr = data.variables['rain_number'][:,-1]
        ni = data.variables['ice_number'][:,-1]
        ns = data.variables['snow_number'][:,-1]

        data.close()

        # if (args.Warm):
        #     state = np.concatenate((qv, qc, qr))
        # else:
        #     state = np.concatenate((qv, qc, qr, qi, qs))

        err_qv[i] = np.sqrt(np.average(np.square(qv-qv_ref)))

        err_qc[i] = np.sqrt(np.average(np.square(qc-qc_ref)))
        err_qr[i] = np.sqrt(np.average(np.square(qr-qr_ref)))

        err_nc[i] = np.sqrt(np.average(np.square(nc-nc_ref)))
        err_nr[i] = np.sqrt(np.average(np.square(nr-nr_ref)))

        if (not args.Warm):
            err_qi[i] = np.sqrt(np.average(np.square(qi-qi_ref)))
            err_qs[i] = np.sqrt(np.average(np.square(qs-qs_ref)))

            err_ni[i] = np.sqrt(np.average(np.square(ni-ni_ref)))
            err_ns[i] = np.sqrt(np.average(np.square(ns-ns_ref)))

        # err_state[i] = np.sqrt(np.average(np.square(state-state_ref)))

    # plot convergence
    # ---------------------------------------------------------------------------
    stepsizes = np.array(args.stepsizes)

    plt.figure(1)
    plt.loglog(stepsizes, err_qv, label='qv', color='blue')
    plt.loglog(stepsizes, err_qc, label='qc', color='green')
    plt.loglog(stepsizes, err_qr, label='qr', color='red')

    if (not args.Warm):
        plt.loglog(stepsizes, err_qi, label='qi', color='cyan')
        plt.loglog(stepsizes, err_qs, label='qs', color='magenta')

    plt.title(args.Label.replace("_"," "))
    plt.xlabel('Physics Step Size')
    plt.ylabel('RMS Error Final Time')
    plt.grid()
    plt.legend(loc='best')
    plt.tight_layout()

    if (args.show):
        plt.show()
    else:
        figname = args.Label+'_convergence_q_finaltime.pdf'
        plt.savefig(figname)
        plt.close()

    # plot convergence
    # ---------------------------------------------------------------------------
    plt.figure(1)
    plt.loglog(args.stepsizes, err_nc, label='nc', color='green')
    plt.loglog(args.stepsizes, err_nr, label='nr', color='red')

    if (not args.Warm):
        plt.loglog(args.stepsizes, err_ni, label='ni', color='cyan')
        plt.loglog(args.stepsizes, err_ns, label='ns', color='magenta')

    plt.title(args.Label.replace("_"," "))
    plt.xlabel('Physics Step Size')
    plt.ylabel('RMS Error Final Time')
    plt.grid()
    plt.legend(loc='best')
    plt.tight_layout()

    if (args.show):
        plt.show()
    else:
        figname = args.Label+'_convergence_n_finaltime.pdf'
        plt.savefig(figname)
        plt.close()

    # estimate order from last three values
    m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_qv[0:3]), 1)
    print "qv: ",m

    m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_qc[0:3]), 1)
    print "qc: ",m

    m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_qr[0:3]), 1)
    print "qr: ",m

    if (not args.Warm):
        m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_qi[0:3]), 1)
        print "qi: ",m
        
        m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_qs[0:3]), 1)   
        print "qs: ",m
        
    # m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_state[0:3]), 1)
    # print "state: ",m

    if (args.WriteArrays):
        with open(args.Label+'_errors_finaltime.txt','w') as fout:
            write_convergence(stepsizes, err_qv, 'qv', fout)
            write_convergence(stepsizes, err_qc, 'qc', fout)
            write_convergence(stepsizes, err_qr, 'qr', fout)
            if (not args.Warm):
                write_convergence(stepsizes, err_qi, 'qi', fout)
                write_convergence(stepsizes, err_qs, 'qs', fout)

            write_convergence(stepsizes, err_nc, 'nc', fout)
            write_convergence(stepsizes, err_nr, 'nr', fout)
            if (not args.Warm):
                write_convergence(stepsizes, err_ni, 'ni', fout)
                write_convergence(stepsizes, err_ns, 'ns', fout)

            # write_convergence(stepsizes, err_state, 'state', fout)

    # all time
    # ---------------------------------------------------------------------------

    # reference data file
    data = Dataset(args.reffile, mode='r')

    # solutions at final time
    qv_ref = np.ravel(data.variables['vapour'][...])

    qc_ref = np.ravel(data.variables['cloud_mass'][...])
    qr_ref = np.ravel(data.variables['rain_mass'][...])
    qi_ref = np.ravel(data.variables['ice_mass'][...])
    qs_ref = np.ravel(data.variables['snow_mass'][...])

    nc_ref = np.ravel(data.variables['cloud_number'][...])
    nr_ref = np.ravel(data.variables['rain_number'][...])
    ni_ref = np.ravel(data.variables['ice_number'][...])
    ns_ref = np.ravel(data.variables['snow_number'][...])
               
    data.close()

    # if (args.Warm):
    #     state_ref = np.concatenate((qv_ref, qc_ref, qr_ref))
    # else:
    #     state_ref = np.concatenate((qv_ref, qc_ref, qr_ref, qi_ref, qs_ref))
        
    # print "Shape: ",np.shape(state_ref)

    err_qv = np.empty([noutfiles]) 

    err_qc = np.empty([noutfiles])  
    err_qr = np.empty([noutfiles])  
    err_qi = np.empty([noutfiles])  
    err_qs = np.empty([noutfiles])  

    err_nc = np.empty([noutfiles])  
    err_nr = np.empty([noutfiles])  
    err_ni = np.empty([noutfiles])  
    err_ns = np.empty([noutfiles])  

    # err_state = np.empty([noutfiles]) 

    for i in range(noutfiles):

        print "Reading:",args.outfiles[i]

        # data file
        data = Dataset(args.outfiles[i], mode='r')

        qv = np.ravel(data.variables['vapour'][...])

        qc = np.ravel(data.variables['cloud_mass'][...])
        qr = np.ravel(data.variables['rain_mass'][...])
        qi = np.ravel(data.variables['ice_mass'][...])
        qs = np.ravel(data.variables['snow_mass'][...])

        nc = np.ravel(data.variables['cloud_number'][...])
        nr = np.ravel(data.variables['rain_number'][...])
        ni = np.ravel(data.variables['ice_number'][...])
        ns = np.ravel(data.variables['snow_number'][...])

        data.close()

        # if (args.Warm):
        #     state = np.concatenate((qv, qc, qr))
        # else:
        #     state = np.concatenate((qv, qc, qr, qi, qs))

        err_qv[i] = np.sqrt(np.average(np.square(qv-qv_ref)))

        err_qc[i] = np.sqrt(np.average(np.square(qc-qc_ref)))
        err_qr[i] = np.sqrt(np.average(np.square(qr-qr_ref)))

        err_nc[i] = np.sqrt(np.average(np.square(nc-nc_ref)))
        err_nr[i] = np.sqrt(np.average(np.square(nr-nr_ref)))

        if (not args.Warm):
            err_qi[i] = np.sqrt(np.average(np.square(qi-qi_ref)))
            err_qs[i] = np.sqrt(np.average(np.square(qs-qs_ref)))

            err_ni[i] = np.sqrt(np.average(np.square(ni-ni_ref)))
            err_ns[i] = np.sqrt(np.average(np.square(ns-ns_ref)))

        # err_state[i] = np.sqrt(np.average(np.square(state-state_ref)))

    # plot convergence
    # ---------------------------------------------------------------------------
    stepsizes = np.array(args.stepsizes)

    plt.figure(1)
    plt.loglog(stepsizes, err_qv, label='qv', color='blue')
    plt.loglog(stepsizes, err_qc, label='qc', color='green')
    plt.loglog(stepsizes, err_qr, label='qr', color='red')

    if (not args.Warm):
        plt.loglog(stepsizes, err_qi, label='qi', color='cyan')
        plt.loglog(stepsizes, err_qs, label='qs', color='magenta')

    plt.title(args.Label.replace("_"," "))
    plt.xlabel('Physics Step Size')
    plt.ylabel('RMS Error All Time')
    plt.grid()
    plt.legend(loc='best')
    plt.tight_layout()

    if (args.show):
        plt.show()
    else:
        figname = args.Label+'_convergence_q_alltime.pdf'
        plt.savefig(figname)
        plt.close()

    # plot convergence
    # ---------------------------------------------------------------------------
    plt.figure(1)
    plt.loglog(args.stepsizes, err_nc, label='nc', color='green')
    plt.loglog(args.stepsizes, err_nr, label='nr', color='red')

    if (not args.Warm):
        plt.loglog(args.stepsizes, err_ni, label='ni', color='cyan')
        plt.loglog(args.stepsizes, err_ns, label='ns', color='magenta')

    plt.title(args.Label.replace("_"," "))
    plt.xlabel('Physics Step Size')
    plt.ylabel('RMS Error All Time')
    plt.grid()
    plt.legend(loc='best')
    plt.tight_layout()

    if (args.show):
        plt.show()
    else:
        figname = args.Label+'_convergence_n_alltime.pdf'
        plt.savefig(figname)
        plt.close()

    # estimate order from last three values
    m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_qv[0:3]), 1)
    print "qv: ",m

    m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_qc[0:3]), 1)
    print "qc: ",m

    m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_qr[0:3]), 1)
    print "qr: ",m

    if (not args.Warm):
        m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_qi[0:3]), 1)
        print "qi: ",m
        
        m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_qs[0:3]), 1)   
        print "qs: ",m
        
    # m,b = np.polyfit(np.log10(stepsizes[0:3]), np.log10(err_state[0:3]), 1)
    # print "state: ",m

    if (args.WriteArrays):
        with open(args.Label+'_errors_alltime.txt','w') as fout:
            write_convergence(stepsizes, err_qv, 'qv', fout)
            write_convergence(stepsizes, err_qc, 'qc', fout)
            write_convergence(stepsizes, err_qr, 'qr', fout)
            if (not args.Warm):
                write_convergence(stepsizes, err_qi, 'qi', fout)
                write_convergence(stepsizes, err_qs, 'qs', fout)

            write_convergence(stepsizes, err_nc, 'nc', fout)
            write_convergence(stepsizes, err_nr, 'nr', fout)
            if (not args.Warm):
                write_convergence(stepsizes, err_ni, 'ni', fout)
                write_convergence(stepsizes, err_ns, 'ns', fout)

# ===============================================================================

def write_convergence(dt, err, label, fout):
    
    import numpy as np

    # check lengths
    print >> fout, '# ',10*'-'
    print >> fout, '# ',label
    print >> fout, '# ',10*'-'
    for i in range(len(dt)):
        if (i == len(dt)-1):
            print >> fout, dt[i], '\t', err[i], '\t', 0.0
        else:
            print >> fout, dt[i], '\t', err[i], '\t', np.log10(err[i]/err[i+1])/np.log10(dt[i]/dt[i+1])
               
# ===============================================================================

if __name__ == "__main__":
    import sys
    sys.exit(main())
    
# EOF




    # plt.loglog(args.stepsizes, err_state, label='state')

    # ordpow = 1.0
    # X = stepsizes[0:5]
    # Y = err_state[0] * (X / X[0])**ordpow

    # # horizontal shift
    # if (args.OrderHShift < 0): 
    #     # left shift
    #     X = X/np.log10(10 + abs(args.OrderHShift))
    # elif (args.OrderHShift > 0):
    #     # right shift
    #     X = X*np.log10(10 + args.OrderHShift)
    
    # # vertical shift
    # Y = args.OrderVShift * Y


    # left shift
    # X = X/np.log10(10 + 10)
               
    # vertical shift
    # Y = args.OrderVShift * Y
        
    # plt.loglog(X, Y,
    #            linestyle='--', 
    #            color='black', 
    #            label=None)
