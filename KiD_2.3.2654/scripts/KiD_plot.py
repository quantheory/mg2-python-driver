#!/usr/bin/env python
# -------------------------------------------------------------------------------
# Programmer(s):  David J. Gardner @ LLNL
# -------------------------------------------------------------------------------
# Script to plot integrated results from KiD NetCDF output files 
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
    import matplotlib.cm     as cm   # colormap
    
    parser = argparse.ArgumentParser(description = 'Plot KiD NetCDF files')
    
    # required arguments
    parser.add_argument('label', type=str,
                        help='Name for saving figure')
    
    parser.add_argument('variable', type=str,
                        help='which variable to plot form the NetCDF file')
    
    parser.add_argument('output_dir',type=str,nargs='+',
                        help='Ouptput directory to plot from')

    # plot types
    parser.add_argument('--timeaverage', dest='timeaverage', action='store_true',
                        help='Plot time averaged quantities')
    
    # plotting options
    parser.add_argument('--stepsizes', dest='dtm', type=int, nargs='+',
                        help='Step sizes for plot legend')

    parser.add_argument('--StopPlot', dest='StopPlot', type=int,
                        help='Stop plotting after set number of files')

    parser.add_argument('--SetColors', type=int, nargs='+', dest='SetColors',
                        help='Set line colors')

    parser.add_argument('--LegendOutside', dest='LegendOutside', action='store_true',
                        help='Place legend outside the plot')
    
    parser.add_argument('--show', dest='show', action='store_true',
                        help='Display plot on the screen rather than writing to file')
    
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help='Turn on debugging ouput')
       
    # parse command line args
    args = parser.parse_args()

    # empty list for output files in output directory
    output_files = []

    # plot line styles
    LineStyles = ['-', '--']
    Lstyle = []

    # initialize directory counter for setting plot properties
    count = 0

    # iterate over output directories
    for outdir in args.output_dir:
    
        # check that directory exist 
        if (not os.path.isdir(outdir)):
            print "ERROR: Directory",outdir,"does not exist."
            sys.exit()
    
        print "Searching for netcdf files in",outdir

        # get list of output files in the outdir
        temp_output_files = FileHelper.list_files(outdir,'*.nc')

        if (len(temp_output_files) == 0):
            print "Warning: No output files found in",outdir
        else:
            # sort output files
            temp_output_files = sorted(temp_output_files, 
                                       key=FileHelper.numericalSort)

            # add output files to list for plotting
            for f in temp_output_files:
                output_files.append(f)
                Lstyle.append(LineStyles[count])

        # update directory counter 
        count += 1 

    # number of output files found
    noutfiles = len(output_files)

    if (noutfiles == 0):
        print "ERROR: No output files found"
        sys.exit()
    else:
        print "Number of files to plot",noutfiles
    
    # line color map, qualitative colors from colorbrewer2.org
    if (args.SetColors):    
        totalcolors = max(args.SetColors)+1
    else:
        totalcolors = noutfiles

    if (totalcolors < 10):
        colors = ['#000000','#e41a1c','#377eb8','#4daf4a','#984ea3',
                  '#ff7f00','#a65628','#f781bf','#999999']
    elif (totalcolors < 14):
        colors = ['#000000','#a6cee3','#1f78b4','#b2df8a','#33a02c',
                  '#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',
                  '#6a3d9a','#ffff99','#b15928']
    else:
        colors = cm.Set1(np.linspace(0, 1, totalcolors))

    # counters for plot properties 
    count  = 0
    ncolor = 0

    # plot output files
    for outfile in output_files:

        print "Reading:",outfile

        # name list file
        nmlfile = outfile[:-3] + '.nml'

        # read step sizes from name list
        dt    = 0.0
        mstep = 0
        # dtm   = 0.0
        
        with open(nmlfile) as fn:
            for line in fn:
                text = line.split()

                # dynamics step size
                if ("DT" in text):
                    dt = float(text[2])

                # steps between physics computation
                if ("MSTEP" in text):
                    mstep = int(text[2][:-1])
        
        # effective physics step size
        # dtm = mstep * dt
        
        # read data file
        data = Dataset(outfile, mode='r')
    
        if (args.timeaverage):
            height   = data.variables['z'][:]
            datavals = data.variables[args.variable][...]

            plotdata = np.average(datavals, axis=1)

            if (args.debug):
                print np.size(height), np.size(plotdata)
                print np.shape(height), np.shape(plotdata)
        
        else:
            time     = data.variables['time'][:]
            plotdata = data.variables[args.variable][:]

            # convert seconds to hours
            time = time / 3600

        # convert units from kg/kg to g/kg
        plotdata = 1000 * plotdata
        
        data.close()

        if (args.debug):
            print "Plotting:",outfile

        if (args.timeaverage):       

            plt.figure
            plt.plot(plotdata, height, 
                     color=colors[ncolor], 
                     linestyle=Lstyle[ncolor], 
                     label=str(args.dtm[count]))

        else:
        
            plt.figure
            plt.plot(time, plotdata, 
                     color=colors[ncolor], 
                     linestyle=Lstyle[ncolor], 
                     label=str(args.dtm[count]))

        # update counters
        count  += 1   
        ncolor += 1

        # stop plotting after set number of files
        if (args.StopPlot):
            if (count >= args.StopPlot):
                break
    
    plt.figure(1)

    if (args.timeaverage):       
        plt.xlabel(args.variable.replace("_"," ")+' [g/kg]')
    else:
        plt.xlabel('hours')

    if (args.timeaverage):       
        plt.ylabel('meters')
    else:
        plt.ylabel(args.variable.replace("_"," ")+' [g/kg]')

    plt.grid(True)

    # (horiz shift, vert shift)
    lgd = plt.legend(loc='best', bbox_to_anchor=(1, 1)) 
    plt.title(args.label.replace("_"," "))

    if (args.show):
        # show plot on screen
        plt.show()
    else:
        # save plot to file
        fname=args.label+'_'+args.variable+'.pdf'

        if (args.LegendOutside):
            plt.savefig(fname, bbox_extra_artists=(lgd,), 
                        bbox_inches='tight')
        else:
            plt.savefig(fname)
        
        try:
            # remove extra whitespace from pdf file
            subprocess.call('pdfcrop '+fname, shell=True)
        except:
            print "Warning: unable to crop pdf file"
        else:
            # crop successful remove uncropped file
            os.remove(fname);
        
    plt.close()

# ===============================================================================

if __name__ == "__main__":
    import sys
    sys.exit(main())

# EOF






