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
    from netCDF4 import Dataset # netcdf reader

    parser = argparse.ArgumentParser(description = 'Print info about NetCDF file')
    
    # required arguments
    parser.add_argument('datafile', type=str,
                        help='NetCDF data file')

    # parse command line args
    args = parser.parse_args()
   
    print "Reading:",args.datafile

    data = Dataset(args.datafile, mode='r')
    
    print
    print '----------------------------------------------------------------------'
    print ' Dimensions'
    print '----------------------------------------------------------------------'
    for keys,values in data.dimensions.items():
        print keys
        print values

    print '----------------------------------------------------------------------'
    print ' Variables '
    print '----------------------------------------------------------------------'
    for keys,values in data.variables.items():
        print keys
        print values

    for v in data.variables.keys():
        print v

# ===============================================================================

if __name__ == "__main__":
    import sys
    sys.exit(main())

# EOF



