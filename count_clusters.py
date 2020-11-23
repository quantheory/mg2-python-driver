#!/usr/bin/env python

import numpy as np
import scipy.linalg as la
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

CLUSTER_FILE_NAME = "/home/santos/Data/MG2_data_collection.10_cluster_labels.0001-01-06-00000.nc"

cfile = nc4.Dataset(CLUSTER_FILE_NAME, 'r')

ncluster = len(cfile.dimensions['ncluster'])

labels = cfile.variables["label"][:].flatten()

for i in range(ncluster):
    counter = 0
    for label in labels:
        if label == i:
            counter += 1
    print(counter, " grid points are in cluster ", i)
