#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import gamma

# Constants
from mg2_constants import rh2o, cpair, latvap, rair, tmelt, tboil, epsilo

br = 0.8
gamma_br_const = gamma(2.5 + 0.5*br)
f1r = 0.78
f2r = 0.308
rv = rh2o
cpp = cpair
xxlv = latvap
rhow = 1000.
shape_coef = rhow * np.pi
r = rair
ar = 841.99667
rhosu = 85000. / (r * tmelt)

def goff_gratch_svp(t):
    return 100 * 10.**(
        -7.90298*(tboil/t-1.) +
        5.02808*np.log10(tboil/t) -
        1.3816e-7*(10.**(11.344*(1.-t/tboil)) - 1.) +
        8.1328e-3*(10.**(-3.49149*(tboil/t-1.)) - 1.) +
        np.log10(1013.246))

def calc_qvl(t, p):
    esl = goff_gratch_svp(t)
    assert p > esl, "Saturation pressure greater than actual pressure!"
    return epsilo * esl / (p - (1.-epsilo)*esl)

# Input quantities
# Note that the "0." below stands for a cloud fraction of 0, which is allowed
# in this particular case due to a quirk of rain evaporation.
t = 278.146057129
p = 78384.5625
qclr = (0.00625352561474 - 0.*calc_qvl(t, p)) / (1. - 0.)
qric = 1.18794696391e-6 / 0.0001
dr = (qric / ((3398.73754883 / 0.0001) * shape_coef)) ** (1./3.)

def evaporation_tendency(t, p, qclr, qric, dr):

    # Derived quantities
    lamr = 1./dr
    nric = lamr**3 * qric / shape_coef
    n0r = nric * lamr
    rho = p / (r*t)
    rhof = (rhosu/rho)**0.54
    arn = ar * rhof
    dv = 8.794e-5 * t**1.81 / p
    mu = 1.496e-6 * t**1.5 / (t + 120.)
    sc = mu / (rho*dv)
    qvl = calc_qvl(t, p)

    ab = 1. + xxlv*xxlv*qvl / (rv*cpp*t*t)
    eps = 2.*np.pi*n0r*rho*dv* \
          (f1r/(lamr*lamr) +
           f2r*np.sqrt(arn*rho/mu)*
           sc**(1./3.)*gamma_br_const/
           (lamr**(2.5+0.5*br)))

    return eps*(qclr-qvl)/ab

npoints = 101
diameters = np.linspace(20.e-6, 500.e-6, npoints)
tendencies = np.zeros((npoints,))
for i in range(npoints):
    tendencies[i] = evaporation_tendency(t, p, qclr, qric, diameters[i])

plt.plot(diameters*1.e6, tendencies*1.e3)
plt.xlabel('Diameter (microns)')
plt.ylabel('Tendency (g/kg/s)')
plt.axis('tight')
plt.savefig('evap_rate.eps')
plt.close()
