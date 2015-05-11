'''
Author: Edward J. Kim
Email: jkim575@illinois.edu

Modified from Fadely's Hierarchical Bayesian star/galaxy classifier at
https://github.com/rossfadely/star-galaxy-classification

Also see Fadely et al. ApJ, 760, 15 (2012).
'''
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
from scipy import interpolate

__all__ = ['regrid_sed', 'integrate_sed']


def regrid_sed(z, p_lam, p_val, fine_length, sed_length, fine_p_lam):
    '''
    Regrid input SED profile to a common, finer wavelength grid that
    both SEDs and filters will share.
    '''
    x = p_lam[: sed_length] * (1.0 + z)
    y = p_val[: sed_length]

    tck = interpolate.splrep(x, y)

    fine_p_val = np.zeros(fine_length)

    mask = (fine_p_lam[: fine_length] / (1.0 + z) < p_lam[0])
    fine_p_val[mask] = 0.0
    fine_p_val[~mask] = interpolate.splev(fine_p_lam[~mask], tck)

    # force throughput > 0
    fine_p_val[fine_p_val < 0.0] = 0.0

    return fine_p_val


def integrate_sed(length, sed_val, fil_lam, fil_val):
    '''
    Calculate the model flux for given filter.
    '''
    h = fil_lam[1] - fil_lam[0]

    w = np.ones(length - 1)
    w[0] = 3.0 / 8.0
    w[length - 2] = 3.0 / 8.0
    w[1] = 7.0 / 6.0
    w[length - 3] = 7.0 / 6.0
    w[2] = 23.0 / 24.0
    w[length - 4] = 23.0 / 24.0

    stop = length - 1
    mod_val = (w * h * fil_lam[:stop] * fil_val[:stop] * sed_val[:stop]).sum()

    return mod_val
