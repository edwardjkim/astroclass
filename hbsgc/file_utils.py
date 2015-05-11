'''
Simple file IO functions for HBSGC.
'''

from __future__ import (print_function, division, unicode_literals)
import os
import numpy as np
from scipy import interpolate
import pyfits

__all__ = ['get_num_files', 'get_file_length', 'read_file', 'regrid_filter',
           'calc_normalization', 'count_data']


def get_num_files(file_location):
    '''
    Takes a text file with a list of paths to SED files.
    Returns the number of SEDs listed in the file.
    '''
    # initialize i for empty files
    if not os.path.isabs(file_location):
        file_location = os.path.join(os.getcwd(), file_location)

    with open(file_location) as f:
        for i, line in enumerate(f):
            pass

    return i + 1


def get_file_length(file_num, file_locations):
    '''
    Takes a text file with a list of paths to SED files,
    and returns the number of lines in the n-th file in the list.
    '''
    n = 0
    with open(file_locations) as list_file:
        for i, list_line in enumerate(list_file):
            if i == file_num:
                if os.path.isabs(list_line):
                    path = os.path.abspath(list_line.strip())
                else:
                    path = os.path.join(os.path.dirname(
                        os.path.abspath(file_locations)), list_line.strip())

                with open(path) as sed_file:
                    for line in sed_file:
                        if line.strip():
                            n += 1
    return n


def read_file(file_num, file_locations):
    '''
    Read the wavelengths and the values of the n-th file listed
    in the SED file.
    Returns a pair of (wavelengths, values) numpy arrays.
    '''
    with open(file_locations) as list_file:
        for i, list_line in enumerate(list_file):
            if i == file_num:
                if os.path.isabs(list_line):
                    path = os.path.abspath(list_line.strip())
                else:
                    path = os.path.join(os.path.dirname(
                        os.path.abspath(file_locations)), list_line.strip())

                plam, pval = np.loadtxt(path, unpack=True)

    return plam, pval


def regrid_filter(plam, pval, slength, length):
    '''
    Regrid input filter profile to a common, finer wavelength grid that
    both SEDs and filters will share.
    Return the rescaled (wavelengths, values) arrays.
    '''
    x = plam[: slength]
    y = pval[: slength]

    start = x[0]
    stop = x[-1]
    step = (x[-1] - x[0]) / (length - 1)

    fine_p_lam = np.arange(start, stop, step)

    tck = interpolate.splrep(x, y)

    fine_p_val = interpolate.splev(fine_p_lam, tck)

    # force throughput > 0
    fine_p_val[fine_p_val < 0.0] = 0.0

    return fine_p_lam, fine_p_val


def calc_normalization(plam, pval, length):
    '''
    Calculate the specified fileter's normalization (zero point)
    in flux density (F_lambda), for AB photometric system.
    '''
    abfnu = 3.631e-20
    cinang = 3.0e18

    h = plam[1] - plam[0]

    w = np.ones(length - 1)
    w[0] = 3.0 / 8.0
    w[length - 2] = 3.0 / 8.0
    w[1] = 7.0 / 6.0
    w[length - 3] = 7.0 / 6.0
    w[-2] = 23.0 / 24.0
    w[length - 4] = 23.0 / 24.0

    norm_val = (w * h / plam * pval * abfnu * cinang).sum()

    return norm_val


def count_data(data_file):
        '''
        Count the number of files in the data file.
        Supports FITS format and Numpy .npy files.
        '''
        if data_file.lower().endswith('.fits'):
            with pyfits.open(data_file) as fits:
                nrow = len(fits[1].data)

        if data_file.lower().endswith('.npy'):
            npy = np.load(data_file)
            nrow = len(npy)

        return nrow
