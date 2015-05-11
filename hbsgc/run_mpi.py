#!/usr/bin/env python

'''
Author: Edward J. Kim
Email: jkim575@illinois.edu

Hierarchical Bayesian Star/Galaxy Classifier using MPI.
Uses emcee.
'''

from __future__ import print_function, division, unicode_literals
import hbsgc
from emcee.utils import MPIPool
import sys
import time


def main():
    '''
    A parallel run.
    '''
    pool = MPIPool(loadbalance=True)

    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    clf = hbsgc.HBSGC(pool=pool)

    # save start time
    clf.last_clock = time.clock()

    clf.filter_calcs()

    clf.data_calcs()

    clf.star_model_calcs()

    # if clf.calc_model_mags:
    #     clf.star_model_mags()

    clf.gal_model_calcs()

    # if clf.calc_model_mags:
    #     clf.gal_model_mags()

    clf.fit_calcs()

    clf.count_tot = 0

    clf.sample()

    clf.save_proba()

    if clf.min_chi2_write:
        clf.save_min_chi2()

    pool.close()

if __name__ == '__main__':

    sys.exit(main())
