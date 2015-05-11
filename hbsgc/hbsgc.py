#! /usr/bin/env python
'''
A Hierarchical Bayesian Star-Galaxy Classifier.

Modified from Ross Fadely's code at

https://github.com/rossfadely/star-galaxy-classification

Also see R. Fadely et al. 2012 ApJ 760 15.
'''

from __future__ import (print_function, division, unicode_literals)
import numpy as np
import pyfits
import os
import sys
import time
import emcee
from functools import partial
from file_utils import (get_num_files, get_file_length, read_file,
                        calc_normalization, regrid_filter, count_data)
from sed_utils import regrid_sed, integrate_sed


class HBSGC(object):
    '''
    A Hierarchical Bayesian Star-Galaxy Classifier object.
    '''
    def __init__(self, pool=None):
        '''
        Constructor.
        '''
        # for MPI, pool is an emcee.utils.MPIPool() object. See emcee doc [1].
        # [1] http://dan.iel.fm/emcee/current/user/advanced
        self.pool = pool
        # Change the following parameters appropriately. See docs directory.
        self.filters_input = 'filters/filters.txt'
        self.gals_sed_input = 'seds/galaxies/galaxies_sed_cfhtls.list'
        self.stars_sed_input = 'seds/stars/stars_sed_all.list'
        self.data_file = 'data/vvds_g_r_0_1.fits'
        self.hyp_out_file = 'results/hyp.npy'
        self.ln_like_out_file = 'results/lnlike.npy'
        self.ln_p_tot_file = 'results/ln_p_tot.txt'
        self.median_proba_out_file = 'results/vvds_g_r_0_1_median.hbc'
        self.mean_proba_out_file = 'results/vvds_g_r_0_1_mean.hbc'
        self.zmin = 0.0
        self.zmax = 2.5
        self.nz = 51
        self.noise_fudge = 0.06
        self.nburn = 50
        self.niter = 100
        # optional parameters.
        self.hyp_in_file = 'results/hyp.npy'
        self.star_mod_mags_file = ''
        self.gal_mod_mags_file = ''
        self.star_chi2_file = ''
        self.gal_chi2_file = ''
        self.flux_unit_factor = 5.0e13
        self.prob_frac = 1.0e-10
        self.p_floor = np.finfo(np.float).tiny
        self.ncstep = 7
        self.nwalkers = 1000
        self.use_hyp_in = False
        self.calc_model_mags = False
        self.min_chi2_write = False

    def filter_calcs(self):
        '''
        Reads in the models seds and filter throughput response curves,
        regrids the filters onto a finer grid for interpolation,
        and calls a routine to calculate the filter flux zero-points.
        Testing indicates details of regridding and interpolation scheme
        is not very important.
        '''
        sed_length = 0

        self.n_filter = get_num_files(self.filters_input)

        if self.n_filter < 2:
            raise ValueError

        print('Found %d Filters' % self.n_filter)

        self._n_star_template = get_num_files(self.stars_sed_input)
        self._n_gal_template = get_num_files(self.gals_sed_input)

        # find the longest SED file among the bunch
        # first the star SEDs
        for i in xrange(self._n_star_template):
            n = get_file_length(i, self.stars_sed_input)
            if n * 2 > sed_length:
                sed_length = n * 2

        # do the same for the galaxies
        for i in xrange(self._n_gal_template):
            n = get_file_length(i, self.gals_sed_input)
            if n * 2 > sed_length:
                sed_length = n * 2

        self.filter_lgth_fine = np.zeros(self.n_filter)
        self.norm = np.zeros(self.n_filter)
        self.filter_lamb_fine = {}
        self.filter_thru_fine = {}

        for i in xrange(self.n_filter):
            n = get_file_length(i, self.filters_input)
            filter_length = n
            regrid_factor = np.round(float(sed_length) / float(n))

            filter_lgth_fine = n * regrid_factor

            filt_lamb, filt_thru = read_file(i, self.filters_input)

            filter_lamb_fine, filter_thru_fine = \
                regrid_filter(filt_lamb, filt_thru,
                              filter_length, filter_lgth_fine)

            norm = calc_normalization(filter_lamb_fine, filter_thru_fine,
                                      filter_lgth_fine)

            print("Filter %ld has (AB) zeropoint flux normalization: %f"
                  % (i, norm))

            self.filter_lgth_fine[i] = filter_lgth_fine
            self.norm[i] = norm
            self.filter_lamb_fine[i] = filter_lamb_fine
            self.filter_thru_fine[i] = filter_thru_fine

    def data_calcs(self):
        '''
        Assign data arrays, call calc_data_vals to calculate flux values.
        '''
        self.n_data = count_data(self.data_file)
        print('\nThere are %ld sources in data file' % self.n_data)

        self.data_flux = np.zeros((self.n_filter, self.n_data))
        self.data_flux_err = np.zeros((self.n_filter, self.n_data))

        for i in xrange(self.n_filter):
            dflux, dflux_err = self.calc_data_vals(i)
            self.data_flux[i] = dflux
            self.data_flux_err[i] = dflux_err

    def calc_data_vals(self, i):
        '''
        Read in the FITS data. FITS file should have magnitudes in the first
        n_filter columns, then the uncertainties, then extra columns
        (not read in).
        Data and uncertainties are converted to flux units, with special
        treatment for missing or undetected data. The latter needs to be
        taylored to the catalog, or vice versa.
        '''
        if self.data_file.endswith('.fits'):

            hdulist = pyfits.open(self.data_file)
            data = hdulist[1].data
            nrows = len(data)

            cols_names = hdulist[1].columns.names
            cols_mag = cols_names[i]
            p_mag = data[cols_mag]
            cols_err = cols_names[i + self.n_filter]
            p_err = data[cols_err]

            hdulist.close()

        if self.data_file.endswith('.npy'):

            npy = np.load(self.data_file)
            nrows = len(npy)

            p_mag = npy[:, i]
            p_err = npy[:, i + self.n_filter]

        flux = np.zeros(nrows)
        flux_err = np.zeros(nrows)

        # missing data
        missing = (p_mag < 0.0)
        flux[missing] = np.power(10.0, -0.4 * 10.0) * self.norm[i]
        flux_err[missing] = np.power(10.0, -0.4 * 10.0) * self.norm[i] * \
            np.log(10.0) * 0.4 * 1000000.0
        # undetected
        undetected = (p_mag > 50.0)
        flux[undetected] = 0.0
        flux_err[undetected] = np.power(10.0, -0.4 * p_err[undetected]) * \
            self.norm[i] * 2.0
        # all ok
        all_ok = np.logical_and(~missing, ~undetected)
        flux[all_ok] = np.power(10.0, -0.4 * p_mag[all_ok]) * self.norm[i]
        flux_err[all_ok] = np.power(10.0, -0.4 * p_mag[all_ok]) * \
            self.norm[i] * np.log(10.0) * 0.4 * p_err[all_ok]
        # add noise model
        flux_err = np.sqrt(np.power(flux_err, 2.0) +
                           np.power(flux * self.noise_fudge, 2.0))

        return flux, flux_err

    def star_model_calcs(self):
        '''
        Read in star SEDs, regrid them, calculate the model flux in each filter
        for each redshift.
        '''
        print('Using %ld star templates' % self._n_star_template)

        model_flux_stars = np.zeros((self.n_filter, self._n_star_template))

        for i in xrange(self._n_star_template):
            sed_length = get_file_length(i, self.stars_sed_input)

            p_lam, p_val = read_file(i, self.stars_sed_input)

            for k in xrange(self.n_filter):
                p_val_fine = regrid_sed(
                    0.0, p_lam, p_val, self.filter_lgth_fine[k],
                    sed_length, self.filter_lamb_fine[k])

                mod_val = integrate_sed(
                    self.filter_lgth_fine[k], p_val_fine,
                    self.filter_lamb_fine[k], self.filter_thru_fine[k])

                model_flux_stars[k][i] = mod_val

        self.model_flux_stars = model_flux_stars

    def gal_model_calcs(self):
        '''
        Read in galaxy SEDs, regrid them, calculate the model flux in each
        filter for each redshift.
        '''
        print('Using %ld gal templates over %ld redshifts'
              % (self._n_gal_template, self.nz))

        self.zstep = (self.zmax - self.zmin) / (self.nz - 1)

        model_flux_gals = np.zeros(
            (self.n_filter, self._n_gal_template * self.nz))

        for i in xrange(self._n_gal_template):

            sed_length = get_file_length(i, self.gals_sed_input)

            p_lam, p_val = read_file(i, self.gals_sed_input)

            for j in xrange(self.nz):
                z = self.zmin + j * self.zstep
                for k in xrange(self.n_filter):
                    p_val_fine = regrid_sed(
                        z, p_lam, p_val, self.filter_lgth_fine[k],
                        sed_length, self.filter_lamb_fine[k])
                    mod_val = integrate_sed(
                        self.filter_lgth_fine[k], p_val_fine,
                        self.filter_lamb_fine[k], self.filter_thru_fine[k])

                    model_flux_gals[k][j + i * self.nz] = mod_val

        self.model_flux_gals = model_flux_gals

    def fit_calcs(self):
        '''
        Fits templates to the data, and records likelihood related info.
        '''
        self.star_min_chi = np.zeros(self.n_data)
        self.gal_min_chi = np.zeros(self.n_data)

        self.coeff_calcs()
        self.array_calcs()

        self._n_star_hyper_parms = self._n_star_template
        self._n_gal_hyper_parms = self._n_gal_template
        self._n_hyper_parms = self._n_star_hyper_parms + \
            self._n_gal_hyper_parms + 2

        self.count_tot = 0

    def coeff_calcs(self):
        '''
        Fits star and galaxy templates, determines the coefficient priors.
        '''
        print('\nCalculating values for coefficient priors...')

        self.star_coeff_mean = np.zeros(self._n_star_template)
        self.star_coeff_var = np.zeros(self._n_star_template)
        star_chi2 = np.zeros(self.n_data)
        star_ln_coeff = np.zeros(self.n_data)
        star_coeffvar = np.zeros(self.n_data)

        for i in xrange(self._n_star_template):
            star_ln_coeff, star_coeffvar, star_chi2 = self.fit_star_template(i)
            self.coeff_mean_var(1, i, star_ln_coeff, star_coeffvar, star_chi2)

        self.gal_coeff_mean = np.zeros(self._n_gal_template * self.nz)
        self.gal_coeff_var = np.zeros(self._n_gal_template * self.nz)
        gal_chi2 = np.zeros(self.n_data)
        gal_ln_coeff = np.zeros(self.n_data)
        gal_coeffvar = np.zeros(self.n_data)

        for i in xrange(self._n_gal_template):
            for j in xrange(self.nz):
                gal_ln_coeff, gal_coeffvar, gal_chi2 = \
                    self.fit_gal_template(i, j)
                self.coeff_mean_var(0.0, i * self.nz + j, gal_ln_coeff,
                                    gal_coeffvar, gal_chi2)

    def fit_star_template(self, i):
        '''
        For a given star template to all the data and record the coefficients,
        their uncertainties, and the associated chi2.
        '''
        mflux = self.model_flux_stars[:, i][:, np.newaxis]
        dflux = self.data_flux * self.flux_unit_factor
        dflux_err = self.data_flux_err * self.flux_unit_factor

        lh = np.sum(mflux * mflux / dflux_err / dflux_err, axis=0)
        rh = np.sum(mflux * dflux / dflux_err / dflux_err, axis=0)

        coeff = rh / lh
        coeff_err = np.sqrt(1.0 / lh)

        coeff[coeff <= 0.0] = 1.0

        chi2 = np.sum((dflux - coeff * mflux) * (dflux - coeff * mflux)
                      / dflux_err / dflux_err, axis=0)

        star_ln_coeff = np.log(coeff)
        star_coeffvar = coeff_err**2.0
        star_chi2 = chi2

        return star_ln_coeff, star_coeffvar, star_chi2

    def fit_gal_template(self, i, j):
        '''
        Fits a given galaxy template to all the data, and record
        the coefficients, their uncertainties, and the associated chi2.
        '''
        mflux = self.model_flux_gals[:, j + i * self.nz][:, np.newaxis]
        dflux = self.data_flux * self.flux_unit_factor
        dflux_err = self.data_flux_err * self.flux_unit_factor

        lh = np.sum(mflux * mflux / dflux_err / dflux_err, axis=0)
        rh = np.sum(mflux * dflux / dflux_err / dflux_err, axis=0)

        coeff = rh / lh
        coeff_err = np.sqrt(1.0 / lh)

        coeff[coeff <= 0.0] = 1.0

        chi2 = np.sum((dflux - coeff * mflux) * (dflux - coeff * mflux)
                      / dflux_err / dflux_err, axis=0)

        gal_ln_coeff = np.log(coeff)
        gal_coeffvar = coeff_err**2.0
        gal_chi2 = chi2

        return gal_ln_coeff, gal_coeffvar, gal_chi2

    def coeff_mean_var(self, flag, i, ln_coeff_in, coeff_var, chi2):
        '''
        Calculate the mean and variance of coefficients produced by fitting
        templates to all the data.
        '''
        weight = np.ones(self.n_data) / coeff_var
        ln_coeff = ln_coeff_in

        if flag == 0:
            self.gal_coeff_mean[i] = np.average(ln_coeff, weights=weight)
            self.gal_coeff_var[i] = np.average(
                (ln_coeff - self.gal_coeff_mean[i])**2.0, weights=weight)
            if self.n_data < 10:
                self.gal_coeff_var[i] = self.gal_coeff_mean[i] * -0.5

            self.gal_coeff_var[i] = np.absolute(self.gal_coeff_var[i])

        if flag == 1:
            self.star_coeff_mean[i] = np.average(ln_coeff, weights=weight)
            self.star_coeff_var[i] = np.average(
                (ln_coeff - self.star_coeff_mean[i])**2.0, weights=weight)
            if self.n_data < 10:
                self.star_coeff_var[i] = self.star_coeff_mean[i] * -0.5

            self.star_coeff_var[i] = np.absolute(self.star_coeff_var[i])

    def array_calcs(self):
        '''
        Class functions to save the marginalized likelihood (over coefficients)
        for all templates above the probability threshold.
        '''
        if self.pool is not None:
            M = self.pool.map
        else:
            M = map

        map_calc_P_F_kS = partial(wrap_calc_P_F_kS, hbsgsep=self)

        star_results = list(M(map_calc_P_F_kS,
                            [i for i in xrange(self.n_data)]))

        # for i in xrange(self.n_data):
        #     self.calc_P_F_kS(i)

        self.star_array = np.array(star_results)

        if self.pool is not None and not self.pool.is_master():
            self.pool.wait()
            self.exit(0)

        map_calc_P_F_kG = partial(wrap_calc_P_F_kG, hbsgsep=self)

        gal_results = list(M(map_calc_P_F_kG,
                           [i for i in xrange(self.n_data)]))

        self.gal_array = np.array(gal_results)

        if self.pool is not None and not self.pool.is_master():
            self.pool.wait()
            self.exit(0)

        # self.gal_array = np.zeros((self.n_data, self._n_gal_template))

        # for i in xrange(self.n_data):
        #     self.calc_P_F_kG(i)

        #if self.pool is not None:
        #    self.pool.close()
        #    self.pool = emcee.utils.MPIPool(loadbalance=True)
        #    if not self.pool.is_master():
        #        self.pool.wait()
        #        sys.exit(0)

    def calc_P_F_kS(self, i):
        '''
        Calculate the likelihood of each star for each template,
        marginalized over the coefficient of the fit.
        '''
        mflux = self.model_flux_stars
        dflux = self.data_flux[:, i] * self.flux_unit_factor
        dflux = dflux[:, np.newaxis]
        dflux_err = self.data_flux_err[:, i] * self.flux_unit_factor
        dflux_err = dflux_err[:, np.newaxis]

        det_dflux_err = np.prod(dflux_err)
        lh = np.sum(mflux * mflux / dflux_err / dflux_err, axis=0)
        rh = np.sum(mflux * dflux / dflux_err / dflux_err, axis=0)

        coeff = rh / lh
        coeff_err = np.sqrt(1.0 / lh)

        chi2 = np.sum((dflux - coeff * mflux) * (dflux - coeff * mflux)
                      / dflux_err / dflux_err, axis=0)

        temp_P_F_KS = np.zeros(self._n_star_template)
        p_F_CkS_max = np.zeros(self._n_star_template)
        p_F_CkS = np.zeros(self._n_star_template)
        p_lnC_kS = np.zeros(self._n_star_template)
        p_C_kS = np.zeros(self._n_star_template)

        for k in xrange(self.ncstep):
            h = coeff_err * 3.0 / (self.ncstep - 1) * 2.0

            if k == 0 | k == self.ncstep - 2:
                w = 3.0 / 8.0
            elif k == 1 | k == self.ncstep - 3:
                w = 7.0 / 6.0
            elif k == 2 | k == self.ncstep - 4:
                w = 23.0 / 24.0
            else:
                w = 1.0

            cval = coeff - h * (self.ncstep - 1) / 2.0 + h * k

            mask = (cval > 0.0) & (chi2 < 1400.0)

            p_F_CkS_max[mask] = 1.0 / np.sqrt(2.0 * np.pi) / det_dflux_err * \
                np.exp(-0.5 * chi2[mask])

            p_F_CkS[mask] = p_F_CkS_max[mask] * np.exp(
                -0.5 * (cval[mask] - coeff[mask]) * (cval[mask] - coeff[mask])
                / coeff_err[mask] / coeff_err[mask])

            p_lnC_kS[mask] = 1.0 / np.sqrt(
                2.0 * np.pi * self.star_coeff_var[mask]) * np.exp(
                    -0.5 * (np.log(cval[mask]) - self.star_coeff_mean[mask])
                    * (np.log(cval[mask]) - self.star_coeff_mean[mask])
                    / self.star_coeff_var[mask])

            p_C_kS[mask] = p_lnC_kS[mask] / cval[mask]

            temp_P_F_KS[mask] += w * h[mask] * p_F_CkS[mask] * p_C_kS[mask]

            p_F_CkS[~mask] = self.p_floor
            p_C_kS[~mask] = self.p_floor
            temp_P_F_KS[~mask] += self.p_floor

        # max_P_F_kS = temp_P_F_KS.max()
        min_chi2 = chi2.min()

        # self.star_array[i] = temp_P_F_KS
        self.star_min_chi[i] = min_chi2

        if i == 0:
            print('Marginalizing over star templates...')

        return temp_P_F_KS

    def calc_P_F_kG(self, i):
        '''
        '''
        mflux = self.model_flux_gals
        dflux = self.data_flux[:, i] * self.flux_unit_factor
        dflux = dflux[:, np.newaxis]
        dflux_err = self.data_flux_err[:, i] * self.flux_unit_factor
        dflux_err = dflux_err[:, np.newaxis]

        det_dflux_err = np.prod(dflux_err)

        lh = np.sum(mflux * mflux / dflux_err / dflux_err, axis=0)
        rh = np.sum(mflux * dflux / dflux_err / dflux_err, axis=0)

        coeff = rh / lh
        coeff_err = np.sqrt(1.0 / lh)

        chi2 = np.sum((dflux - mflux * coeff) * (dflux - coeff * mflux)
                      / dflux_err / dflux_err, axis=0)

        temp_p_F_kG = np.zeros(self._n_gal_template)

        p_F_CkzG_max = np.zeros(self._n_gal_template * self.nz)
        p_F_CkzG = np.zeros(self._n_gal_template * self.nz)
        p_lnC_kzG = np.zeros(self._n_gal_template * self.nz)
        p_C_kzG = np.zeros(self._n_gal_template * self.nz)
        p_F_kzG = np.zeros(self._n_gal_template * self.nz)

        for k in xrange(self.ncstep):
            h = coeff_err * 3.0 / (self.ncstep - 1) * 2.0

            if k == 0 or k == self.ncstep - 2:
                w = 3.0 / 8.0
            elif k == 1 or k == self.ncstep - 3:
                w = 7.0 / 6.0
            elif k == 2 or k == self.ncstep - 4:
                w = 23.0 / 24.0
            else:
                w = 1.0

            cval = coeff - h * (self.ncstep - 1) / 2.0 + h * k

            mask = (cval > 0.0) & (chi2 < 1400.0)

            p_F_CkzG_max[mask] = 1.0 / np.sqrt(2.0 * np.pi) / det_dflux_err * \
                np.exp(-0.5 * chi2[mask])
            p_F_CkzG[mask] = p_F_CkzG_max[mask] * np.exp(
                -0.5 * (cval[mask] - coeff[mask])
                * (cval[mask] - coeff[mask])
                / coeff_err[mask]
                / coeff_err[mask])
            p_lnC_kzG[mask] = 1.0 / np.sqrt(
                2.0 * np.pi * self.gal_coeff_var[mask]) * np.exp(
                    -0.5 * (np.log(cval[mask]) - self.gal_coeff_mean[mask])
                    * (np.log(cval[mask]) - self.gal_coeff_mean[mask])
                    / self.gal_coeff_var[mask])
            p_C_kzG[mask] = p_lnC_kzG[mask] / cval[mask]
            p_F_kzG[mask] += w * h[mask] * p_F_CkzG[mask] * p_C_kzG[mask]

            p_F_CkzG[~mask] = self.p_floor
            p_C_kzG[~mask] = self.p_floor
            p_F_kzG[~mask] += self.p_floor

        temp_p_F_kG = np.reshape(p_F_kzG / self.nz,
                                 (self._n_gal_template, self.nz))
        temp_p_F_kG = np.sum(temp_p_F_kG, axis=1)

        max_p_F_kG = temp_p_F_kG.max()
        min_chi2 = chi2.min()

        mask = (temp_p_F_kG / max_p_F_kG > self.prob_frac)

        # self.gal_array[i] = temp_p_F_kG

        self.gal_min_chi[i] = min_chi2

        if i == 0:
            print('Marginalizing over galaxy templates...')

        return temp_p_F_kG

    def sample(self):
        '''
        '''
        self.count = 0

        p0 = [np.zeros(self._n_hyper_parms) for i in xrange(self.nwalkers)]

        for i in xrange(self.nwalkers):
            p0[i][: self._n_star_hyper_parms] = np.random.dirichlet(
                np.ones(self._n_star_hyper_parms))

            start = self._n_star_hyper_parms
            end = self._n_star_hyper_parms + self._n_gal_hyper_parms
            p0[i][start: end] = np.random.dirichlet(
                np.ones(self._n_gal_hyper_parms))

            p0[i][self._n_hyper_parms - 2:] = np.random.dirichlet(np.ones(2))

        if self.pool is not None:
            print('\nStarting Emcee sampler with '
                  '{0} cores...\n'.format(self.pool.size))
        else:
            print('\nStarting Emcee sampler...\n')

        sampler = emcee.EnsembleSampler(self.nwalkers,
                                        self._n_hyper_parms,
                                        wrap_lnlike,
                                        args=[self],
                                        pool=self.pool)

        # burn-in step
        pos, prob, state = sampler.run_mcmc(p0, self.nburn)
        print('Finished burn in...')
        sampler.reset()
        print('Starting MCMC...')
        pos, prob, state = sampler.run_mcmc(pos, self.niter, rstate0=state)

        # the mean acceptance fraction should be between 0.25 and 0.5
        # if everything worked OK.
        print("\nMean acceptance fraction: "
              "{0:.3f}".format(np.mean(sampler.acceptance_fraction)))

        self.flatchain = sampler.flatchain
        self.write_hyp_pars()
        del self.flatchain

        if self.pool is not None and not self.pool.is_master():
            self.pool.wait()
            sys.exit(0)

    def _get_hyp_parms(self, x):
        '''
        '''
        hyp_parms = np.zeros(self._n_hyper_parms, dtype=np.float128)

        # read in current value of hyperparameters from optimizer
        # star hyperparameters
        hyp_parms_stars = np.absolute(x[: self._n_star_hyper_parms])
        hyp_parms_stars = hyp_parms_stars / hyp_parms_stars.sum()
        hyp_parms[: self._n_star_hyper_parms] = hyp_parms_stars

        # galaxy hyperparameters
        start = self._n_star_hyper_parms
        stop = self._n_star_hyper_parms + self._n_gal_hyper_parms
        hyp_parms_gals = np.absolute(x[start: stop])
        hyp_parms_gals = hyp_parms_gals / hyp_parms_gals.sum()
        hyp_parms[start: stop] = hyp_parms_gals

        start = self._n_hyper_parms - 2
        stop = self._n_hyper_parms
        hyp_parms_weights = np.absolute(x[start: stop])
        hyp_parms_weights = hyp_parms_weights / hyp_parms_weights.sum()
        hyp_parms[start: stop] = hyp_parms_weights

        return hyp_parms

    def save_proba(self):
        '''
        '''
        if self.pool is not None:
            M = self.pool.map
        else:
            M = map

        calc_ith_ln_like = partial(wrap_calc_ith_ln_like, hbsgsep=self)

        print('\nCalculating S/G probabilities...')
        results = list(M(calc_ith_ln_like,
                       [i for i in xrange(self.n_data)]))
        results = np.array(results)

        #mean = np.mean(results, axis=1)
        median = results[:, 0]
        lower = results[:, 1]
        upper = results[:, 2]

        fname, ext = os.path.splitext(self.ln_like_out_file)

        #np.save(fname + '_mean' + ext, mean)
        np.save(fname + '_median' + ext, median)
        np.save(fname + '_84th_percentile' + ext, upper)
        np.save(fname + '_16th_percentile' + ext, lower)

        print('\nMedian: {0:d} stars, {1:d} galaxies at '
              'probability cut = 0.5'.format((median >= 0).sum(),
                                      (median < 0).sum()))
        #print('Mean: {0:d} stars, {1:d} galaxies at '
        #      'probability cut = 0.5'.format((mean >= 0).sum(),
        #                              (mean < 0).sum()))

        proba_median, proba_lower, proba_upper = map(
            lambda x: np.exp(x) / (1.0 + np.exp(x)),
            [median, lower, upper])

        # savetxt doesn't accpet unicode
        np.savetxt(self.median_proba_out_file,
                   np.array((proba_median, proba_lower, proba_upper)).T,
                   delimiter=',',
                   fmt=('%1.4f ' * 3).encode('ascii'))
        #np.savetxt(self.mean_proba_out_file, proba_mean,
        #           fmt='%1.4f'.encode('ascii'))
        #print('\nSaving to ' + self.ln_like_out_file)
        #np.save(self.ln_like_out_file, results)

        return results

    def calc_ith_ln_like(self, i):
        '''
        '''
        flatchain = np.load(self.hyp_out_file)
        h = np.absolute(flatchain)
        hyp_parms_stars = h[:, : self._n_star_hyper_parms]
        hyp_parms_stars = hyp_parms_stars / hyp_parms_stars.sum()
        hyp_parms_gals = h[:, self._n_star_hyper_parms:
                           self._n_star_hyper_parms +
                           self._n_gal_hyper_parms]
        hyp_parms_gals = hyp_parms_gals / hyp_parms_gals.sum()
        hyp_parms_weights = h[:, self._n_hyper_parms - 2:
                                      self._n_hyper_parms]
        hyp_parms_weights = hyp_parms_weights / hyp_parms_weights.sum()

        p_F_S = np.dot(hyp_parms_stars, self.star_array.T[:, i])
        p_F_G = np.dot(hyp_parms_gals, self.gal_array.T[:, i])

        p_S = p_F_S * hyp_parms_weights[:, 0]
        p_G = p_F_G * hyp_parms_weights[:, 1]
        p_S[p_S <= 0] = self.p_floor
        p_G[p_G <= 0] = self.p_floor

        ln_like = np.log(p_S) - np.log(p_G)

        #ln_like = np.load(self.ln_like_out_file)
        # uncertainties based on the 16th, 50th, and 84th percentiles
        median = np.percentile(ln_like, 50, axis=0)
        lower = np.percentile(ln_like, 16, axis=0)
        upper = np.percentile(ln_like, 84, axis=0)
        # mean is in separte file
        #mean = np.mean(ln_like, axis=0)

        if i % 1000 == 0:
            print('Saving {}th data point...'.format(i))
        
        return median, lower, upper

        # np.save(self.ln_like_out_file + '.' + str(i // self.nchunk) +
        # '.npy', ln_like)

    def loglikelihood(self, x):
        '''
        '''
        hyp_parms = self._get_hyp_parms(x)
        hyp_parms_stars = hyp_parms[: self._n_star_hyper_parms]
        hyp_parms_gals = hyp_parms[self._n_star_hyper_parms:
                                   self._n_star_hyper_parms +
                                   self._n_gal_hyper_parms]
        hyp_parms_weights = hyp_parms[self._n_star_hyper_parms - 2:
                                      self._n_star_hyper_parms]

        # calculate P_F_S
        p_F_S = np.dot(hyp_parms_stars, self.star_array.T)
        p_F_G = np.dot(hyp_parms_gals, self.gal_array.T)

        p_F_G[p_F_G <= 0] = self.p_floor
        p_F_S[p_F_S <= 0] = self.p_floor

        # ln_p_tot_G = np.log(p_F_G).sum()
        # ln_p_tot_S = np.log(p_F_S).sum()
        ln_p_tot = np.log(p_F_S * hyp_parms_weights[0] +
                          p_F_G * hyp_parms_weights[1]).sum()

        # record starting point for the record
        # if self.count_tot == 0:
        #     self.ln_p_tot_start = ln_p_tot
        #     self.old_ln_p_tot = ln_p_tot

        # if ln_p_tot > self.old_ln_p_tot:
        #     self.write_hyp_pars(hyp_parms)
        #     self.old_ln_p_tot = ln_p_tot

        # doesn't work properly if run in parallel
        if self.pool is None and self.count_tot % self.nwalkers == 0:
            print("{0:d} iterations, {1:.0f} seconds".format(
                self.count_tot // self.nwalkers,
                time.clock() - self.last_clock))

        self.count += 1
        self.count_tot += 1

        # make ln_p_tot positive for scipy minimizer
        # ln_p_tot = np.absolute(ln_p_tot)

        return ln_p_tot

    def read_hyp_pars(self):
        '''
        Read hyperparameters from an npy file.
        '''
        hyp_parms = np.load(self.hyp_in_file)

        return hyp_parms

    def write_hyp_pars(self):
        '''
        Save hyperparamaters in an npy file.
        '''
        file_path = self.hyp_out_file
        dir_path = os.path.dirname(os.path.abspath(file_path))

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        np.save(file_path, self.flatchain)


# The following wrap_* functions are necessary because only global functions
# can be pickled while methods cannot be pickled.
def wrap_lnlike(x, hbsgsep):
    return hbsgsep.loglikelihood(x)


def wrap_calc_P_F_kS(i, hbsgsep):
    return hbsgsep.calc_P_F_kS(i)


def wrap_calc_P_F_kG(i, hbsgsep):
    return hbsgsep.calc_P_F_kG(i)


def wrap_calc_ith_ln_like(i, hbsgsep):
    return hbsgsep.calc_ith_ln_like(i)


def main():
    '''
    A serial run. For parallized version using MPI, use run_mpi.py
    '''
    clf = HBSGC()

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

    if clf.min_chi2_write:
        clf.save_min_chi2()


if __name__ == '__main__':
    sys.exit(main())
