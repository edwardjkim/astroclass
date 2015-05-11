'''
Author: Edward J. Kim
Email: jkim575@ilinois.edu
'''
from __future__ import print_function, division, unicode_literals
import numpy as np
import os
from mlz.ml_codes import SOMZ

class BMC(object):
    '''
    Bayesian Model Combination (Monteith et al. 2011).
    A Hybrid Ensemble Approach to Star-galaxy Classification (Kim, Brunner & Carrasco Kind 2015).
    '''
    def __init__(self, q=3, ncomb=1200, nburn=200, ntop=10, icore=0):
        '''
        Constructor.
        
        Parameters
        ----------
        q: The number of combinations drawn each time from Dirichlet distribution.
        ncomb: The total number of ensembles.
        nburn: The number of iterations in the burn-in step.
        ntop: The side length of a rectangular SOM.
        icore: The number of cores (0 unless parallelized).
        '''
        self.q = q
        self.ncomb = ncomb
        self.nburn = nburn
        self.ntop = ntop
        self.icore = icore
        self.pfloor = 1e-300

    def get_SOM_cells(self, topology='grid', ntop=10, iterations=200,
                      periodic='no', iproc=0):
        '''
        Creates a SOM representation.
        Only supports rectangular topology.
        '''
        if iproc is not 0:
            print('Creating {}th SOM map...'.format(iproc))

        if topology == 'grid':
            self.n_cell = ntop**2

        # Calls the SOMZ mode
        self.M = SOMZ.SelfMap(self.X_train, self.Y_train,
                         topology=topology, Ntop=ntop,
                         iterations=iterations, periodic=periodic)
        # creates a map
        self.M.create_mapF()
        # evaluates it with the Y entered, or anyoher desired colum
        self.M.evaluate_map()

        cells_train = np.zeros(len(self.X_train))

        for i in range(len(self.X_train)):
            cells_train[i] = self.M.som_best_cell(self.X_train[i])[0]

        np.save(self.file_train_cells, cells_train)

    def get_log_likelihood(self, label, array, weight):
        '''
        The likelihood function.
        
        Parameters
        ----------
        label: A numpy array. The ground truth.
        array: The (probabilistic) outputs of base classifiers.  
        weight: The weights of base classifiers.
        
        Returns
        -------
        A numpy array.
        '''
        y = label
        x = np.sum(array * weight, axis=1)
        error = np.abs(y - x)
        
        # Cromwell's rule:
        # I beseech you, in the bowels of Christ,
        # think it possible that you may be mistaken.
        # and also the fact that log of 0 diverges.
        error[error == 0] = self.pfloor
        error[error == 1] = 1 - self.pfloor
        lnlike = np.log(1 - error).sum()

        return lnlike

    def fit(self, attr, X, Y, create_cells=True):
        '''
        Calculates weights in each SOM cell using the cross-vaildation data.
        
        Parameters
        ----------
        attr: Input attributes, eg magnitudes or colors.
        X: Probability estimates from cross-validation data.
        Y: Truth values of training data.
        create_cells: If false, use previous SOM map.
        
        Returns
        -------
        None. Use predict_proba() to create outputs.
        '''
        if create_cells:
            self.make_cells = True
            self.X_train = attr
            self.Y_train = Y
            self.folder_cells = 'cells'
            self.file_train_cells = 'cells/cv_cells.{}.npy'.format(self.icore)
            self.file_test_cells = 'cells/test_cells.{}.npy'.format(self.icore) 
            if not os.path.exists(self.folder_cells):
                os.mkdir(self.folder_cells)

            self.get_SOM_cells(ntop=self.ntop, iproc=self.icore)
        else:
            self.make_cells = False

        nxrows, nxcols = X.shape
        self.m = nxcols

        train_som_cell = np.load(self.file_train_cells)

        post_all = np.zeros(nxrows)

        # create folders
        folder_weights = 'weights'
        if not os.path.exists(folder_weights):
            os.mkdir(folder_weights)

        if not os.path.exists('results'):
            os.mkdir('results')

        for c in range(self.n_cell):

            icell = np.where(train_som_cell == c)[0]

            weight_all = np.zeros(((self.ncomb - self.nburn) * self.q, self.m))

            if len(icell) > 0:

                alpha = np.ones(self.m)
                log_p_comb = np.zeros(self.q)

                for i in range(self.ncomb):
                    weight = np.random.dirichlet(alpha, self.q)

                    for j in range(self.q):
                        w = self.get_log_likelihood(Y[icell],
                                                    X[icell],
                                                    weight[j])
                        log_p_comb[j] = w
                        
                        if i >= self.nburn:
                            weight_all[(i - self.nburn) * self.q + j] = weight[j]

                    best_weight = weight[log_p_comb.argmax()]
                    alpha += best_weight

                a = np.dot(X[icell], weight_all.T)

                post = np.sum(a, axis=1) / ((self.ncomb - self.nburn) * self.q)

                post_all[icell] = post

            np.save('weights/bmc_weights.{0}.{1}.npy'.format(self.icore, c),
                    weight_all)

        # save OOB probabilities
        np.save('results/bmc_oob.{0}.npy'.format(self.icore), post_all)

        del weight_all
        del post_all
        
        print('Done creating SOM map.')

    def predict_proba(self, attr, X):
        '''
        Make predictions of test data.
        
        Parameters
        ----------
        attr: Attributes of test data, eg magnitudes or colors.
        X: Probability estimates of base classifiers.
        
        Returns
        -------
        A numpy array.
        '''
        test_som_cell = np.zeros(len(attr))

        for i in range(len(attr)):
            test_som_cell[i] = self.M.som_best_cell(attr[i])[0]

        post_sum = np.zeros(len(X))

        for i in range(self.n_cell):
            icell = np.where(test_som_cell == i)[0]
            
            if len(icell) > 0:
                w = np.load('weights/bmc_weights.{0}.{1}.npy'.format(
                    self.icore, i))
                a = np.dot(X[icell], w.T)
                post = np.sum(a, axis=1) / len(w)
                post_sum[icell] = post
                
        return post_sum
