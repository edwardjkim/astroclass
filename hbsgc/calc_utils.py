'''
Matrix calculation utils for HBSGC.
'''

import numpy as np


def dot(left_array, right_array, ncols=1):
    '''
    Performs matrix dot product column by column to avoid memory error.
    '''
    m, left_inner = left_array.shape
    right_inner, n = right_array.shape

    if not left_inner == right_inner:
        raise ValueError('shapes are not alighed')

    out = np.zeros((m, n))
    for i in xrange(0, n, ncols):
        out[:, i] = np.dot(left_array, right_array[:, i])

    return out
