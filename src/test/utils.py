import numpy as np
from numba import jit

def spectral_analysis(S):
    """Computes normalized and commulative spectral analysis for a given matrix with the SVD"""
    u,sigma,vt = np.linalg.svd(S)
    sigma_total = np.sum(sigma)
    normed_sigma = sigma/sigma_total
    cum_sigma = np.cumsum(normed_sigma)
    return cum_sigma,normed_sigma

@jit
def calc_vector_rel_err(v_est,v_true):
    """computes the relative error of v_est to v_true"""
    return np.linalg.norm(v_est-v_true)/np.linalg.norm(v_true)

@jit
def calc_vector_rmse(yhat,y):
    return np.sqrt(np.nanmean((yhat-y)**2))
@jit
def calc_mae(yhat,y):
    return np.nanmean(np.abs(yhat-y))

@jit
def calc_max_error(yhat,y):
    abs_err = np.abs(yhat-y)
    return np.max(abs_err)