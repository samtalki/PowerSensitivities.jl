import numpy as np
from numba import jit
from measurements import constrained_linear_measurement_operator

def make_sens_ts(dvdp,dvdq,n=274):
    """Make timeseries of sensitivity matrices"""
    assert dvdp.shape[0] == dvdq.shape[0]
    m_tot = dvdp.shape[0] #total measurments
    m = int(m_tot/n) #Timeseries interval
    svp,svq = [],[] #timeseries list of svp and svq matrices
    for t in range(m):
        svp.append(dvdp[t*n:(t+1)*n,:])
        svq.append(dvdq[t*n:(t+1)*n,:])
    return {'svp':svp,
            'svq':svq}

def make_deviations(data):
    """Make finite differences"""
    p,q,v = data
    (dp,dq,dv) = [np.diff(x) for x in (p,q,v)]
    dx = np.vstack((dp,dq))
    return (dx,dv)

def make_S_tilde(svp,svq):
    """Make wide S_tilde matrix"""
    return np.vstack((svp.T,svq.T)).T

def make_S_0(S_tilde,pct_obs):
    """Make initial observed matrix"""
    O,_ = constrained_linear_measurement_operator(S_tilde,pct_obs)
    O,S_0 = np.asarray(O),np.multiply(np.asarray(O),S_tilde)
    return O,S_0

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