import numpy as np
from numba import jit
import pandas as pd
from utils import make_deviations,make_S_tilde,make_sens_ts

@jit
def calc_complex_estimate(sens,data,K):
    """
    Given sensitivities, perturbations, and power factor matrice(s) K, compute dp_hat and dq_hat
    """
    dx,dv = data
    svp,svq = sens
    
    #Num buses
    n,_ = dv[0].shape
    
    #Decompose dx into the true dp and dq
    dp,dq = dx[:n,:],dx[n:,:]

    #Estimate of dp and dq
    dp_hat,dq_hat = np.zeros_like(dp),np.zeros_like(dq)

    for t,(svp_t,svq_t,dv_t) in enumerate(zip(svp,svq,dv.T)):
        S_dag_t = svp_t + svq_t@K #Make the implicit sensitivity matrix
        dp_hat[:,t] = np.linalg.inv(S_dag_t)@dv_t #Compute the active power estimate
        dq_hat[:,t] = K@dp_hat[:,t] #Compute the reactive power estimate
    
    
    return dp_hat,dq_hat