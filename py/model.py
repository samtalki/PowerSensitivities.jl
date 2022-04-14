import numpy as np
from sklearn.linear_model import ElasticNetCV,MultiTaskElasticNetCV

random_state=42

def dx_estimator(S_tilde,dv,l1_ratio = [.1, .5, .7, .9, .95, .99, 1],max_iter=5e3):
    """Estimate complex power given a vector of voltage magnitude measurements delta v and a wide S tilde matrix"""
    if dv.ndim>1:
        reg = MultiTaskElasticNetCV(l1_ratio=l1_ratio,
                                    max_iter=max_iter,
                                    n_jobs=5,random_state=random_state)
    else:
        reg = ElasticNetCV(l1_ratio=l1_ratio,
                           max_iter=max_iter,
                           random_state=random_state)
    return reg


