import numpy as np

def linear_measurement_operator(S,pct_obs):
    """Returns the linear measurmeent OPERATOR, not the partiall observed matrix."""
    m,n = S.shape
     # pick which entries we observe uniformly at random
    O = np.random.binomial(1,p=pct_obs,size=(m,n))
    return O

def constrained_linear_measurement_operator(S,pct_obs):
    """Returns a linear measurement operator subject to the constraint that we observe at least one measurement from each column"""
    m,n = S.shape
    #pick which entries we observe uniformly at random
    O = np.random.binomial(1,p=pct_obs,size=(m,n))
    new_entries =0 #count the new entries
    for O_row in O:
        if 1 not in O_row: #if there is not a one in a O_row, add to it randomly.
            idx = np.random.randint(0,high=len(O_row))
            O_row[idx] = 1
            new_entries+=1
    for O_col in O.T:
        if 1 not in O_col:
            idx = np.random.randint(0,high=len(O_col))
            O_col[idx] = 1
            new_entries+=1
    pct_obs = ((m*n)*pct_obs + new_entries)/(m*n)
    return O,pct_obs
            

def linear_col_measurement_operator(S,pct_obs):
    """Generates an MxN incomplete matrix where pct_unobs COLUMNS are not seen"""
    m,n = S.shape
    operator_shape = np.ones((m,n))
    #pick which columns we observe uniformly at random
    O = np.random.binomial(1,p=pct_obs,size=(n))
    return np.multiply(operator_shape,O)


def linear_measurements(S,pct_obs):
    """Generates an MxN incomplete matrix where pct_unobs ENTRIES are not seen"""
    O = linear_measurement_operator(S,pct_obs)
    return np.multiply(S,O)

def linear_col_measurements(S,pct_obs):
    """Generates an MxN incomplete matrix where pct_unobs COLUMNS are not seen"""
    O = linear_col_measurement_operator(S,pct_obs)
    return np.multiply(S,O)
