import numpy as np
import cvxpy as cp
from measurements import constrained_linear_measurement_operator

def predict(S_hat,dx):
    """Predict voltage magnitude perturbations from complex power"""
    assert S_hat.shape[1] == dx.shape[0]
    dv_hat = S_hat @ dx
    return dv_hat

def calc_rel_err(M_hat,M):
    """L2 norm relative error"""
    return (np.linalg.norm(M_hat-M,))/(np.linalg.norm(M,))

def K(pf):
    """Compute the K matrix given bus power factors"""
    k = lambda pf_ : np.sqrt(1-pf_**2)/pf_
    return np.diag([k(pf_i) for pf_i in pf])

    
def mat_rec_problem(S,S_0,dx,dv,lamb,o,delta):
    """
    Make the matrix recovery problem
    Params:
        S: sensitivity matrix decision variable
        S_0: Inital sensitivity matrix. Optionally warm start, partially filled sensitivity matrix
        dx: complex power perturbations
        dv: voltage magnitude perturbations
        o: measurement operator matrix    
        delta: penalty to match this warm start.
    """
    known_values_idx = np.where(o==1)
    known_values = S_0[known_values_idx] #known values
    #Matrix recovery loss function with nuclear norm regularization
    obj = cp.Minimize(cp.sum_squares(S@dx - dv) + lamb*cp.norm(S,"nuc"))
    #Constraint on straying too far from known values
    cons = [cp.norm(S[known_values_idx] - known_values) <= delta] 
    prob = cp.Problem(objective=obj,constraints=cons)
    return prob

def mat_rec_solution(S_0,dx,dv,lamb,o,delta=1e-3,verbose=True):
    """
    Compute the sensitivity matrix recovery solution
    Params:
        S_0: Initialized sensitivity matrix
        data: (dx,dv) complex power/voltage magnitude deviations
        o: constrained measurement operator
        lambd: nuclear norm penalty
        delta: radius of the masked l2 ball constraint
    Returns:
        (S_fit,prob)
    """
    S = cp.Variable(S_0.shape)
    prob = mat_rec_problem(S,S_0,dx,dv,lamb,o,delta)
    prob.solve(verbose=verbose)
    return S,prob


def mat_rec_problem_implicit(S,S_0,K,dx,dv,lamb,o,delta1,delta2):
    """
    Make the matrix recovery problem with implicit power factor constraints
    Params:
        S: sensitivity matrix decision variable
        S_0: Inital sensitivity matrix. Optionally warm start, partially filled sensitivity matrix
        dx: complex power perturbations
        dv: voltage magnitude perturbations
        o: measurement operator matrix    
        delta: penalty to match this warm start.
    """
    (n,_) = S.shape
    known_values_idx = np.where(o==1)
    known_values = S_0[known_values_idx] #known values
    #Matrix recovery loss function with nuclear norm regularization
    obj = cp.Minimize(cp.sum_squares(S@dx - dv) + lamb*cp.norm(S,"nuc"))
    #Constraint on straying too far from known values
    cons = [
        cp.norm(S[known_values_idx] - known_values) <= delta1,
        cp.norm((S[:,:n] + S[:,n:]@K)@dx[:n,:] - dv) <= delta2 #implicit function theorem constraint
    ] 
    prob = cp.Problem(objective=obj,constraints=cons)
    return prob

