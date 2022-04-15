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


def mat_rec_solution(S_0,dx,dv,lamb,o,delta=1e-3,verbose=True):
    """
    Compute the sensitivity matrix recovery solution with a single nuclear norm penalty applied to both the active and reactive power matrices equally.
    Params:
        S_0: Initialized sensitivity matrix
        (dx,dv) complex power/voltage magnitude deviations
        lamb: nuclear norm penalty
        o: constrained measurement operator
        delta: radius of the masked l2 ball constraint
    Returns:
        (S_fit,prob)
    """
    S = cp.Variable(S_0.shape)
    prob = mat_rec_problem(S,S_0,dx,dv,lamb,o,delta)
    prob.solve(verbose=verbose,requires_grad=True)
    return S,prob

def mat_rec_problem(S,S_0,dx,dv,lamb,o,delta):
    """
    Make the matrix recovery problem with a single nuclear norm penalty applied to both the active and reactive power matrices equally.
    Params:
        S: wide sensitivity matrix decision variable
        S_0: Inital sensitivity matrix. Optionally warm start, partially filled sensitivity matrix
        dx: complex power perturbations
        dv: voltage magnitude perturbations
        o: measurement operator matrix    
        lamb: Nuclear norm penalty on the active and reactive power sensitivity matrices
        delta: penalty to match this warm start.
    """
    (n,tn) = S_0.shape #Wide Stilde matrix shape
    known_values_idx = np.where(o==1)
    known_values = S_0[known_values_idx] #known values
    
    #Unpack dx into (dp,dq)
    dp,dq = dx[:n,:], dx[n:,:]
    
    #Matrix recovery loss function with nuclear norm regularization
    obj = cp.Minimize(cp.sum_squares((S[:,:n]@dp + S[:,n:]@dq) - dv) + lamb*cp.norm(S[:,:n],"nuc") + lamb*cp.norm(S[:,n:],"nuc"))
    
    #Constraints
    cons = [
        (cp.norm(S[known_values_idx] - known_values) <= delta),  #Constraint on straying too far from known values
        (S[:,:n] >> 0), #The voltage-active power sensitivity matrix must be positive semidefinite.
        (S[:,n:] >> 0) #The voltage-reactive power sensitivity matrix must be positive semidefinite.
    ] 
    
    #Return the problem
    prob = cp.Problem(objective=obj,constraints=cons)
    return prob

def mat_rec_problem_implicit(S,S_0,K,dx,dv,lamb,o,delta):
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
    (n,tn) = S_0.shape #Wide Stilde matrix shape
    known_values_idx = np.where(o==1)
    known_values = S_0[known_values_idx] #known values
    
    #Unpack dx into (dp,dq)
    dp,dq = dx[:n,:], dx[n:,:]
    
    #Matrix recovery loss function with nuclear norm regularization
    obj = cp.Minimize(
        cp.sum_squares(S@dx - dv) + lamb*cp.norm(S,"nuc") + cp.sum_squares((S[:,:n] + S[:,n:]@K)@dx[:n,:] - dv)
    )
    
    #Constraint on straying too far from known values
    cons = [cp.norm(S[known_values_idx] - known_values) <= delta] 
    
    #Return the problem
    prob = cp.Problem(objective=obj,constraints=cons)
    return prob


    
def mat_rec_problem_rect_nuc(S,S_0,dx,dv,p_lamb,q_lambd,o,delta):
    """
    Make sensitivity matrix recovery problem 
    With individual nuclear norm penalties on the active and reactive power matrices.
    Params:
        S: wide sensitivity matrix decision variable
        S_0: Inital sensitivity matrix. Optionally warm start, partially filled sensitivity matrix
        dx: complex power perturbations
        dv: voltage magnitude perturbations
        o: measurement operator matrix    
        p_lambd: Nuclear norm penalty on the active power sensitivity matrix
        q_lambd: Nuclear norm penalty on the reactive power sensitivity matrix
        delta: penalty to match this warm start.
    """
    (n,tn) = S_0.shape #Wide Stilde matrix shape
    known_values_idx = np.where(o==1)
    known_values = S_0[known_values_idx] #known values
    
    #Unpack dx into (dp,dq)
    dp,dq = dx[:n,:], dx[n:,:]
    
    #Matrix recovery loss function with nuclear norm regularization
    obj = cp.Minimize(cp.sum_squares((S[:,:n]@dp + S[:,n:]@dq) - dv) + p_lamb*cp.norm(S[:,:n],"nuc") + q_lambd*cp.norm(S[:,n:],"nuc"))
    cons = [
        (cp.norm(S[known_values_idx] - known_values) <= delta),  #Constraint on straying too far from known values
        (S[:,:n] >> 0), #The voltage-active power sensitivity matrix must be positive semidefinite.
        (S[:,n:] >> 0), #The voltage-reactive power sensitivity matrix must be positive semidefinite.
        (cp.diag(S[:,n:]) >= cp.diag(S[:,:n]))
    ] 
    prob = cp.Problem(objective=obj,constraints=cons)
    return prob


