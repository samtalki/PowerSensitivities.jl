from matrix_sensitivities import predict,mat_rec_solution,mat_rec_problem
from cvxpylayers.jax import CvxpyLayer
import cvxpy as cp
import jax.numpy as jnp
import jax


"""
Returns the gradient of the matrix recovery solution with respect to the matrix recovery parameters
"""
def diff_matrix_recovery_solution(problem,S,lambd,gamma):
    assert problem.is_dpp()
    layer = CvxpyLayer(problem, parameters=[lambd, gamma], variables=[S])
    key = jax.random.PRNGKey(0)
    key, k1, k2 = jax.random.split(key, 3)
    A_jax = jax.random.normal(k1, shape=(m, n))
    b_jax = jax.random.normal(k2, shape=(m,))

    solution, = layer(A_jax, b_jax)

    # compute the gradient of the summed solution with respect to A, b
    dcvxpylayer = jax.grad(lambda A, b: sum(layer(A, b)[0]), argnums=[0, 1])
    gradA, gradb = dcvxpylayer(A_jax, b_jax)



def matrix_recovery_loss(S,dx,dv,params):
    (lambd,delta) = params #hyperparams
    return jnp.linalg.norm(S@dx-dv,ord='fro') + lambd*jnp.linalg.norm(S,ord='nuc')
