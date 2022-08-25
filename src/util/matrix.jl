const DEFAULT_COND_TOL = 1e6
using LinearAlgebra

"""
Checks if a matrix M is positive definite
"""
ispd(M) = all([real(eig)>0 for eig in eigvals(M)])

"""
Checks if a matrix M is negative definite
"""
isnd(M) = all([real(eig)<0 for eig in eigvals(M)])
isnsd(M,ϵ=1e-12) = all([real(eig)<=ϵ for eig in eigvals(M)])

"""
Checks if a matrix M is invertible
"""
# function isinvertible(M::Matrix; cond_tol::Number=DEFAULT_COND_TOL)
#     return issquare(M) && cond(M) < cond_tol
# end

function isinvertible(M::Matrix)
    return 0 ∉ eigvals(M)
end
"""
Returns distance between M and M transpose
"""
symmetricdiff(M) = norm(M-transpose(M))

"""
Check if symmetric part of a matrix is negative definite
"""
symmetric_part_nd(M) = isnd(0.5.*(M + transpose(M)))
symmetric_part_nsd(M) = isnsd(M + transpose(M./2))

"""
Check if symmetric part of a matrix is positive definite
"""
symmetric_part_pd(M) = ispd(0.5*(M +transpose(M)))
