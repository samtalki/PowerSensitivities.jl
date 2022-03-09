
module PowerSensitivities

using PowerModels,OPFLearn
using JuMP, Ipopt 
using TimeSeries
using LinearAlgebra
using Gadfly
using Convex 
using Mosek
import SparseArrays

include("jacobian_matrix.jl")
include("sensitivities.jl")
include("matrix_completion.jl")

"""
Jacobian info
"""
struct JacobianInvertibility
    dpdth::Bool
    dqdth::Bool
    J::Bool
end

"""
Test the positive definite angle submatrices of the Jacobian
"""
function test_pd_angle_matrices(network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/")
    names = readdir(network_data_path)
    paths = readdir(network_data_path,join=true)
    results = Dict()
    for (name,path) in zip(names,paths)
        data = make_basic_network(parse_file(path))
        J = try
            calc_jacobian_matrix(data,1)
        catch
            println("Jacobian cannot be computed for "*name)
            continue
        end
        info = Dict(
            "dpdth_pd" => ispd(Matrix(J.dpdth)),
            "dqdth_nd" => isnd(Matrix(J.dqdth)),
            "dpdth+dqdth_invertible" => isinvertible(Matrix(J.dpdth).+Matrix(J.dqdth)),
            "J_invertible" => isinvertible(Matrix(J.matrix)),
            "dpdth_norm" => norm(J.dpdth),
            "dqdth_norm" => norm(J.dqdth)
        )
        results[name] = info
    end
    return results
end

"""
Checks if a matrix M is positive definite
"""
ispd(M) = all([real(eig)>0 for eig in eigvals(M)])

"""
Checks if a matrix M is negative definite
"""
isnd(M) = all([real(eig)<0 for eig in eigvals(M)])

isinvertible(x,ϵ=1e-12) = applicable(inv, x) && norm(I(size(x)[1]) - inv(x)*x) < ϵ

end
