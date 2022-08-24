#using PowerSensitivities: VoltageSensitivityMatrix,PowerFlowJacobian
# include("jacobian_matrix.jl")
# include("../sens/voltage.jl")
# using LinearAlgebra
using LinearAlgebra:cond

"""
Spectral analysis struct
    Parameters:
        - matrix: The matrix to be analyzed
        - svd: The SVD object for `matrix`
        - scuml: the cummulative distribution of the singular values of the matrix
        - snorm: the normalized distribution of the singular values of the matrix
"""
struct SpectralAnalysis
    matrix::Union{AbstractMatrix,PowerFlowJacobian,VoltageSensitivityMatrix}
    svd::SVD
    scuml::AbstractArray
    snorm::AbstractArray
end


"""
Given an arbitrary matrix A∈R(mxn), analyze the normalized and cummulative singular values of A
"""
function calc_spectral_analysis(A::AbstractMatrix)
    A_svd = svd(A)
    Σ = A_svd.S
    Σnorm = Σ./sum(Σ)
    Σcuml = cumsum(Σ)
    return SpectralAnalysis(A,A_svd,Σcuml,Σnorm)
end

"""
Given a network data dict and a type, construct a SpectralAnalysis of the voltage sensitivity matrix (default) or the power flow jacobian.
Params:
    -net: PowerModels data dict
    -matrix_constructor::Function -  a function to construct the jac/sens matrix to be analyzed.
"""
function calc_spectral_analysis(net::Dict{String,<:Any})
    S = calc_voltage_sensitivity_matrix(net) 
    Svp,Svq = S.vp,S.vq
    Svp_spec,Svq_spec = calc_spectral_analysis(Svp),calc_spectral_analysis(Svq)
    return Dict(
        "svp"=>Svp_spec,
        "svq"=>Svq_spec
        )
end

function calc_condition_number(net::Dict)
    S = calc_voltage_sensitivity_matrix(net)
    Svp,Svq = S.vp,S.vq
    κ_vp,κ_vq = con(Svp),con(Svq)
    return [κ_vp,κ_vq]
end

function calc_condition_number(net::Dict,sel_bus_types=[1])
    S = calc_voltage_sensitivity_matrix(net,sel_bus_types)
    Svp,Svq = S.vp,S.vq
    κ_vp,κ_vq = con(Svp),con(Svq)
    return [κ_vp,κ_vq]
end 

function calc_condition_number(file::String)
    net = make_basic_network(parse_file(path))
    return calc_condition_number(net)
end 



# function calc_spectral_analysis(J::PowerFlowJacobian)
#     Jpv,Jqv = J.pv, J.qv
#     Jpv_spec,Jqv_spec = calc_spectral_analysis(Jpv),calc_spectral_analysis(Jqv)
#     return [Jpv_spec,Jqv_spec]
# end




