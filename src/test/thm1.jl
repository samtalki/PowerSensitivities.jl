###
#Algorithms for testing Theorem 1.
###
include("../PowerSensitivities.jl")
include("../util/matrix.jl")
using PowerModels
using LinearAlgebra
import .PowerSensitivities

"""
Given a network data dict, tests theorem 1 conditions
"""
function test_thm1(network::Dict,sel_bus_types=[1],ϵ=1e-6)
    #Compute the indeces that will be considered
    bad_idx = PowerSensitivities.calc_bad_idx(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]
    #Compute the matrices of interest
    K = PowerSensitivities.calc_K_matrix(network)[study_idx,study_idx]
    ∂p∂θ = PowerSensitivities.calc_pth_jacobian(network)[study_idx,study_idx]
    ∂q∂θ = PowerSensitivities.calc_qth_jacobian(network)[study_idx,study_idx]
    #Check the sizes
    @assert size(∂p∂θ,1) == length(study_idx) && size(∂p∂θ,2) == length(study_idx)
    @assert size(∂q∂θ,1) == length(study_idx) && size(∂q∂θ,2) == length(study_idx)
    @assert size(K,1) == length(study_idx) && size(K,2) == length(study_idx)
    #Compute the theorem quantities
    k_max = maximum(diag(K))
    ΔK = k_max*I - K
    M = k_max*∂p∂θ - ∂q∂θ
    return opnorm(inv(M)*ΔK*∂p∂θ)
end


"""
Given a network data dict,
Calculate RHS of inequality "Δk_max" with the option drop_bad_idx.
If drop_bad_idx==True, then the buses that meet the conditions in calc_bad_idx are discarded.
"""
function calc_Δk_bound(network::Dict{String,<:Any},sel_bus_types=[1],ϵ=1e-6)
    #Compute the indeces that will be considered
    bad_idx = PowerSensitivities.calc_bad_idx(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]
    #Compute the matrices of interest
    K = PowerSensitivities.calc_K_matrix(network)[study_idx,study_idx]
    ∂p∂θ = PowerSensitivities.calc_pth_jacobian(network)[study_idx,study_idx]
    ∂q∂θ = PowerSensitivities.calc_qth_jacobian(network)[study_idx,study_idx]
    k_max = maximum(diag(K))
    M = k_max*∂p∂θ - ∂q∂θ
    #Check the sizes
    @assert size(∂p∂θ,1) == length(study_idx) && size(∂p∂θ,2) == length(study_idx)
    @assert size(∂q∂θ,1) == length(study_idx) && size(∂q∂θ,2) == length(study_idx)
    @assert size(K,1) == length(study_idx) && size(K,2) == length(study_idx)
    #Compute the bound on Δk
    Δk_bound = opnorm(inv(M))^(-1)*opnorm(∂p∂θ)^(-1)
    return Δk_bound
end

"""
Given a network data dict,
Calculate LHS of inequality "Δk" with the option drop_bad_idx.
If drop_bad_idx==True, then the buses that meet the conditions in calc_bad_idx are discarded.
"""
function calc_Δk_observed(network::Dict{String,<:Any},sel_bus_types=[1],ϵ=1e-6)
    #Compute the indeces that will be considered
    bad_idx = PowerSensitivities.calc_bad_idx(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]
    #Compute the matrix K and the bus power factors
    K = PowerSensitivities.calc_K_matrix(network)[study_idx,study_idx]
    pf = PowerSensitivities.calc_basic_power_factor(network)[study_idx] #Compute bus power factors
    #Check the sizes are consistent
    @assert size(K,1) == length(study_idx) && size(K,2) == length(study_idx)
    @assert length(pf) == length(study_idx) #Check sizes are consistent
    k_diag = diag(K)
    k_max,k_min = maximum(k_diag),minimum(k_diag)
    pf_max,pf_min = maximum(pf),minimum(pf)
    Δk_obs = k_max-k_min#abs(maximum(k_diag)-minimum(k_diag)) #old
    Δpf_obs = maximum(pf) - minimum(pf) #abs(maximum(pf) - minimum(pf)) #old
    return Dict(
        "Δk_obs" => Δk_obs,
        "Δpf_obs" => Δpf_obs,
        "k_max" => k_max,
        "k_min" => k_min,
        "pf_max" => pf_max,
        "pf_min" => pf_min
    )
end

"""
Given a network data dict,
Compute the maximum difference between nodal power factors so that complex power injections can be modeled from voltage magnitudes
"""
function calc_thm1_data(network::Dict{String,<:Any},sel_bus_types=[1],ϵ=1e-6)
    #Compute the indeces that will be considered
    bad_idx = PowerSensitivities.calc_bad_idx(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]
    #Compute the matrices of interest
    K = PowerSensitivities.calc_K_matrix(network)[study_idx,study_idx]
    pf = PowerSensitivities.calc_basic_power_factor(network)[study_idx]
    ∂p∂θ = PowerSensitivities.calc_pth_jacobian(network)[study_idx,study_idx]
    ∂q∂θ = PowerSensitivities.calc_qth_jacobian(network)[study_idx,study_idx]
    k_max = maximum(diag(K))
    M = k_max*∂p∂θ - ∂q∂θ
    ΔK = k_max*I - K
    #Check the sizes
    @assert size(∂p∂θ,1) == length(study_idx) && size(∂p∂θ,2) == length(study_idx)
    @assert size(∂q∂θ,1) == length(study_idx) && size(∂q∂θ,2) == length(study_idx)
    @assert size(K,1) == length(study_idx) && size(K,2) == length(study_idx)
    return Dict(
        "holds" => opnorm(inv(M)*ΔK*∂p∂θ)<1,
        "M" => M,
        "K" => K,
        "ΔK" => ΔK,
        "k_max" => k_max,
        "spth" => ∂p∂θ,
        "sqth" => ∂q∂θ,
        "pf" => pf,
        "study_idx"=> study_idx)
end
