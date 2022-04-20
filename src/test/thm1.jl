###
#Algorithms for testing Theorem 1.
###
include("../PowerSensitivities.jl")
include("../util/matrix.jl")
using PowerModels
using LinearAlgebra
import .PowerSensitivities


### Functions to test the basic conditions of Theorem 1

"""
Given a network data dict, tests theorem 1 conditions with automatically chosen study indeces.
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
Given a network data dict and manually chosen study indeces, tests theorem 1 conditions
"""
function test_thm1_study_idx(network::Dict,study_idx::AbstractVector)
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

### Functions to compute the Δk upper bound

"""
Given a network data dict,
Calculate RHS of inequality "Δk_max" with automatically chosen study indeces.
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
Given a network data dict, and manually chosen study indeces,
Calculate RHS of inequality "Δk_max"
"""
function calc_Δk_bound_study_idx(network::Dict{String,<:Any},study_idx::AbstractVector)
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

###
# Functions to compute the Δk observed from the power factor data.
### 


"""
Given a network data dict,
Calculate LHS of inequality "Δk" with automatically chosen study indeces.
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
Given a network data dict, and manually chosen study indeces,
Calculate LHS of inequality "Δk" 
"""
function calc_Δk_observed_study_idx(network::Dict{String,<:Any},study_idx::AbstractVector)
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


###
# Functions to compute the overall summary data for Theorem 1.
### 

"""
Given a network data dict,
Calculate the overall Theorem 1 data with automatically chosen study indeces.
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
    #Compute the net bus injection vector
    bus_injection = calc_basic_bus_injection(network)[study_idx]
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
        "bus_injection"=>bus_injection,
        "study_idx"=> study_idx)
end

"""
Given a network data dict, and manually chosen study indeces,
Calculate the overall Theorem 1 data.
"""
function calc_thm1_data_study_idx(network::Dict{String,<:Any},study_idx::AbstractVector)
    #Compute the matrices of interest
    K = PowerSensitivities.calc_K_matrix(network)[study_idx,study_idx]
    pf = PowerSensitivities.calc_basic_power_factor(network)[study_idx]
    ∂p∂θ = PowerSensitivities.calc_pth_jacobian(network)[study_idx,study_idx]
    ∂q∂θ = PowerSensitivities.calc_qth_jacobian(network)[study_idx,study_idx]
    k_max = maximum(diag(K))
    M = k_max*∂p∂θ - ∂q∂θ
    ΔK = k_max*I - K
    #Compute the net bus injection vector
    bus_injection = calc_basic_bus_injection(network)[study_idx]
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
        "bus_injection"=>bus_injection,
        "study_idx"=> study_idx)
end

