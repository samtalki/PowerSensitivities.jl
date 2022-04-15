###
#Test Theorem 1 at the AC power flow solution.
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
    Δk_bound =  opnorm(inv(M))^(-1)*opnorm(∂p∂θ)^(-1)
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
    Δk_obs = k_max-k_min#abs(maximum(k_diag)-minimum(k_diag))
    Δpf_obs = maximum(pf) - minimum(pf)
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






#####################################################################################
#Radial cases
#Test case path and parameters
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
names = readdir(network_data_path);
paths = readdir(network_data_path,join=true);

#Results dicts
thm1_results = Dict()
thm1_weak_rhs = Dict()
thm1_weak_lhs = Dict()
thm1_data = Dict()
for (name,path) in zip(names,paths)
    network =  try
        make_basic_network(parse_file(path));
    catch
        println("PM cannot parse "*name)
        continue
    end
    if PowerSensitivities.is_radial(network)
        try  ### Compute the AC power flow solution first!
            compute_ac_pf!(network)  
        catch
            println("AC Power Flow solution failed for: ",name)
            continue
        end
        try
            thm1_results[name] = test_thm1(network)
            thm1_weak_rhs[name] = calc_Δk_bound(network)
            thm1_weak_lhs[name] = calc_Δk_observed(network)
            thm1_data[name] = calc_thm1_data(network)
        catch
            thm1_results[name] = nothing
            thm1_weak_rhs[name] = nothing
            thm1_weak_lhs[name] = nothing
            thm1_data[name] = nothing
        end
    end
end

#Test for non-radial cases
####
allow_mesh = true #Whether to allow meshed test cases/require test feeder is radial
#Results dicts
mesh_thm1_results = Dict()
mesh_thm1_results_pq_pv = Dict()
mesh_thm1_weak_rhs = Dict()
mesh_thm1_weak_lhs = Dict()
mesh_thm1_data = Dict()
for (name,path) in zip(names,paths)
    network =  try
        make_basic_network(parse_file(path));
    catch
        println("PM cannot parse "*name)
        continue
    end
    if !PowerSensitivities.is_radial(network) 
        try ### Compute the AC power flow solution first!
            compute_ac_pf!(network) 
        catch
            println("AC Power Flow solution failed for: ",name)
            continue
        end
        try
            mesh_thm1_results[name] = test_thm1(network)
            mesh_thm1_results_pq_pv[name] = test_thm1(network,[1,2])
            mesh_thm1_weak_rhs[name] = calc_Δk_bound(network)
            mesh_thm1_weak_lhs[name] = calc_Δk_observed(network)
            mesh_thm1_data[name] = calc_thm1_data(network)
        catch
            continue
        end
        
    end
end


#Check where theorem 1 holds
radial_thm1_holds = Dict()
for (name,data) in thm1_data
    radial_thm1_holds[name] = try 
        data["holds"]
    catch
        false
    end
end

mesh_thm1_holds = Dict()
for (name,data) in mesh_thm1_data
    mesh_thm1_holds[name] = try 
        data["holds"]
    catch
        false
    end
end

n_radial_cases_holds = sum([1 for (name,holds) in radial_thm1_holds if holds])
n_mesh_cases_holds = sum([1 for (name,holds) in mesh_thm1_holds if holds])