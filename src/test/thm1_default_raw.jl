###
#Test Theorem 1 at the default operating point reflected in the network data with no preprocessing steps.
###
include("../PowerSensitivities.jl")
include("../util/matrix.jl")
using PowerModels
using LinearAlgebra
using Statistics: mean
import .PowerSensitivities


"""
Given a network data dict, with MINIMAL PREPROCESSING,
calculate the indeces do not satisfy the assumptions of Theorem 1
"""
function calc_bad_idx_minimal(data::Dict{String,<:Any},ϵ=1e-6)
    s = calc_basic_bus_injection(data);
    pnet,qnet,pf = real.(s),imag.(s),PowerSensitivities.calc_basic_power_factor(data) 
    bad_idx = [] #Array of indeces with zero pnet or zero MVA injections to be discarded
    for (i,pf_i) in enumerate(pf)
        #Ignore buses with zero power injection
        if abs(pnet[i]) ≤ ϵ && abs(qnet[i]) ≤ ϵ 
            push!(bad_idx,i)
        #Ignore buses with zero power factor
        elseif abs(pf_i) <= 1e-5 || pf_i == NaN || pnet[i] == 0.0
            push!(bad_idx,i) #K[i,i] = 0
        # #Ignore buses with capacitive injections
        # elseif qnet[i] > 0 
        #    push!(bad_idx,i)
        else
            continue
        end
    end
    return bad_idx
end


"""
Given a network data dict, with MINIMAL PREPROCESSING,
Compute K matrix where K = diag(√(1-pf_i^2)/pf_i)
"""
function calc_K_matrix_minimal(network::Dict{String,<:Any},ϵ=1e-6)
    bad_idx = calc_bad_idx_minimal(network,ϵ) #Array of indeces with zero p or zero MVA injections to be discarded
    s = calc_basic_bus_injection(network);
    pnet,qnet,pf = real.(s),imag.(s),PowerSensitivities.calc_basic_power_factor(network) 
    n = length(pf)
    K = zeros((n,n))
    for (i,(pf_i,qnet_i)) in enumerate(zip(pf,qnet))
        K[i,i] = PowerSensitivities.k(pf_i,qnet_i) #<--- Note that we used the signed version here!
    end
    k_mean = mean(diag(K)[[i for i in 1:n if i ∉ bad_idx]]) #Replace the bad indeces with the mean of the other k entries
    for i in bad_idx
        K[i,i] = k_mean ## Replace k entry at bad idx with the mean of other ks.
    end
    return K
end



### Function to test the basic conditions of Theorem 1
"""
Given a network data dict, tests theorem 1 conditions with automatically chosen study indeces.
"""
function test_thm1(network::Dict,sel_bus_types=[1],ϵ=1e-6)
  
    #RAW/MINIMAL METHOD with  bus type selection and LIMITED preprocesing allowing both inductive/capacitive injections
    bad_idx = calc_bad_idx_minimal(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]

    #Compute the matrices of interest
    K = calc_K_matrix_minimal(network)[study_idx,study_idx]
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

### Function to compute the Δk upper bound
"""
Given a network data dict,
Calculate RHS of inequality "Δk_max" with automatically chosen study indeces.
"""
function calc_Δk_bound(network::Dict{String,<:Any},sel_bus_types=[1],ϵ=1e-6)
   
    #RAW/MINIMAL METHOD with  bus type selection and LIMITED preprocesing allowing both inductive/capacitive injections
    bad_idx = calc_bad_idx_minimal(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]

    #Compute the matrices of interest
    K = calc_K_matrix_minimal(network)[study_idx,study_idx]
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

# Function to compute the Δk observed from the power factor data.
"""
Given a network data dict,
Calculate LHS of inequality "Δk" with automatically chosen study indeces.
"""
function calc_Δk_observed(network::Dict{String,<:Any},sel_bus_types=[1],ϵ=1e-6)

    #RAW/MINIMAL METHOD with  bus type selection and LIMITED preprocesing allowing both inductive/capacitive injections
    bad_idx = calc_bad_idx_minimal(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]

    #Compute the matrix K and the bus power factors
    K = calc_K_matrix_minimal(network)[study_idx,study_idx]
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
# Function to compute the overall summary data for Theorem 1.
### 
"""
Given a network data dict,
Calculate the overall Theorem 1 data with automatically chosen study indeces.
"""
function calc_thm1_data(network::Dict{String,<:Any},sel_bus_types=[1],ϵ=1e-6)
    #Compute the indeces that will be considered
    
    #RAW/MINIMAL METHOD with  bus type selection and LIMITED preprocesing allowing both inductive/capacitive injections
    bad_idx = calc_bad_idx_minimal(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]
    
    
    #Compute the matrices of interest
    K = calc_K_matrix_minimal(network)[study_idx,study_idx]
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


#####################################################################################
#Radial cases
#Test case path and parameters
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
names = readdir(network_data_path);
paths = readdir(network_data_path,join=true);

#Results dicts for automatic indexing.
#PQ bus only, unsuitable indexes removed.
radial_thm1_results = Dict()
radial_thm1_results_pq_pv = Dict()
radial_thm1_weak_rhs = Dict()
radial_thm1_weak_lhs = Dict()
radial_thm1_data = Dict()
radial_network_dicts = Dict() #Save PowerModels networks
for (name,path) in zip(names,paths)
    network =  try
        make_basic_network(parse_file(path));
    catch
        println("PM cannot parse "*name)
        continue
    end
    if PowerSensitivities.is_radial(network)
        try  ### Compute the AC power flow solution first?
            #compute_ac_pf!(network)  
        catch
            println("AC Power Flow solution failed for: ",name)
            continue
        end
        try
            radial_thm1_results[name] = test_thm1(network)
            radial_thm1_results_pq_pv[name] = test_thm1(network,[1,2])
            radial_thm1_weak_rhs[name] = calc_Δk_bound(network)
            radial_thm1_weak_lhs[name] = calc_Δk_observed(network)
            radial_thm1_data[name] = calc_thm1_data(network)
            radial_network_dicts[name] = network
        catch
           radial_thm1_results[name] = nothing
           radial_thm1_weak_rhs[name] = nothing
           radial_thm1_weak_lhs[name] = nothing
           radial_thm1_data[name] = nothing
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
mesh_network_dicts = Dict() #Save PowerModels networks
for (name,path) in zip(names,paths)
    network =  try
        make_basic_network(parse_file(path));
    catch
        println("PM cannot parse "*name)
        continue
    end
    if !PowerSensitivities.is_radial(network) 
        try ### Compute the AC power flow solution first?
            #compute_ac_pf!(network) 
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
            mesh_network_dicts[name] = network
        catch
            continue
        end
        
    end
end


#Check where theorem 1 holds
radial_thm1_holds = Dict()
for (name,data) in radial_thm1_data
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

radial_n_cases_holds = sum([1 for (name,holds) in radial_thm1_holds if holds])
mesh_n_cases_holds = sum([1 for (name,holds) in mesh_thm1_holds if holds])