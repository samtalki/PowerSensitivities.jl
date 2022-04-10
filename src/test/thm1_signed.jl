#Test theorem 1 while allowing non-inducitve injections
include("../PowerSensitivities.jl")
include("../util/matrix.jl")
using PowerModels
using LinearAlgebra
using Statistics
import .PowerSensitivities

calc_basic_power_factor(network::Dict{String,<:Any}) =  [cos(θi) for θi in angle.(calc_basic_bus_injection(network))]

"""
Given a network data dict, calculate the "bad indeces" that do not satisfy the assumptions of Theorem 1
"""
function calc_bad_idx(network::Dict{String,<:Any},ϵ=1e-6)
    s = calc_basic_bus_injection(network);
    p,q,pf = real.(s),imag.(s),calc_basic_power_factor(network) 
    bad_idx = [] #Array of indeces with zero p or zero MVA injections to be discarded
    for (i,pf_i) in enumerate(pf)
        #Ignore buses with zero power injection
        if abs(p[i]) ≤ ϵ && abs(q[i]) ≤ ϵ 
            push!(bad_idx,i)
        #Ignore buses with zero power factor
        elseif abs(pf_i) <= 1e-5 || pf_i == NaN || p[i] == 0.0
            push!(bad_idx,i) #K[i,i] = 0
        # #Ignore buses with capacitive injections
        # elseif q[i] > 0 
        #    push!(bad_idx,i)
        else
            continue
        end
    end
    return bad_idx
end

"""
Given a real power factor pf return the implicit function theorem representation of pf
"""
function k(pf::Real)
    #@assert pf<= 1 && pf > 0 "Power factors must be between (0,1]"
    return sqrt(1-pf^2)/pf
end

"""
Given a real power factor pf and the reactive power injection return the SIGNED implicit function theorem representation of pf
This is used to make entries of the K matrix.
"""
function k(pf::Real,q::Real)
    #@assert pf<= 1 && pf > 0 "Power factors must be between (0,1]"
    return sign(q)*(sqrt(1-pf^2)/pf)
end

"""
Compute K matrix where K = diag(√(1-pf_i^2)/pf_i)
"""
function calc_K_matrix(network::Dict{String,<:Any},ϵ=1e-6)
    bad_idx = calc_bad_idx(network) #Array of indeces with zero p or zero MVA injections to be discarded
    pf = calc_basic_power_factor(network) 
    n = length(pf)
    K = zeros((n,n))
    q = imag.(calc_basic_bus_injection(network))
    for (i,pf_i) in enumerate(pf)
        q_i = q[i]
        K[i,i] = k(pf_i,q_i) #k(pf_i,q[i])#sqrt(1-pf_i^2)/pf_i
    end
    k_mean = mean(diag(K)[[i for i in 1:n if i ∉ bad_idx]]) #Replace the bad indeces with the mean of the other k entries
    for i in bad_idx
        K[i,i] = k_mean ## Replace k entry at bad idx with the mean of other ks.
    end
    return K
end


"""
Given a network data dict, tests theorem 1 conditions
"""
function test_thm1(network::Dict,sel_bus_types=[1,2,3],ϵ=1e-6)
    #Compute the indeces that will be considered
    bad_idx = calc_bad_idx_no_q_constraints(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]
    #Compute the matrices of interest
    K = calc_K_matrix(network)[study_idx,study_idx]
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
If drop_bad_idx==True, then the buses that meet the conditions in calc_bad_idx_no_q_constraints are discarded.
"""
function calc_Δk_bound(network::Dict{String,<:Any},sel_bus_types=[1,2,3],ϵ=1e-6)
    #Compute the indeces that will be considered
    bad_idx = calc_bad_idx_no_q_constraints(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]
    #Compute the matrices of interest
    K = calc_K_matrix(network)[study_idx,study_idx]
    ∂p∂θ = PowerSensitivities.calc_pth_jacobian(network)[study_idx,study_idx]
    ∂q∂θ = PowerSensitivities.calc_qth_jacobian(network)[study_idx,study_idx]
    k_max = maximum(diag(K))
    M = k_max*∂p∂θ - ∂q∂θ
    #Check the sizes
    @assert size(∂p∂θ,1) == length(study_idx) && size(∂p∂θ,2) == length(study_idx)
    @assert size(∂q∂θ,1) == length(study_idx) && size(∂q∂θ,2) == length(study_idx)
    @assert size(K,1) == length(study_idx) && size(K,2) == length(study_idx)
    #Compute the bound on Δk
    Δk_bound = try
        opnorm(inv(M))^(-1)*opnorm(∂p∂θ)^(-1)
    catch
        Δk_bound = nothing
    end
    return Δk_bound
end

"""
Given a network data dict,
Calculate LHS of inequality "Δk" with the option drop_bad_idx.
If drop_bad_idx==True, then the buses that meet the conditions in calc_bad_idx_no_q_constraints are discarded.
"""
function calc_Δk_observed(network::Dict{String,<:Any},sel_bus_types=[1,2,3],ϵ=1e-6)
    #Compute the indeces that will be considered
    bad_idx = calc_bad_idx_no_q_constraints(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]
    #Compute the matrix K and the bus power factors
    K = calc_K_matrix(network)[study_idx,study_idx]
    pf = calc_basic_power_factor(network)[study_idx]
    #Check the sizes are consistent
    @assert size(K,1) == length(study_idx) && size(K,2) == length(study_idx)
    @assert length(pf) == length(study_idx) #Check sizes are consistent
    Δpf_obs = try
        abs(maximum(pf) - minimum(pf))
    catch
        Δpf_obs = nothing
    end
    Δk_obs = try 
        k_diag = diag(K)
        abs(maximum(k_diag)-minimum(k_diag))
    catch
        Δk_obs = nothing
    end
    return Dict(
        "Δk_obs" => Δk_obs,
        "Δpf_obs" => Δpf_obs
    )
end


"""
Given a network data dict,
Compute the maximum difference between nodal power factors so that complex power injections can be modeled from voltage magnitudes
"""
function calc_thm1_data(network::Dict{String,<:Any},sel_bus_types=[1,2,3],ϵ=1e-6)
    #Compute the indeces that will be considered
    bad_idx = calc_bad_idx_no_q_constraints(network,ϵ)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]
    #Compute the matrices of interest
    K = calc_K_matrix(network)[study_idx,study_idx]
    pf = calc_basic_power_factor(network)[study_idx]
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


begin
    #Test case path and parameters
    #network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/" #Folder with meshed and radial systems
    network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
    allow_mesh = false #Whether to allow meshed test cases/require test feeder is radial

    names = readdir(network_data_path,sort=false);
    paths = readdir(network_data_path,join=true,sort=false);

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
        if PowerSensitivities.is_radial(network) || allow_mesh
            thm1_results[name] = test_thm1(network)
            thm1_weak_rhs[name] = calc_Δk_bound(network)
            thm1_weak_lhs[name] = calc_Δk_observed(network)
            thm1_data[name] = calc_thm1_data(network)
        end
    end
    thm1_holds = Dict()
    for (name,data) in thm1_data
        thm1_holds[name] = data["holds"]
    end
    n_cases_holds = sum([1 for (name,holds) in thm1_holds if holds])
end


#Non-radial cases
begin
    #Results dicts
    mesh_thm1_results = Dict()
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
            mesh_thm1_results[name] = test_thm1(network)
            mesh_thm1_weak_rhs[name] = calc_Δk_bound(network)
            mesh_thm1_weak_lhs[name] = calc_Δk_observed(network)
            mesh_thm1_data[name] = calc_thm1_data(network)
        end
    end
end
