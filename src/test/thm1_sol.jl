###
#Test Theorem 1 at the AC power flow solution.
###
include("../PowerSensitivities.jl")
include("../test/thm1.jl")
include("../util/matrix.jl")
using PowerModels
using LinearAlgebra
import .PowerSensitivities


#####################################################################################
#Radial cases
#Test case path and parameters
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
names = readdir(network_data_path);
paths = readdir(network_data_path,join=true);

#Results dicts for automatic indexing.
#PQ bus only, unsuitable indexes removed.
radial_thm1_results = Dict()
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
        try  ### Compute the AC power flow solution first!
            compute_ac_pf!(network)  
        catch
            println("AC Power Flow solution failed for: ",name)
            continue
        end
        try
            radial_thm1_results[name] = test_thm1(network)
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

#Mesh cases
#Test case path and parameters
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/pm_matpower/" #Folder with radial-only systems
names = readdir(network_data_path);
paths = readdir(network_data_path,join=true);
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