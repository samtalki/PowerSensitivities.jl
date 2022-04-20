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