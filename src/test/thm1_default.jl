###
#Test Theorem 1 at the default operating point reflected in the network data.
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
#network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/" #Folder with meshed and radial systems
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
allow_mesh = false #Whether to allow meshed test cases/require test feeder is radial
names = readdir(network_data_path);
paths = readdir(network_data_path,join=true);
begin    
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
mesh_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/" #Folder with meshed and radial systems
names = readdir(mesh_data_path);
paths = readdir(mesh_data_path,join=true);
begin
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
end
