## Calculate the minimum power factor
include("../PowerSensitivities.jl")
include("../util/matrix.jl")
using PowerModels: parse_file,make_basic_network,calc_basic_jacobian_matrix,calc_basic_bus_injection,calc_basic_bus_voltage
import .PowerSensitivities

#Test case path and parameters
#network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/" #Folder with meshed and radial systems
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
allow_mesh = false #Whether to allow meshed test cases/require test feeder is radial

"""
Given:
    a network data dict under study,
    a chosen maximum power factor pf_max,
    and the SELECTED BUS TYPES under study
Compute:
    the minimum power factor pf_min such that the complex power injections are observable
"""
function calc_pf_min(network::Dict{String,<:Any},sel_bus_types=[1,2],pf_max::Real=1)
    bad_idx = PowerSensitivities.calc_bad_idx(network,sel_bus_types)
    ∂p∂θ = PowerSensitivities.calc_pth_jacobian(network,sel_bus_types)[bad_idx,bad_idx]
    ∂q∂θ = PowerSensitivities.calc_qth_jacobian(network,sel_bus_types)[bad_idx,bad_idx]
    K = PowerSensitivities.calc_K_matrix(network,sel_bus_types)[bad_idx,bad_idx]
    k_max = maximum(diag(K))
    M = PowerSensitivities.calc_M_matrix(k_max,∂p∂θ,∂q∂θ)
    if pf_max==1
        pf_min = PowerSensitivities.kinv(
            opnorm(inv(M))^(-1)*opnorm(∂p∂θ)^(-1)
        )
    else
        pf_min = PowerSensitivities.kinv(
            PowerSensitivities.k(pf_max) + opnorm(inv(M))^(-1)*opnorm(∂p∂θ)^(-1)
        )
    end
    return pf_min
end

"""
Given:
    a network data dict under study,
    and the selected bus types under study
Compute:
    the maximum value of the K matrix k_max, which is used to compute:
    the minimum power factor pf_min such that the complex power injections are observable
"""
function calc_pf_min(network::Dict{String,<:Any},sel_bus_types=[1,2])
    M = PowerSensitivities.calc_M_matrix(network,sel_bus_types) #Compute the M matrix
    ∂p∂θ = PowerSensitivities.calc_pth_jacobian(network,sel_bus_types) #Compute the power-angle jacobian
    k_max = PowerSensitivities.calc_k_max(network,sel_bus_types) #Compute the actual observed k_max
    pf_min = PowerSensitivities.kinv(
        k_max + opnorm(inv(M))^(-1)*opnorm(∂p∂θ)^(-1)
    )
    return pf_min
end

"""
Given chosen bus types under study and folder of test cases under study,
Compute the minimum power factor implied by Theorem 1 of Talkington and Turizo et al. for the feeder models in network_data_path
"""
function test_pf_min(sel_bus_types=[1,2],network_data_path=network_data_path)
    names = readdir(network_data_path);
    paths = readdir(network_data_path,join=true);
    results = Dict()
    for (name,path) in zip(names,paths)
        network =  try
            make_basic_network(parse_file(path));
        catch
            println("PM cannot parse "*name)
            continue
        end
        if PowerSensitivities.is_radial(network) || allow_mesh
            results[name] = try
                calc_pf_min(network,sel_bus_types)
            catch
                results[name] = nothing
            end
        end
    end
    return results
end

"""
Given chosen bus types under study, a maximum power factor, and folder of test cases under study,
Compute the minimum power factor implied by Theorem 1 of Talkington and Turizo et al. for the feeder models in network_data_path
"""
function test_pf_min(sel_bus_types=[1,2],pf_max::Real=1,network_data_path=network_data_path)
    names = readdir(network_data_path);
    paths = readdir(network_data_path,join=true);
    results = Dict()
    for (name,path) in zip(names,paths)
        network =  try
            make_basic_network(parse_file(path));
        catch
            println("PM cannot parse "*name)
            continue
        end
        if PowerSensitivities.is_radial(network) || allow_mesh
            results[name] = try
                calc_pf_min(network,sel_bus_types,pf_max)
            catch
                results[name] = nothing
            end
        end
    end
    return results
end

#Compute the minimum bus power factors
#Actual kmax of operating point
pf_min_pq = test_pf_min([1])
pf_min_pq_pv = test_pf_min([1,2])
#Maximum power factor of unity
pf_min_unity_pq = test_pf_min([1],1.)
pf_min_unity_pq_pv = test_pf_min([1,2],1.)
