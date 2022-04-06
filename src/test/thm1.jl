include("../PowerSensitivities.jl")
include("../util/matrix.jl")
using PowerModels: parse_file,make_basic_network,calc_basic_jacobian_matrix,calc_basic_bus_injection,calc_basic_bus_voltage
import .PowerSensitivities
using LinearAlgebra
using JuMP
using Ipopt
using Gadfly

#Test case path
#network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/" #Folder with meshed and radial systems
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
allow_mesh = false #Whether to allow meshed test cases/require test feeder is radial

"""
Given a network data dict,
Calculate RHS of inequality "Δk_max" with the option drop_bad_idx.
If drop_bad_idx==True, then the buses that meet the conditions in calc_bad_idx are discarded.
"""
function calc_rhs(network::Dict{String,<:Any},sel_bus_types=[1,2],drop_bad_idx=true)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    n = length(idx_sel_bus_types)
    if drop_bad_idx
        bad_idx = PowerSensitivities.calc_bad_idx(network,sel_bus_types)
        good_idx = [i for i in 1:n if i ∉ bad_idx]
    else
        good_idx = [i for i in 1:n]
    end
    ∂p∂θ = PowerSensitivities.calc_pth_jacobian(network,sel_bus_types)[good_idx,good_idx]
    M = PowerSensitivities.calc_M_matrix(network,sel_bus_types)[good_idx,good_idx]
    Δk_max = try
        (1/(opnorm(inv(Matrix(M)),2)))*(1/(opnorm(Matrix(∂p∂θ),2)))
    catch
        Δk_max = nothing
    end
    return Δk_max
end

"""
Given a network data dict,
Calculate LHS of inequality "Δk" with the option drop_bad_idx.
If drop_bad_idx==True, then the buses that meet the conditions in calc_bad_idx are discarded.
"""
function calc_lhs(network::Dict{String,<:Any},sel_bus_types=[1,2],drop_bad_idx=true)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    n = length(idx_sel_bus_types)
    if drop_bad_idx
        bad_idx = PowerSensitivities.calc_bad_idx(network,sel_bus_types)
        good_idx = [i for i in 1:n if i ∉ bad_idx]
    else
        good_idx = [i for i in 1:n]
    end
    k(pf::Real) = sqrt(1-pf^2)/pf
    pf = PowerSensitivities.calc_basic_power_factor(network,sel_bus_types)[good_idx]
    Δpf = try
        abs(maximum(pf) - minimum(pf))
    catch
        Δpf = nothing
    end
    Δk = try 
        abs(k(minimum(abs.(pf))) - k(maximum(abs.(pf))))
    catch
        Δk = nothing
    end
    return Dict(
        "Δk" => Δk,
        "Δpf" => Δpf
    )
end

"""
Test LHS for all feeder models in network_data_path
"""
function test_rhs(sel_bus_types=[1,2],network_data_path=network_data_path)
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
                calc_rhs(network,sel_bus_types)
            catch
                continue
            end
        end
    end
    return results
end

"""
Test LHS for all feeder models in network_data_path
"""
function test_lhs(sel_bus_types=[1,2],network_data_path=network_data_path)
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
                calc_lhs(network,sel_bus_types)
            catch
                continue
            end
        end
    end
    return results
end

"""
Test theorem 1 for all feeder models in network_data_path
"""
function test_thm1(sel_bus_types=[1,2],network_data_path=network_data_path)
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
                PowerSensitivities.calc_vmag_condition(network,sel_bus_types)
            catch
                continue
            end
        end
    end
    return results
end

function thm1_satisfied(lhs,rhs)
    res = Dict()
    for (n, lhs_n) in lhs
        rhs_n = rhs[n]
        res[n] = lhs_n["Δk"]<rhs_n
    end
    return res
end


#Test Theorem 1 bounds
thm1_pq = test_thm1([1]);
thm1_pq_pv = test_thm1([1,2]);

#Thm1 inequality LHS/RHS
lhs_pq = test_lhs([1]);
rhs_pq = test_rhs([1]);
lhs_pq_pv = test_lhs([1,2]);
rhs_pq_pv = test_rhs([1,2]);
thm1_hold_pq = thm1_satisfied(lhs_pq,rhs_pq)
thm1_hold_pq_pv = thm1_satisfied(lhs_pq_pv,rhs_pq_pv)

#Test for the maximum power factor distances
delta_pf_max_pq,delta_pf_max_pq_pv = Dict(),Dict()
for (name,d) in thm1_pq
    delta_pf_max_pq[name] = d["Δpf_max"]
end
for (name,d) in thm1_pq_pv
    delta_pf_max_pq_pv[name] = d["Δpf_max"]
end

# Auxillary data for testing
case3 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case3.m"))
case5 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case5.m"))
case24 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case24.m"))
case14 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case14.m"))
case18 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/radial_test/case18.m"))
J_c3 = PowerSensitivities.calc_jacobian_matrix(case3)
J_c5 = PowerSensitivities.calc_jacobian_matrix(case5)
J_c24_full = calc_basic_jacobian_matrix(case24)
J_c24_blocks = PowerSensitivities.calc_jacobian_matrix(case24)
J_c24_blocks_pq_pv = PowerSensitivities.calc_jacobian_matrix(case24,[1,2])
J_c24_blocks_pq = PowerSensitivities.calc_jacobian_matrix(case24,[1])





#####Plotting functions

# p1 = spy(J_c5.pth,Guide.xlabel("Bus Index k"),Guide.ylabel("Bus Index i"),Guide.title("∂pᵢ/∂θₖ"))
# p2 = spy(J_c5.qth,Guide.xlabel("Bus Index k"),Guide.ylabel("Bus Index i"),Guide.title("∂qᵢ/∂θₖ"))
# title(hstack(p1,p2),"IEEE Case 5 Voltage Phase Angle Jacobians")

# dpi=100
# plot_network(case3; 
#     filename="case3.pdf",
#     label_nodes=true,
#     label_edge=true,
#     node_size_limits=[15,20],
#     edge_width_limits=[3, 4],
#     fontsize=12,
#     plot_size=(floor(Int64,3.5*dpi),floor(Int64,3.5/1.61828*dpi)),
#     plot_dpi=dpi);

# graph = plot_network(case5;
#     node_size_limits=[10, 15],
#     edge_width_limits=[2, 3],
#     label_nodes=true,
#     fontsize=10,
#     plot_size=(600,600),
#     plot_dpi=100,
#     filename="case5.pdf");
