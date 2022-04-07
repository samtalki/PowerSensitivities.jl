include("../PowerSensitivities.jl")
include("../util/matrix.jl")
using PowerModels: parse_file,make_basic_network,calc_basic_jacobian_matrix,calc_basic_bus_injection,calc_basic_bus_voltage
import .PowerSensitivities
using LinearAlgebra
using JuMP
using Ipopt
using Gadfly

#Test case path and parameters
#network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/" #Folder with meshed and radial systems
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
allow_mesh = false #Whether to allow meshed test cases/require test feeder is radial


### Calculation functions

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
    n_selected = length(good_idx)
    ∂p∂θ = PowerSensitivities.calc_pth_jacobian(network,sel_bus_types)[good_idx,good_idx]
    M = PowerSensitivities.calc_M_matrix(network,sel_bus_types)[good_idx,good_idx]
    @assert size(M,1) == n_selected && size(M,2) == n_selected && size(∂p∂θ,1) == n_selected &&  size(∂p∂θ,2) == n_selected #Check the lengths are consistent
    Δk_max = try
        opnorm(inv(M))^(-1)*opnorm(∂p∂θ)^(-1)
        #(1/(opnorm(inv(Matrix(M)),2)))*(1/(opnorm(Matrix(∂p∂θ),2)))
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
    k(pf::Real) = sqrt(1-pf^2)/pf
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    n = length(idx_sel_bus_types)
    if drop_bad_idx
        bad_idx = PowerSensitivities.calc_bad_idx(network,sel_bus_types)
        good_idx = [i for i in 1:n if i ∉ bad_idx]
    else
        good_idx = [i for i in 1:n]
    end
    n_selected = length(good_idx)
    pf = PowerSensitivities.calc_basic_power_factor(network,sel_bus_types)[good_idx]
    @assert length(pf) == n_selected #Check sizes are consistent
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
Given:
    a network data dict under study,
    a chosen maximum power factor pf_max,
    and the SELECTED BUS TYPES under study
Compute:
    the minimum power factor pf_min such that the complex power injections are observable
"""
function calc_pf_min(network::Dict{String,<:Any},pf_max::Real=1,sel_bus_types=[1,2])
    bad_idx = PowerSensitivities.calc_bad_idx(network,sel_bus_types)
    M = PowerSensitivities.calc_M_matrix(network,sel_bus_types)[bad_idx,bad]
    ∂p∂θ = PowerSensitivities.calc_pth_jacobian(network,sel_bus_types)
    if pf_max==1
        pf_min = PowerSensitivities.kinv(
            opnorm(inv(M))^(-1)*opnorm(∂p∂θ)^(-1)
        )
    else
        pf_min = PowerSensitivities. kinv(
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


### Testing functions

"""
Given chosen bus types under study and folder of test cases under study,
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
Given chosen bus types under study and folder of test cases under study,
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
Given chosen bus types under study and folder of test cases under study,
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
                throw(Exception)
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
                continue
            end
        end
    end
    return results
end


"""
Given a lhs (Δk=k_max - k_min) and a RHS (1/opnorm(inv(M)))(1/opnorm(∂p∂θ)),
Check if Theorem 1 from Talkington, Turizo, et al. is satisfied.
    TODO: Add a second implementation that takes a network data dict.
"""
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
#Check if thm1 holds 
thm1_hold_pq = thm1_satisfied(lhs_pq,rhs_pq)
thm1_hold_pq_pv = thm1_satisfied(lhs_pq_pv,rhs_pq_pv)


#Compute the minimum bus power factors
#Actual kmax of operating point
pf_min_pq = test_pf_min([1])
pf_min_pq_pv = test_pf_min([1,2])
#Maximum power factor of unity
pf_min_unity_pq = test_pf_min([1],1)
pf_min_unity_pq_pv = test_pf_min([1,2],1)


# #Test for the maximum power factor distances
# delta_pf_max_pq,delta_pf_max_pq_pv = Dict(),Dict()
# for (name,d) in thm1_pq
#     delta_pf_max_pq[name] = d["Δpf_max"]
# end
# for (name,d) in thm1_pq_pv
#     delta_pf_max_pq_pv[name] = d["Δpf_max"]
# end

# Auxillary data for testing
case5 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case5.m"));
case14 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case14.m"));
case18 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/radial_test/case18.m"));
J_c5 = PowerSensitivities.calc_jacobian_matrix(case5);
J_c24_full = calc_basic_jacobian_matrix(case24);
J_c24_blocks = PowerSensitivities.calc_jacobian_matrix(case24);
J_c24_blocks_pq_pv = PowerSensitivities.calc_jacobian_matrix(case24,[1,2]);
J_c24_blocks_pq = PowerSensitivities.calc_jacobian_matrix(case24,[1]);





# Plotting functions

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
