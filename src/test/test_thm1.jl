include("../PowerSensitivities.jl")
import .PowerSensitivities
#using PowerModels
#using PowerModelsAnalytics
using LinearAlgebra
using JuMP
using Ipopt
using Gadfly
PowerModels.silence()


"""
Checks if a matrix M is positive definite
"""
ispd(M) = all([real(eig)>0 for eig in eigvals(M)])

"""
Checks if a matrix M is negative definite
"""
isnd(M) = all([real(eig)<0 for eig in eigvals(M)])
isnsd(M,ϵ=1e-12) = all([real(eig)<=ϵ for eig in eigvals(M)])

"""
Checks if a matrix M is invertible
"""
isinvertible(x) = applicable(inv, x)

"""
Returns distance between M and M transpose
"""
symmetricdiff(M) = norm(M-transpose(M))

"""
Check if symmetric part of a matrix is negative definite
"""
symmetric_part_nd(M) = isnd(0.5.*(M + transpose(M)))
symmetric_part_nsd(M) = isnsd(M + transpose(M./2))

"""
Check if symmetric part of a matrix is positive definite
"""
symmetric_part_pd(M) = ispd(0.5*(M +transpose(M)))

#Test case path
#network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/"
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/"
not_require_radial = false #whether to require test feeder is radial

"""
Test assumption 1 for all feeder models in folder network_data_path
"""
function test_assumption1(sel_bus_types=[1,2],network_data_path=network_data_path)
	names = readdir(network_data_path);
    paths = readdir(network_data_path,join=true);
    results, network_dicts, full_J_matrices = Dict(),Dict(),Dict()
    for (name,path) in zip(names,paths)
        data = try
            make_basic_network(parse_file(path));
        catch
            println("PM cannot parse "*name)
            continue
        end
        if PowerSensitivities.is_radial(data) || not_require_radial
            J = try
                PowerSensitivities.calc_jacobian_matrix(data,sel_bus_types);
            catch
                println("Jacobian cannot be computed for "*name)
                continue
            end
            Y = Matrix(calc_basic_admittance_matrix(data));
            Jmat = Matrix(J.matrix)
            spth = Matrix(J.pth)
            sqth = Matrix(J.qth)
            info = Dict(
                "n_slack_node" => length(findall(d -> d["bus_type"] ==3, data["bus"])),
                "slack_nodes" => findall(d -> d["bus_type"] ==3, data["bus"]),
                "dpdth_pd" => ispd(spth),
                "y_symmetric" => issymmetric(Y),
                "dpdth_symmetric" => issymmetric(spth),
                "J_nonsingular" => isinvertible(Jmat),
                "sym_dpdth_pd" => symmetric_part_pd(spth),
                "sym_dqdth_nsd" => symmetric_part_nsd(sqth))
            network_dicts[name] = data
            full_J_matrices[name] = calc_basic_jacobian_matrix(data)
            results[name] = info
        else
            continue
        end
    end
    return results,full_J_matrices,network_dicts
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
        if PowerSensitivities.is_radial(network) || not_require_radial
            results[name] = PowerSensitivities.calc_vmag_condition(network,sel_bus_types)
        end
    end
    return results
end

#Test assumption 1
results_all,J_all,network_dicts = test_assumption1([1,2,3]);
results_pq_pv,J_pq_pv,network_dicts = test_assumption1([1,2]);
results_pq,J_pq,network_dicts = test_assumption1([1]);


#Test Theorem 1 bounds
thm1_pq = test_thm1([1]);
thm1_pq_pv = test_thm1([1,2]);

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

p1 = spy(J_c5.pth,Guide.xlabel("Bus Index k"),Guide.ylabel("Bus Index i"),Guide.title("∂pᵢ/∂θₖ"))
p2 = spy(J_c5.qth,Guide.xlabel("Bus Index k"),Guide.ylabel("Bus Index i"),Guide.title("∂qᵢ/∂θₖ"))
title(hstack(p1,p2),"IEEE Case 5 Voltage Phase Angle Jacobians")

dpi=100
plot_network(case3; 
    filename="case3.pdf",
    label_nodes=true,
    label_edge=true,
    node_size_limits=[15,20],
    edge_width_limits=[3, 4],
    fontsize=12,
    plot_size=(floor(Int64,3.5*dpi),floor(Int64,3.5/1.61828*dpi)),
    plot_dpi=dpi);

graph = plot_network(case5;
    node_size_limits=[10, 15],
    edge_width_limits=[2, 3],
    label_nodes=true,
    fontsize=10,
    plot_size=(600,600),
    plot_dpi=100,
    filename="case5.pdf");
