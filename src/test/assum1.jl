include("../util/matrix.jl")
using PowerModels: parse_file,make_basic_network,calc_basic_jacobian_matrix
#Test case path
#network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/"
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/"
allow_mesh = false #Whether to allow meshed test cases/require test feeder is radial

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
        if PowerSensitivities.is_radial(data) || allow_mesh
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


#Test assumption 1
assum1_all,J_all,network_dicts = test_assumption1([1,2,3]);
assum1_pq_pv,J_pq_pv,network_dicts = test_assumption1([1,2]);
assum1_pq,J_pq,network_dicts = test_assumption1([1]);
