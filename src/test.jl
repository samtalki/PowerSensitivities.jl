include("PowerSensitivities.jl")
using .PowerSensitivities
using PowerModels
using LinearAlgebra


"""
Checks if a matrix M is positive definite
"""
ispd(M) = all([real(eig)>0 for eig in eigvals(M)])

"""
Checks if a matrix M is negative definite
"""
isnd(M) = all([real(eig)<0 for eig in eigvals(M)])

"""
Checks if a matrix M is invertible
"""
isinvertible(x,ϵ=1e-12) = applicable(inv, x) && norm(I(size(x)[1]) - inv(x)*x) < ϵ


"""
Returns distance between M and M transpose
"""
symmetricdiff(M) = norm(M-transpose(M))


"""
Test assumption 1
"""
function test_assumption1(sel_bus_types=[1,2],network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/")
	names = readdir(network_data_path)
    paths = readdir(network_data_path,join=true)
    results = Dict()
	full_J_matrix = Dict()
    for (name,path) in zip(names,paths)
        data = make_basic_network(parse_file(path))
		J = try
            calc_jacobian_matrix(data,sel_bus_types)
        catch
            println("Jacobian cannot be computed for "*name)
            continue
        end
		Y = Matrix(calc_basic_admittance_matrix(data))
		Jmat = Matrix(J.matrix)
		spth = Matrix(J.spth)
        
		info = Dict(
			"n_slack_node" => length(findall(d -> d["bus_type"] ==3, data["bus"])),
			"slack_nodes" => findall(d -> d["bus_type"] ==3, data["bus"]),
			"norm(spth-spth')" => symmetricdiff(spth),
			"rel(spth-spth')" => symmetricdiff(spth)/norm(spth),
            "dpdth_pd" => ispd(spth),
			"y_symmetric" => issymmetric(Y),
			"dpdth_symmetric" => issymmetric(spth),
            "J_nonsingular" => isinvertible(Jmat)
            #"dpdth_norm" => norm(J.spth),
            #"dqdth_norm" => norm(J.sqth)
        )
		full_J_matrix[name] = calc_basic_jacobian_matrix(data)
        results[name] = info
    end
    return results,full_J_matrix
end



case24 = "/home/sam/github/PowerSensitivities.jl/data/matpower/case24.m"
case24 = make_basic_network(parse_file(case24))
J_case_24_baseline = calc_basic_jacobian_matrix(case24)

results_all,J_all = test_assumption1([1,2,3])
results_pq_pv,J_pq_pq = test_assumption1([1,2])
results_pq,J_pq = test_assumption1([1])


# #Load network data
# data = parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case14.m")
# data = make_basic_network(data)
# Y = calc_admittance_matrix(data)
# J = calc_basic_jacobian_matrix(data)
# vm,q = abs.(calc_basic_bus_voltage(data)),imag(calc_basic_bus_injection(data))
# println(length(vm),length(q))


function make_timeseries(data,N=100)
	outputs = ["p_gen", "q_gen","vm_gen","vm_bus","va_bus"]
	train_data = create_samples(data,N,output_vars=outputs)
	return train_data
end


function make_ami_dataset(data,N=100,bus_type="1")
	timeseries = make_timeseries(data,N)
    num_bus = length(data["bus"])
	load_bus = [data["load"]["$i"]["load_bus"] for i in 1:length(data["load"])]
	gen_bus = [data["gen"][string(i)]["gen_bus"] for i in 1:length(data["gen"])]
	idx_sel_bus_types = get_idx_bus_types(data,bus_type)
	num_type = length(idx_sel_bus_types)
	P,Q,V = zeros((N,num_bus)),zeros((N,num_bus)),zeros((N,num_bus))
	for (bus_idx,load_bus_idx,gen_bus_idx) in zip(1:num_bus,1:num_load)
		P[:,k] = dataset["inputs"]["pd"][:,k] - dataset["outputs"]["p_gen"][:,bus_idx]
		#Q[:,k] = 
		#V[:,k] = 
	end
	return ami_dataset
end