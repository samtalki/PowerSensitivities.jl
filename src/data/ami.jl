"""
Given a network data dict, find the bus indeces that match `sel_bus_types`.
"""
# function get_idx_bus_types(data,sel_bus_types=1)
# 	num_bus = length(data["bus"])
# 	bus_types = [data["bus"][string(i)]["bus_type"] for i in 1:num_bus]
#     idx_sel_bus_types = findall(bus_idx-> bus_idx ∈ sel_bus_types,bus_types)
#     return idx_sel_bus_types
# end

function make_timeseries_dataset(network,N=100)
	#outputs = ["p_gen", "q_gen","vm_gen","vm_bus","va_bus"]
	train_data = create_samples(network,N)
	return train_data
end

"""
Given a network data dict and an OPFLearn timeseries dict, compute the true Jacobians at each timestep
"""
function calc_jacobian_timeseries(network,dataset)
	network = copy(network)
	pd,qd = dataset["inputs"]["pd"],dataset["inputs"]["qd"]
	n = length(network["bus"])
	m = size(dataset["inputs"]["qd"])[1]
	J_timeseries = zeros((m,2*n,2*n)) #Tensor timeseries of Jacobian matrices
	for (t,(pd_t,qd_t)) in enumerate(zip(eachrow(pd),eachrow(qd)))
		xd_t = [pd_t;qd_t] #2nx1 vector of active and reactive power demand
		set_network_load(network,xd_t)
		J_timeseries[t,:,:] = calc_basic_jacobian_matrix(network)
	end
	return J_timeseries
end

"""
Given a network data dict, generate an AMI dataset corresponding to net active/reactive power injections and voltage magnitudes at all PQ buses.
"""
function make_ami_dataset(network::Dict{String,<:Any},sel_bus_types=1,M=100)
	bus_type_idxs = calc_bus_idx_of_type(network,sel_bus_types)
	#outputs = ["p_gen", "q_gen","vm_gen","vm_bus","va_bus"]
	ts = create_samples(network,M) #Make time series with OPF learn
	num_bus = length(network["bus"])
	#load_bus_idx = [network["load"]["$i"]["load_bus"] for i in 1:length(network["load"])]
	#gen_bus_idx = [network["gen"][string(i)]["gen_bus"] for i in 1:length(network["gen"])]
	num_type = length(bus_type_idxs)
	pnet,qnet,vm = zeros((M,num_bus)),zeros((M,num_bus)),ts["outputs"]["vm_bus"]
	for (load_idx,load) in network["load"]
		load_bus_idx = load["load_bus"]
		load_idx = load["index"]
		if load_bus_idx ∈ bus_type_idxs
			pnet[:,load_bus_idx] -= ts["inputs"]["pd"][:,load_idx]
			qnet[:,load_bus_idx] -= ts["inputs"]["qd"][:,load_idx]
		end
	end
	for (gen_idx,gen) in network["gen"]
		gen_bus_idx = gen["gen_bus"]
		gen_idx = gen["index"]
		if gen_bus_idx ∈ bus_type_idxs
			pnet[:,gen_bus_idx] += ts["outputs"]["p_gen"][:,gen_idx]
			qnet[:,gen_bus_idx] += ts["outputs"]["q_gen"][:,gen_idx]
		end
	end
	return Dict(
		"pnet" => pnet[:,bus_type_idxs], 
		"qnet" => qnet[:,bus_type_idxs], 
		"vm" => vm[:,bus_type_idxs])
end


"""
Compute finite differences of time series data
"""
function calc_finite_differences(dataset::Dict)
	pnet,qnet,vm = dataset["pnet"],dataset["qnet"],dataset["vm"]
	dp,dq,dv = -diff(pnet,dims=1),-diff(qnet,dims=1),diff(vm,dims=1)
	return Dict("dp" => dp, "dq" => dq, "dv" => dv)
end

