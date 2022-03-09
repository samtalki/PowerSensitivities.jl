"""
Given a network data dict, find the bus indeces that match `sel_bus_types`.
"""
function get_idx_bus_types(data,sel_bus_types=1)
	num_bus = length(data["bus"])
	bus_types = [data["bus"][string(i)]["bus_type"] for i in 1:num_bus]
    idx_sel_bus_types = findall(bus_idx-> bus_idx ∈ sel_bus_types,bus_types)
    return idx_sel_bus_types
end

function make_timeseries_dataset(network_data,N=100)
	outputs = ["p_gen", "q_gen","vm_gen","vm_bus","va_bus"]
	train_data = create_samples(network_data,N,output_vars=outputs)
	return train_data
end

function make_ami_dataset(network_data,dataset,bus_type="1")
	num_bus = length(network_data["bus"])
	num_load = length(network_data["load"])
	num_gen = length(network_data["gen"])
	gen_bus = [network_data["gen"][string(i)]["gen_bus"] for i in 1:length(network_data["gen"])]
	idx_sel_bus_types = get_idx_bus_types(network_data,bus_type)
	num_type = length(idx_sel_bus_types)
	P,Q,V = zeros((N,num_bus)),zeros((N,num_bus)),zeros((N,num_bus))
	for (bus_idx,load_idx) in zip(1:num_bus,1:num_load)
		P[:,k] = dataset["inputs"]["pd"][:,k] - dataset["outputs"]["p_gen"][:,bus_idx]
		#Q[:,k] = 
		#V[:,k] = 
	end
	return ami_dataset