include("PowerSensitivities.jl")
using .PowerSensitivities


#Load network data
data = parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case14.m")
data = make_basic_network(data)
Y = calc_admittance_matrix(data)
J = calc_basic_jacobian_matrix(data)
vm,q = abs.(calc_basic_bus_voltage(data)),imag(calc_basic_bus_injection(data))
println(length(vm),length(q))


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