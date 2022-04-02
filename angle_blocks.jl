### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# ╔═╡ 3a02936c-9fcc-11ec-2bee-270c8626402c
begin
	using PowerModels
	using Ipopt
	using OPFLearn
	import SparseArrays
	using JuMP
	using Statistics
	using LinearAlgebra
	using DataFrames
	using Gadfly
end

# ╔═╡ b8d2050e-9bdb-43cc-8c2b-22c80a3b236a
pwd()

# ╔═╡ dfb0816c-aeda-4124-957a-f439fe8dde28
function plot_matrix(M)
	set_default_plot_size(3.5inch, 3.5/1.61828inch)
	spy(M)
end

# ╔═╡ 436a43d6-4d70-4eba-9f29-b198f4ba9ebe
"""
Stores data related to a Jacobian Matrix.  Only supports
sparse matrices.

* `idx_to_bus` - a mapping from 1-to-n bus idx values to data model bus ids
* `bus_to_idx` - a mapping from data model bus ids to 1-to-n bus idx values
* `matrix` - the sparse Jacobian matrix values
* `spth` - the sparse active power-angle sensitivity submatrix values
* `sqth` - the sparse reactive power-angle sensitivity submatrix values
* `spv` - the sparse active power-voltage magnitude sensitivity submatrix values
* `sqv` - the sparse reactive power-voltage magnitude sensitivity submatrix values
"""
struct JacobianMatrix{T}
    idx_to_bus::Vector{Int}
    bus_to_idx::Dict{Int,Int}
    matrix::SparseArrays.SparseMatrixCSC{T,Int}
    spth::SparseArrays.SparseMatrixCSC{T,Int}
    sqth::SparseArrays.SparseMatrixCSC{T,Int}
    spv::SparseArrays.SparseMatrixCSC{T,Int}
    sqv::SparseArrays.SparseMatrixCSC{T,Int}
end

# ╔═╡ cc904d3d-d426-4da0-81ac-052fb021d262
pq_buses(network) = filter(p->p.second["bus_type"]∈[1],network["bus"])

# ╔═╡ a8fdadc0-a662-4767-b9b9-5860d2b4f510
"""
Given a network network dict, find the bus indeces that match bus_types.
"""
function get_idx_bus_types(network,bus_types)
    bus_of_type = filter(d->d.second["bus_type"]∈bus_types,network["bus"])
    return [d.second["index"] for d in bus_of_type]
end

# ╔═╡ cf94627d-6cba-4720-a07a-fe685e5a31d0
get_bus_of_type(network,bus_types=[1]) = filter(d->d.second["bus_type"]∈bus_types,network["bus"])

# ╔═╡ 672e2a92-c90a-43e2-a143-cfb1ea4f6c0b
"""
Calculate power flow Jacobian submatrix corresponding to specified bus_type
"""
function calc_jacobian_matrix(data::Dict{String,<:Any},sel_bus_types=1)
    num_bus = length(data["bus"])
    Y = calc_admittance_matrix(data)
    J = calc_basic_jacobian_matrix(data)
    idx_to_bus,bus_to_idx = Y.idx_to_bus,Y.bus_to_idx
    #TODO: Add bus_to_idx and idx_to_bus filter
    #idx_to_bus = idx_to_bus[bus_types==bus_type]
    #bus_to_idx = filter( d_i -> d_i==bus_type, bus_to_idx)
    bus_types = [data["bus"][string(i)]["bus_type"] for i in 1:num_bus]
    idx_sel_bus_types = findall(bus_idx-> bus_idx ∈ sel_bus_types,bus_types)
    J_idx_sel_bus_types = [idx_sel_bus_types; idx_sel_bus_types .+ num_bus] #Get indeces from all blocks
    num_sel_bus_type = length(idx_sel_bus_types)
    J = J[J_idx_sel_bus_types,J_idx_sel_bus_types]
    spth,sqth = J[1:num_sel_bus_type,1:num_sel_bus_type],J[num_sel_bus_type+1:end,1:num_sel_bus_type] #Angle submatrices
    spv,sqv = J[1:num_sel_bus_type,num_sel_bus_type+1:end],J[num_sel_bus_type+1:end,num_sel_bus_type+1:end] #Voltage magnitude submatrices
    return JacobianMatrix(idx_to_bus,bus_to_idx,J,spth,sqth,spv,sqv)
end

# ╔═╡ 8c1e987a-8124-4a46-a89c-ef31aeba4f32
"""
Calculate ∂p/∂θ block of the power flow Jacobian given voltage magnitudes vm, net reactive injections qnet and block ∂q/∂v.
"""
function calc_spth_jacobian_block(sqv,vm,qnet)
    n = length(vm)
    spth = zeros((n,n))
    for (i,q_i) in enumerate(qnet)
        for (j,v_j) in enumerate(vm)
            if i==j
                spth[i,j] = v_j*sqv[i,j] - 2*q_i
            else
                spth[i,j] = v_j*sqv[i,j]
            end
        end
    end
    return spth
end

# ╔═╡ 9a66a455-aeb0-4e18-b50b-5945227b3e63
function calc_spth_jacobian_block(data::Dict{String,<:Any},sel_bus_types=1)
    idx_bus_types = get_idx_bus_types(data,sel_bus_types)
    J = calc_jacobian_matrix(data,sel_bus_types)
    vm = abs.(calc_basic_bus_voltage(data))[idx_bus_types]
	q = imag(calc_basic_bus_injection(data))[idx_bus_types]
    return calc_spth_jacobian_block(J.sqv,vm,q)
end

# ╔═╡ b4236c4e-cd36-4d77-900b-c4538ae4dc59
function calc_sqth_jacobian_block(spv,vm,pnet)
    n = length(vm)
    sqth = zeros((n,n))
    for (i,p_i) in enumerate(pnet)
        for (j,v_j) in enumerate(vm)
            if i==j
                sqth[i,j] = -v_j*spv[i,j] + 2*p_i
            else
                sqth[i,j] = -v_j*spv[i,j]
            end
        end
    end
    return sqth
end

# ╔═╡ 738f2979-c6a2-4c0d-9339-d7143164b511
"""
Given a network data dict, calculate the `∂q/∂θ` block of the power flow Jacobian. 
"""
function calc_sqth_jacobian_block(network_data::Dict{String,<:Any},sel_bus_types=1)
    idx_bus_types = get_idx_bus_types(network_data,sel_bus_types)
    J = calc_jacobian_matrix(network_data,sel_bus_types)
    vm,p = abs.(calc_basic_bus_voltage(network_data))[idx_bus_types],real(calc_basic_bus_injection(network_data))[idx_bus_types]
    return calc_sqth_jacobian_block(J.spv,vm,p)
end

# ╔═╡ 24193ad5-8e28-4351-bcc7-ebe902b9b8a7
begin
	#Load network data
	data = parse_file("data/matpower/case14.m")
	data = make_basic_network(data)
	Y = calc_admittance_matrix(data)
	J = calc_jacobian_matrix(data,1)
	vm,q = abs.(calc_basic_bus_voltage(data)),imag(calc_basic_bus_injection(data))
end

# ╔═╡ 88ad817d-5292-44f7-ba65-8e36163d4780
idx = get_idx_bus_types(data,[1])

# ╔═╡ b390b38c-e5ae-4f4b-bfdb-ed54dcad1582
function find_bus_of_type(network,sel_bus_types)
    bus_of_type = Dict()
    for (key,bus) in network["bus"]
        type = bus["bus_type"]
        if(type ∈ sel_bus_types)
            push!(bus_of_type, key => bus)
        end
    end
    return bus_of_type
end

# ╔═╡ 251026c2-d8ba-444b-94e8-10b6eb8c7afe
idx_bus_types = get_idx_bus_types(data,1)

# ╔═╡ b76bc7d2-3d31-41a4-8ece-e07ef187f99c


# ╔═╡ 1805e082-677b-4901-80db-2dca2a873b2d
for (i,bus) in data["bus"]
	println("========================")
	println(i)
	println(bus["bus_i"])
	println(bus["index"])
end

# ╔═╡ 034c3789-dc83-4bf3-8045-5ac677bfb2b7
sqth = calc_sqth_jacobian_block(data)

# ╔═╡ 28d59501-e7e8-4169-9517-c040f3b37d5e
est_sqth = J.sqth

# ╔═╡ c8a29be0-8279-428d-aafc-988884df8b24
spth = calc_spth_jacobian_block(data)

# ╔═╡ 2d566319-c015-4615-ade7-2e9e4a60d28e
est_spth = J.spth

# ╔═╡ a557566f-220f-4636-afd7-6ace1a85a7ad
begin
	p1 = spy(sqth,Guide.xlabel("Bus Index"),Guide.ylabel("Bus Index"))
	p2 = spy(est_sqth,Guide.xlabel("Bus Index"),Guide.ylabel("Bus Index"))
	title(hstack(p1,p2),"Estimating Angle Submatrices")
end

# ╔═╡ fa3b25ab-d8a3-4797-b8ad-241c3951c7f8


# ╔═╡ f5708a08-b41c-49a7-880a-be63090ff033
println(norm(sqth-est_sqth)/norm(sqth)*100,"\n",norm(spth-est_spth)/norm(spth)*100)

# ╔═╡ 87909e3e-ae63-4894-ae64-62d36bbd6053
"""
Given a network data dict, generate an AMI dataset corresponding to net active/reactive power injections and voltage magnitudes at all PQ buses.
"""
function make_ami_dataset(data::Dict{String,<:Any},sel_bus_types=1,M=100)
	bus_type_idxs = get_idx_bus_types(data,sel_bus_types)
	outputs = ["p_gen", "q_gen","vm_gen","vm_bus","va_bus"]
	ts = create_samples(data,M,output_vars=outputs) #Make time series with OPF learn
	num_bus = length(data["bus"])
	#load_bus_idx = [data["load"]["$i"]["load_bus"] for i in 1:length(data["load"])]
	#gen_bus_idx = [data["gen"][string(i)]["gen_bus"] for i in 1:length(data["gen"])]
	num_type = length(bus_type_idxs)
	pnet,qnet,vm = zeros((M,num_bus)),zeros((M,num_bus)),ts["outputs"]["vm_bus"]
	for (load_idx,load) in data["load"]
		load_bus_idx = load["load_bus"]
		load_idx = load["index"]
		if load_bus_idx ∈ bus_type_idxs
			pnet[:,load_bus_idx] += ts["inputs"]["pd"][:,load_idx]
			qnet[:,load_bus_idx] += ts["inputs"]["qd"][:,load_idx]
		end
	end
	for (gen_idx,gen) in data["gen"]
		gen_bus_idx = gen["gen_bus"]
		gen_idx = gen["index"]
		if gen_bus_idx ∈ bus_type_idxs
			pnet[:,gen_bus_idx] -= ts["outputs"]["p_gen"][:,gen_idx]
			qnet[:,gen_bus_idx] -= ts["outputs"]["q_gen"][:,gen_idx]
		end
	end
	return Dict(
		"pnet" => pnet[:,bus_type_idxs], 
		"qnet" => qnet[:,bus_type_idxs], 
		"vm" => vm[:,bus_type_idxs])
end

# ╔═╡ 5a16335f-5209-48a1-8587-5e6f31fa58d9
dataset = make_ami_dataset(data)

# ╔═╡ 6c67f202-686b-40a5-9965-a2f82ef00c81
#Test index functionality of powermodels
begin
	Ytest = calc_admittance_matrix(data)
	println(Ytest)
	for (load_idx,load) in data["load"]
		load_bus_idx = load["load_bus"]
		load_idx = load["index"]
		bus_type = data["bus"][string(load_bus_idx)]
		println((load_bus_idx,load_idx))
		println((bus_type))
	end
end


# ╔═╡ c61ffa45-fe3c-4d68-a942-f54098b453e8
PowerModels.calc

# ╔═╡ 8595f5e2-1b6e-4819-a4f0-cb693161da81
"""
Compute finite differences of time series data
"""
function get_finite_diferences(dataset::Dict)
	pnet,qnet,vm = dataset["pnet"],dataset["qnet"],dataset["vm"]
	dp,dq,dv = -diff(pnet,dims=1),-diff(qnet,dims=1),diff(vm,dims=1)
	return Dict("dp" => dp, "dq" => dq, "dv" => dv)
end

# ╔═╡ 535bfc2a-dd28-4d42-8cf3-9bc8f3ae94db
diff_data = get_finite_diferences(dataset)

# ╔═╡ 95e9551e-172b-4a57-8734-71972334e589
est_power_voltage_sens(Δx,Δv,lambd) = inv(Δv'*Δv + lambd*I(size(Δv)[2]))*Δv'*Δx

# ╔═╡ a3f59dd4-b500-484b-96ed-15d220ba3216
"""
Given dictionary of deviations, estimate power-to-voltage sensitivities
"""
function est_power_voltage_sens(diff_data::Dict,λ=1e-15)
	Δp,Δq,Δv = diff_data["dp"],diff_data["dq"],diff_data["dv"]
	spv = est_power_voltage_sens(Δp,Δv,λ)
	sqv = est_power_voltage_sens(Δq,Δv,λ)
	return spv,sqv
end

# ╔═╡ 3e146902-1c43-4d87-85c5-5a6577007583
"""
Estimates the power-to-angle sensitivities

Params:
Δp - MxN matrix of deviations of active power
Δq - MxN matrix of deviations of reactive power
Δv - MxN matrix of deviations of bus voltage magnitudes
λ - ℓ2 regularization

Returns:
Estimate of ∂p/∂θ,∂q/∂θ

"""
function est_power_angle_sens(Δp,Δq,Δv,vm,p,q,λ=1e-15)
	spv,sqv = est_power_voltage_sens(Δp,Δq,Δv,λ)
	spth,sqth = calc_spth_jacobian_block(sqv,vm,q),calc_sqth_jacobian_block(spv,vm,p)
	return spth,sqth
end

# ╔═╡ 1d1bf5c5-a01b-45af-ba66-ca816900b9ab
"""
Given dictionary of deviations, estimate the power-to-angle sensitivities

Params:
*`Δp` - MxN matrix of deviations of active power
Δq - MxN matrix of deviations of reactive power
Δv - MxN matrix of deviations of bus voltage magnitudes
λ - ℓ2 regularization

Returns:
Estimate of ∂p/∂θ,∂q/∂θ

"""
function est_power_angle_sens(diff_data::Dict,network_data,λ=1e-15,bus_types::Int=1)
	Δp,Δq,Δv = diff_data["dp"],diff_data["dq"],diff_data["dv"]
	idx_bus_types = get_idx_bus_types(network_data,bus_types)
	vm = abs.(calc_basic_bus_voltage(network_data))[idx_bus_types]
	p = real(calc_basic_bus_injection(network_data))[idx_bus_types]
	q = imag(calc_basic_bus_injection(network_data))[idx_bus_types]
	spv,sqv = est_power_voltage_sens(diff_data,λ)
	spth,sqth = calc_spth_jacobian_block(sqv,vm,q),calc_sqth_jacobian_block(spv,vm,p)
	return spth,sqth
end

# ╔═╡ 479c7dae-62ba-404a-b14b-0dfa44cddb2f
hat_spth,hat_sqth = est_power_angle_sens(diff_data,data)

# ╔═╡ b91d73dc-b94a-45f2-ab5b-58ec6336030e
calc_spth_jacobian_block(data)

# ╔═╡ 155d88c3-d3d2-46fb-aae4-9e5fc6a9e94e
hat_spth

# ╔═╡ b19a5dc3-8ab3-49b2-9ce2-7b4a7548601f
hat_sqth

# ╔═╡ c0779644-7e11-4b4b-aac2-28db5f135a00
pq_buses(data)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
OPFLearn = "426fd8f4-da10-4833-b320-d1e833f8a26f"
PowerModels = "c36e90e8-916a-50a6-bd94-075b64ef4655"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
DataFrames = "~1.3.2"
Gadfly = "~1.3.4"
Ipopt = "~0.7.0"
JuMP = "~0.21.10"
OPFLearn = "~0.1.1"
PowerModels = "~0.18.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.ASL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6252039f98492252f9e47c312c8ffda0e3b9e78d"
uuid = "ae81ac8f-d209-56e5-92de-9978fef736f9"
version = "0.1.3+0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "9f8186bc19cd1c129d367cb667215517cc03e144"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "5.0.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "3b60064cb48efe986179359e08ffb568a6d510a2"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.3"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "9a2695195199f4f20b94898c8a8ac72609e165a4"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.3"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.CoupledFields]]
deps = ["LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "6c9671364c68c1158ac2524ac881536195b7e7bc"
uuid = "7ad07ef1-bdf2-5661-9d2b-286fd4296dac"
version = "0.2.0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ae02104e835f219b8930c7664b8012c93475c340"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.2"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "dd933c4ef7b4c270aacd4eb88fa64c147492acf0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.10.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "9d3c0c762d4666db9187f363a76b47f7346e673b"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.49"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "90b158083179a6ccbce2c7eb1446d5bf9d7ae571"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.7"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "505876577b5481e50d089c1c68899dfb6faebc62"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.6"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "0dbc5b9683245f905993b51d2814202d75b34f1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.1"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "56956d1e4c1221000b7781104c58c34019792951"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Gadfly]]
deps = ["Base64", "CategoricalArrays", "Colors", "Compose", "Contour", "CoupledFields", "DataAPI", "DataStructures", "Dates", "Distributions", "DocStringExtensions", "Hexagons", "IndirectArrays", "IterTools", "JSON", "Juno", "KernelDensity", "LinearAlgebra", "Loess", "Measures", "Printf", "REPL", "Random", "Requires", "Showoff", "Statistics"]
git-tree-sha1 = "13b402ae74c0558a83c02daa2f3314ddb2d515d3"
uuid = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
version = "1.3.4"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.Hexagons]]
deps = ["Test"]
git-tree-sha1 = "de4a6f9e7c4710ced6838ca906f81905f7385fd6"
uuid = "a1b4810d-1bce-5fbd-ac56-80944d57a21f"
version = "0.2.0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "65e4589030ef3c44d3b90bdc5aac462b4bb05567"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.8"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.InfrastructureModels]]
deps = ["JuMP", "MathOptInterface", "Memento"]
git-tree-sha1 = "9ef28829a626a0d99c78344ac42c6b44b1e6c270"
uuid = "2030c09a-7f63-5d83-885d-db604e0e9cc0"
version = "0.6.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b15fc0a95c564ca2e0a7ae12c1f095ca848ceb31"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.5"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.Ipopt]]
deps = ["BinaryProvider", "Ipopt_jll", "Libdl", "LinearAlgebra", "MathOptInterface", "MathProgBase"]
git-tree-sha1 = "380786b4929b8d18d76e909c6b2eca355b7c3bd6"
uuid = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
version = "0.7.0"

[[deps.Ipopt_jll]]
deps = ["ASL_jll", "Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "MUMPS_seq_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "82124f27743f2802c23fcb05febc517d0b15d86e"
uuid = "9cc047cb-c261-5740-88fc-0cf96f7bdcc7"
version = "3.13.4+2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JSONSchema]]
deps = ["HTTP", "JSON", "URIs"]
git-tree-sha1 = "2f49f7f86762a0fbbeef84912265a1ae61c4ef80"
uuid = "7d188eb4-7ad8-530c-ae41-71a32a6d4692"
version = "0.3.4"

[[deps.JuMP]]
deps = ["Calculus", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MathOptInterface", "MutableArithmetics", "NaNMath", "Printf", "Random", "SparseArrays", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "4358b7cbf2db36596bdbbe3becc6b9d87e4eb8f5"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "0.21.10"

[[deps.Juno]]
deps = ["Base64", "Logging", "Media", "Profile"]
git-tree-sha1 = "07cb43290a840908a771552911a6274bc6c072c7"
uuid = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
version = "0.8.4"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "46efcea75c890e5d820e670516dc156689851722"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.5.4"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "3f7cb7157ef860c637f3f4929c8ed5d9716933c6"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.7"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.METIS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1d31872bb9c5e7ec1f618e8c4a56c8b0d9bddc7e"
uuid = "d00139f3-1899-568f-a2f0-47f597d42d70"
version = "5.1.1+0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[deps.MUMPS_seq_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "METIS_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "1a11a84b2af5feb5a62a820574804056cdc59c39"
uuid = "d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d"
version = "5.2.1+4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "JSON", "JSONSchema", "LinearAlgebra", "MutableArithmetics", "OrderedCollections", "SparseArrays", "Test", "Unicode"]
git-tree-sha1 = "575644e3c05b258250bb599e57cf73bbf1062901"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "0.9.22"

[[deps.MathProgBase]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9abbe463a1e9fc507f12a69e7f29346c2cdc472c"
uuid = "fdba3010-5040-5b88-9595-932c9decdf73"
version = "0.7.8"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Media]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "75a54abd10709c01f1b86b84ec225d26e840ed58"
uuid = "e89f7d12-3494-54d1-8411-f7d8b9ae1f27"
version = "0.5.0"

[[deps.Memento]]
deps = ["Dates", "Distributed", "JSON", "Serialization", "Sockets", "Syslogs", "Test", "TimeZones", "UUIDs"]
git-tree-sha1 = "19650888f97362a2ae6c84f0f5f6cda84c30ac38"
uuid = "f28f55f0-a522-5efc-85c2-fe41dfb9b2d9"
version = "1.2.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "8d9496b2339095901106961f44718920732616bb"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.22"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OPFLearn]]
deps = ["Dates", "DelimitedFiles", "Distributed", "InfrastructureModels", "Ipopt", "JSON", "JuMP", "LinearAlgebra", "MathOptInterface", "OnlineStats", "OnlineStatsBase", "PowerModels", "Random", "StatsBase"]
git-tree-sha1 = "8355fc191601a20e90adc1e1bea4894da9e445dd"
uuid = "426fd8f4-da10-4833-b320-d1e833f8a26f"
version = "0.1.1"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.OnlineStats]]
deps = ["AbstractTrees", "Dates", "LinearAlgebra", "OnlineStatsBase", "OrderedCollections", "Random", "RecipesBase", "Statistics", "StatsBase"]
git-tree-sha1 = "78feae582915781a235912de43280f98e6962a78"
uuid = "a15396b6-48d5-5d58-9928-6d29437db91e"
version = "1.5.13"

[[deps.OnlineStatsBase]]
deps = ["AbstractTrees", "Dates", "LinearAlgebra", "OrderedCollections", "Statistics", "StatsBase"]
git-tree-sha1 = "287bd0f7ee1cc2a73f08057a7a6fcfe0c23fe4b0"
uuid = "925886fa-5bf2-5e8e-b522-a9147a512338"
version = "1.4.9"

[[deps.OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c6c2ed4b7acd2137b878eb96c68e63b76199d0f"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.17+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "7e2166042d1698b6072352c74cfd1fca2a968253"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.6"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.PowerModels]]
deps = ["InfrastructureModels", "JSON", "JuMP", "LinearAlgebra", "MathOptInterface", "Memento", "NLsolve", "SparseArrays"]
git-tree-sha1 = "6c604c9c2bc232f420c428e8903e7fcb5a7df0fb"
uuid = "c36e90e8-916a-50a6-bd94-075b64ef4655"
version = "0.18.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "de893592a221142f3db370f48290e3a2ef39998f"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.4"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "cbf21db885f478e4bd73b286af6e67d1beeebe4c"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.4"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "87e9954dfa33fd145694e42337bdd3d5b07021a6"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.6.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "25405d7016a47cf2bd6cd91e66f4de437fd54a07"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.16"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.Syslogs]]
deps = ["Printf", "Sockets"]
git-tree-sha1 = "46badfcc7c6e74535cc7d833a91f4ac4f805f86d"
uuid = "cea106d9-e007-5e6c-ad93-58fe2094e9c4"
version = "0.3.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "2d4b6de8676b34525ac518de36006dc2e89c7e2e"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.7.2"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═3a02936c-9fcc-11ec-2bee-270c8626402c
# ╠═b8d2050e-9bdb-43cc-8c2b-22c80a3b236a
# ╠═dfb0816c-aeda-4124-957a-f439fe8dde28
# ╠═436a43d6-4d70-4eba-9f29-b198f4ba9ebe
# ╠═cc904d3d-d426-4da0-81ac-052fb021d262
# ╠═a8fdadc0-a662-4767-b9b9-5860d2b4f510
# ╠═cf94627d-6cba-4720-a07a-fe685e5a31d0
# ╠═88ad817d-5292-44f7-ba65-8e36163d4780
# ╠═672e2a92-c90a-43e2-a143-cfb1ea4f6c0b
# ╠═8c1e987a-8124-4a46-a89c-ef31aeba4f32
# ╠═9a66a455-aeb0-4e18-b50b-5945227b3e63
# ╠═b4236c4e-cd36-4d77-900b-c4538ae4dc59
# ╠═738f2979-c6a2-4c0d-9339-d7143164b511
# ╠═24193ad5-8e28-4351-bcc7-ebe902b9b8a7
# ╠═b390b38c-e5ae-4f4b-bfdb-ed54dcad1582
# ╠═251026c2-d8ba-444b-94e8-10b6eb8c7afe
# ╠═b76bc7d2-3d31-41a4-8ece-e07ef187f99c
# ╠═1805e082-677b-4901-80db-2dca2a873b2d
# ╠═034c3789-dc83-4bf3-8045-5ac677bfb2b7
# ╠═28d59501-e7e8-4169-9517-c040f3b37d5e
# ╠═c8a29be0-8279-428d-aafc-988884df8b24
# ╠═2d566319-c015-4615-ade7-2e9e4a60d28e
# ╠═a557566f-220f-4636-afd7-6ace1a85a7ad
# ╠═fa3b25ab-d8a3-4797-b8ad-241c3951c7f8
# ╠═f5708a08-b41c-49a7-880a-be63090ff033
# ╠═87909e3e-ae63-4894-ae64-62d36bbd6053
# ╠═5a16335f-5209-48a1-8587-5e6f31fa58d9
# ╠═6c67f202-686b-40a5-9965-a2f82ef00c81
# ╠═c61ffa45-fe3c-4d68-a942-f54098b453e8
# ╠═8595f5e2-1b6e-4819-a4f0-cb693161da81
# ╠═535bfc2a-dd28-4d42-8cf3-9bc8f3ae94db
# ╠═95e9551e-172b-4a57-8734-71972334e589
# ╠═a3f59dd4-b500-484b-96ed-15d220ba3216
# ╠═3e146902-1c43-4d87-85c5-5a6577007583
# ╠═1d1bf5c5-a01b-45af-ba66-ca816900b9ab
# ╠═479c7dae-62ba-404a-b14b-0dfa44cddb2f
# ╠═b91d73dc-b94a-45f2-ab5b-58ec6336030e
# ╠═155d88c3-d3d2-46fb-aae4-9e5fc6a9e94e
# ╠═b19a5dc3-8ab3-49b2-9ce2-7b4a7548601f
# ╠═c0779644-7e11-4b4b-aac2-28db5f135a00
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
