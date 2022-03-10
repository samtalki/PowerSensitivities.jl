using PowerModels
using SparseArrays


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


"""
Given a network network dict, find the bus indeces that have bus_types in sel_bus_types .
"""
function calc_bus_idx_of_type(network,sel_bus_types)
    bus_ordered = get_bus_ordered(network)
    bus_types = [bus["bus_type"] for bus in bus_ordered]
    idx_sel_bus_types = findall(bus_idx-> bus_idx ∈ sel_bus_types,bus_types)
    return idx_sel_bus_types
end

"""
Given a network network dict, return the bus dictionaries that have bus_types in sel_bus_types.
"""
function get_bus_of_type(network,sel_bus_types)
    bus_of_type = filter(d->d.second["bus_type"]∈bus_types,network["bus"])
    return bus_of_type
end

"""
Get ordered vector of buses
"""
function get_bus_ordered(network)
    b = [bus for (i,bus) in network["bus"] if bus["bus_type"] != 4]
    bus_ordered = sort(b, by=(x) -> x["index"])
    return bus_ordered
end


"network should be a PowerModels network network model; only supports networks with exactly one reference bus"
function calc_jacobian_matrix(network::Dict{String,<:Any})
    num_bus = length(network["bus"])
    Y = calc_admittance_matrix(network)
    J = calc_basic_jacobian_matrix(network)
    idx_to_bus,bus_to_idx = Y.idx_to_bus,Y.bus_to_idx
    spth,sqth = J[1:num_bus,1:num_bus],J[num_bus+1:end,1:num_bus] #Angle submatrices
    spv,sqv = J[1:num_bus,num_bus+1:end],J[num_bus+1:end,num_bus+1:end] #Voltage magnitude submatrices
    return JacobianMatrix(idx_to_bus,bus_to_idx,J,spth,sqth,spv,sqv)
end

"""
Calculate power flow Jacobian submatrix corresponding to specified bus_type
"""
function calc_jacobian_matrix(network::Dict{String,<:Any},sel_bus_types=[1,2])
    num_bus = length(network["bus"])
    Y = calc_admittance_matrix(network)
    idx_to_bus,bus_to_idx = Y.idx_to_bus,Y.bus_to_idx
    J = calc_basic_jacobian_matrix(network)
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    J_idx_sel_bus_types = [idx_sel_bus_types; idx_sel_bus_types .+ num_bus] #Shift up matrix indeces to cover all blocks
    num_sel_bus_type = length(idx_sel_bus_types)
    J = J[J_idx_sel_bus_types,J_idx_sel_bus_types]
    spth,sqth = J[1:num_sel_bus_type,1:num_sel_bus_type],J[num_sel_bus_type+1:end,1:num_sel_bus_type] #Angle submatrices
    spv,sqv = J[1:num_sel_bus_type,num_sel_bus_type+1:end],J[num_sel_bus_type+1:end,num_sel_bus_type+1:end] #Voltage magnitude submatrices
    return JacobianMatrix(idx_to_bus,bus_to_idx,J,spth,sqth,spv,sqv)
end

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

"""
Given network data dict, calculate the ∂p/∂θ block of the power flow Jacobian.
"""
function calc_spth_jacobian_block(network::Dict{String,<:Any},bus_types=1)
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    J = calc_jacobian_matrix(network,sel_bus_types)
    vm = abs.(calc_basic_bus_voltage(network))[idx_sel_bus_types]
	q = imag(calc_basic_bus_injection(network))[idx_sel_bus_types]
    return calc_spth_jacobian_block(J.sqv,vm,q)
end


"""
Calculate the ∂q/∂θ Jacobian block given vm, the ∂p/∂v Jacobian block, and pnet
"""
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

"""
Given a network data dict, calculate the `∂q/∂θ` block of the power flow Jacobian. 
"""
function calc_sqth_jacobian_block(network::Dict{String,<:Any},sel_bus_types=1)
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    J = calc_jacobian_matrix(network,sel_bus_types)
    vm,p = abs.(calc_basic_bus_voltage(network))[idx_sel_bus_types],real(calc_basic_bus_injection(network))[idx_sel_bus_types]
    return calc_sqth_jacobian_block(J.spv,vm,p)
end




"""
Computes the Jacobian submatrices
H = dP/dδ
L = dQ/dV
The derivatives are ordered according to the bus number, 
H[1, 1] is ∂P_1 / ∂δ_1
H[1, 2] is ∂P_1 / ∂δ_2
It does not exclude PV buses rows and columns
"""
function calc_basic_decoupled_jacobian_matrices(data::Dict{String,<:Any})
    if !get(data, "basic_network", false)
        Memento.warn(_LOGGER, "calc_basic_decoupled_jacobian_matrices requires basic network data and given data may be incompatible. make_basic_network can be used to transform data into the appropriate form.")
    end
    num_bus = length(data["bus"])
    v = calc_basic_bus_voltage(data)
    vm, va = abs.(v), angle.(v)
    Y = calc_basic_admittance_matrix(data)
    neighbors = [Set{Int}([i]) for i in 1:num_bus]
    I, J, V = findnz(Y)
    for nz in eachindex(V)
        push!(neighbors[I[nz]], J[nz])
        push!(neighbors[J[nz]], I[nz])
    end
    H0_I = Int64[]; H0_J = Int64[]; H0_V = Float64[];
    L0_I = Int64[]; L0_J = Int64[]; L0_V = Float64[];
    for i in 1:num_bus
        for j in neighbors[i]
            push!(H0_I, i); push!(H0_J, j); push!(H0_V, 0.0)
            push!(L0_I, i); push!(L0_J, j); push!(L0_V, 0.0)
        end
    end
    H = sparse(H0_I, H0_J, H0_V)
    L = sparse(L0_I, L0_J, L0_V)
    for i in 1:num_bus
        for j in neighbors[i]
            bus_type = data["bus"]["$(j)"]["bus_type"]
            if bus_type == 1
                if i == j
                    y_ii = Y[i,i]
                    H[i, j] =                      + vm[i] * sum( -real(Y[i,k]) * vm[k] * sin(va[i] - va[k]) + imag(Y[i,k]) * vm[k] * cos(va[i] - va[k]) for k in neighbors[i] if k != i )
                    L[i, j] = - 2*imag(y_ii)*vm[i] +         sum(  real(Y[i,k]) * vm[k] * sin(va[i] - va[k]) - imag(Y[i,k]) * vm[k] * cos(va[i] - va[k]) for k in neighbors[i] if k != i )
                else
                    y_ij = Y[i,j]
                    H[i, j] = vm[i] * vm[j] * (  real(y_ij) * sin(va[i] - va[j]) - imag(y_ij) * cos(va[i] - va[j]) )
                    L[i, j] =         vm[i] * (  real(y_ij) * sin(va[i] - va[j]) - imag(y_ij) * cos(va[i] - va[j]) )
                end
            elseif bus_type == 2
                if i == j
                    y_ii = Y[i,i]
                    H[i, j] = vm[i] * sum( -real(Y[i,k]) * vm[k] * sin(va[i] - va[k]) + imag(Y[i,k]) * vm[k] * cos(va[i] - va[k]) for k in neighbors[i] if k != i )
                    L[i, j] = 1.0
                else
                    y_ij = Y[i,j]
                    H[i, j] = vm[i] * vm[j] * (  real(y_ij) * sin(va[i] - va[j]) - imag(y_ij) * cos(va[i] - va[j]) )
                    L[i, j] = 0.0
                end
            elseif bus_type == 3
                if i == j
                    H[i, j] = 1.0
                    L[i, j] = 1.0
                end
            else
                @assert false
            end
        end
    end
    return H, L
end


