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
Given a network network dict, find the bus indeces that match sel_bus_types.
"""
function get_idx_bus_types(network,sel_bus_types)
    bus_types = [network["bus"][string(i)]["bus_type"] for i in 1:num_bus]
    idx_sel_bus_types = findall(bus_idx-> bus_idx ∈ sel_bus_types,bus_types)
    return idx_sel_bus_types
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
function calc_jacobian_matrix(network::Dict{String,<:Any},sel_bus_types=1)
    num_bus = length(network["bus"])
    Y = calc_admittance_matrix(network)
    J = calc_basic_jacobian_matrix(network)
    idx_to_bus,bus_to_idx = Y.idx_to_bus,Y.bus_to_idx
    #TODO: Add bus_to_idx and idx_to_bus filter
    #idx_to_bus = idx_to_bus[bus_types==bus_type]
    #bus_to_idx = filter( d_i -> d_i==bus_type, bus_to_idx)
    bus_types = [network["bus"][string(i)]["bus_type"] for i in 1:num_bus]
    idx_sel_bus_types = findall(bus_idx-> bus_idx ∈ sel_bus_types,bus_types)
    J_idx_sel_bus_types = [idx_sel_bus_types; idx_sel_bus_types .+ num_bus] #Get indeces from all blocks
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
function calc_spth_jacobian_block(network::Dict{String,<:Any},sel_bus_types=1)
    idx_bus_types = get_idx_bus_types(network,sel_bus_types)
    J = calc_jacobian_matrix(network,sel_bus_types)
    vm = abs.(calc_basic_bus_voltage(network))[idx_bus_types]
	q = imag(calc_basic_bus_injection(network))[idx_bus_types]
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
    idx_bus_types = get_idx_bus_types(network,sel_bus_types)
    J = calc_jacobian_matrix(network,sel_bus_types)
    vm,p = abs.(calc_basic_bus_voltage(network))[idx_bus_types],real(calc_basic_bus_injection(network))[idx_bus_types]
    return calc_sqth_jacobian_block(J.spv,vm,p)
end

