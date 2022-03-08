using PowerModels
using SparseArrays

"""
Stores data related to a Jacobian Matrix.  Only supports
sparse matrices.

* `idx_to_bus` - a mapping from 1-to-n bus idx values to data model bus ids
* `bus_to_idx` - a mapping from data model bus ids to 1-to-n bus idx values
* `matrix` - the sparse Jacobian matrix values
* `dpdth` - the sparse active power-angle sensitivity submatrix values
* `dqdth` - the sparse reactive power-angle sensitivity submatrix values
* `dpdv` - the sparve active power-voltage magnitude sensitivity submatrix values
* `dqdv` - the sparve reactive power-voltage magnitude sensitivity submatrix values
"""
struct JacobianMatrix{T}
    idx_to_bus::Vector{Int}
    bus_to_idx::Dict{Int,Int}
    matrix::SparseArrays.SparseMatrixCSC{T,Int}
    dpdth::SparseArrays.SparseMatrixCSC{T,Int}
    dqdth::SparseArrays.SparseMatrixCSC{T,Int}
    dpdv::SparseArrays.SparseMatrixCSC{T,Int}
    dqdv::SparseArrays.SparseMatrixCSC{T,Int}
end

"data should be a PowerModels network data model; only supports networks with exactly one reference bus"
function calc_jacobian_matrix(data::Dict{String,<:Any})
    num_bus = length(data["bus"])
    Y = calc_admittance_matrix(data)
    J = calc_basic_jacobian_matrix(data)
    idx_to_bus,bus_to_idx = Y.idx_to_bus,Y.bus_to_idx
    dpdth,dqdth = J[1:num_bus,1:num_bus],J[num_bus+1:end,1:num_bus] #Angle submatrices
    dpdv,dqdv = J[1:num_bus,num_bus+1:end],J[num_bus+1:end,num_bus+1:end] #Voltage magnitude submatrices
    return JacobianMatrix(idx_to_bus,bus_to_idx,J,dpdth,dqdth,dpdv,dqdv)
end

"""
Calculate power flow Jacobian submatrix corresponding to specified bus_type
"""
function calc_jacobian_matrix(data::Dict{String,<:Any},sel_bus_types)
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
    dpdth,dqdth = J[1:num_sel_bus_type,1:num_sel_bus_type],J[num_sel_bus_type+1:end,1:num_sel_bus_type] #Angle submatrices
    dpdv,dqdv = J[1:num_sel_bus_type,num_sel_bus_type+1:end],J[num_sel_bus_type+1:end,num_sel_bus_type+1:end] #Voltage magnitude submatrices
    return JacobianMatrix(idx_to_bus,bus_to_idx,J,dpdth,dqdth,dpdv,dqdv)
end

"""
Calculate ∂p/∂θ block of the power flow Jacobian given voltage magnitudes vm, net reactive injections qnet and block ∂q/∂v.
"""
function calc_dpdth_jacobian_block(dqdv,vm,qnet)
    n = length(vm)
    dpdth = zeros((n,n))
    for (i,q_i) in enumerate(qnet)
        for (j,v_j) in enumerate(vm)
            if i==j
                dpdth[i,j] = v_j*dqdv[i,j] - 2*q_i
            else
                dpdth[i,j] = v_j*dqdv[i,j]
            end
        end
    end
end

"""
Calculate ∂q/∂θ given vm, ∂p/∂v, and pnet
"""
function calc_dqdth_jacobian_block(dpdv,vm,pnet)
    n = length(vm)
    dpdth = zeros((n,n))
    for (i,p_i) in enumerate(pnet)
        for (j,v_j) in enumerate(vm)
            if i==j
                dqdth[i,j] = -v_j*dpdv[i,j] - 2*p_i
            else
                dqdth[i,j] = -v_j*dpdv[i,j]
            end
        end
    end
end
