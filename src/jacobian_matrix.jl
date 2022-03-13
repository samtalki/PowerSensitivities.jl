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
function calc_jacobian_matrix(network::Dict{String,<:Any},sel_bus_types)
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
given a basic network data dict, returns a sparse real valued Jacobian matrix
of the ac power flow problem.  The power variables are ordered by p and then q
while voltage values are ordered by voltage angle and then voltage magnitude.
"""
function calc_basic_jacobian_matrix(data::Dict{String,<:Any})
    if !get(data, "basic_network", false)
        Memento.warn(_LOGGER, "calc_basic_jacobian_matrix requires basic network data and given data may be incompatible. make_basic_network can be used to transform data into the appropriate form.")
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
    J0_I = Int[]
    J0_J = Int[]
    J0_V = Float64[]
    for i in 1:num_bus
        f_i_r = i
        f_i_i = i + num_bus
        for j in neighbors[i]
            x_j_fst = j + num_bus
            x_j_snd = j
            push!(J0_I, f_i_r); push!(J0_J, x_j_fst); push!(J0_V, 0.0)
            push!(J0_I, f_i_r); push!(J0_J, x_j_snd); push!(J0_V, 0.0)
            push!(J0_I, f_i_i); push!(J0_J, x_j_fst); push!(J0_V, 0.0)
            push!(J0_I, f_i_i); push!(J0_J, x_j_snd); push!(J0_V, 0.0)
        end
    end
    J = sparse(J0_I, J0_J, J0_V)
    for i in 1:num_bus
        i1 = i
        i2 = i + num_bus
        for j in neighbors[i]
            j1 = j
            j2 = j + num_bus
            bus_type = data["bus"]["$(j)"]["bus_type"]
            if bus_type == 1
                if i == j
                    y_ii = Y[i,i]
                    J[i1, j1] =                      + vm[i] * sum( -real(Y[i,k]) * vm[k] * sin(va[i] - va[k]) + imag(Y[i,k]) * vm[k] * cos(va[i] - va[k]) for k in neighbors[i] if k != i )
                    J[i1, j2] = + 2*real(y_ii)*vm[i] +         sum(  real(Y[i,k]) * vm[k] * cos(va[i] - va[k]) + imag(Y[i,k]) * vm[k] * sin(va[i] - va[k]) for k in neighbors[i] if k != i )
                    J[i2, j1] =                        vm[i] * sum(  real(Y[i,k]) * vm[k] * cos(va[i] - va[k]) + imag(Y[i,k]) * vm[k] * sin(va[i] - va[k]) for k in neighbors[i] if k != i )
                    J[i2, j2] = - 2*imag(y_ii)*vm[i] +         sum(  real(Y[i,k]) * vm[k] * sin(va[i] - va[k]) - imag(Y[i,k]) * vm[k] * cos(va[i] - va[k]) for k in neighbors[i] if k != i )
                else
                    y_ij = Y[i,j]
                    J[i1, j1] = vm[i] * vm[j] * (  real(y_ij) * sin(va[i] - va[j]) - imag(y_ij) * cos(va[i] - va[j]) )
                    J[i1, j2] =         vm[i] * (  real(y_ij) * cos(va[i] - va[j]) + imag(y_ij) * sin(va[i] - va[j]) )
                    J[i2, j1] = vm[i] * vm[j] * ( -real(y_ij) * cos(va[i] - va[j]) - imag(y_ij) * sin(va[i] - va[j]) )
                    J[i2, j2] =         vm[i] * (  real(y_ij) * sin(va[i] - va[j]) - imag(y_ij) * cos(va[i] - va[j]) )
                end
            elseif bus_type == 2
                if i == j
                    J[i1, j1] = vm[i] * sum( -real(Y[i,k]) * vm[k] * sin(va[i] - va[k]) + imag(Y[i,k]) * vm[k] * cos(va[i] - va[k]) for k in neighbors[i] if k != i )
                    J[i1, j2] = 0.0
                    J[i2, j1] = vm[i] * sum(  real(Y[i,k]) * vm[k] * cos(va[i] - va[k]) + imag(Y[i,k]) * vm[k] * sin(va[i] - va[k]) for k in neighbors[i] if k != i )
                    J[i2, j2] = 1.0
                else
                    y_ij = Y[i,j]
                    J[i1, j1] = vm[i] * vm[j] * (  real(y_ij) * sin(va[i] - va[j]) - imag(y_ij) * cos(va[i] - va[j]) )
                    J[i1, j2] = 0.0
                    J[i2, j1] = vm[i] * vm[j] * ( -real(y_ij) * cos(va[i] - va[j]) - imag(y_ij) * sin(va[i] - va[j]) )
                    J[i2, j2] = 0.0
                end
            elseif bus_type == 3
                if i == j
                    J[i1, j1] = 1.0
                    J[i1, j2] = 0.0
                    J[i2, j1] = 0.0
                    J[i2, j2] = 1.0
                end
            else
                @assert false
            end
        end
    end
    return J
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

