import PowerModels as PM
import LinearAlgebra as LA
"""
An AC Power Flow Solver from scratch. Only compatible with basic networks.
"""
function compute_basic_ac_pf!(data::Dict{String, Any})
    Y = PM.calc_basic_admittance_matrix(data)
    tol = 1e-4
    itr_max = 20
    itr = 0

    while itr < itr_max

        # STEP 1: Compute mismatch and check convergence
        V = calc_differentiable_bus_voltage(data)
        S = calc_differentiable_bus_injection(data)
        Si = V .* conj(Y * V)
        Δp, Δq = real(S - Si), imag(S - Si)
        if LinearAlgebra.normInf([Δp; Δq]) < tol
            break
        end

        # STEP 2 and 3: Compute the jacobian and update step
        J = calc_differentiable_jacobian_matrix(data)
        vph = J \ [Δp; Δq]

        # STEP 4 and 5
        # update voltage variables
        data = update_voltage_state!(data,vph)
        # update power variables
        data = update_injection_state!(data,[Δp; Δq])

        # update iteration counter
        itr += 1
    end
    if itr == itr_max
        @assert false "Max iteration limit"
    end
end

"""

"""
function calc_basic_mismatch(data::Dict{String,Any})

end

"""
given a basic network data dict, returns a complex valued vector of bus voltage
values in rectangular coordinates as they appear in the network data.
"""
function calc_differentiable_bus_voltage(data::Dict{String,<:Any})


    b = [bus for (i,bus) in data["bus"] if bus["bus_type"] != 4]
    bus_ordered = sort(b, by=(x) -> x["index"])

    return [bus["vm"]*cos(bus["va"]) + bus["vm"]*sin(bus["va"])im for bus in bus_ordered]
end

"""
computes the power injection of each bus in the network, with a focus on the
needs of Power Flow solvers.
excludes voltage-dependent components (e.g. shunts), these should be addressed
as needed by the calling functions.  note that voltage dependent components are
resolved during an AC Power Flow solve and are not static.
data should be a PowerModels network data model
"""
function calc_differentiable_bus_injection(data::Dict{String,<:Any})
    bus_values = Dict(bus["index"] => Dict() for (i,bus) in data["bus"])

    for (i,bus) in data["bus"]
        bvals = bus_values[bus["index"]]
        bvals["vm"] = bus["vm"]

        bvals["pd"] = 0.0
        bvals["qd"] = 0.0

        bvals["ps"] = 0.0
        bvals["qs"] = 0.0

        bvals["pg"] = 0.0
        bvals["qg"] = 0.0
    end

    for (i,load) in data["load"]
        if load["status"] != 0
            bvals = bus_values[load["load_bus"]]
            bvals["pd"] += load["pd"]
            bvals["qd"] += load["qd"]
        end
    end

    for (i,storage) in data["storage"]
        if storage["status"] != 0
            bvals = bus_values[storage["storage_bus"]]
            bvals["ps"] += storage["ps"]
            bvals["qs"] += storage["qs"]
        end
    end

    for (i,gen) in data["gen"]
        if gen["gen_status"] != 0
            bvals = bus_values[gen["gen_bus"]]
            bvals["pg"] += gen["pg"]
            bvals["qg"] += gen["qg"]
        end
    end

    Δps = Dict()
    Δqs = Dict()
    for (i,bus) in data["bus"]
        if bus["bus_type"] != 4
            bvals = bus_values[bus["index"]]
            Δp = bvals["pg"] - bvals["ps"] - bvals["pd"]
            Δq = bvals["qg"] - bvals["qs"] - bvals["qd"]
        else
            Δp = NaN
            Δq = NaN
        end

        Δps[bus["index"]] = Δp
        Δqs[bus["index"]] = Δq
    end
    bi_dict =  (Δps, Δqs)
    return [bi_dict[1][i] + bi_dict[2][i]im for i in 1:length(data["bus"])]
end


"""
Given a network data dict, calculate the number of generators per bus
"""
function _calc_gen_per_bus(data::Dict{String,<:Any})
    # Count the number of generators per bus
    gen_per_bus = Dict()
    for (i, gen) in data["gen"]
        bus_i = gen["gen_bus"]
        gen_per_bus[bus_i] = get(gen_per_bus, bus_i, 0) + 1
        # Update set point in PV buses
        if data["bus"]["$bus_i"]["bus_type"] == 2
            data["bus"]["$bus_i"]["vm"] = gen["vg"]
        end
    end
    return gen_per_bus
end


"""
Given a network data dict `data`, and vector of voltage phase angles and voltage magnitudes `vph`, 
Do: Update the network voltage variables
"""
function update_voltage_state!(data::Dict{String,<:Any},vph::AbstractVector)
    n_bus = length(data["bus"])
    #@assert length(vph) === 2*n_bus "Voltage phasor solution not equal to number of buses"
    for i in 1:n_bus
        bus_type = data["bus"]["$(i)"]["bus_type"]
        if bus_type == 1
            data["bus"]["$(i)"]["va"] = data["bus"]["$(i)"]["va"] + vph[i]
            data["bus"]["$(i)"]["vm"] = data["bus"]["$(i)"]["vm"] + vph[i+n_bus] * data["bus"]["$(i)"]["vm"] 
        end
        if bus_type == 2
            data["bus"]["$(i)"]["va"] = data["bus"]["$(i)"]["va"] + vph[i]
        end
    end
    return data
end


"""
Given a network data dict `data`, and vector of changes in active and reactive power injections `x`, 
Update the network power variable -- net injection state
"""
function update_injection_state!(data::Dict{String,<:Any},x::AbstractVector)
    n_bus = length(data["bus"])
    gen_num = length(data["gen"])
    gen_per_bus = _calc_gen_per_bus(data) #number of generators per bus
    #@assert length(x) === 2*n_bus "Injection state solution not equal to number of buses"
    #Separate the injection state vector
    Δp,Δq = x[1:n_bus-1], x[n_bus:end]
    for i in 1:gen_num
        bus_i = data["gen"]["$i"]["gen_bus"]
        bus_type = data["bus"]["$bus_i"]["bus_type"]
        num_gens = gen_per_bus[bus_i]
        if bus_type == 2
            data["gen"]["$i"]["qg"] = data["gen"]["$i"]["qg"] - Δq[bus_i] / num_gens # TODO it is ok for multiples gens in same bus?
        else bus_type == 3
            data["gen"]["$i"]["qg"] = data["gen"]["$i"]["qg"] - Δq[bus_i] / num_gens
            data["gen"]["$i"]["pg"] = data["gen"]["$i"]["pg"] - Δp[bus_i] / num_gens
        end
    end
    return data
end


"""
given a basic network data dict, returns a sparse real valued Jacobian matrix
of the ac power flow problem.  The power variables are ordered by p and then q
while voltage values are ordered by voltage angle and then voltage magnitude.
"""
function calc_differentiable_jacobian_matrix(data::Dict{String,<:Any})
    num_bus = length(data["bus"])
    v = calc_differentiable_bus_voltage(data)
    vm, va = abs.(v), angle.(v)
    Y = PM.calc_basic_admittance_matrix(data)
    neighbors = [Set{Int}([i]) for i in 1:num_bus]
    I, J, V = findnz(Y)
    for nz in eachindex(V)
        push!(neighbors[I[nz]], J[nz])
        push!(neighbors[J[nz]], I[nz])
    end
    J0_I = []
    J0_J = []
    J0_V = []
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
