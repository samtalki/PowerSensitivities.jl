# BSD 3-Clause License

# Copyright (c) 2022, Samuel Talkington
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import PowerModels as PM
import LinearAlgebra as LA


"""
Given complex voltage and powers, and an admittance matrix Y, calculate the power flow mismatch.
"""
function calc_basic_mismatch(V,S,Y)
    # STEP 1: Compute mismatch and check convergence
    Si = V .* conj(Y * V)
    Δp, Δq = real(S - Si), imag(S - Si)
    Δx = [Δp ; Δq]
    return Δx
end


"""
The Newton-Raphson Power Flow Algorithm from scratch. Only compatible with basic networks.
Credit ccoffrin, et al.
"""
function compute_basic_ac_pf!(data::Dict{String, Any})
    Y = PM.calc_basic_admittance_matrix(data)
    tol = 1e-4
    itr_max = 20
    itr = 0

    while itr < itr_max

        # STEP 1: Compute mismatch and check convergence
        V = PM.calc_basic_bus_voltage(data) #∈C
        S = PM.calc_basic_bus_injection(data) #∈C
        Δx = calc_basic_mismatch(V,S,Y)
        if LinearAlgebra.normInf(Δx) < tol
            break
        end
        
        # STEP 2 and 3: Compute the jacobian and update step
        J = PM.calc_basic_jacobian_matrix(data)
        vph = J \ Δx

        # STEP 4 and 5
        # update voltage variables
        data = update_voltage_state!(data,vph)
        # update power variables
        data = update_injection_state!(data,Δx)

        # update iteration counter
        itr += 1
    end
    if itr == itr_max
        @assert false "Max iteration limit"
    end
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

