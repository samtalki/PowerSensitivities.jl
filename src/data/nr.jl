# Newton-raphson power flow with mutable injection and voltage states to enable automatic differentiation.
using ForwardDiff: derivative,jacobian
using PowerModels: compute_ac_pf

abstract type SensitivityModel end

"""
Type that contains global sensitivity models for the network
"""
struct VoltageSensitivityModel <: SensitivityModel
    data::Dict{String,<:Any}
    ∂v::Function
end

function VoltageSensitivityModel(data::Dict{String,<:Any})
    voltage_response = x -> compute_voltage_response!(data,x)
    ∂v = x -> ForwardDiff.jacobian(voltage_response,x)
    return VoltageSensitivityModel(data,∂v)
end

"""
Given a network data dict, and an active-reactive power injection state change x,
Calculate the voltage response of the network
"""
function compute_voltage_response!(data::Dict{String,<:Any},x)
    updated_data = update_injection_state!(data,x)
    sol = compute_ac_pf(updated_data)
    return calc_basic_bus_voltage(sol)
end


"""
Given a network data dict, calculate the number of generators per bus
"""
function calc_gen_per_bus(data::Dict{String,<:Any})
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
function update_voltage_state!(data::Dict{String,<:Any},vph)
    n_bus = length(data["bus"])
    @assert length(vph) === 2*n_bus "Voltage phasor solution not equal to number of buses"
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
function update_injection_state!(data::Dict{String,<:Any},x)
    n_bus = length(data["bus"])
    gen_num = length(data["gen"])
    gen_per_bus = calc_gen_per_bus(data) #number of generators per bus
    @assert length(x) === 2*n_bus "Injection state solution not equal to number of buses"
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
