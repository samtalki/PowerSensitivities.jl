# Newton-raphson power flow with mutable injection and voltage states to enable automatic differentiation.

import ForwardDiff
using LinearAlgebra
using SparseArrays
using PowerModels: parse_file,make_basic_network
include("/home/sam/github/PowerSensitivities.jl/src/prob/nr_basic.jl")



abstract type SensitivityModel end

"""
Type that contains global sensitivity models for the network
"""
struct VoltageSensitivityModel <: SensitivityModel
    data
    ∂v::Function
end

function VoltageSensitivityModel(data)
    voltage_response = x -> compute_voltage_response!(data,x)
    ∂v = x::AbstractArray -> ForwardDiff.jacobian(voltage_response,x)
    return VoltageSensitivityModel(data,∂v)
end

"""
Given a network data dict, and an active-reactive power injection state change x,
Calculate the voltage response of the network
"""
function compute_voltage_response!(data,x)
    updated_data = update_injection_state!(data,x)
    compute_basic_ac_pf!(updated_data)
    return calc_rect_bus_voltage(updated_data)
end

"""
given a basic network data dict, returns a complex valued vector of bus voltage
values in rectangular coordinates as they appear in the network data.
"""
function calc_rect_bus_voltage(data)
    b = [bus for (i,bus) in data["bus"] if bus["bus_type"] != 4]
    bus_ordered = sort(b, by=(x) -> x["index"])
    va,vm = [bus["va"] for bus in bus_ordered],[bus["vm"] for bus in bus_ordered]
    return [va ; vm]
end


case5 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case5.m"))
Y = PM.calc_basic_admittance_matrix(case5)
sens = VoltageSensitivityModel(case5)
x = randn(10)
dvdx = sens.∂v(x)