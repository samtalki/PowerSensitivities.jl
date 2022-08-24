include("../test/thm1.jl")
include("../util/matrix.jl")
include("phobs.jl")
import PowerModels as PM
import PowerSensitivities as PS
using LinearAlgebra

#Initialize a run counter
test_run_counter = 1

"""
Given a network model dict, get the buses with leading pf
Asusmed PQ
"""
function get_leading_pf_idx(net::Dict)
    ξ = PS.calc_basic_power_factor(net)
    s = PM.calc_basic_bus_injection(net)
    @assert length(s) == length(ξ)
    n_bus = length(s)
    cap_bus = []
    for (i,(ξᵢ,sᵢ)) in enumerate(zip(ξ,s))
        bus_type = net["bus"][String(i)]["type"]
        if imag(sᵢ)>0
            push!(cap_bus,i)
        end
    end
    return cap_bus
end

function get_leading_pf_idx(net::Dict,sel_bus_types::Union{AbstractArray,AbstractSet})
    ξ = PS.calc_basic_power_factor(net)
    s = PM.calc_basic_bus_injection(net)
    @assert length(s) == length(ξ)
    n_bus = length(s)
    cap_bus = []
    for (i,(ξᵢ,sᵢ)) in enumerate(zip(ξ,s))
        bus_type = net["bus"][String(i)]["type"]
        if imag(sᵢ)>0 && i ∈ sel_bus_types
            push!(cap_bus,i)
        end
    end
    return cap_bus
end

begin


    rts_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/RTS_GMLC.m"
    twok_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/ACTIVSg2000.m"
    tenk_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/case_ACTIVSg10k.m"

    paths = [rts_path, twok_path, tenk_path]
    names = ["rts-gmlc", "2k-active", "10k-active"]

    thm1_res_pq,thm1_res_pqpv = Dict(),Dict()
    thm2_res_pq,thm2_res_pqpv = Dict(),Dict()
    thm1_data_pq,thm2_data_pq = Dict(),Dict()
    thm1_data_pqpv,thm2_data_pqpv = Dict(),Dict()
    for (path,name) in zip(paths,names)
    
        #Solve the AC power flow equations
        net = PM.make_basic_network(PM.parse_file(path))
        PM.compute_ac_pf!(net)

        pq_leading_idx = get_leading_pf_idx(net,[1])
        pq_pv_leading_idx = get_leading_pf_idx(net,[1,2])
        
    end
end