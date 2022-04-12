##### Test the jacobian solutions and the voltage sensitivities
using PowerModels
using ForwardDiff

"""
Given a vector of injection states for active and reactive power numbered bus 1,...,n,
    compute a dictionary of new results
"""
function make_injection_state(net::Dict,p_inj::Vector,q_inj::Vector)
    n_bus = length(net["bus"])
    n_gen = length(net["gen"])
    n_load = length(net["load"])
    @assert n_bus == length(p_inj) == length(q_inj)
    new_data = Dict()
    for (bus_i,(p_i,q_i)) in enumerate(zip(p_inj,q_inj))
        new_data["bus"][string(bus_i)][""]
    end
    sol = compute_ac_pf(net)["solution"]

end
"""
Given a basic network data dict and an injection state, compute the voltage response
"""
function calc_bus_voltage_response(net::Dict,inj_state::Vector)

end

begin
    
    net = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case5.m"))
    
    #Compute the values of the matrix
    J_base = calc_basic_jacobian_matrix(net)
    s_base = calc_basic_bus_injection(net)
    v_base = calc_basic_bus_injection(net)
    
    #Solve the AC power flow equations
    compute_ac_pf!(net)
    J_sol = calc_basic_jacobian_matrix(net)
    s_sol = calc_basic_bus_injection(net)
    v_sol = calc_basic_bus_voltage(net)

end