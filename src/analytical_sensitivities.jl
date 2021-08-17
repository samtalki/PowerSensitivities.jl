using PowerSystems
using LinearAlgebra
using DifferentialEquations

mutable struct NetworkSensitivities
    system_model::System
    pq_buses::AbstractArray
    slack_buses::AbstractArray
    ybus::Matrix
    v0::AbstractArray
end


function NetworkSensitivities(system_model::System)
    pq_buses = get_pq_buses(system_model)
    slack_buses = get_slack_buses(system_model)
    ybus = Ybus(system_model)
    v0 = voltages()
    return AnalyitcSensitivities(system_model,pq_buses,slack_buses,ybus)
end

function real_complex_voltage_sens(ℓ::Int64,model::NetworkSensitivities)
    v0 = zeros(length(model.pq_buses))
    v_span = [i for i in 0.8:0.01:1.2] #go between 0.8 and 1.2 pu
    params = (model.ybus,ℓ)
    problem = ODEProblem(dvph_dp_l!,v0,v_span,params)
    sol = solve(problem)
    return sol
end

function imag_complex_voltage_sens(n::Int64)

end

#ΔV_ == Conjugate of sensitivities
function dvph_dp_l!(ΔV_,v,p::Tuple{Matrix,Int64},p_inj) 
    ybus = p[1]
    ℓ = p[2]
    n = length(v)
    for i in 1:n  
        ΔV_[i] = ((i==ℓ) - v[ℓ]*ybus[ℓ,:]'*conj(ΔV_))/(ybus[ℓ,:]'*v) 
    end
end

function dvph_dq_l!(ΔV_,v,p::Tuple{Matrix,Int64},p_inj)
    ybus = p[1]
    ℓ = p[2]
    n = length(v)
    for i in 1:n
        ΔV_[i] = (-im*(i==ℓ) - v[ℓ]*ybus[ℓ,:]'*conj(ΔV_))/(ybus[ℓ,:]'*v)
    end
end

function dvdq_system()

end


function get_pq_buses(sys::System)
    buses = collect(get_components(Bus,sys))
    PQ_buses = [bus for bus in buses if get_bustype(bus) == BusTypes.PQ]
    return PQ_buses
end

function get_slack_buses(system::System)
    buses = collect(get_components(Bus,sys))
    return [bus for bus in buses if get_bustype(bus) == BusTypes.SLACK]
end