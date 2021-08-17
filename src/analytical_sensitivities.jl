using PowerSystems
using LinearAlgebra
using DifferentialEquations

mutable struct LinearSensitivityModel
    system_model::System
    pq_buses::AbstractArray
    slack_buses::AbstractArray
    ybus::AbstractArray
    vph_0::AbstractArray
end


function LinearSensitivityModel(system_model::System)
    pq_buses = get_pq_buses(system_model)
    slack_buses = get_slack_buses(system_model)
    ybus = Ybus(system_model).data
    vph_0 = get_voltage_phasors_rectangular(system_model)
    return LinearSensitivityModel(system_model,pq_buses,slack_buses,ybus,vph_0)
end

function vph_p_sens(ℓ::Int64,model::LinearSensitivityModel)
    vph_0 = model.vph_0
    v_span = (0.8,1.2) #go between 0.8 and 1.2 pu
    params = (model.ybus,ℓ)
    problem = ODEProblem(dvph_dp_l!,vph_0,v_span,params)
    sol = solve(problem)
    return sol
end

function vph_q_sens(ℓ::Int64,model::LinearSensitivityModel)
    vph_0 = model.vph_0
    v_span = (0.8,1.2) #go between 0.8 and 1.2 pu
    params = (model.ybus,ℓ)
    problem = ODEProblem(dvph_dq_l!!,vph_0,v_span,params)
    sol = solve(problem)
    return sol
end

#ΔV_ == Conjugate of sensitivities
function dvph_dp_l!(ΔV_,v,p::Tuple,p_inj) 
    ybus = p[1]
    ℓ = p[2]
    n = length(v)
    for i in 1:n  
        # println("Ybus Row Size: ",size(ybus[ℓ,:]))
        # println("ΔV_ Size: ",size(ΔV_))
        # println("Ybus Row: ",ybus[ℓ,:])
        # println("ΔV_: ",ΔV_)
        # println("conj(ΔV_)",conj(ΔV_))
        # println("v shape: ",size(v))
        # println("v: ",v)
        # println("Numerator: (i==ℓ) - v[ℓ]*ybus[ℓ,:]'*conj(ΔV_): ",(i==ℓ) - v[ℓ]*ybus[ℓ,:]'*conj(ΔV_))
        # println("Denom: ybus[ℓ,:]'*v: ",ybus[ℓ,:]'*v)
        ΔV_[i] = ((i==ℓ) - v[ℓ]*ybus[ℓ,:]'*conj(ΔV_))/(ybus[ℓ,:]'*v) 
        #ΔV_[i] = ((i==ℓ) - v[ℓ]*sum(ybus[ℓ,:].*conj(ΔV_)))/(ybus[ℓ,:]'*v) 
    end
end

function dvph_dq_l!(ΔV_,v,p::Tuple{Matrix,Int64},p_inj)
    ybus = p[1]
    ℓ = p[2]
    n = length(v)
    for i in 1:n
        ΔV_[i] = (-im*(i==ℓ) - v[ℓ]*sum(ybus[ℓ,:].*conj(ΔV_)))/(ybus[ℓ,:]'*v)
    end
end

function dvdq_system()

end


function get_pq_buses(sys::System)
    buses = collect(get_components(Bus,sys))
    PQ_buses = [bus for bus in buses if get_bustype(bus) == BusTypes.PQ]
    return PQ_buses
end

function get_slack_buses(sys::System)
    buses = collect(get_components(Bus,sys))
    return [bus for bus in buses if get_bustype(bus) == BusTypes.SLACK]
end