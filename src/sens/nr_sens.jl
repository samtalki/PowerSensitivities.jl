# Newton-raphson power flow with mutable injection and voltage states to enable automatic differentiation.

using LinearAlgebra, SparseArrays
import ForwardDiff
import PowerModels as PM

"""
Type that comprises a standard Newton-Raphson sensitivity model for a power network
    data:: PowerModels network data
    Y:: Admittance matrix 
    J:: A function J: v_ph -> [ ∂p/∂θ  ∂q/∂θ ; ∂p/∂v ∂q/∂v](v_ph) that evaluates the power flow jacobian given a vector of voltage phasors.
"""
struct SensitivityModel 
    data::Dict{String, Any} #PowerModels network
    Y::SparseMatrixCSC #Network admittance matrix
    J::Function #The  power flow Jacobian found through automatic differentiatiation
end

"""
Constructor for the basic network sensitivty model.
    data:: PowerModels network
"""
function SensitivityModel(data::Dict{String, Any})
    s_inj = PM.calc_basic_bus_injection(data)
    Y = PM.calc_basic_admittance_matrix(data)
    mismatch = v_ph::AbstractArray -> calc_mismatch(v_ph,s_inj,Y)
    J = v_ph::AbstractArray -> ForwardDiff.jacobian(mismatch,v_ph) #Note: v_ph = [θ ; vmag]
    return SensitivityModel(data,Y,J)
end

"""
Given complex voltage and powers, and an admittance matrix Y, calculate the power flow mismatch.
"""
function calc_mismatch(v_ph,s,Y)
    n_bus = size(Y,1)
    # Convert phasor to rectangular
    function calc_rect_bus_voltage(v_ph::AbstractArray)
        θ,vm = v_ph[1:n_bus],v_ph[n_bus+1:end]
        [vmᵢ*(cos(θᵢ) + sin(θᵢ)*im) for (θᵢ,vmᵢ) in zip(θ,vm)]
    end
    v_rect = calc_rect_bus_voltage(v_ph)
    si = v_rect .* conj(Y * v_rect) #compute the injection
    Δp, Δq = real(s - si), imag(s - si)
    Δx = [Δp ; Δq]
    return Δx
end

"""
1st definition: 
given a PowerModels network, returns a [vangle; vmag] vector of bus voltage phasor quantities
as they appear in the network data.
"""
function calc_phasor_bus_voltage(data::Dict{String, Any})
    b = [bus for (i,bus) in data["bus"] if bus["bus_type"] != 4]
    bus_ordered = sort(b, by=(x) -> x["index"])
    θ,vm = [bus["va"] for bus in bus_ordered],[bus["vm"] for bus in bus_ordered]
    return [θ ; vm]
end

"""
2nd definition (Multiple dispatch!): 
Given vector of rectangular complex voltages [e+jf], calculate the phasor form of the voltages [vangle; vmag]
"""
function calc_phasor_bus_voltage(v_rect::AbstractArray)
    θ,vm = [angle(v_i) for v_i in v_rect],[abs(v_i) for v_i in v_rect]
    return [θ ; vm]
end

#Load a basic test case
case5 = PM.make_basic_network(PM.parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case5.m"))

#Solve the AC Power flow equations
#PM.compute_ac_pf!(case5)

#Compute the phasor voltages at the solution [θ ; vmag]
v_ph_sol = calc_phasor_bus_voltage(case5)

#Calculate the analytical Jacobian using the power flow equations and the SensitivityModel version. Check they are the same.
model = SensitivityModel(case5) #Create a sensitivity model
J_analytic = PM.calc_basic_jacobian_matrix(case5) #Equations-based jacobian
J_automatic = model.J(v_ph_sol)
@assert norm(J_analytic-J_automatic) < 1e-3