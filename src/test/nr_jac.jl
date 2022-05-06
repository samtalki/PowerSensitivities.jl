##### Test the jacobian solutions and the voltage sensitivities
using ForwardDiff
using LinearAlgebra
import PowerSensitivities as PS
include("/home/sam/github/PowerSensitivities.jl/src/sens/nr_sens.jl")

"""
Given a Jacobian J and vector of sel_bus_idx, return a sliced Jacobian where all 4 blocks are sliced according to sel_bus_idx.
"""
function slice_jacobian(J::AbstractMatrix,sel_bus_idx::AbstractArray)
    num_bus = Int(size(J,1)/2) #number of buses
    J_sel_idx = [sel_bus_idx; sel_bus_idx .+ num_bus] #Shift up matrix indeces to cover all blocks
    return J[J_sel_idx, J_sel_idx]
end

#Load a basic test case and calculate solution data
data = PM.make_basic_network(PM.parse_file("/home/sam/github/PowerSensitivities.jl/data/radial_test/case4_dist.m"))
J_base = PM.calc_basic_jacobian_matrix(data)  #Calculate the analytical Jacobian using the power flow equations 
s_base = PM.calc_basic_bus_injection(data)
v_base = PM.calc_basic_bus_injection(data) #Rectangular bus voltages v_r + j v_i at basecase
v_ph_base = calc_phasor_bus_voltage(data)  #Compute the phasor voltages  [θ ; vmag] at basecase

#Solve the AC Power flow solution and calculate solution data
PM.compute_ac_pf!(data)
J_sol = PM.calc_basic_jacobian_matrix(data) #Calculate the analytical Jacobian using the power flow equations 
s_sol = PM.calc_basic_bus_injection(data)
v_sol = PM.calc_basic_bus_voltage(data) #Rectangular bus voltages v_r + j v_i at the solution
v_ph_sol = calc_phasor_bus_voltage(data) #Compute the phasor voltages  [θ ; vmag] at the solution

#Get the PQ and PQ+PV bus indeces
pq_idx = PS.calc_bus_idx_of_type(data,[1])
pq_pv_idx = PS.calc_bus_idx_of_type(data,[1,2])

#Calculate the PFJ with SensitivityModel
model = SensitivityModel(data) #Create a sensitivity model
J_auto = model.J(v_ph_sol)

# Check they are the same on the specified bus types (PQ only)

@assert norm(slice_jacobian(J_sol,pq_idx)-slice_jacobian(J_auto,pq_idx)) < 1e-3