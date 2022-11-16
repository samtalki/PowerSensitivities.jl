include("../PowerSensitivities.jl")
import .PowerSensitivities
using PowerModels
using LinearAlgebra,Plots,JuMP,SCS

#3 bus test case
net = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/pm_matpower/case3.m"))

#Admittance matrix
Y = calc_basic_admittance_matrix(net)

#injections
s0 = calc_basic_bus_injection(net)
p0,q0 = real.(s),imag.(s)

#voltages
vph0 = calc_basic_bus_voltage(net)
vm0,va0 = abs.(vph),angle.(vph)

#Jacobian
J0 = calc_basic_jacobian_matrix(net)

#Solve AC power flow
compute_ac_pf!(net)

#injections
s = calc_basic_bus_injection(net)
p,q = real.(s),imag.(s)

#voltages
vph = calc_basic_bus_voltage(net)
vm,va = abs.(vph),angle.(vph)

#Jacobian
J = calc_basic_jacobian_matrix(net)

#Flat start
θ_2,θ_3,V_2 = 0,0,1
#Initial guess for state
x0 = [θ_2,θ_3,V_2]
