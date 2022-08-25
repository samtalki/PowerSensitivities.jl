include("/home/sam/github/PowerSensitivities.jl/src/sens/analytical.jl")

data = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/radial_test/case4_dist.m"))
compute_ac_pf!(data)
Y = calc_basic_admittance_matrix(data)
v = calc_basic_bus_voltage(data)
vm,va = abs.(v),angle.(v)

∂vf = calc_vf_sensitivities(data) #phasor sensitivities
∂v = calc_vmag_sensitivities(data) #magnitude sensitivities

