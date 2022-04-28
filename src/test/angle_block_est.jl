### Case the angle block recovery problem for a PowerModels.jl Radial Distribution Network.
using PowerModels
using LinearAlgebra
using PowerSensitivities
using Plots, LaTeXStrings

function test_∂θ_computation(data,sel_bus_types=[1],solve_ac_opf=true)
	"""
	Given a basic network data dict, test the ∂pθ and ∂qθ recovery problems.
	"""
	#Solve/do not solve AC opf
	if solve_ac_opf
		compute_ac_pf!(data)
	end

	#Calculate bus indeces of type selected
	sel_idx = calc_bus_idx_of_type(data,sel_bus_types)

	#Get state vectors
	s,v = calc_basic_bus_injection(data),calc_basic_bus_voltage(data)
	p,q = real.(s),imag.(s)
	vmag = abs.(v)

	#Compute the Jacobian and pull out submatrices
	J = calc_jacobian_matrix(data)
	∂pθ,∂qθ,∂pv,∂qv = J.pth,J.qth,J.pv,J.qv

	#Calculate the angle blocks from the submatrices
	hat_∂pθ,hat_∂qθ = calc_pth_jacobian(∂qv,vmag,p),calc_qth_jacobian(∂pv,vmag,q)

	return Dict(
		"pth" => ∂pθ[sel_idx,sel_idx],
		"qth" => ∂qθ[sel_idx,sel_idx],
		"hat_pth" => hat_∂pθ[sel_idx,sel_idx],
		"hat_qth" => hat_∂qθ[sel_idx,sel_idx],
		"diff_pth" => abs.(∂pθ - hat_∂pθ)[sel_idx,sel_idx],
		"diff_qth" => abs.(∂qθ - hat_∂qθ)[sel_idx,sel_idx]
		)
end

#--- Experiment paths and networks under study
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/case4_dist.m"
figure_path = "/home/sam/github/PowerSensitivities.jl/figures/spring_22/"
sel_bus_types = [1] #Selected bus types

#--- Start experiment
data = make_basic_network(parse_file(network_data_path))
results = test_∂θ_computation(data,sel_bus_types)

#--- Make plots
xs = ["2" "3"]
ys = ["2" "3"]

hp = heatmap(
	results["pth"],
	#title=L"\frac{\partial \boldsymbol{v}}{\partial \boldysmbol{\theta}}"
	)
hq = heatmap(
	results["qth"],
	#title=L"\frac{\partial \boldsymbol{v}}{\partial \boldsymbol{\theta}}"
	)
hat_hp = heatmap(results["hat_pth"])
hat_hq = heatmap(results["hat_qth"])
fig = plot(hp,hq,hat_hp,hat_hq)

#savefig(figure_path*"phase_angle_blocks.pdf")

