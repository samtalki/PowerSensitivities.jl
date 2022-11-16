### Plot the angle block recovery problem for a PowerModels.jl Network.
### Plot the gershdics for a PowerModels.jl Network.
using PowerModels
using LinearAlgebra
using Plots, LaTeXStrings,ColorSchemes
include("../PowerSensitivities.jl")
import .PowerSensitivities

function test_∂θ_computation(data,sel_bus_types=[1],solve_ac_opf=true)
	"""
	Given a basic network data dict, test the ∂pθ and ∂qθ recovery problems.
	"""
	#Solve/do not solve AC opf
	if solve_ac_opf
		compute_ac_pf!(data)
	end

	#Calculate bus indeces of type selected
	sel_idx = PowerSensitivities.calc_bus_idx_of_type(data,sel_bus_types)

	#Get state vectors
	s,v = calc_basic_bus_injection(data),calc_basic_bus_voltage(data)
	p,q = real.(s),imag.(s)
	vmag = abs.(v)

	#Compute the Jacobian and pull out submatrices
	J = PowerSensitivities.calc_jacobian_matrix(data)
	∂pθ,∂qθ,∂pv,∂qv = J.pth,J.qth,J.pv,J.qv

	#Calculate the angle blocks from the submatrices
	hat_∂pθ,hat_∂qθ = PowerSensitivities.calc_pth_jacobian(∂qv,vmag,p),PowerSensitivities.calc_qth_jacobian(∂pv,vmag,q)

	return Dict(
		"pth" => Matrix(∂pθ[sel_idx,sel_idx]),
		"qth" => Matrix(∂qθ[sel_idx,sel_idx]),
		"hat_pth" => Matrix(hat_∂pθ[sel_idx,sel_idx]),
		"hat_qth" => Matrix(hat_∂qθ[sel_idx,sel_idx]),
		"diff_pth" => Matrix(abs.(∂pθ - hat_∂pθ)[sel_idx,sel_idx]),
		"diff_qth" => Matrix(abs.(∂qθ - hat_∂qθ)[sel_idx,sel_idx])
		)
end

#--- Experiment paths and networks under study
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/pm_matpower/case14.m"
figure_path = "/home/sam/github/PowerSensitivities.jl/figures/fall_22/"
sel_bus_types = [1] #Selected bus types

#--- Start experiment
data = make_basic_network(parse_file(network_data_path))
results = test_∂θ_computation(data,sel_bus_types)

#--- Make plots
xs = ["2" "3"]
ys = ["2" "3"]

hp = heatmap(
	results["pth"],
	title=L"$\frac{\partial p}{\partial \theta}$ (At power flow sol.)",
	c=:curl,
	ylabel="Bus index",
	)
hq = heatmap(
	results["qth"],
	title=L"$\frac{\partial q}{\partial \theta}$ (At power flow sol.)",
	c=:curl,
	
	)
hat_hp = heatmap(
	results["hat_pth"],
	title=L"$\frac{\partial p}{\partial \theta}$ ($\theta$-free expression)",
	c=:curl,
	xlabel="Bus index",
	ylabel="Bus index",
	)
hat_hq = heatmap(
	results["hat_qth"],
	title=L"$\frac{\partial q}{\partial \theta}$ ($\theta$-free expression)",
	c=:curl,
	xlabel="Bus index",
)
fig = plot(
	hp,hq,hat_hp,hat_hq,
	titlefontsize=11,
	tickfontsize=6,
	labelfontsize=10,
	grid=true,
	cb=:best
)

savefig(fig,figure_path*"phase_angle_blocks_case14.pdf")

discs = PowerSensitivities.plot_gershdisc(data)
savefig(discs,figure_path*"discs_case14.pdf")

