### Case the angle block recovery problem for a PowerModels.jl Radial Distribution Network.


function test_∂θ_computation(data,solve_ac_opf=true)
	"""
	Given a basic network data dict, test the ∂pθ and ∂qθ recovery problems.
	"""
	#Solve/do not solve AC opf
	if solve_ac_opf
		compute_ac_pf!(data)
	end

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
		"pth" => ∂pθ,
		"qth" => ∂qθ,
		"hat_pth" => hat_∂pθ,
		"hat_qth" => hat_∂qθ,
		"diff_pth" => abs.(∂pθ - hat_∂pθ),
		"diff_qth" => abs.(∂qθ - hat_∂qθ)
		)
end




