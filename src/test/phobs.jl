###
#Test Th2 at the AC power flow solution.
###
include("../PowerSensitivities.jl")
include("../test/thm1.jl")
include("../util/matrix.jl")
using PowerModels
using LinearAlgebra
using Plots
using LaTeXStrings
import .PowerSensitivities

# function is_pos_def(A::AbstractMatrix)
# 	#Hermitian part test for positive definitineness of real/complex matrices
# 	Aherm = (1/2)*(A + conj(transpose(A))
# 		)
# 	return eigmin(Aherm) > 0
# end


struct ObservabilityData
	eigmin_A::Union{Real,Complex}
	eigmin_B::Union{Real,Complex}
	A::AbstractMatrix
	B::AbstractMatrix
	K::AbstractMatrix
	vp::AbstractMatrix
	vq::AbstractMatrix
	observable::Bool
end

function calc_thm2_data(data::Dict,sel_bus_types=[1])
	K = PowerSensitivities.calc_K_matrix(data)
	pf = PowerSensitivities.calc_basic_power_factor(data)
	S = PowerSensitivities.calc_voltage_sensitivity_matrix(data)
	#dimensionality setup
	bus_idx_sel_type = PowerSensitivities.calc_bus_idx_of_type(data,sel_bus_types)
	n_bus_sel_type = length(bus_idx_sel_type)
	n_bus = length(pf)

	@assert size(K) == size(S.vp) == size(S.vq) "Mismatch in sensitivity dimensions"
	@assert n_bus == size(K,1) == size(K,2)
	
	#Truncate the test matrices to choose only the selected bus types
	pf = pf[bus_idx_sel_type]
	K = K[bus_idx_sel_type,bus_idx_sel_type]
	Svp,Svq = S.vp[bus_idx_sel_type,bus_idx_sel_type],S.vq[bus_idx_sel_type,bus_idx_sel_type]

	#Test values
	K_min,K_max = minimum(diag(K)),maximum(diag(K))
	
	#Test matrices
	A = Svp + Svq*K
	B = Svp*inv(K) + Svq
	is_observable = isinvertible(A) || isinvertible(B)
	#is_observable = eigmin(A)>0 || eigmin(B)>0#Test
	return Dict(
		"A" => A,
		"B" => B,
		"is_observable" => is_observable,
		"K" => K,
		"pf" => pf,
		"S" => S,
	)
end

function test_observability(data::Dict{String,<:Any},sel_bus_types=[1])
	K = PowerSensitivities.calc_K_matrix(data)
	pf = PowerSensitivities.calc_basic_power_factor(data)
	S = PowerSensitivities.calc_voltage_sensitivity_matrix(data)
	
	#dimensionality setup
	bus_idx_sel_type = PowerSensitivities.calc_bus_idx_of_type(data,sel_bus_types)
	n_bus_sel_type = length(bus_idx_sel_type)
	n_bus = length(pf)

	@assert size(K) == size(S.vp) == size(S.vq) "Mismatch in sensitivity dimensions"
	@assert n_bus == size(K,1) == size(K,2)
	
	#Truncate the test matrices to choose only the selected bus types
	pf = pf[bus_idx_sel_type]
	K = K[bus_idx_sel_type,bus_idx_sel_type]
	Svp,Svq = S.vp[bus_idx_sel_type,bus_idx_sel_type],S.vq[bus_idx_sel_type,bus_idx_sel_type]

	#Test values
	K_min,K_max = try 
		minimum(diag(K)),maximum(diag(K))
	catch
		nothing,nothing
	end
	#Test matrices
	A = Svp + Svq*K
	B = Svp*inv(K) + Svq
	is_observable = try
		eigmin(A)>0 || eigmin(B)>0#Test
	catch
		is_observable = false
	end
	#is_observable = is_pos_def(A) || is_pos_def(B)
	if is_observable
		return Dict(
			"eigmin_A" => eigmin(A),
			"eigmin_B" => eigmin(B),
			"eigs_A" => eigvals(A),
			"eigs_B" => eigvals(B),
			"A"=>A,
			"B"=>B,
			"K"=>K,
			"K_invertible"=>isinvertible(K),
			"svp"=>Svp,
			"svq"=>Svq,
			"is_observable"=>is_observable,
			"pf_min"=> minimum(pf),
			"pf_max"=> maximum(pf),
			"k_min" => K_min,
			"k_max" => K_max,
			)
	else
		return nothing
	end
end

function calc_ac_sol_data!(data::Dict{String,<:Any})
	compute_ac_pf!(data)
	Y = calc_basic_admittance_matrix(data)
	s = calc_basic_bus_injection(data)
	v = calc_basic_bus_voltage(data)
	p,q = real.(s),imag.(s)
	J = PowerSensitivities.calc_jacobian_matrix(data)
	S = PowerSensitivities.calc_voltage_sensitivity_matrix(data)
	return Dict(
		"Y" => Y,
		"p" => p,
		"q" => q,
		"J" => J,
		"S" => S
		)
end

"""
Given a network data dict, plot the eigenvalues of the phobs matrices
"""
function plot_eigvals(data::Dict{String,<:Any})
	K = PowerSensitivities.calc_K_matrix(data)
	S = PowerSensitivities.calc_voltage_sensitivity_matrix(data)
	ξ = PowerSensitivities.calc_basic_power_factor(data)
	
	#Phaseless observability matrices
	A = S.vp + S.vq*K
	B = S.vp*inv(K) + S.vq

	pA = plot_eigvals(A)
	title!(L"$\operatorname{eigs}\left(\frac{\partial v}{\partial p} + \frac{\partial v}{\partial q} K \right)$")
	pB = plot_eigvals(B)
	title!(L"$\operatorname{eigs}\left(\frac{\partial v}{\partial p} K^{-1} + \frac{\partial v}{\partial q} \right)$")
	p = plot(pA,pB)
	return p
end

function plot_eigvals(d::ObservabilityData)
	λA,λB = eigvals(d.A),eigvals(d.B)
	plot_A = plot_eigvals(d.A)
	plot_B = plot_eigvals(d.B)
	p = plot(plot_A,plot_B)
	return p 
end

function plot_eigvals!(data::Dict{String,<:Any})
	compute_ac_pf!(data)
	return plot_eigvals(data)
end

function plot_eigvals(A::AbstractMatrix;sorted=true)
	if sorted
		λ = sort(
			real.(
				eigvals(
					Matrix(A)
					)
				)
			)
	else 
		λ = real.(eigvals(Matrix(A)))
	end
	return plot(λ,xlabel="Bus Index",ylabel="Eigenvalue")
end


"""
Given dict of results dicts and key(s), print the value of each dict at that key
"""
function sample_results(key::String,results_dicts_dict::Dict)
	for (name,d) in results_dicts_dict
		try
			println(name*": "*String(d[key]))
		catch
			continue
		end
	end
end
