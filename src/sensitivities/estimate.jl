struct VoltageSensitivities
	spv::Matrix
	sqv::Matrix
end

"""
Estimate voltage magnitude sensitivities

Params:
Δp - MxN matrix of deviations of active power
Δq - MxN matrix of deviations of reactive power
Δv - MxN matrix of deviations of bus voltage magnitudes
λ - ℓ2 regularization

Returns ∂v/∂p,∂v/∂q
"""
function est_voltage_power_sens(Δp,Δq,Δv,λ=nothing)
	m,n = size(Δp)
	svp = inv(Δp'*Δp + λ*eye(n))*Δp'*Δv
	svq = inv(Δq'*Δq + λ*eye(n))*Δq'*Δv
	return svp,svq
end

est_voltage_power_sens(x,Δv,λ) = sxv = inv(x'*x + λ*eye(n))*x'*Δv

function est_voltage_power_sens(x,Δv,λ=nothing)
	m,n = size(x)
	sxv = inv(x'*x + λ*eye(n))*x'*Δv
	return sxv
end

"""
Estimate power-to-voltage magnitude sensitivities

Params:
Δp - MxN matrix of deviations of active power
Δq - MxN matrix of deviations of reactive power
Δv - MxN matrix of deviations of bus voltage magnitudes
λ - ℓ2 regularization

Returns estimate of ∂p/∂v,∂q/∂v
"""
function est_power_voltage_sens(Δp,Δq,Δv,λ=nothing)
	m,n = size(Δp)
	spv = inv(Δv'*Δv + λ*eye(n))*Δv'*Δp
	sqv = inv(Δv'*Δv + λ*eye(n))*Δv'*Δq
	return spv,sqv
end

est_power_voltage_sens(x,Δv,λ=nothing) = inv(Δv'*Δv + λ*eye(size(Δv)[2]))*Δv'*x


"""
Estimates the power-to-angle sensitivities

Params:
Δp - MxN matrix of deviations of active power
Δq - MxN matrix of deviations of reactive power
Δv - MxN matrix of deviations of bus voltage magnitudes
λ - ℓ2 regularization

Returns:
Estimate of ∂p/∂θ,∂q/∂θ

"""
function est_power_angle_sens(Δp,Δq,Δv,vm,p,q,λ=nothing)
	spv,sqv = est_power_voltage_sens(Δp,Δq,Δv,λ)
	spth,sqth = calc_spth_jacobian_block(sqv,vm,q),calc_sqth_jacobian_block(spv,vm,p)
	return spth,sqth
end



"""
Given dictionary of deviations, estimate power-to-voltage sensitivities
"""
function est_power_voltage_sens(diff_data::Dict,λ=1e-15)
	one_vec = ones(size(diff_data["dp"])[1])
	Δp,Δq,Δv = [diff_data["dp"] one_vec],[diff_data["dq"] one_vec],[diff_data["dv"] one_vec]
	spv = est_power_voltage_sens(Δp,Δv,λ)
	sqv = est_power_voltage_sens(Δq,Δv,λ)
	return spv,sqv
end