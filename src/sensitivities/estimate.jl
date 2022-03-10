struct AMI
	p::Matrix
	q::Matrix
	v::Matrix
end

struct VoltageSensitivities
	spv::Matrix
	sqv::Matrix
	svp::Matrix
	svq::Matrix
end

"""
Estimate voltage magnitude sensitivities

Params:
Δp - MxN matrix of deviations of active power
Δq - MxN matrix of deviations of reactive power
Δv - MxN matrix of deviations of bus voltage magnitudes
lambd - ℓ2 regularization

Returns ∂v/∂p,∂v/∂q
"""
function est_voltage_power_sens(Δp,Δq,Δv,lambd=nothing)
	m,n = size(Δp)
	svp = inv(Δp'*Δp + lambd*eye(n))*Δp'*Δv
	svq = inv(Δq'*Δq + lambd*eye(n))*Δq'*Δv
	return svp,svq
end

est_voltage_power_sens(x,Δv,lambd) = sxv = inv(x'*x + lambd*eye(n))*x'*Δv

function est_voltage_power_sens(x,Δv,lambd=nothing)
	m,n = size(x)
	sxv = inv(x'*x + lambd*eye(n))*x'*Δv
	return sxv
end

"""
Estimate power-to-voltage magnitude sensitivities

Params:
Δp - MxN matrix of deviations of active power
Δq - MxN matrix of deviations of reactive power
Δv - MxN matrix of deviations of bus voltage magnitudes
lambd - ℓ2 regularization

Returns estimate of ∂p/∂v,∂q/∂v
"""
function est_power_voltage_sens(Δp,Δq,Δv,lambd=nothing)
	m,n = size(Δp)
	spv = inv(Δv'*Δv + lambd*eye(n))*Δv'*Δp
	sqv = inv(Δv'*Δv + lambd*eye(n))*Δv'*Δq
	return spv,sqv
end

est_power_voltage_sens(x,Δv,lambd=nothing) = inv(Δv'*Δv + lambd*eye(size(Δv)[2]))*Δv'*x

"""
Compute finite differences of time series data
"""
function get_finite_diferences(pd,pg,qd,qg,vm)
	pnet,qnet = pd-pg,qd-qg
	dp,dq,dv = diff(pnet),diff(qnet),diff(vm)
	return Dict("dp" => dp, "dq" => dq, "dv" => dv)
end

"""
Estimates the power-to-angle sensitivities

Params:
Δp - MxN matrix of deviations of active power
Δq - MxN matrix of deviations of reactive power
Δv - MxN matrix of deviations of bus voltage magnitudes
lambd - ℓ2 regularization

Returns:
Estimate of ∂p/∂θ,∂q/∂θ

"""
function est_power_angle_sens(Δp,Δq,Δv,vm,p,q,lambd=nothing)
	spv,sqv = est_power_voltage_sens(Δp,Δq,Δv,lambd)
	spth,sqth = calc_spth_jacobian_block(sqv,vm,q),calc_sqth_jacobian_block(spv,vm,p)
	return spth,sqth
end



"""
Given dictionary of deviations, estimate power-to-voltage sensitivities
"""
function est_power_voltage_sens(diff_data::Dict,λ=1e-15)
	Δp,Δq,Δv = diff_data["dp"],diff_data["dq"],diff_data["dv"]
	spv = est_power_voltage_sens(Δp,Δv,λ)
	sqv = est_power_voltage_sens(Δq,Δv,λ)
	return spv,sqv
end