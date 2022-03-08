#using Flux

#loss(S,Δx,Δv) = norm(S*Δx - Δv)^2
#∇l(S,Δx,Δv) = gradient(loss)[1]


"""
Estimate voltage magnitude sensitivities

Params:
Δp - MxN matrix of deviations of active power
Δq - MxN matrix of deviations of reactive power
Δv - MxN matrix of deviations of bus voltage magnitudes
lambd - ℓ2 regularization
"""
function calc_v_sens(Δp,Δq,Δv,lambd=nothing)
	m,n = size(Δp)
	spv = inv(Δp'*Δp + lambd*eye(n))*Δp'*Δv
	sqv = inv(Δq'*Δq + lambd*eye(n))*Δq'*Δv
	return spv,sqv
end

calc_v_sens(x,Δv,lambd) = sxv = inv(x'*x + lambd*eye(n))*x'*Δv

function calc_v_sens(x,Δv,lambd=nothing)
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
"""
function calc_pq_sens(Δp,Δq,Δv,lambd=nothing)
	m,n = size(Δp)
	svp = inv(Δv'*Δv + lambd*eye(n))*Δv'*Δp
	svq = inv(Δv'*Δv + lambd*eye(n))*Δv'*Δq
	return svp,svq
end

calc_pq_sens(x,Δv,lambd=nothing) = inv(Δv'*Δv + lambd*eye(size(Δv)[2]))*Δv'*x

