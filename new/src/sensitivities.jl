"""
Estimate voltage magnitude sensitivities

Params:
p - MxN matrix of deviations of active power
q - MxN matrix of deviations of reactive power
v - MxN matrix of deviations of bus voltage magnitudes
lambd - ℓ2 regularization
"""
function calc_v_sens(p,q,v,lambd=nothing)
	m,n = size(p)
	spv = inv(p'*p + lambd*eye(n))*p'*v
	sqv = inv(q'*q + lambd*eye(n))*q'*v
	return spv,sqv
end

calc_v_sens(x,v,lambd) = sxv = inv(x'*x + lambd*eye(n))*x'*v

function calc_v_sens(x,v,lambd=nothing)
	m,n = size(x)
	sxv = inv(x'*x + lambd*eye(n))*x'*v
	return sxv
end

"""
Estimate power-to-voltage magnitude sensitivities

Params:
p - MxN matrix of deviations of active power
q - MxN matrix of deviations of reactive power
v - MxN matrix of deviations of bus voltage magnitudes
lambd - ℓ2 regularization
"""
function calc_pq_sens(p,q,v,lambd=nothing)
	m,n = size(p)
	svp = inv(v'*v + lambd*eye(n))*v'*p
	svq = inv(v'*v + lambd*eye(n))*v'*q
	return svp,svq
end

calc_pq_sens(x,v,lambd=nothing) = inv(v'*v + lambd*eye(size(v)[2]))*v'*x

