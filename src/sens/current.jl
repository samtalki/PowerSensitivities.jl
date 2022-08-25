
#--- Current Sensitivities
abstract type CurrentSensitivities <: Sensitivities end
"""
Branch current flow single branch
"""
function l(Yij::Complex,vi::Complex,vj::Complex)
	return Yij*(vi-vj)
end

"""
Branch current flow matrix
"""
function l(Y::AbstractMatrix,v::AbstractArray)
	n = length(v)
	L = zeros((n,n))
	for (i,yᵢ) in enumerate(eachrow(Y))
		for (j,yᵢⱼ) in enumerate(yᵢ)
			L[i,j] = Y[i,j]*(v[i]-v[j])
		end
	end
	return L
end

function calc_current_sensitivities(Y,∂vf∂p,∂vf∂q)
	n = size(∂vf∂p,1)
	∂lp,∂lq= zeros((n,n)),zeros((n,n))
	for j in 1:n #for each injection location
		for (i,(dvf_i_dp,dvf_i_dq)) in enumerate(zip(eachrow(∂vf∂p),eachrow(∂vf∂p)))
			for (l,(dvf_i_dp_l,dvf_i_dq_l)) in enumerate(zip(dvf_i_dp,dvf_i_dq))	
				dvf_j_dp_l,dvf_j_dq_l = 
				∂lp[i,j] = Y[i,j]*(dvf_i_dp[i] - dvf_idq_l*im)
				∂lY[i,j]*(v[i]-v[j])
			end
		end
		return L
	end
	
	function calc_current_sensitivities(Y,∂vf∂p,∂vf∂q)
		n = size(∂vf∂p,1)
		∂cp,∂cq= zeros((n,n)),zeros((n,n))
		for j in 1:n #for each injection location
			for (i,(dvf_i_dp,dvf_i_dq)) in enumerate(zip(eachrow(∂vf∂p),eachrow(∂vf∂p)))
				#for (l,(dq[i,j] = Y[i,j]
				#end
			end
		end
	end
end
