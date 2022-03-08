
"""
Calculate ∂p/∂θ block of the power flow Jacobian given Δv, ∂q/∂v, and Δq
"""
function calc_dp_dth(dqdv,Δv,Δq)
    n = length(Δv)
    dpdth = zeros((n,n))
    for (i,q_i) in enumerate(Δq)
        for (j,v_j) in enumerate(Δv)
            if i==j
                dpdth[i,j] = v_j*dqdv[i,j] - 2*q_i
            else
                dpdth[i,j] = v_j*dqdv[i,j]
            end
        end
    end
end

"""
Calculate ∂q/∂θ given Δv, ∂p/∂v, and Δp
"""
function calc_dq_dth(dpdv,Δv,Δp)
    n = length(Δv)
    dpdth = zeros((n,n))
    for (i,p_i) in enumerate(Δp)
        for (j,v_j) in enumerate(Δv)
            if i==j
                dqdth[i,j] = -v_j*dpdv[i,j] - 2*p_i
            else
                dqdth[i,j] = -v_j*dpdv[i,j]
            end
        end
    end
end
