
"""
Calculate ∂p/∂θ given v, ∂q/∂v, and q_net
"""
function calc_dp_dth(dqdv,v,q)
    n = length(v)
    dpdth = zeros((n,n))
    for (i,q_i) in enumerate(q)
        for (j,v_j) in enumerate(v)
            if i==j
                dpdth[i,j] = v_j*dqdv[i,j] - 2*q_i
            else
                dpdth[i,j] = v_j*dqdv[i,j]
            end
        end
    end
end

"""
Calculate ∂q/∂θ given v, ∂p/∂v, and p_net
"""
function calc_dq_dth(dpdv,v,p)
    n = length(v)
    dpdth = zeros((n,n))
    for (i,p_i) in enumerate(p)
        for (j,v_j) in enumerate(v)
            if i==j
                dqdth[i,j] = -v_j*dpdv[i,j] - 2*p_i
            else
                dqdth[i,j] = -v_j*dpdv[i,j]
            end
        end
    end
end
