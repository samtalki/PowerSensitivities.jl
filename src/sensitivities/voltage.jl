"""
Given a network data dict, calculate the ∂p/∂v plock of the power flow Jacobian.
Uses the voltage and power injections as depicted in the basic network.
"""
function calc_pv_jacobian(network::Dict{})
    v = calc_basic_bus_voltage(network)
    vm,va = abs.(v),angle.(v)
    n = length(v)
    p = real(calc_basic_bus_injection(network))
    Y = calc_basic_admittance_matrix(network)
    G,B = real(Y),imag(Y)
    ∂p∂v = zeros((n,n))
    for (i,q_i) in enumerate(q)
        for j(v_j) in enumerate(vm)
            if i==j
                ∂p∂v[i,j] = q_i/v_j - B[i,j]*v_j
            else
                ∂p∂v[i,j] = v_i*(G[i,j]*sin(va[i]-va[j])-B[i,j]*cos(va[i]-va[j]))
            end
        end
    end
    return ∂p∂v
    
end
"""
Given a network data dict, calculate ∂q/∂v block of the power flow Jacobian.
Uses the voltage and power injections as depicted in the basic network.
"""
function calc_qv_jacobian(network::Dict{String,<:Any})
    v = calc_basic_bus_voltage(network)
    vm,va = abs.(v),angle.(v)
    n = length(v)
    q = imag(calc_basic_bus_injection(network))
    Y = calc_basic_admittance_matrix(network)
    G,B = real(Y),imag(Y)
    ∂q∂v = zeros((n,n))
    for (i,q_i) in enumerate(q)
        for j(v_j) in enumerate(vm)
            if i==j
                ∂q∂v[i,j] = q_i/v_j - B[i,j]*v_j
            else
                ∂q∂v[i,j] = v_i*(G[i,j]*sin(va[i]-va[j])-B[i,j]*cos(va[i]-va[j]))
            end
        end
    end
    return ∂q∂v
end

"""
Given voltage magnitudes and power injections, and  
"""