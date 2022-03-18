"""
Given a network data dict, and sel_bus_types ⊂{1,2,3}:
Calculate power flow Jacobian block ∂p/∂v  
"""
function calc_pv_jacobian(network::Dict{String,<:Any},sel_bus_types::Union{Vector,Set})
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    ∂p∂v = calc_pv_jacobian(network)
    return ∂p∂v[idx_sel_bus_types,idx_sel_bus_types]
end

"""
Given a network data dict, 
Calculate the ∂p/∂v plock of the power flow Jacobian.
Usees the voltage and power injections as depicted in the basic network.
"""
function calc_pv_jacobian(network::Dict{String,<:Any})
    v = calc_basic_bus_voltage(network)
    vm,va = abs.(v),angle.(v)
    n = length(v)
    p = real.(calc_basic_bus_injection(network))
    Y = calc_basic_admittance_matrix(network)
    G,B = real.(Y),imag.(Y)
    ∂p∂v = zeros((n,n))
    for (i,p_i) in enumerate(p)
        v_i = vm[i]
        for (k,v_k) in enumerate(vm)
            if i==k
                ∂p∂v[i,i] = (p_i/v_k) + G[i,k]*v_k
            else
                θik = va[i] - va[k]
                ∂p∂v[i,k] = v_i*( G[i,k]*cos(θik)+ B[i,k]*sin(θik) )
            end
        end
    end
    return ∂p∂v
end

"""
Given a network data dict, and sel_bus_types ⊂{1,2,3}:
Calculate power flow Jacobian block ∂p/∂v  
"""
function calc_qv_jacobian(network::Dict{String,<:Any},sel_bus_types::Union{Vector,Set})
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    ∂q∂v = calc_qv_jacobian(network)
    return ∂q∂v[idx_sel_bus_types,idx_sel_bus_types]
end


"""
Given a network data dict, 
Calculate ∂q/∂v block of the power flow Jacobian.
Uses the voltage and power injections as depicted in the basic network.
"""
function calc_qv_jacobian(network::Dict{String,<:Any})
    v = calc_basic_bus_voltage(network)
    vm,va = abs.(v),angle.(v)
    n = length(v)
    q = imag.(calc_basic_bus_injection(network))
    Y = calc_basic_admittance_matrix(network)
    G,B = real.(Y),imag.(Y)
    ∂q∂v = zeros((n,n))
    for (i,q_i) in enumerate(q)
        v_i = vm[i]
        for (k,v_k) in enumerate(vm)
            if i==k
                ∂q∂v[i,i] = q_i/v_i - B[i,i]*v_i
            else
                θik = va[i]-va[k]
                ∂q∂v[i,k] = v_i*(G[i,k]*sin(θik)-B[i,k]*cos(θik))
            end
        end
    end
    return ∂q∂v
end

"""
Given voltage magnitudes, power injections, and the ∂  
"""