"""
Calculate ∂p/∂θ block of the power flow Jacobian given voltage magnitudes vm, net reactive injections qnet and block ∂q/∂v.
"""
function calc_spth_jacobian_block(∂q∂v,vm,qnet)
    n = length(vm)
    ∂p∂θ = zeros((n,n))
    for (i,q_i) in enumerate(qnet)
        for (j,v_j) in enumerate(vm)
            if i==j
                ∂p∂θ[i,j] = v_j*∂q∂v[i,j] - 2*q_i
            else
                ∂p∂θ[i,j] = v_j*∂q∂v[i,j]
            end
        end
    end
    return ∂p∂θ
end

"""
Given network data dict, calculate the ∂p/∂θ block of the power flow Jacobian.
"""
function calc_spth_jacobian_block(network::Dict{String,<:Any},bus_types=1)
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    J = calc_jacobian_matrix(network,sel_bus_types)
    vm = abs.(calc_basic_bus_voltage(network))[idx_sel_bus_types]
	q = imag(calc_basic_bus_injection(network))[idx_sel_bus_types]
    return calc_spth_jacobian_block(J.sqv,vm,q)
end


"""
Given a vector of voltage magnitudes vm, the ∂p/∂v Jacobian block, and a vector of real power injections pm, 
Calculate the ∂q/∂θ Jacobian block 
"""
function calc_qth_jacobian(∂p∂v,vm,pnet)
    n = length(vm)
    ∂q∂θ = zeros((n,n))
    for (i,p_i) in enumerate(pnet)
        for (j,v_j) in enumerate(vm)
            if i==j
                ∂q∂θ[i,j] = -v_j*∂p∂v[i,j] + 2*p_i
            else
                ∂q∂θ[i,j] = -v_j*∂p∂v[i,j]
            end
        end
    end
    return ∂q∂θ
end

"""
Given a network data dict, and sel_bus_types ⊂{1,2,3}
Calculate the `∂q/∂θ` block of the power flow Jacobian. 
"""
function calc_qth_jacobian(network::Dict{String,<:Any},sel_bus_types)
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    J = calc_jacobian_matrix(network,sel_bus_types)
    vm,p = abs.(calc_basic_bus_voltage(network))[idx_sel_bus_types],real(calc_basic_bus_injection(network))[idx_sel_bus_types]
    return calc_qth_jacobian(J.spv,vm,p)
end
