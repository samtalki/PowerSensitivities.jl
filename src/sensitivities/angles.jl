"""
Given network data dict and sel_bus_types ⊂{1,2,3}, calculate the ∂p/∂θ block of the power flow Jacobian.
"""
function calc_pth_jacobian(network::Dict{String,<:Any},sel_bus_types::Union{Vector,Set})
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    J = calc_jacobian_matrix(network,sel_bus_types)
    vm = abs.(calc_basic_bus_voltage(network))[idx_sel_bus_types]
	q = imag(calc_basic_bus_injection(network))[idx_sel_bus_types]
    return calc_pth_jacobian(J.sqv,vm,q)
end

"""
Given network data dict, calculate power flow jacobian block ∂p/∂θ
"""
function calc_pth_jacobian(network::Dict{String,<:Any})
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
end

"""
Calculate ∂p/∂θ block of the power flow Jacobian given voltage magnitudes vm, net reactive injections qnet and block ∂q/∂v.
"""
function calc_pth_jacobian(∂q∂v::AbstractMatrix,vm::AbstractVector,q::AbstractVector)
    n = length(vm)
    ∂p∂θ = zeros((n,n))
    for (i,q_i) in enumerate(q)
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
Given a network data dict, and sel_bus_types ⊂{1,2,3}, calculate power flow Jacobian block ∂q/∂θ  
"""
function calc_qth_jacobian(network::Dict{String,<:Any},sel_bus_types::Union{Vector,Set})
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    J = calc_jacobian_matrix(network,sel_bus_types)
    vm,p = abs.(calc_basic_bus_voltage(network))[idx_sel_bus_types],real(calc_basic_bus_injection(network))[idx_sel_bus_types]
    return calc_qth_jacobian(J.spv,vm,p)
end

"""
Given network data dict, calculate power flow jacobian blcok ∂q/∂θ
"""
function calc_qth_jacobian(network::Dict{String,<:Any})
    v = calc_basic_bus_voltage(network)
    vm,va = abs.(v),angle.(v)
    n = length(v)
    p = real(calc_basic_bus_injection(network))
    Y = calc_basic_admittance_matrix(network)
    G,B = real(Y),imag(Y)
    ∂q∂θ = zeros((n,n))
    for (i,p_i) in enumerate(p)
        v_i = vm[i]
        for (j,v_j) in enumerate(vm)
            if i==j
                ∂q∂θ[i,j] = p_i - G[i,j]*v_j^2
            else
                θij = va[i]-va[j]
                ∂q∂θ[i,j] = -v_j*v_i*(G[i,j]*cos(θij) + B[i,j]*sin(θij))
            end
        end
    end
    return ∂q∂θ
end

"""
Given a vector of voltage magnitudes vm, vector of real power injections pm, and the ∂p/∂v Jacobian block calculate the ∂q/∂θ Jacobian block 
"""
function calc_qth_jacobian(∂p∂v::AbstractMatrix,vm::AbstractVector,p::AbstractVector)
    n = length(vm)
    ∂q∂θ = zeros((n,n))
    for (i,p_i) in enumerate(p)
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

