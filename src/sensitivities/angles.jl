"""
Given network data dict and sel_bus_types ⊂{1,2,3}:
Calculate the ∂p/∂θ block of the power flow Jacobian.
"""
function calc_pth_jacobian(network::Dict{String,<:Any},sel_bus_types::Union{Vector,Set})
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    ∂p∂θ = calc_pth_jacobian(network)
    return ∂p∂θ[idx_sel_bus_types,idx_sel_bus_types]
end

"""
Given network data dict:
Calculate power flow jacobian block ∂p/∂θ
"""
function calc_pth_jacobian(network::Dict{String,<:Any})
    q = imag.(calc_basic_bus_injection(network))
    v = calc_basic_bus_voltage(network)
    vm,va = abs.(v),angle.(v)
    n = length(vm)
    Y = calc_basic_admittance_matrix(network)
    G,B = real.(Y),imag.(Y)
    ∂p∂θ = zeros((n,n))
    for (i,q_i) in enumerate(q)
        v_i = vm[i]
        for (k,v_k) in enumerate(vm)
            if i==k
                ∂p∂θ[i,k] = -q_i - B[i,k]*v_i^2
            else
                θik = va[i]-va[k]
                ∂p∂θ[i,k] = v_i*v_k*(G[i,k]*sin(θik) - B[i,k]*cos(θik))
            end
        end
    end
    return ∂p∂θ
end

"""
Given voltage magnitudes vm, net reactive injections qnet and block ∂q/∂v:
Calculate ∂p/∂θ block of the power flow Jacobian
"""
function calc_pth_jacobian(∂q∂v::AbstractMatrix,vm::AbstractVector,q::AbstractVector)
    n = length(vm)
    ∂p∂θ = zeros((n,n))
    for (i,q_i) in enumerate(q)
        for (k,v_k) in enumerate(vm)
            if i==k
                ∂p∂θ[i,i] = v_k*∂q∂v[i,i] - 2*q_i
            else
                ∂p∂θ[i,k] = v_k*∂q∂v[i,k]
            end
        end
    end
    return ∂p∂θ
end

"""
Given a network data dict, and sel_bus_types ⊂{1,2,3}:
Calculate power flow Jacobian block ∂q/∂θ  
"""
function calc_qth_jacobian(network::Dict{String,<:Any},sel_bus_types::Union{Vector,Set})
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    ∂q∂θ = calc_qth_jacobian(network)
    return ∂q∂θ[idx_sel_bus_types,idx_sel_bus_types]
end

"""
Given network data dict:
Calculate power flow jacobian blcok ∂q/∂θ
"""
function calc_qth_jacobian(network::Dict{String,<:Any})
    v = calc_basic_bus_voltage(network)
    vm,va = abs.(v),angle.(v)
    n = length(v)
    p = real(calc_basic_bus_injection(network))
    Y = calc_basic_admittance_matrix(network)
    G,B = real.(Y),imag.(Y)
    ∂q∂θ = zeros((n,n))
    for (i,p_i) in enumerate(p)
        v_i = vm[i]
        for (k,v_k) in enumerate(vm)
            if i==k
                ∂q∂θ[i,i] = p_i - G[i,i]*v_i^2
            else
                θik = va[i]-va[k]
                ∂q∂θ[i,k] = -v_i*v_k*(G[i,k]*cos(θik) + B[i,k]*sin(θik))
            end
        end
    end
    return ∂q∂θ
end

"""
Given a vector of voltage magnitudes vm, vector of real power injections pm, and the ∂p/∂v Jacobian block:
Calculate the ∂q/∂θ Jacobian block 
"""
function calc_qth_jacobian(∂p∂v::AbstractMatrix,vm::AbstractVector,p::AbstractVector)
    n = length(vm)
    ∂q∂θ = zeros((n,n))
    for (i,p_i) in enumerate(p)
        for (k,v_k) in enumerate(vm)
            if i==k
                ∂q∂θ[i,k] = -v_k*∂p∂v[i,k] + 2*p_i
            else
                ∂q∂θ[i,k] = -v_k*∂p∂v[i,k]
            end
        end
    end
    return ∂q∂θ
end


