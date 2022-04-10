"""
Given a network data dict, calculate the nodal power factors
"""
calc_basic_power_factor(network::Dict{String,<:Any}) =  [abs(cos(θi)) for θi in angle.(calc_basic_bus_injection(network))]
#calc_basic_power_factor(network::Dict{String,<:Any}) =  cos.(angle.(calc_basic_bus_injection(network)))

"""
Given a real power factor pf return the implicit function theorem representation of pf
"""
function k(pf::Real)
    #@assert pf<= 1 && pf > 0 "Power factors must be between (0,1]"
    return sqrt(1-pf^2)/pf
end

"""
Given a real power factor pf and the reactive power injection return the SIGNED implicit function theorem representation of pf
This is used to make entries of the K matrix.
"""
function k(pf::Real,q::Real)
    #@assert pf<= 1 && pf > 0 "Power factors must be between (0,1]"
    return -sign(q)*abs(sqrt(1-pf^2)/pf)
end

"""
Given a network data dict return a vector of implicit function theorem representations for all node power factors
Makes entries of the K matrix.
"""
function k(network::Dict{String,<:Any})
    pf = calc_basic_power_factor(network)
    return [k(pf_i) for pf_i in pf]
end

"""
Given a real valued β, compute inverse function of k
"""
function kinv(β::Real)
    #@assert β<= 1 && β > 0 "Arugment of k^{-1} must be between (0,1]"
    if β<=0
        return nothing
    else
        return sqrt(1/(k(β)^2 + 1))
    end
end

"""
Given a network data dict, compute inverse function of k
"""
function kinv(network::Dict{String,<:Any})
    pf = calc_basic_power_factor(network)
    return [k_inv(pf_i) for pf_i in pf]
end

"""
Compute K matrix where K = diag(√(1-pf_i^2)/pf_i)
"""
function calc_K_matrix(network::Dict{String,<:Any},ϵ=1e-6)
    bad_idx = calc_bad_idx(network) #Array of indeces with zero p or zero MVA injections to be discarded
    pf = calc_basic_power_factor(network) 
    n = length(pf)
    K = zeros((n,n))
    for (i,pf_i) in enumerate(pf)
        K[i,i] = k(pf_i) #k(pf_i,q[i])#sqrt(1-pf_i^2)/pf_i
    end
    k_mean = mean(diag(K)[[i for i in 1:n if i ∉ bad_idx]]) #Replace the bad indeces with the mean of the other k entries
    for i in bad_idx
        K[i,i] = k_mean ## Replace k entry at bad idx with the mean of other ks.
    end
    return K
end

"""
Given network data dict, calculate the maximum or minimum values of the K matrix
"""
calc_k_max(network::Dict{String,<:Any}) = maximum(diag(calc_K_matrix(network)))
calc_k_min(network::Dict{String,<:Any}) = minimum(diag(calc_K_matrix(network)))

"""
Given a network data dict, calculate the Δk=k_max-k_min value for the current operating point
"""
function calc_delta_k(network::Dict{String,<:Any},ϵ=1e-6)
    K = calc_K_matrix(network)
    Δk = maximum(diag(K)) - minimum(diag(K))
    return Δk
end

function calc_delta_K_matrix(network::Dict{String,<:Any},ϵ=1e-6)
    K = calc_K_matrix(network)
    k_max = calc_k_max(network)
    return k_max*I - K
end

"""
Given network data dict and sel_bus_types, calculate the M Matrix defined in Theorem 1
"""
function calc_M_matrix(network::Dict{String,<:Any})
    ∂p∂θ = calc_pth_jacobian(network)
    ∂q∂θ = calc_qth_jacobian(network)
    k_max = calc_k_max(network)
    M = k_max*∂p∂θ - ∂q∂θ
    return M
end

"""
Given k_max, ∂p∂θ, and ∂q∂θ, calculate the M Matrix defined in Theorem 1
"""
calc_M_matrix(k_max::Real,∂p∂θ::AbstractMatrix,∂q∂θ::AbstractMatrix) = k_max*∂p∂θ - ∂q∂θ
    



