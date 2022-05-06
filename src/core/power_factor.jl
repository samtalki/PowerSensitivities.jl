"""
Given a network data dict, calculate the nodal power factors
"""
calc_basic_power_factor(network::Dict{String,<:Any}) =  [abs(cos(θi)) for θi in angle.(calc_basic_bus_injection(network))]
calc_basic_power_factor(s::AbstractArray) =  [abs(cos(θi)) for θi in angle.(s)]

"""
Implicit function theorem representation of reactive power. This is used to make entries of the K matrix.
    Params:
    pf: Power factor(s)
    q: Reactive power injections(s) used for signed representation.
"""
k(pf::Real) =  pf>0 ? sqrt(1-pf^2)/pf : NaN
k(pf::Real,q::Real) = pf>0 ? -sign(q)*abs(sqrt(1-pf^2)/pf) : NaN #SIGNED implicit representation of reactive power
k(network::Dict{String,<:Any}) = [k(pf_i) for pf_i in calc_basic_power_factor(network)] #Given a network data dict return a vector of implicit representations of reactive power for all node power factors
k(network::Dict{String,<:Any},q::Vector) = [k(pf_i,q_i) for (pf_i,q_i) in zip(calc_basic_power_factor(network),q)] #Given a network data dict return a vector of implicit representations for all node power factors

"""
Implicit function theorem representation of reactive power.
"""
kinv(β::Real) = β>0 ? sqrt(1/(k(β)^2 + 1)) : NaN
kinv(network::Dict{String,<:Any}) = [k_inv(pf) for pf in calc_basic_power_factor(network)] #Given a network data dict, compute vector of inverse representations

"""
Compute K matrix where K = diag(√(1-pf_i^2)/pf_i)
"""
function calc_K_matrix(network::Dict{String,<:Any},ϵ=1e-6)
    bad_idx = calc_bad_idx(network,ϵ) #Array of indeces with zero p or zero MVA injections to be discarded
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
calc_k_max(network::Dict{String,<:Any},ϵ=1e-6) = maximum(diag(calc_K_matrix(network,ϵ)))
calc_k_min(network::Dict{String,<:Any},ϵ=1e-6) = minimum(diag(calc_K_matrix(network,ϵ)))

"""
Given a network data dict, calculate the Δk=k_max-k_min value for the current operating point
"""
function calc_delta_k(network::Dict{String,<:Any},ϵ=1e-6)
    K = calc_K_matrix(network,ϵ)
    Δk = maximum(diag(K)) - minimum(diag(K))
    return Δk
end

function calc_delta_K_matrix(network::Dict{String,<:Any},ϵ=1e-6)
    K = calc_K_matrix(network,ϵ)
    k_max = calc_k_max(network)
    return k_max*I - K
end

"""
Given network data dict and sel_bus_types, calculate the M Matrix defined in Theorem 1
"""
function calc_M_matrix(network::Dict{String,<:Any},ϵ=1e-6)
    k_max = calc_k_max(network,ϵ)
    ∂p∂θ = calc_pth_jacobian(network)
    ∂q∂θ = calc_qth_jacobian(network)
    M = k_max*∂p∂θ - ∂q∂θ
    return M
end

"""
Given k_max, ∂p∂θ, and ∂q∂θ, calculate the M Matrix defined in Theorem 1
"""
calc_M_matrix(k_max::Real,∂p∂θ::AbstractMatrix,∂q∂θ::AbstractMatrix) = k_max*∂p∂θ - ∂q∂θ
    



