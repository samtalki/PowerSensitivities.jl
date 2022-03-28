"""
Given a power factor return the implicit function theorem representation
"""
function k(pf)
    if pf==0
        return nothing
    else
        return sqrt(1-pf^2)/pf
    end
end

function kinv(pf)
    if pf==0
        return nothing
    else
        return sqrt(1/(k(pf)^2 + 1))
    end
end

function k(network::Dict{String,<:Any},sel_bus_types=[1,2])
    pf = calc_basic_power_factor(network,sel_bus_types)
    return [k(pf_i) for pf_i in pf]
end

function kinv(network::Dict{String,<:Any},sel_bus_types=[1,2])
    pf = calc_basic_power_factor(network,sel_bus_types)
    return [k_inv(pf_i) for pf_i in pf]
end

"""
Given a network data dict, return the minimum power factor to satisfy Theorem 1 as a function of the maximum power factor allowable.
"""
function pf_min(network::Dict{String,<:Any},sel_bus_types=[1,2])

end


"""
Given a network data dict, calculate the nodal power factors
"""
function calc_basic_power_factor(network::Dict{String,<:Any},sel_bus_types=[1,2])
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    s = calc_basic_bus_injection(network)[idx_sel_bus_types]
    θ = angle.(s)
    pf = [abs(cos(θi)) for θi in θ]
    #p,q = real.(s),imag.(s)
    #pf = [p_i/(sqrt(p_i^2+q_i^2)) for (p_i,q_i) in zip(p,q)]
    return pf
end

################################### TODO: Fix [idx_sel_bus_types]
"""
Given a network data dict, calculate the Δk value for the current operating point
"""
function calc_delta_k(network::Dict{String,<:Any},sel_bus_types=[1,2],ϵ=1e-3)
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    pf = calc_basic_power_factor(network,sel_bus_types)[idx_sel_bus_types]
    s = calc_basic_bus_injection(network)[idx_sel_bus_types]
    p,q = real.(s),imag.(s)
    pf = [pf_i for (i,pf_i) in enumerate(pf) if abs(p[i])≥ϵ && abs(q[i]≥ϵ)]
    try 
        abs(k(maximum(pf)) - k(minimum(pf)))
    catch
        Δk = nothing
    end
    return Δk
end

"""
Calculate the M Matrix for theorem 1
"""
function calc_M_matrix(network::Dict{String,<:Any},sel_bus_types=[1,2])
    ∂p∂θ = calc_pth_jacobian(network,sel_bus_types)
    ∂q∂θ = calc_qth_jacobian(network,sel_bus_types)
    k_max = maximum(k(network,sel_bus_types))
    #k_max = calc_k_max(network,sel_bus_types)
    M = k_max*∂p∂θ - ∂q∂θ
    return M
end


"""
Compute K matrix where K = diag(√(1-pf_i^2)/pf_i)
"""
function calc_K_matrix(network::Dict{String,<:Any},sel_bus_types=[1,2],ϵ=1e-3)
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    s = calc_basic_bus_injection(network)[idx_sel_bus_types];
    p,q = real.(s),imag.(s)
    #pf = [cos(atan(q_i/p_i)) for (p_i,q_i) in zip(p,q)]
    pf = calc_basic_power_factor(network,sel_bus_types)
    n = length(pf)
    K = zeros((n,n))
    zero_inj_indeces = [] #Array of indeces with zero p and q injections
    for (i,pf_i) in enumerate(pf)
        if abs(p[i]) ≤ ϵ && abs(q[i]) ≤ ϵ
            push!(zero_inj_indeces,i)
        elseif abs(pf_i) <= 1e-5 || pf_i == NaN || p[i] == 0.0
            K[i,i] = 0
        else
            K[i,i] = sqrt(1-pf_i^2)/pf_i
        end
    end
    k_mean = mean(diag(K))
    for i in zero_inj_indeces
        K[i,i] = k_mean ######################################## Store indeces where this occurs and replace the entry with the mean.
    end
    return K
end



"""
Given network data dict, calculate the maximum value of the K matrix
"""
function calc_k_max(network::Dict{String,<:Any},sel_bus_types=[1,2])
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types);
    s = calc_basic_bus_injection(network)[idx_sel_bus_types];
    p,q = real.(s),imag.(s)
    θ = angle.(s)
    n = length(s)
    k_max = 0
    for (i,(p_i,q_i)) in enumerate(zip(p,q))
        if abs(p_i)<1e-3
            k_ii = 0
        elseif abs(q_i)<1e-3
            k_ii = 0
        else
            #pf_i = p_i/sqrt(p_i^2+q_i^2)
            θi = θ[i]
            pf_i = cos(θi)
            k_ii = abs(sqrt(1-pf_i^2)/pf_i)
        end
        if k_ii>k_max
            k_max = k_ii     
        end
    end
    return k_max
end


"""
Given a maximum difference Δk_max between the Q-P sensitivities,
Calculate the maximum power factor distance
"""
function calc_max_pf_distance(Δk_max)
    model = Model(Ipopt.Optimizer);
    @variable(model, 1 >= pf_min >= 0);
    @variable(model, 1>= pf_max >= 0);
    @objective(model, Max, pf_max - pf_min);
    @constraint(model, pf_max>=pf_min);
    @NLconstraint(model, (sqrt(1-pf_min^2)/pf_min) - (sqrt(1-pf_max^2)/pf_max) <= Δk_max  );
    optimize!(model);
    pf_min,pf_max = value(pf_min),value(pf_max);
    Δpf_max = pf_max - pf_min;
    return Δpf_max
end

"""
Given a network data dict,
Compute the maximum difference between nodal power factors so that complex power injections can be modeled from voltage magnitudes
"""
function calc_vmag_condition(network::Dict{String,<:Any},sel_bus_types=[1,2])
    ∂p∂θ = calc_pth_jacobian(network,sel_bus_types)
    M = calc_M_matrix(network,sel_bus_types)
    pf = calc_basic_power_factor(network,sel_bus_types)
    k_max = calc_k_max(network,sel_bus_types)
    M_nonsingular = true
    Δk_max = try
        (1/(opnorm(inv(M))))*(1/(opnorm(∂p∂θ)))
    catch
        Δk_max = nothing
        M_nonsingular = false
    end
    Δk = calc_delta_k(network,sel_bus_types)
    Δpf_max = calc_max_pf_distance(Δk_max)
    return Dict(
        "Δk" => Δk,
        "M_nonsingular" => M_nonsingular,
        "Δk_max" => Δk_max,
        "Δpf_max" => Δpf_max,
        "k_max" => k_max,
        "M" => M,
        "pf" => pf)
end