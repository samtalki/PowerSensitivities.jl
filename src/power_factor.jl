"""
Given a network data dict, calculate the nodal power factors
"""
function calc_basic_power_factor(network::Dict{String,<:Any},sel_bus_types=[1,2])
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    s = calc_basic_bus_injection(network)[idx_sel_bus_types]
    θ = angle.(s)
    #p,q = real.(s),imag.(s)
    #pf = [p_i/(sqrt(p_i^2+q_i^2)) for (p_i,q_i) in zip(p,q)]
    pf = [cos(θi) for θi in θ]
    return pf
end

#Calculate power factor without bus truncation
calc_basic_power_factor(network::Dict{String,<:Any}) =  [cos(θi) for θi in angle.(calc_basic_bus_injection(network))]

"""
Given a real power factor pf return the implicit function theorem representation of pf
"""
function k(pf::Real)
    if pf==0
        return nothing
    else
        return sqrt(1-pf^2)/pf
    end
end

"""
Given a real power factor pf and the reactive power injection return the SIGNED implicit function theorem representation of pf
This is used to make entries of the K matrix.
"""
function k(pf::Real,q::Real)
    if pf==0
        return nothing
    else
        return sign(q)*abs(sqrt(1-pf^2)/pf)
    end
end

"""
Given a network data dict return a vector of implicit function theorem representations for all node power factors
Makes entries of the K matrix.
"""
function k(network::Dict{String,<:Any},sel_bus_types=[1,2])
    pf = calc_basic_power_factor(network,sel_bus_types)
    return [k(pf_i) for pf_i in pf]
end

"""
Given a real valued β, compute inverse function of k
"""
function kinv(β::Real)
    if β≤0
        return nothing
    else
        return sqrt(1/(k(β)^2 + 1))
    end
end

"""
Given a network data dict, compute inverse function of k
"""
function kinv(network::Dict{String,<:Any},sel_bus_types=[1,2])
    pf = calc_basic_power_factor(network,sel_bus_types)
    return [k_inv(pf_i) for pf_i in pf]
end



"""
Compute K matrix where K = diag(√(1-pf_i^2)/pf_i)
"""
function calc_K_matrix(network::Dict{String,<:Any},sel_bus_types=[1,2],ϵ=1e-6)
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    s = calc_basic_bus_injection(network)[idx_sel_bus_types];
    p,q,pf = real.(s),imag.(s),calc_basic_power_factor(network,sel_bus_types) 
    n = length(pf)
    K = zeros((n,n))
    bad_idx = [] #Array of indeces with zero p or zero MVA injections to be discarded
    for (i,pf_i) in enumerate(pf)
        if abs(p[i]) ≤ ϵ && abs(q[i]) ≤ ϵ #If there is no apparent power injection, it doesn't make sense
            push!(bad_idx,i)
        elseif abs(pf_i) <= 1e-5 || pf_i == NaN || p[i] == 0.0 #If there is no real power injection, it doesn't make sense
            push!(bad_idx,i) #K[i,i] = 0
        elseif q[i] > 0 #Ignore capacitive injections
            push!(bad_idx,i)
        else
            K[i,i] = k(pf_i)#k(pf_i,q[i])#sqrt(1-pf_i^2)/pf_i
        end
    end
    k_mean = mean(diag(K)[[i for i in 1:n if i ∉ bad_idx]]) #Replace the bad indeces with the mean of the other k entries
    for i in bad_idx
        K[i,i] = k_mean ## Replace k entry at bad idx with the mean of other ks.
    end
    return K
end

function calc_delta_K_matrix(network::Dict{String,<:Any},sel_bus_types=[1,2],ϵ=1e-6)
    K = calc_K_matrix(network,sel_bus_types)
    k_max = maximum(K)
    return k_max .- K
end

"""
Given network data dict, calculate the maximum value of the K matrix
"""
calc_k_max(network::Dict{String,<:Any},sel_bus_types=[1,2]) = maximum(calc_K_matrix(network,sel_bus_types))
calc_k_min(network::Dict{String,<:Any},sel_bus_types=[1,2]) = minimum(calc_K_matrix(network,sel_bus_types))



"""
Given a network data dict, calculate the Δk=k_max-k_min value for the current operating point
"""
function calc_delta_k(network::Dict{String,<:Any},sel_bus_types=[1,2],ϵ=1e-6)
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    pf = calc_basic_power_factor(network,sel_bus_types)
    s = calc_basic_bus_injection(network)[idx_sel_bus_types]
    p,q = real.(s),imag.(s)
    pf = [pf_i for (i,pf_i) in enumerate(pf) if abs(p[i])≥ϵ || abs(q[i]≥ϵ)]
    Δk = try 
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
    #k_max = maximum(k(network,sel_bus_types))
    k_max = calc_k_max(network,sel_bus_types)
    M = k_max*∂p∂θ - ∂q∂θ
    return M
end


"""
Given a network data dict,
Compute the maximum difference between nodal power factors so that complex power injections can be modeled from voltage magnitudes
"""
function calc_vmag_condition(network::Dict{String,<:Any},sel_bus_types=[1,2])
    ∂p∂θ = calc_pth_jacobian(network,sel_bus_types)
    M = try
        calc_M_matrix(network,sel_bus_types)
    catch
        throw(ArgumentError("No valid bus"))
    end
    pf = calc_basic_power_factor(network,sel_bus_types)
    k_max = calc_k_max(network,sel_bus_types)
    M_nonsingular = true
    Δk_max = try
        (1/(opnorm(inv(Matrix(M)),2)))*(1/(opnorm(Matrix(∂p∂θ),2)));
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


#Deprecated
# function calc_k_max(network::Dict{String,<:Any},sel_bus_types=[1,2])
#     idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types);
#     s = calc_basic_bus_injection(network)[idx_sel_bus_types];
#     p,q = real.(s),imag.(s)
#     θ = angle.(s)
#     k_max = 0
#     for (i,(p_i,q_i)) in enumerate(zip(p,q))
#         if abs(p_i)<1e-3
#             k_ii = 0
#         elseif abs(q_i)<1e-3
#             k_ii = 0
#         else
#             #pf_i = p_i/sqrt(p_i^2+q_i^2)
#             pf_i = cos(θ[i])
#             k_ii = k(pf_i)
#         end
#         if k_ii>k_max
#             k_max = k_ii     
#         end
#     end
#     return k_max
# end

 #p,q = real.(s),imag.(s)
    #pf = [p_i/(sqrt(p_i^2+q_i^2)) for (p_i,q_i) in zip(p,q)]