include("PowerSensitivities.jl")
using .PowerSensitivities
using PowerModels
using LinearAlgebra
using JuMP
using Ipopt

"""
Checks if a matrix M is positive definite
"""
ispd(M) = all([real(eig)>0 for eig in eigvals(M)])

"""
Checks if a matrix M is negative definite
"""
isnd(M) = all([real(eig)<0 for eig in eigvals(M)])
isnsd(M,ϵ=1e-12) = all([real(eig)<=ϵ for eig in eigvals(M)])

"""
Checks if a matrix M is invertible
"""
isinvertible(x) = applicable(inv, x)

"""
Returns distance between M and M transpose
"""
symmetricdiff(M) = norm(M-transpose(M))

"""
Check if symmetric part of a matrix is negative definite
"""
symmetric_part_nd(M) = isnd(0.5.*(M + transpose(M)))
symmetric_part_nsd(M) = isnsd(M + transpose(M./2))

"""
Check if symmetric part of a matrix is positive definite
"""
symmetric_part_pd(M) = ispd(0.5*(M +transpose(M)))


"""
Compute K matrix where K = diag(√(1-pf_i^2)/pf_i)
"""
function calc_K_matrix(network,sel_bus_types=[1,2])
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    s = calc_basic_bus_injection(network)[idx_sel_bus_types];
    p,q = real.(s),imag.(s)
    pf = [cos(atan(q_i/p_i)) for (p_i,q_i) in zip(p,q)]
    n = length(pf)
    K = zeros((n,n))
    for (i,pf_i) in enumerate(pf)
        if(abs(pf_i) <= 1e-5 || pf_i == NaN || p[i] == 0.0)
            K[i,i] = 0
        else
            K[i,i] = sqrt(1-pf_i^2)/pf_i
        end
    end
    return K
end

"""
Given network data dict, calculate the maximum difference between the nodal power factors
"""
function calc_k_max(network,sel_bus_types=[1,2])
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
        if(k_ii>k_max)
            k_max = k_ii     
        end
    end
    return k_max
end

"""
Given a network data dict, calculate the nodal power factors
"""
function calc_basic_power_factor(network,sel_bus_types=[1,2])
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    s = calc_basic_bus_injection(network)[idx_sel_bus_types]
    p,q = real.(s),imag.(s)
    pf = [p_i/(sqrt(p_i^2+q_i^2)) for (p_i,q_i) in zip(p,q)]
    return pf
end

"""
Calculate the M Matrix for theorem 1
"""
function calc_M_matrix(network,sel_bus_types=[1,2])
    ∂p∂θ = calc_pth_jacobian(network,sel_bus_types)
    ∂q∂θ = calc_qth_jacobian(network,sel_bus_types)
    k_max = calc_k_max(network,sel_bus_types)
    M = k_max*∂p∂θ - ∂q∂θ
    return M
end

"""
Calculate the maximum power factor distance given a maximum difference Δk_max between the Q-P sensitivities 
"""
function calc_max_pf_distance(Δk_max)
    model = Model(Ipopt.Optimizer)
    @variable(model, 1 >= pf_min >= 0)
    @variable(model, 1>= pf_max >= 0)
    @objective(model, Max, pf_max - pf_min)
    @NLconstraint(model, (sqrt(1-pf_min^2)/pf_min) - (sqrt(1-pf_max^2)/pf_max) <= Δk_max  )
    optimize!(model)
    pf_min,pf_max = value(pf_min),value(pf_max)
    Δpf_max = pf_max - pf_min
    return Δpf_max
end

"""
Given a network data dict,
Calculate the maximum power factor distance
"""
function calc_max_pf_distance(network::Dict{String,<:Any})
    a=1+2
end

"""
Given a network data dict,
Compute the maximum difference between nodal power factors so that complex power injections can be modeled from voltage magnitudes
"""
function calc_thm1_condition(network,sel_bus_types=[1,2])
    ∂p∂θ = calc_pth_jacobian(network,sel_bus_types)
    M = calc_M_matrix(network,sel_bus_types)
    pf = calc_basic_power_factor(network,sel_bus_types)
    k_max = calc_k_max(network,sel_bus_types)
    M_nonsingular = true
    Δk_max =try
        (1/(norm(inv(M))))*(1/(norm(∂p∂θ)))
    catch
        Δk_max = nothing
        M_nonsingular = false
    end
    Δpf_max = calc_max_pf_distance(Δk_max)
    return Dict(
        "M_nonsingular" => M_nonsingular,
        "Δk_max" => Δk_max,
        "Δpf_max" => Δpf_max,
        "k_max" => k_max,
        "M" => M,
        "pf" => pf
    )
end

function test_thm1(sel_bus_types=[1,2],network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/")
    names = readdir(network_data_path);
    paths = readdir(network_data_path,join=true);
    results = Dict()
    for (name,path) in zip(names,paths)
        network = make_basic_network(parse_file(path));
        results[name] = calc_thm1_condition(network,sel_bus_types)
    end
    return results
end


"""
Test assumption 1
"""
function test_assumption1(sel_bus_types=[1,2],network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/")
	names = readdir(network_data_path);
    paths = readdir(network_data_path,join=true);
    results = Dict()
	full_J_matrices = Dict()
    for (name,path) in zip(names,paths)
        data = make_basic_network(parse_file(path));
		J = try
            calc_jacobian_matrix(data,sel_bus_types);
        catch
            println("Jacobian cannot be computed for "*name)
            continue
        end
		Y = Matrix(calc_basic_admittance_matrix(data))
		Jmat = Matrix(J.matrix)
		spth = Matrix(J.spth)
        sqth = Matrix(J.sqth)
		info = Dict(
			"n_slack_node" => length(findall(d -> d["bus_type"] ==3, data["bus"])),
			"slack_nodes" => findall(d -> d["bus_type"] ==3, data["bus"]),
            "dpdth_pd" => ispd(spth),
			"y_symmetric" => issymmetric(Y),
			"dpdth_symmetric" => issymmetric(spth),
            "J_nonsingular" => isinvertible(Jmat),
            "sym_dpdth_pd" => symmetric_part_pd(spth),
            "sym_dqdth_nsd" => symmetric_part_nsd(sqth),
            #"norm(spth-spth')" => symmetricdiff(spth),
			#"rel(spth-spth')" => symmetricdiff(spth)/norm(spth)
        )
		full_J_matrices[name] = calc_basic_jacobian_matrix(data)
        results[name] = info
    end
    return results,full_J_matrices
end

#Test assumption 1
results_all,J_all = test_assumption1([1,2,3])
results_pq_pv,J_pq_pv = test_assumption1([1,2])
results_pq,J_pq = test_assumption1([1])


#Test Theorem 1 bounds
thm1_pq = test_thm1([1])
thm1_pq_pv = test_thm1([1,2])

#Get the maximum power factor distances
delta_pf_max_pq,delta_pf_max_pq_pv = Dict(),Dict()
for (name,d) in thm1_pq
    delta_pf_max_pq[name] = d["Δpf_max"]
end
for (name,d) in thm1_pq_pv
    delta_pf_max_pq_pv[name] = d["Δpf_max"]
end

# Auxillary data for testing
case24 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case24.m"))
case14 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case14.m"))
J_c24_full = calc_basic_jacobian_matrix(case24)
J_c24_blocks = calc_jacobian_matrix(case24)
J_c24_blocks_pq_pv = calc_jacobian_matrix(case24,[1,2])
J_c24_blocks_pq = calc_jacobian_matrix(case24,[1])
