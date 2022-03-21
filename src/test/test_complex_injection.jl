struct Sensitivities
    vp::AbstractMatrix
    vq::AbstractMatrix
    qp::AbstractMatrix
end

"""
Given a sensitivity struct calculcate an estimate for S†:
    S† = ∂v/∂p + ∂v/∂q * ∂q/∂p 
"""
est_sdag(s::Sensitivities) = s.vp + s.vq*s.qp

"""Solve a least squares problem"""
est_lsq(X,y,λ) = inv(X'*X + λ*I(size(X)[2]))*X'*y

"""
Given dictionary of deviations, estimate power-to-voltage sensitivities
"""
function est_sensitivities(diff_data::Dict,λ=1e-15)
	one_vec = ones(size(diff_data["dp"])[1])
	Δp,Δq,Δv = [diff_data["dp"] one_vec],[diff_data["dq"] one_vec],diff_data["dv"]
	svp = est_lsq(Δp,Δv,λ)
	svq = est_lsq(Δq,Δv,λ)
    sqp = est_lsq(Δp,Δq[:,1],λ)
	return Sensitivities(svp,svq,sqp)
end

"""
Given a network data dict, calculate the S† matrix.
"""
function calc_sdag_matrix(network,sel_bus_types)
    #Compute the analytical S̃ matrix
    sel_bus_idxs = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    J = PowerSensitivities.calc_jacobian_matrix(network,sel_bus_types)
    invJ = inv(Matrix(J.matrix)) #Calculate inverse of jacobian
    n = length(sel_bus_idxs)
    #Calculate S_tilde, SPV, and SQV
    S̃ = invJ[n+1:end,1:end] 
    SPV = S̃[:,1:n]
    SQV = S̃[:,n+1:end]
    #Calculate the S† = (Spv + SqvK) matrix
    K = PowerSensitivities.calc_K_matrix(network,sel_bus_types)
    Sdag = SPV + SQV*K
    return Sdag
end


"""
Given an AMI finite difference dataset, and estimated sensitivities, 
compute the performance of the inverse problems
"""
function est_deviations(ami_diff::Dict,s::Sensitivities)
    dp,dq,dv = ami_diff["dp"],ami_diff["dq"],ami_diff["dv"]
    s_dag = s.vp + s.vq*s.qp
    hat_dp,hat_dq,hat_dv = zeros(size(dp)),zeros(size(dq)),zeros(size(dv))
    for (t,(dp_t,dq_t,dv_t)) in enumerate(zip(eachrow(dp),eachrow(dq),eachrow(dv)))
        hat_dv[t,:] = s.vp*dp_t + s.vq*dq_t
        hat_dp[t,:] = inv(s_dag)*dv_t
        hat_dq[t,:] = s.qp*hat_dp_t
    end
    return Dict("dp"=>hat_dp,"dq"=>hat_dq,"dv"=>hat_dv)
end


"""
Given a file of test cases, calculate a test of the compelx power injection estimation. 
Does these things:
1. Generates the true S† matrices,
2. Calculates an AMI dataset
3. Computes the estimated S† matrix
4. Calculates the estimated Δx at each timestep from the AMI dataset.
"""
function test_complex_injection_est(diff_datasets,sel_bus_types=[1,2],network_data_path=network_data_path)
    names = readdir(network_data_path);
    paths = readdir(network_data_path,join=true);
    results = Dict()
    for (name,path) in zip(names,paths)
        network =  try
            make_basic_network(parse_file(path));
        catch
            println("PM cannot parse "*name)
            continue
        end
        if PowerSensitivities.is_radial(network)
            sdag_true = calc_sdag_matrix(network,sel_bus_types)
            ami_diff = diff_datasets[name] 
            s = est_sensitivities(ami_diff)
            deviations = est_deviations(ami_diff,s)
            results[name] = Dict(
                "sdag_true"=>sdag_true,
                "hat_sdag"=>est_sdag(s),
                "hat_s"=>s,
                "deviations"=>deviations)
        end
    end
    return results
end

"""
Given path of network files, make ami and ami_diff datasets
"""
function make_ami_datasets(network_data_path=network_data_path,sel_bus_types=[1])
    names = readdir(network_data_path);
    paths = readdir(network_data_path,join=true);
    ami_datasets,diff_datasets = Dict(),Dict()
    for (name,path) in zip(names,paths)
        network =  try
            make_basic_network(parse_file(path));
        catch
            println("PM cannot parse "*name)
            continue
        end
        if PowerSensitivities.is_radial(network)
            ami_dataset = try
                PowerSensitivities.make_ami_dataset(network,sel_bus_types)
                ami_diff = PowerSensitivities.calc_finite_differences(ami_dataset) 
                ami_datasets[name] = ami_dataset
                diff_datasets[name] = ami_diff
            catch
                continue
            end
        end
    end
    return ami_datasets,diff_datasets
end

#Test the complex power injection estimation
#complex_est_results_pq = test_complex_injection_est([1])
#complex_est_results_pq_pv = test_complex_injection_est([1,2])
