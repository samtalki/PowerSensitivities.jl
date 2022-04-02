""" 
Sets the load values in the net to the given load profile
specified in a 2*n_load x 1 array of n_load real then n_load
reactive load powers.
"""
function set_network_load(network::Dict{String,<:Any}, new_load; scale_load=true)
    if scale_load
		new_load = new_load * 1 / network["baseMVA"]
	end
    num_loads = length(network["load"])
    for load in values(network["load"])
        load_index = load["index"]
        load["pd"] = float(new_load[load_index][1])
        load["qd"] = float(new_load[load_index + num_loads][1])
	end
end

"""
Given a network data dict, calculate the "bad indeces" that have 0 injection
"""
function calc_bad_idx(network::Dict{String,<:Any},sel_bus_types=[1,2],ϵ=1e-6)
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    s = calc_basic_bus_injection(network)[idx_sel_bus_types];
    p,q,pf = real.(s),imag.(s),calc_basic_power_factor(network,sel_bus_types) 
    bad_idx = [] #Array of indeces with zero p or zero MVA injections to be discarded
    for (i,pf_i) in enumerate(pf)
        if abs(p[i]) ≤ ϵ && abs(q[i]) ≤ ϵ #If there is no apparent power injection, it doesn't make sense
            push!(bad_idx,i)
        elseif abs(pf_i) <= 1e-5 || pf_i == NaN || p[i] == 0.0 #If there is no real power injection, it doesn't make sense
            push!(bad_idx,i) #K[i,i] = 0
        #elseif q[i] < 0 #Capacitor banks
        #    push!(bad_idx,i)
        else
            continue
        end
    end
    return bad_idx
end

"""
Given a network data dict, check if the network is radial
"""
function is_radial(network::Dict{String,<:Any})
    
    function upper_off_diag(A::AbstractMatrix)
        [A[i] for i in CartesianIndices(A) if i[1] > i[2]]
    end 
    
    function lower_off_diag(A::AbstractMatrix)
        [A[i] for i in CartesianIndices(A) if i[1]<i[2]]
    end

    #Get admittance matrix and the sum of nonzero upper and lower off diagonal elements
    Y = calc_basic_admittance_matrix(network)
    nz_upper = sum([1 for y_ij in upper_off_diag(Y) if y_ij != 0])
    nz_lower = sum([1 for y_ij in lower_off_diag(Y) if y_ij != 0])
    n = size(Y)[1]
    if sum(nz_upper)>n-1 || sum(nz_lower) >n-1
        return false
    else
        return true
    end
end

