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

