"""
Given a network network dict, find the bus indeces that have bus_types in sel_bus_types .
"""
function calc_bus_idx_of_type(network,sel_bus_types)
    bus_ordered = get_bus_ordered(network)
    bus_types = [bus["bus_type"] for bus in bus_ordered]
    idx_sel_bus_types = findall(bus_idx-> bus_idx ∈ sel_bus_types,bus_types)
    return idx_sel_bus_types
end

"""
Given a network network dict, return the bus dictionaries that have bus_types in sel_bus_types.
"""
function get_bus_of_type(network,sel_bus_types)
    bus_of_type = filter(d->d.second["bus_type"]sel_bus_types,network["bus"])
    return bus_of_type
end

"""
Get ordered vector of buses
"""
function get_bus_ordered(network)
    b = [bus for (i,bus) in network["bus"] if bus["bus_type"] != 4]
    bus_ordered = sort(b, by=(x) -> x["index"])
    return bus_ordered
end

function get_bus_ordered(network,sel_bus_types)
    b = [bus for (i,bus) in network["bus"] if bus["bus_type"]∈sel_bus_types]
    bus_ordered = sort(b, by=(x) -> x["index"])
    return bus_ordered
end

"""
Given a network data dict, calculate the "bad indeces" that do not satisfy the assumptions of Theorem 1
"""
function calc_bad_idx(network::Dict{String,<:Any},ϵ=1e-6)
    s = calc_basic_bus_injection(network);
    p,q,pf = real.(s),imag.(s),calc_basic_power_factor(network) 
    bad_idx = [] #Array of indeces with zero p or zero MVA injections to be discarded
    for (i,pf_i) in enumerate(pf)
        #Ignore buses with zero power injection
        if abs(p[i]) ≤ ϵ && abs(q[i]) ≤ ϵ 
            push!(bad_idx,i)
        #Ignore buses with zero power factor
        elseif abs(pf_i) <= 1e-5 || pf_i == NaN || p[i] == 0.0
            push!(bad_idx,i) #K[i,i] = 0
        #Ignore buses with capacitive injections
        elseif q[i] > 0 
           push!(bad_idx,i)
        else
            continue
        end
    end
    return bad_idx
end

"""
Given a network data dict and !SELECTED BUS TYPES!, calculate the "bad indeces" that do not satisfy the assumptions of Theorem 1
"""
function calc_bad_idx(network::Dict{String,<:Any},sel_bus_types=[1,2],ϵ=1e-6)
    idx_sel_bus_types = calc_bus_idx_of_type(network,sel_bus_types)
    n_selected = length(idx_sel_bus_types) #Number of selected buses
    s = calc_basic_bus_injection(network)[idx_sel_bus_types];
    p,q,pf = real.(s),imag.(s),calc_basic_power_factor(network,sel_bus_types) 
    @assert length(s) == n_selected && length(pf) == n_selected #Check the lengths are consistent
    bad_idx = [] #Array of indeces with zero p or zero MVA injections to be discarded
    for i in 1:n_selected
        pf_i,p_i,q_i = pf[i],p[i],q[i]
        #Ignore buses with zero power injection
        if abs(p_i) ≤ ϵ && abs(q_i) ≤ ϵ 
            push!(bad_idx,i)
        #Ignore buses with zero power factor
        elseif abs(pf_i) <= 1e-5 || pf_i == NaN || p_i == 0.0
            push!(bad_idx,i) #K[i,i] = 0
        #Ignore buses with capacitive injections
        elseif q_i > 0 
           push!(bad_idx,i)
        else
            continue
        end
    end
    return bad_idx
end