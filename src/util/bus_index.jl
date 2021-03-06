"""
Given a network data dict, calculate the "bad indeces" that do not satisfy the assumptions of Theorem 1
"""
function calc_bad_idx(data::Dict{String,<:Any},ϵ=1e-6)
    s = calc_basic_bus_injection(data);
    pnet,qnet,pf = real.(s),imag.(s),calc_basic_power_factor(data) 
    bad_idx = [] #Array of indeces with zero pnet or zero MVA injections to be discarded
    for (i,pf_i) in enumerate(pf)
        #Ignore buses with zero power injection
        if abs(pnet[i]) ≤ ϵ && abs(qnet[i]) ≤ ϵ 
            push!(bad_idx,i)
        #Ignore buses with zero power factor
        elseif abs(pf_i) <= 1e-5 || pf_i == NaN || pnet[i] == 0.0
            push!(bad_idx,i) #K[i,i] = 0
        #Ignore buses with capacitive injections
        elseif qnet[i] > 0 
           push!(bad_idx,i)
        else
            continue
        end
    end
    return bad_idx
end


"""
Given a network data dict, calculate the bus idxs that meet the following conditions:
    1. Only the selected bus types
    2. Nonzero apparent power injections (with tolerance ϵ)
    3. Nonzero real power injection
    4. Inductive injection 
"""
function calc_study_idx(data::Dict{String,<:Any},sel_bus_types=[1,2],ϵ=1e-6)
    n_bus = length(data["bus"])
    #First, compute the index of bus types
    idx_sel_bus_types = calc_bus_idx_of_type(data,sel_bus_types)
    #Compute the good idx amoung all buses.
    idx_good_bus = [i for i in 1:n_bus if i ∉ calc_bad_idx(data,ϵ)]
    #Compute the indeces that are of the selected type and satify all of the assumptions.
    return [i for i in 1:n_bus if i ∈ idx_good_bus || i ∈ idx_sel_bus_types]
end

"""
Given a network data dict, calculate the indeces that are net-inductive
"""
function calc_net_inductive_idx(data::Dict{String,<:Any},ϵ=1e-6)
    s = calc_basic_bus_injection(data);
    pnet,qnet = real.(s),imag.(s)
    inductive_idx = [] #Array of indeces with inductive net injections
    for (i,qᵢ) in enumerate(qnet)
        if qᵢ < 0 
           push!(inductive_idx,i)
        else
            continue
        end
    end
    return inductive_idx
end

"""
Given a network data dict, calculate the indeces that are net-inductive
"""
function calc_net_capacitive_idx(data::Dict{String,<:Any},ϵ=1e-6)
    s = calc_basic_bus_injection(data);
    pnet,qnet = real.(s),imag.(s)
    capacitive_idx = [] #Array of indeces with inductive net injections
    for (i,qᵢ) in enumerate(qnet)
        if qᵢ > 0 
           push!(capacitive_idx,i)
        else
            continue
        end
    end
    return capacitive_idx
end


"""
Given a network network dict, find the bus indeces that have bus_types in sel_bus_types .
"""
function calc_bus_idx_of_type(data,sel_bus_types)
    bus_ordered = get_bus_ordered(data)
    bus_types = [bus["bus_type"] for bus in bus_ordered]
    idx_sel_bus_types = findall(bus_idx-> bus_idx ∈ sel_bus_types,bus_types)
    return idx_sel_bus_types
end

"""
Given a network network dict, return the bus dictionaries that have bus_types in sel_bus_types.
"""
function get_bus_of_type(data,sel_bus_types)
    bus_of_type = filter(d->d.second["bus_type"]sel_bus_types,data["bus"])
    return bus_of_type
end

"""
Get ordered vector of buses
"""
function get_bus_ordered(data)
    b = [bus for (i,bus) in data["bus"] if bus["bus_type"] != 4]
    bus_ordered = sort(b, by=(x) -> x["index"])
    return bus_ordered
end

function get_bus_ordered(data,sel_bus_types)
    b = [bus for (i,bus) in data["bus"] if bus["bus_type"]∈sel_bus_types]
    bus_ordered = sort(b, by=(x) -> x["index"])
    return bus_ordered
end

