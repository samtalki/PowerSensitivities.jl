using LinearAlgebra
using PowerModels
using PowerSensitivities

function calc_condition_number(net::Dict)
    S = calc_voltage_sensitivity_matrix(net)
    Svp,Svq = S.vp,S.vq
    κ_vp,κ_vq = cond(Svp),cond(Svq)
    return [κ_vp,κ_vq]
end

function calc_condition_number(net::Dict,sel_bus_types=[1])
    S = calc_voltage_sensitivity_matrix(net,sel_bus_types)
    Svp,Svq = S.vp,S.vq
    κ_vp,κ_vq = cond(Svp),cond(Svq)
    return [κ_vp,κ_vq]
end 

function calc_condition_number(file::String)
    net = make_basic_network(parse_file(path))
    return calc_condition_number(net)
end 

function calc_condition_number(files::Vector{String})
    conds = Dict()
    for file in files
        net = make_basic_network(parse_file(file))
        name = net["name"]
        conds[name] = calc_condition_number(net)
    end
    return conds
end 
