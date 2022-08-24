using LinearAlgebra:cond
using PowerModels:parse_file,make_basic_network
using PowerSensitivities:calc_voltage_sensitivity_matrix

function calc_condition_number(net::Dict)
    S = calc_voltage_sensitivity_matrix(net)
    Svp,Svq = S.vp,S.vq
    κ_vp,κ_vq = con(Svp),con(Svq)
    return [κ_vp,κ_vq]

function calc_condition_number()