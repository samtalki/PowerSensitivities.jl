module PowerSensitivities

#Imports
using Base
using PowerSystems
using PowerModels
using Ipopt
using DifferentialEquations
using LinearAlgebra
using DataFrames
using TimeSeries 
using Dates 
using GlobalSensitivity
using QuasiMonteCarlo


#Exports
export LinearSensitivityModel
export get_voltage_phasors_rectangular,get_voltage_magnitudes
export vph_p_sens,vph_q_sens
export RealPerturb,ReactivePerturb
export vph_p_sens,vph_q_sens
export system_model_path,system_model,PQ_buses,set_system_model,check_system_model #system model interface
export calc_basic_jacobian_matrix,calc_basic_decoupled_jacobian_matrices #Analytical AC power flow jacobian

#Structs 
#Real power purturbation
mutable struct RealPerturb
    inj_state::AbstractArray
end

#Reactive power perturbation
mutable struct ReactivePerturb
    inj_state::AbstractArray
end

abstract type SensitivityModel end
#General sensitivity model
mutable struct LinearSensitivityModel <: SensitivityModel
    system_model::System
    pq_buses::AbstractArray
    slack_buses::AbstractArray
    ybus::AbstractArray
    vph_0::AbstractArray
    powerflow_results

    function LinearSensitivityModel(system_model::System)
        pq_buses = get_pq_buses(system_model)
        slack_buses = get_slack_buses(system_model)
        ybus = Ybus(system_model).data
        results = solve_powerflow(system_model)
        vph_0 = get_voltage_phasors_rectangular(results)
        return new(system_model,pq_buses,slack_buses,ybus,vph_0,results)
    end
end


#Module scripts
include("injections.jl")
include("sensitivities.jl")
include("analytical_sensitivities.jl")
include("jacobian.jl")

#System model global variables
global system_model = nothing
global PQ_buses = nothing
global system_model_path = nothing

function set_system_model(system_model::System)
    """Sets the global system model parameters"""
    if(system_model !== nothing)
        println("Setting global system model")
        global system_model = model
    end
end

function set_system_model(data_folder_path::String="matpower/case14.m")
    base_dir = PowerSystems.download(PowerSystems.TestData; branch = "master");
    if !@isdefined(system_model) || (system_model===nothing)
        global system_model_path = joinpath(base_dir, data_folder_path)
        global system_model = System(system_model_path);
        global PQ_buses = get_PQ_buses(system_model);
    end
end

function check_system_model(modelx::System)
    if(system_model !== nothing) && (modelx==system_model)
        println("Global System model PASSED")
    else
        set_system_model(modelx)
    end
end

end