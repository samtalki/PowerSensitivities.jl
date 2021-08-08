using PowerSystems
using TimeSeries
using Dates
using DataFrames
using LinearAlgebra
#using PowerSensitivities : inj_range_matrix
#using PowerSensitivities : sensitivity_global
#include("sensitivities.jl")

function get_voltages(system_model::System=nothing)
    results = solve_powerflow(system_model)
    return results["bus_results"][!,"Vm"]
end

function get_PQ_buses(sys::System)
    buses = collect(get_components(Bus,sys))
    PQ_buses = [bus for bus in buses if get_bustype(bus) == BusTypes.PQ]
    return PQ_buses
end

function get_voltages(bus_results)
    return bus_results["Vm"]
end

function place_inj(P_inj::Float64,Q_inj::Float64,bus::Bus,name,system_model::System)
    load = PowerLoad(
        name = name,
        available = true,
        bus = bus,
        model = LoadModels.ConstantPower,
        active_power = P_inj,
        reactive_power = Q_inj,
        base_power = get_base_power(system_model),
        max_active_power = P_inj,
        max_reactive_power = Q_inj,
        #services = Vector{Service}
       # dynamic_injector::Union{Nothing, DynamicInjection}
       # ext::Dict{String, Any}
       # time_series_container::InfrastructureSystems.TimeSeriesContainer
       # internal::InfrastructureSystemsInternal
    )
    add_component!(system_model,load)
    return load
end

function place_inj(P_inj::Float64,Q_inj::Float64,inj_bus::Int,name::String,system_model::System)
    load = PowerLoad(
        name = name,
        available = true,
        bus = get_component(Bus,system_model, string(inj_bus)),
        model = LoadModels.ConstantPower,
        active_power = P_inj,
        reactive_power = Q_inj,
        base_power = get_base_power(system_model),
        max_active_power = P_inj,
        max_reactive_power = Q_inj,
    )
    add_component!(system_model,load)
end

function remove_inj(load::PowerLoad,system_model::System)
    remove_component!(system_model,load)
end

function get_inj_range_voltages(inj_bus::Bus,range_P_inj::Vector,range_Q_inj::Vector,system_model)
    """Gets the voltage mags for each measured bus for each value of the injection in P and Q inj_range at inj_bus"""
    voltages_per_inj = []
    for (idx,(P_inj,Q_inj)) in enumerate(zip(range_P_inj,range_Q_inj))
        println(string(idx),typeof(string(idx)))
        load = place_inj(P_inj,Q_inj,inj_bus,string(idx),system_model)
        voltages = get_voltages(system_model)
        remove_inj(load,system_model)
        push!(voltages_per_inj,voltages)
    end
    return voltages_per_inj
end

function get_P_inj_voltages(inj_bus::Int,P_inj::Int64,system_model::System)
    """Gets the vector of voltage magnitudes for an injection on bus num inj_bus of value P_inj"""
    Q_inj = 0
    load = place_inj(P_inj,Q_inj,inj_bus,"Inj "+string(inj_bus),system_model)
    voltages = get_voltages(system_model)
    remove_inj(load,system_model)
    return voltages
end


function get_P_voltage_matrix(inj_matrix)
    voltage_matrix = zeros(size(inj_matrix))
    for (i,row) in enumerate(eachrow(inj_matrix))
        voltage_matrix[i,!] = get_P_inj_voltages(i,row[i],sys)
    end
    return voltage_matrix
end



function inj_range_matrix(inj_min::Float64,inj_max::Float64,n_injections)
    #A = Array{Tuple{Float64,Float64}}(undef,n_injections,n_injections)
    A = map(tuple, (0.0 for i=1:n_injections,j=1:n_injections),zeros(n_injections,n_injections))

    #!fill(A,(0,0))
    #A = Array{Tuple}(nothing,n_injections,n_injections)
    for i=1:n_injections,j=1:n_injections
        if i==j
            print(i,j)
            A[i,j] = (inj_min,inj_max) 
        end
    end
    return A
    #Diagonal([(inj_min,inj_max) for i in 1:n_injections])
end

function sensitivity_global(inj_range_matrix,f_v=get_inj_voltages,method=Sobol,order=2,N=1000)
    """Numeric global sensitivities for a single inj_bus"""
    sens = gsa(f_v,method,inj_range_matrix; N, batch=false,order=order)
    return sens
end


base_dir = PowerSystems.download(PowerSystems.TestData; branch = "master");
sys = System(joinpath(base_dir, "matpower/case14.m"));
@time begin
    vmag_0 = get_voltages(sys);
    PQ_buses = get_PQ_buses(sys);
    inj_bus = PQ_buses[1];
    range_P_inj = [i for i in 0:0.001:1.0];
    range_Q_inj = [0.0 for i in 0:0.001:1.0];
    voltages_per_inj = get_inj_range_voltages(inj_bus,range_P_inj,range_Q_inj,sys);
end

#GSA Experiment (Sobol Indeces)
@time begin
    P_inj_min = 0.0001
    P_inj_max = 0.5
    P_inj_params = inj_range_matrix(P_inj_min,P_inj_max,length(PQ_buses))
    #gsa_results = DataFrame(Inj_Bus = PQ_buses, inj_range = inj_range)
end

