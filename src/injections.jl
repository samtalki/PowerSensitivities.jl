using PowerSystems
using TimeSeries
using Dates

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
    #load = PowerLoad(available=true,bus=bus,active_power=P_inj,reactive_power=Q_inj,name=name)
    add_component!(system_model,load)
    return load
end

function remove_inj(load::PowerLoad,system_model::System)
    remove_component!(system_model,load)
end

function get_inj_voltages(inj_bus::Bus,range_P_inj::Vector,range_Q_inj::Vector,system_model)
    """Gets the voltage mags for each measured bus for each value of the injection in inj_range at inj_bus"""
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

base_dir = PowerSystems.download(PowerSystems.TestData; branch = "master");
sys = System(joinpath(base_dir, "matpower/case14.m"));
vmag_0 = get_voltages(sys);
PQ_buses = get_PQ_buses(sys);
inj_bus = PQ_buses[1];
range_P_inj = [i for i in 0:0.01:0.9];
range_Q_inj = [0.0 for i in 0:0.01:0.9];
voltages_per_inj = get_inj_voltages(inj_bus,range_P_inj,range_Q_inj,sys);