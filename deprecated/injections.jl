function get_inj_voltages(injection::RealPerturb;system::System=system_model)
    """
    Gets a vector of voltage measurements for a real power injection state vector ∈ R^L
    """
    check_system_model(system)
    for (inj_bus,inj) in enumerate(injection.inj_state)
        if(inj == 0)
            continue
        elseif inj != 0
            println("P Inj: ",string(inj),"On bus: ",string(inj_bus))
            voltages = get_p_inj_voltages(inj_bus,inj,system_model=system)
        end
    end
    return voltages
end


function get_inj_voltages(injection::ReactivePerturb;system::System=system_model)
    """
    Gets a vector of voltage measurements for a real power injection state vector ∈ R^L
    """
    check_system_model(system)
    for (inj_bus,inj) in enumerate(injection.inj_state)
        if(inj == 0)
            continue
        elseif inj != 0
            println("Q Inj: ",string(inj),"On bus: ",string(inj_bus))
            voltages = get_q_inj_voltages(inj_bus,inj,system)
        end
    end
    return voltages
end

function get_PQ_buses(sys::System)
    buses = collect(get_components(Bus,sys))
    PQ_buses = [bus for bus in buses if get_bustype(bus) == BusTypes.PQ]
    return PQ_buses
end

#Get voltage phasors from all observable buses in system
function get_voltage_phasors_rectangular(system_model::System=system_model)
    results = solve_powerflow(system_model)
    Vm = results["bus_results"][!,"Vm"]
    θ = results["bus_results"][!,"θ"]
    results["bus_results"][!,"Vph"] = get_rectangular_from_phasor(Vm,θ)
    return results["bus_results"][!,"Vph"]
end

function get_voltage_phasors_rectangular(results::Dict)
    Vm = results["bus_results"][!,"Vm"]
    θ = results["bus_results"][!,"θ"]
    results["bus_results"][!,"Vph"] = get_rectangular_from_phasor(Vm,θ)
    return results["bus_results"][!,"Vph"]
end

#Get voltage magnitudes from all observable buses in system
function get_voltage_magnitudes(system_model::System=system_model)
    results = solve_powerflow(system_model)
    return results["bus_results"][!,"Vm"]
end

#Get voltage magnitudes from all observable buses in system
function get_voltage_magnitudes(bus_results::DataFrame)
    return bus_results[!,"Vm"]
end

#Convert phasors to rectangular
function get_rectangular_from_phasor(Vm,θ)
    return Vm*(cos(θ)+im*sin(θ))
end

function get_rectangular_from_phasor(Vm::AbstractArray,θ::AbstractArray)
    n = length(Vm)
    return [Vm[i]*(cos(θ[i])+im*sin(θ[i])) for i in 1:n]
end




function place_inj(P_inj::Float64,Q_inj::Float64,bus::Bus,name::String;system::System=system_model)
    check_system_model(system)
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


function place_inj(P_inj::Float64,Q_inj::Float64,inj_bus::Int,name::String;system::System=system_model)
    check_system_model(system)
    load = PowerLoad(
        name = name,
        available = true,
        bus = collect(get_components(Bus,system_model))[inj_bus],
        model = LoadModels.ConstantPower,
        active_power = P_inj,
        reactive_power = Q_inj,
        base_power = get_base_power(system_model),
        max_active_power = P_inj,
        max_reactive_power = Q_inj,
    )
    add_component!(system_model,load)
    return load
end

function remove_inj(load::PowerLoad,system_model::System)
    remove_component!(system_model,load)
end




function get_p_inj_voltages(inj_bus::Number,P_inj::Number;system_model::System=system_model)
    """Gets the vector of voltage magnitudes for an injection on bus num inj_bus of value P_inj"""
    Q_inj = 0.0
    load = place_inj(P_inj,Q_inj,inj_bus,string(inj_bus),system=system_model)
    voltages = get_voltage_magnitudes(system_model)
    remove_inj(load,system_model)
    return voltages
end

function get_p_inj_voltages_by_bus(p_inj_by_bus::Vector{Float64};system_model=system_model)
    check_system_model(system_model)
    q_inj=0.0
    n_measurement_columns = length(p_inj_by_bus)
    loads = []
    println(p_inj_by_bus)
    for (inj_bus,p_inj) in enumerate(p_inj_by_bus)
        println(inj_bus,p_inj)
        load = place_inj(p_inj,q_inj,inj_bus,string(inj_bus),system=system_model)
        push!(loads,load)
    end
    voltages = get_voltage_magnitudes(system_model)
    for load in loads
        remove_inj(load,system_model)
    end
    return voltages
end




function get_q_inj_voltages(inj_bus::Int,Q_inj::Int64,system_model::System)
    """Gets the vector of voltage magnitudes for an injection on bus num inj_bus of value P_inj"""
    P_inj = 0
    load = place_inj(P_inj,Q_inj,inj_bus,string(inj_bus),system=system_model)
    voltages = get_voltage_magnitudes(system_model)
    remove_inj(load,system_model)
    return voltages
end



function get_inj_range_voltages(inj_bus::Bus,range_P_inj::Vector,range_Q_inj::Vector,system_model)
    """Gets the voltage mags for each measured bus for each value of the injection in P and Q inj_range at inj_bus"""
    voltages_per_inj = []
    for (idx,(P_inj,Q_inj)) in enumerate(zip(range_P_inj,range_Q_inj))
        println(string(idx),typeof(string(idx)))
        load = place_inj(P_inj,Q_inj,inj_bus,string(idx),system=system_model)
        voltages = get_voltage_magnitudes(system_model)
        remove_inj(load,system_model)
        push!(voltages_per_inj,voltages)
    end
    return voltages_per_inj
end

function get_P_voltage_matrix(inj_matrix)
    voltage_matrix = zeros(size(inj_matrix))
    for (i,row) in enumerate(eachrow(inj_matrix))
        voltage_matrix[i,!] = get_P_inj_voltages(i,row[i],sys)
    end
    return voltage_matrix
end



function inj_range_matrix(inj_min,inj_max,n_injections)
    println("entered inj_range_matrix")
    A = map(tuple, (0.0 for i=1:n_injections,j=1:n_injections),zeros(n_injections,n_injections))
    for i=1:n_injections,j=1:n_injections
        if i==j
            #println(i,j)
            A[i,j] = (inj_min,inj_max) 
        end
    end
    return A
    #Diagonal([(inj_min,inj_max) for i in 1:n_injections])
end


# run_1 = false
# if run_1
#     @time begin
#         vmag_0 = get_voltage_magnitudes(sys);        
#         inj_bus = PQ_buses[1];
#         range_P_inj = [i for i in 0:0.001:1.0];
#         range_Q_inj = [0.0 for i in 0:0.001:1.0];
#         voltages_per_inj = get_inj_range_voltages(inj_bus,range_P_inj,range_Q_inj,sys);
#     end
# end

#GSA Experiment (Sobol Indeces)


# P_inj_min = [0.1, 0.1, 0.1]
# P_inj_max = [0.5,0.5,0.5]
# N=1000
# sampler = SobolSample()
# A,B = QuasiMonteCarlo.generate_design_matrices(N,P_inj_min,P_inj_max,sampler)

# println(P_inj_min,P_inj_max)
# P_inj_params = inj_range_matrix(P_inj_min,P_inj_max,length(PQ_buses))
# sens = sensitivity_global(P_inj_params)
# #gsa_results = DataFrame(Inj_Bus = PQ_buses, inj_range = inj_range)



