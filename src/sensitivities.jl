using Core: Vector
using PowerSystems
using TimeSeries
using Dates
using ForwardDiff

DATA_DIR = download(PowerSystems.UtilsData.TestData, folder = pwd())

function parse_network_model(DATA_DIR,system="matpower/case5_re.m")
    system_data = System(joinpath(DATA_DIR,system))
    return system_data
end

function injection_voltages(inj_bus,sensors::Vector,inj_range::Vector,system_data::System)
    """Gets the voltage measurements at sensors for a specified inj_range at inj_bus"""
    1+2
end

function sensitivity(sensors::Vector,inj_range::Vector)
    """Computes the gradient of the nodal voltages w.r.t. the inj_range"""
    dvdx = ForwardDiff.gradient(f,x) #g = ∇f
end

function sens_mat_perturb(sensors,injections,system_data::System)
    """Makes a perturb-and-observe sensitivity matrix for a model"""
    M = size(sensors)
    L = size(injections)
    S = zeros(M,L)
    for (sensor,injection) in zip(sensors,injections)
        voltages_inj = injection_voltages()
    end
end

function sens_mat_analytical(sensors,injections,system_data::System)
    """Analytical sensitivity matrix, like Christakou et al"""
    ybus = Ybus(system_data)
end

function diffeq_sensitivities(sensors,injections,system_data::System)
    """Solves with DifferentialEquations.jl"""
    1+2
end


function installed_capacity(system::System; technology::Type{T} = Generator) where T <: Generator
    installed_capacity = 0.0
    for g in get_components(T, system)
        installed_capacity += get_max_active_power(g)
    end
    return installed_capacity
end