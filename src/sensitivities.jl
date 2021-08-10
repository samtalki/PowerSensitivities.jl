using Dates: include
using Core: Vector
using PowerSystems
using TimeSeries
using Dates
using ForwardDiff
using DiffEqSensitivity
include("injections.jl")

DATA_DIR = download(PowerSystems.UtilsData.TestData, folder = pwd())

function parse_network_model(DATA_DIR,system="matpower/case5_re.m")
    system_data = System(joinpath(DATA_DIR,system))
    return system_data
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
        1+2
        #voltages_inj = injection_voltages()
    end
end

function sensitivity_mat_analytic(sensors,injections,system_data::System)
    """Analytical sensitivity matrix, like Christakou et al"""
    ybus = Ybus(system_data)
end

function sensitivity_global(inj_range_matrix,f_v=get_inj_range_voltages,method=Sobol,order=2,N=1000)
    """Numeric global sensitivities for a single inj_bus"""
    sens = gsa(f_v,method,inj_range_matrix; N, batch=false,order=order)
    return sens
end



function sensitivity_diffeq(sensors,injections,system_data::System)
    """Solves with DifferentialEquations.jl"""
    1+2
end

function sensitivity_perturb()

end

function installed_capacity(system::System; technology::Type{T} = Generator) where T <: Generator
    installed_capacity = 0.0
    for g in get_components(T, system)
        installed_capacity += get_max_active_power(g)
    end
    return installed_capacity
end

function inj_range_matrix(inj_min::Float64,inj_max::Float64,n_injections::Int)
    return diagonal([(inj_min,inj_max) for i in range(length(n_injections))])
end

