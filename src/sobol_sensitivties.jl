using GlobalSensitivity
using QuasiMonteCarlo
include("injections.jl")

function set_system_model(system_model::System)
    """Sets the global system model parameters"""
    if(system_model !== nothing)
        global system_model = model
    end
end

function make_sobol_sens_matrix(param_range::Vector,f::Function;
    order::Vector=[0,1,2],n_boot::Int64=100,conf_int::Float64=0.95)
    """
    Construncts an estimate for the global sensitivity matrix of the model
        Parameters:
        ---
        param_range::Vector ∈ R^N with injections ranges of real or reactive power for buses i=1,…,N
        f::Function Gets the voltage magnitudes for a given injection at bus i
        
    """
    if (system_model === nothing)
        println("No global system model, first call set_system_model(system_model::System)")
    end
    method = Sobol{
        order:order,
        nboot:n_boot,
        conf_int:conf_int
    }
    res = gsa(f,method,param_range; N, batch=false)
    return res
end

function make_vmag_sens_matrix()
    
end

function make_sens_matrix()

end

function make_design_matrices(n,inj_min::Vector,inj_max::Vector)
    A,B = QuasiMonteCarlo.generate_design_matrices(n,inj_min,inj_max,sampler,n)
    return A,B
end



