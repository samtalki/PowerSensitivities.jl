using GlobalSensitivity
using QuasiMonteCarlo
using PowerSystems

function sobol_sens_matrix(model::System,n_order,n_boot,conf_int=0.95)
    """Construncts an estimate for the global sensitivity matrix of the model"""
    sampler = SobolSample()
    A,B = QuasiMonteCarlo.generate_design_matrices(n,lb,ub,sampler)
    res = gsa(get_voltages,Sobol(order=[0,1,2]),A,B)
    return res

end

function make_vmag_sens_matrix()
    
end

function make_sens_matrix()
end

function make_design_matrices(n,inj_min::Vector,inj_max::Vector)
    A,B = QuasiMonteCarlo.generate_design_matrices(n,inj_min,inj_max,sampler,n)

end

