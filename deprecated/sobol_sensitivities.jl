function make_sobol_sens_matrix(inj_range::Tuple,L_nodes::Int64,inj_type::String;
    N_simulations::Int64=1000,order::Vector=[0,1,2],n_boot::Int64=100,conf_int::Float64=0.95)
    #Matrix whose rows are a vector with L entries, with only one nonzero element containing tuples for the bus to be measured
    param_range_matrix = inj_range_matrix(inj_range[0],inj_range[1],L_nodes)
    for row in eachrow(param_range_matrix)
        param_range = Vector(row)
        if(occursin("p",inj_type) || occursin("P",inj_type))
            results = make_sobol_sens_vector(param_range,get_P_inj_voltages,N_simulations;
            order=order,n_boot=n_boot,conf_int=conf_int)
        elseif occursin("q",inj_type) || occursin("Q",inj_type)
            results = make_sobol_sens_vector(param_range,get_Q_inj_voltages,N_simulations;
            order=order,n_boot=n_boot,conf_int=conf_int)
        end
    end
    return sensitivity_matrix
end

function make_sobol_sens_vector(param_range::Vector,f::Function,N_simulations::Int64;
    order=[0,1,2],n_boot=5,conf_int=0.95)
    """
    Construncts an estimate for the global sensitivity matrix of the model
        Parameters:
        ---
        param_range::Vector ∈ R^L with injections ranges of real or reactive power for buses i=1,…,N
        f::Function Gets the voltage magnitudes for a given injection at bus i
        
    """
    if (system_model === nothing)
        println("No global system model, first call set_system_model(system_model::System)")
    end
    method = Sobol(order=order,nboot=n_boot,conf_level=conf_int)
    res = gsa(f,method,param_range; N=N_simulations, batch=false)
    return res
end


function make_sens_matrix()

end

function make_design_matrices(n,inj_min::Vector,inj_max::Vector)
    A,B = QuasiMonteCarlo.generate_design_matrices(n,inj_min,inj_max,sampler,n)
    return A,B
end



function sensitivity_global(inj_range_matrix,f_v=get_inj_range_voltages,
    method=Sobol(),order=2,N=1000,batch=true,nboot=100,conf_int=0.95)
    """Numeric global sensitivities for a single inj_bus"""
    sens = GlobalSensitivity.gsa(f_v,method,inj_range_matrix; N, batch=batch,order=order,nboot=nboot,conf_int=conf_int)
    return sens
end




