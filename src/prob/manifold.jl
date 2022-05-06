using JuMP
using LinearAlgebra
using Ipopt

"""
JuMP model for phase retrieval 
"""
function phase_retrieval(; verbose=true)
    model = Model(Ipopt.Optimizer)
    @variable(model,θ[1:n_bus]) #The phase angle difference between voltage and current
    @variable(model,p[1:n_bus])
    @variable(model,q[1:n_bus])
    @objective(model,Min,
        
    )
end 

function jacobian_retrieval(; verbose=true)
    @variable(model, J[1:n_bus,1:n_bus])
end