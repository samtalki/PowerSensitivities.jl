using Base: String
using Convex 
using Mosek
using JuMP
using GLPK
using LinearAlgebra


solver = MosekSolver(LOG=1)

function mat_complete(X,write_location::String)
    model = Model(GLPK.Optimizer)
    set_time_limit_sec(model,60**2)
    n_measurements,n_locations = size(X)
    @variable(model,x[1:n_measurements,1:n_locations])
    return model 
end

function cons_mat_complete(X,bounds,write_location::String)
    model = mat_complete(X,write_location)
    
end

function add_component_to_model(model::JuMP.Model)
    x = model[:x]
end

function matrix_completion(X,Y)
    X = Convex.Variable(size(X))
    problem = minimize(nuclearnorm(X))
    problem.constraints += X[obsidx] == Y[obsidx]
    @time solve!(problem,solver)
end


model = Model(GLPK.Optimizer)
