using Convex: norm_fro
using Base: String
using Convex 
using Mosek
using JuMP
using GLPK
using LinearAlgebra


solver = MosekSolver(LOG=1)

function mat_complete_model(X,write_location::String)
    model = Model(GLPK.Optimizer)
    set_time_limit_sec(model,60^2)
    n_measurements,n_locations = size(X)
    @variable(model,x[1:n_measurements,1:n_locations])
    return model 
end

function cons_mat_complete_model(X,δ,β,g,write_location::String)
    model = mat_complete_model(X,write_location)
    @objective(model,Min,nuclearnorm(X)) #Minimize the nuclear norm
    @constraint(model,data_con,norm(X_Ψ-M_Ψ)<=δ) #Least-squares sense constraint
    @constraint(model,phys_con,norm(g(x))<=β) #System physics constraints
    return model
end
 

function solve_model(model,verbose=True)
    optimize!(model)
    if verbose
        solution_summary(model)
    else
        print(model)
    end
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
