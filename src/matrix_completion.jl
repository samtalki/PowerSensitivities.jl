

"""
Given an incomplete sensitivity matrix S, return a matrix completion model
"""
function mat_complete_model(S_0,Δv,Δx)
    model = Model(Ipopt.Optimizer)
    m_measurements,n_injections = size(S_0)
    @variable(model,S[1:m_measurements,1:n_injections])
    @objective(model,Min,norm(S - S_0))
    return model 
end

function cons_mat_complete_model(S_0,Δv,Δx,δ,β,g)
    model = mat_complete_model(X,Δv,Δx)
    @objective(model,Min,nuclearnorm(X)) #Minimize the nuclear norm
    @constraint(model,data_con,norm(X_Ψ-M_Ψ)<=δ) #Least-squares sense constraint
    @constraint(model,phys_con,norm(g(x))<=β) #System physics constraints
    return model
end

