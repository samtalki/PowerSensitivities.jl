# BSD 3-Clause License

# Copyright (c) 2022, Samuel Talkington
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

using JuMP
using LinearAlgebra
import PowerModels as PM
import PowerSensitivities as PS
import SCS
import Ipopt
include("/home/sam/github/PowerSensitivities.jl/src/sens/nr_sens.jl")



function variable_jacobian(model,data::Dict)   
    n_bus = length(data["bus"])
    @variable(model, J[1:2*n_bus,1:2*n_bus])
    @variable(model,∂pv[1:n_bus,1:n_bus])
    @variable(model,∂qv[1:n_bus,1:n_bus])
    @variable(model,∂pθ[1:n_bus,1:n_bus])
    @variable(model,∂qθ[1:n_bus,1:n_bus])
end

"""
Jacobian symmetry constraints
"""
function constraint_jacobian_physics(model,data::Dict)
    vph = PM.calc_basic_bus_voltage(data)
    s = PM.calc_basic_bus_injection(data)
    p,q = real.(s),imag.(s)
    vm,θ = abs.(vph),angle.(vph)

    offdiag2(A::AbstractMatrix) = [A[ξ] for ξ in CartesianIndices(A) if ξ[1] ≠ ξ[2]] 
    
    @constraint(model,
        tr(∂pθ) == tr(diagm(vm)*∂qv - 2*diagm(q))
    )
    @constraint(model,
        tr(∂qθ) == tr(-1 .*diagm(vm)*∂pv + 2*diagm(p))
    )
    @constraint(model,
        offdiag2(∂pθ) == offdiag2(diagm(vm)*∂qv)
    )
    @constraint(model,
        offdiag2(∂qθ) == -1*offdiag2(diagm(vm)*∂pv)
    )
end

function model_est_bus_voltage_phase(data::Dict)
    compute_ac_pf!(data)
    
    s = PM.calc_basic_bus_injection(data)
    v = PM.calc_basic_bus_voltage(data)
    vmag,θ_true = abs.(v),angle.(v)
    p,q = real.(s),imag.(s)
    signs_q = [-1*sign(q_i) for q_i in q]
    signs  = [signs_q ; signs_q]
    #Get Jacobian 
    J = PS.calc_jacobian_matrix(data)

    #Make problem data
    n = length(p)
    A = [J.pv ; J.qv]
    B = [J.pth ; J.qth]
    x = [p ; q]
    Apinv = pinv(Matrix(A))
    diagX = Diagonal(x)

  

    #Find v∈C^n s.t. 
    model = Model(Ipopt.Optimizer)
    @variable(model,θ[1:n])
    @variable(model,ξ[1:2*n]) #Phase to be learned
    @variable(model,t)
    @constraint(model,[t;ξ] in SecondOrderCone())
    @constraint(model,t==n)
    @objective(model,Min,
        (1/2)*transpose((A*vmag + B*θ) - Diagonal(x)*ξ)*((A*vmag + B*θ) - Diagonal(x)*ξ)
    )
   

    # for i in 1:n
    #     @NLconstraint(model,
    #         k[i] == (1/pf[i])*sqrt(1 - pf[i]^2)
    #     )
    # end
    # for i in n+1:2*n
    #     @NLconstraint(model,
    #         k[i] == pf[i]/(sqrt(1-pf[i]^2))
    #     )
    # end
    return model
end

function est_bus_voltage_phase!(A::Matrix{Complex},b::Vector{Real})
    model = make_phase_retrieval_model(A::Matrix{Complex},b::Vector{Real})
    optimize!(model)
    X_val = value.(X)
    return calc_closest_rank_r(X_val,1)[:,1]
end

"""
Find the closest rank-R approximate matrix of A
"""
function calc_closest_rank_r(A::Matrix,r::Integer)
    (m,n) = size(A)
    U,Σ,V = svd(A)
    for (i,s_i) in enumerate(Σ)
        if i > r 
            Σ[i] = 0
        end
    end
    return U * Diagonal(Σ) * V' 
end
"""
WIP
"""
function model_phase_jacobian_retrieval(data::Dict; verbose=true)
    model = Model(SCS.Optimizer)
    
    constraint_jacobian_physics(model,data)
    #The phase angle difference between voltage and current
    @variable(model,θ[1:n_bus]) 
    #Real and imaginary parts of complex voltage e,f s.t. v = e + jf
    @variable(model,e[1:n_bus])
    @variable(model,f[1:n_bus])
    #Real and imaginary parts of complex power p,q s.t. s = p + jq
    @variable(model,p[1:n_bus])
    @variable(model,q[1:n_bus])


    #Constraint on the measurements
    @constraint(model,abs.(b) == y)
    @objective(model,Min,x)
    return model
end 

# """
# Given network data dict return whole jacobian estimated from magnitudes
# """
# function model_jacobian_retrieval(data::Dict; verbose=true)
    
# end


