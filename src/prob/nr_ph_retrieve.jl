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



"""
Given a basic network data dict, construct the physics-informed phase retrieval problem.
"""
function nr_phase_retrieval(network::Dict)
    
    compute_ac_pf!(network)

    #Voltage and power injections
    vphasor = calc_basic_bus_voltage(network)
    vm,va = abs.(vphasor),angle.(vphasor)
    
    s = calc_basic_bus_injection(network)
    p,q = real.(s),imag.(s)


    J = calc_jacobian_matrix(network)
    Spth,Sqth,Spv,Sqv = J.pth,J.qth,J.pv,J.qv
    

    

end


"""
Given complex voltage and powers, and an admittance matrix Y, calculate the power flow mismatch.
"""
function calc_mismatch(v_ph,s,Y)
    n_bus = size(Y,1)
    # Convert phasor to rectangular
    function calc_rect_bus_voltage(v_ph::AbstractArray)
        θ,vm = v_ph[1:n_bus],v_ph[n_bus+1:end]
        [vmᵢ*(cos(θᵢ) + sin(θᵢ)*im) for (θᵢ,vmᵢ) in zip(θ,vm)]
    end
    v_rect = calc_rect_bus_voltage(v_ph)
    si = v_rect .* conj(Y * v_rect) #compute the injection
    Δp, Δq = real(s - si), imag(s - si)
    Δx = [Δp ; Δq]
    return Δx
end

"""
1st definition: 
given a PowerModels network, returns a [vangle; vmag] vector of bus voltage phasor quantities
as they appear in the network data.
"""
function calc_phasor_bus_voltage(data::Dict{String, Any})
    b = [bus for (i,bus) in data["bus"] if bus["bus_type"] != 4]
    bus_ordered = sort(b, by=(x) -> x["index"])
    θ,vm = [bus["va"] for bus in bus_ordered],[bus["vm"] for bus in bus_ordered]
    return [θ ; vm]
end

"""
2nd definition (Multiple dispatch!): 
Given vector of rectangular complex voltages [e+jf], calculate the phasor form of the voltages [vangle; vmag]
"""
function calc_phasor_bus_voltage(v_rect::AbstractArray)
    θ,vm = [angle(v_i) for v_i in v_rect],[abs(v_i) for v_i in v_rect]
    return [θ ; vm]
end

function variable_jacobian(model::Model,data::Dict)   
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


"""
Given a PowerModels basic network dict, 
construct a JuMP model to estimate the bus voltage phase angles 
"""
function model_est_bus_voltage_phase(data::Dict)
    
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


