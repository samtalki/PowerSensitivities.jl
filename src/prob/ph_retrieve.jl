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
import SCS


"""
Given A - Matrix{Complex,(N,N)}
and b - Vector{Real,M}
    Find x ∈ Vector{Complex,N}
    such that |Ax| = b elementwise
"""
struct PhaseRetrievalModel
    A::AbstractMatrix
    b::Vector{Real}
    M::AbstractMatrix
    X::AbstractMatrix
    f_obj::Real
    success::Bool
    model::Model
end


function make_phase_retrieval_model(A::Matrix{Complex},b::Vector{Real})
    model = Model(SCS.Optimizer)

    #Measurement dimension check
    (mA,nA),(mb,nb) = size(A),size(b)
    @assert mA == mb "Measurements b (M x 1) and problem data A (M x N) dimensionality mismatch in M"
    m,n = mA,nA

    #Construct M matrix (fixed phase representation)
    I = Matrix(1.0 * I, m, m)
    M = Diagonal(b)*(I - A * pinv(A))*Diagonal(b)
    Mr,Mi = real.(M),imag.(M)
    n = size(A,1)

    #Make variables and model
    @variable(model,Xr[1:n,1:n])
    @variable(model,Xi[1:n,1:n])
    @constraint(model, [Xr Xi; -Xi Xr]>=0,PSDCone())
    @objective(model,Min,trc([Mr Mi; -Mi Mr]*[Xr Xi; -Xi Xr]))
    return model
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
Construct a rank one matrix X= u * conj(tranpose(u))
"""
function variable_rank_one_matrix(n_bus::Int)
    @variable(model,X[1:n_bus,1:n_bus])
    @constraint(model, X>=0,PSDCone())
end




"""
gerchbergsaxton(Ain, Aout, kw...)
    reconstruct complex signals from their amplitude values based on the Gerchberg-Saxton algorithm.

    References
    1. R. W. Gerchberg and W. O. Saxton, "A practical algorithm for the determination of the phase from image and diffraction plane pictures," Optik 35, 237 (1972)
"""
function gerchbergsaxton(Ain, Aout; 
        maxiter=100,
        projector=(fft, ifft),
        errorfunc=msd,
        tol=eps(),
        init=randn(ComplexF64, size(Ain))
    )
    Cin = init
    @. Cin = Ain * exp(im * angle(Cin))
    Cout = similar(Aout)
    for i in 1:maxiter
        Cout = projector[1](Cin)
        @. Cout = Aout * exp(im * angle(Cout))
        Cin = projector[2](Cout)
        if errorfunc(abs.(Cin), Ain) < tol
            break
        end
        @. Cin = Ain * exp(im * angle(Cin))
    end
    (Cin, Cout)
end
