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


# --- Topology type checkers: methods for assessing the topology of a distribution network
"""
Given a network data dict, check if the network is radial
"""
function is_radial(network::Dict{String,<:Any})
    #Get admittance matrix and the sum of nonzero upper and lower off diagonal elements
    Y = calc_basic_admittance_matrix(network)
    return is_radial(Y)
end
"""
Given a nodal admittance matrix Y, see if the network is radial
"""
function is_radial(Y::AbstractMatrix) 
    n = size(Y)[1]
    #Upper and lower off-diagonal elements
    U(A::AbstractMatrix) = [A[i] for i in CartesianIndices(A) if i[1]>i[2]]
    L(A::AbstractMatrix) = [A[i] for i in CartesianIndices(A) if i[1]<i[2]]
    
    #Get the nonzero upper and lower off diagonal elements
    nz_upper = [1 for y_ij in U(Y) if y_ij != 0]
    nz_lower = [1 for y_ij in L(Y) if y_ij != 0]
    
    return !(sum(nz_upper)>n-1 || sum(nz_lower) >n-1)
end

is_mesh(x::Union{Dict{String,<:Any},AbstractMatrix}) = not(is_radial(x))



