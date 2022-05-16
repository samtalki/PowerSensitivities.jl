#Model based analytical voltage sensitivity analysis toolkit
#Based on:
#"Efficient Computation of Sensitivity Coefficients of Node Voltages and Line Currents in Unbalanced Radial Electrical Distribution Networks", 
#K. Christakou, et al., IEEE Transactions on Power Systems, 2013.

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

using PowerModels
using LinearAlgebra

abstract type Sensitivities end
abstract type VoltageSensitivities <: Sensitivities end
abstract type TapChangeSensitivities <: VoltageSensitivities end

struct VoltagePhasorSensitivities <: VoltageSensitivities
	p::AbstractMatrix #Voltage phasor - active power sensitivities
	q::AbstractMatrix	#voltage phasor- reactive power sensitiviies
end
struct VoltageMagnitudeSensitivities <: VoltageSensitivities
	p::AbstractMatrix
	q::AbstractMatrix
end

#--- Voltage phasor sensitivities

function calc_vf_sensitivities(data::Dict{String,<:Any})
	compute_ac_pf!(data)
	Y = calc_basic_admittance_matrix(data)
	vf = calc_basic_bus_voltage(data)
	return calc_vf_sensitivites(Y,vf)
end
function calc_vf_sensitivites(Y::AbstractArray,vf::AbstractArray)
	n = size(Y,1)
	B = [diagm(ones(n)) zeros((n,n)); zeros((n,n)) -1*diagm(ones(n))*im] 
	A =  [Y*diagm(conj.(vf)) Y*diagm(vf); Y*diagm(vf) Y*diagm(conj.(vf))] 
	#B = X A
	X = B * inv(A)
	Xre = X[:,1:n]
	#non-conj phasor voltage sensitivities
	∂vf∂p, ∂vf∂q = Xre[1:n,1:n],Xre[n+1:end,n+1:end] 
	return VoltagePhasorSensitivities(∂vf∂p,∂vf∂q) 
end

#--- Voltage magnitude sensitivities

function calc_vmag_sensitivities(data::Dict{String,<:Any})
	compute_ac_pf!(data)
	Y = calc_basic_admittance_matrix(data)
	vf = calc_basic_bus_voltage(data)
	∂vf = calc_vf_sensitivites(Y,vf)
	∂v = calc_vmag_sensitivities(∂vf)
	return calc_vmag_sensitivities(Y,vf)
end
function calc_vmag_sensitivities(∂vf::VoltagePhasorSensitivities)
	∂vf∂p,∂vf∂q = ∂vf.p,∂vf.q
	@assert size(∂vf∂p,1) == size(∂vf∂q,1) "Mistmatch in phasor voltage sensitivity dimensionality"
	n = size(∂vf∂p,1)
	∂vp,∂vq = zeros((n,n)), zeros((n,n)) #Voltage magnitude sensitivities
	for (i,(dvfp_i,dvfq_i)) in enumerate(zip(eachrow(∂vf∂p),eachrow(∂vf∂q)))
		∂vp[i,:] = 1/(abs(vf[i])) .* real.(conj(vf[i]) .* dvfp_i)
		∂vq[i,:] = 1/(abs(vf[i])) .* real.(conj(vf[i]) .* dvfq_i)
	end
	return VoltageMagnitudeSensitivities(∂vp,∂vq)
end



# --- Tap changing transformer sensitivities

function calc_tap_change_sensitivities(data::Dict{String,<:Any})
	compute_ac_pf!(data)
	vph = calc_basic_bus_voltage(data)
	vm,va = abs.(vph),angle.(vph)
end