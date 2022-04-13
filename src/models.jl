include("PowerSensitivities.jl")
import .PowerSensitivities
using PowerModels
using LinearAlgebra
using Flux 
using Flux.Losses: mse
using Flux.Optimise: update!
using Convex


norm_nuc(X::Matrix) = sum(svdvals(X))

"""
Nuclear-norm regularized matrix completion loss
"""
function matrec_loss(Δx,Δv)
    Δv̂ = S*Δx
    return Flux.Losses.mse(Δv̂,Δv) + λ*norm_nuc(S)
end

c5 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case5.m"))
dataset = PowerSensitivities.make_ami_dataset(c5,[1],500)


θ = params(S,λ)

grads = gradient(() -> matrec_loss(Δx,Δv),θ)
