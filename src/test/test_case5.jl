### Test the results for case 5
include("/home/sam/github/PowerSensitivities.jl/src/PowerSensitivities.jl")
import .PowerSensitivities
using PowerModels
using PowerModelsAnalytics
using LinearAlgebra
using Ipopt
using JuMP
using Gadfly



c5 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case5.m"))
Jc5 = calc_basic_jacobian_matrix(case5)

#Make timeseries dataset
dataset = make_ami_dataset(c5,1,100)
diff_data = calc_finite_differences(dataset)

#Estimate sensitivities
λ=0 #ℓ2 regularization
S_0 = PowerSensitivities.est_power_voltage_sens(diff_data,λ)