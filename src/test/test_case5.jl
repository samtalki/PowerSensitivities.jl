### Test the results for case 5
include("/home/sam/github/PowerSensitivities.jl/src/PowerSensitivities.jl")
using PowerModels
using PowerModelsAnalytics
using LinearAlgebra
using Ipopt
using JuMP
using Gadfly
import .PowerSensitivities
const ps = PowerSensitivities


#Load Case 5
c5 = make_basic_network(parse_file("/home/sam/github/PowerSensitivities.jl/data/matpower/case5.m"))

#Make timeseries dataset
ami_dataset = make_ami_dataset(c5)
diff_data = calc_finite_differences(dataset)
pnet,qnet,vm = diff_data["pnet"],diff_data[""]
pnet,qnet,vm = diff_data["pnet"],diff_data[""]

#Calculate full Jacobian
J = ps.calc_jacobian_matrix(case5)

#Calculate Jacobian for PQ buses only
J_pq = ps.calc_jacobian_matrix(case5,[1]) #1 -- pq buses only
S = inv(Matrix(J_pq.matrix)) #get the sensitivity matrix
dpdv_true,dqdv_true = S[2,1],S[2,2]

#Estimate sensitivities
λ=0 #ℓ2 regularization
pf = ps.calc_basic_power_factor(c5,[1])
spv,sqv = ps.est_power_voltage_sens(diff_data,λ)