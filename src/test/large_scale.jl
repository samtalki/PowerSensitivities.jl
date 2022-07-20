#Large scale test(s) of voltage sensitivity estimation conditions.
import PowerModels as PM
import PowerSensitivities as PS
import PowerModelsAnalytics as  PMA
using Plots
include("../test/thm1.jl")
include("../util/matrix.jl")
include("phobs.jl")

#Initialize a run counter
test_run_counter = 1
plots = []

rts_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/RTS_GMLC.m"
twok_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/ACTIVSg2000.m"
tenk_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/case_ACTIVSg10k.m"
east_us_70k_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/ACTIVSg70k/case_ACTIVSg70k.m"
paths = [rts_path, twok_path, tenk_path] #east_us_70k_path]
names = ["rts-gmlc", "2k-active", "10k-active"]#, "70k-east-us"]

thm1_res_pq,thm1_res_pqpv = Dict(),Dict()
thm2_res_pq,thm2_res_pqpv = Dict(),Dict()
thm1_data_pq,thm2_data_pq = Dict(),Dict()
thm1_data_pqpv,thm2_data_pqpv = Dict(),Dict()
for (path,name) in zip(paths,names)

    #Solve the AC power flow equations
    net = PM.make_basic_network(PM.parse_file(path))
    PM.compute_ac_pf!(net)

    #Compute the jacobians
    J_pq,J_pq_pq = PS.calc_jacobian_matrix(net,[1]),PS.calc_jacobian_matrix(net,[1,2])
    
    #Compute Theorem data
    thm1_res_pq[name],thm1_res_pqpv[name] = test_thm1(net),test_thm1(net,[1,2])
    thm2_res_pq[name],thm2_res_pqpv[name] = test_observability(net),test_observability(net,[1,2])
    thm1_data_pq[name],thm1_data_pqpv[name] = calc_thm1_data(net),calc_thm1_data(net,[1,2])
    thm2_data_pq[name],thm2_data_pqpv[name] = calc_thm2_data(net),calc_thm2_data(net,[1,2])

    #On the first run of the session, plot the network models.
    if test_run_counter == 1
        save_path = "/home/sam/github/PowerSensitivities.jl/figures/spring_22/"
        plot = PMA.plot_network(net,filename=save_path*name*".pdf")
        push!(plots,plot)
    end
end

test_run_counter+=1
