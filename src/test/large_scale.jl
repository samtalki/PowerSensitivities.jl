#Large scale test(s) of voltage sensitivity estimation conditions.
import PowerModels as PM
import PowerSensitivities as PS
import PowerModelsAnalytics as  PMA
using Plots
include("../test/thm1.jl")
include("../util/matrix.jl")
include("phobs.jl")

#--- Initalize Parameters
test_run_counter = 1 #Run counter
include_academic = true #Boolean on whether to include the academic radial cases or not
#plots = [] # plots

network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
rts_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/RTS_GMLC.m"
texas_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/ACTIVSg2000.m"
tenk_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/case_ACTIVSg10k.m"
#east_us_70k_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/ACTIVSg70k/case_ACTIVSg70k.m"

# paths and names
large_scale_paths = [rts_path, texas_path, tenk_path] #east_us_70k_path]
large_scale_names = ["rts-gmlc", "texas", "tenk"]#, "70k-east-us"]
default_names = readdir(network_data_path);
default_paths = readdir(network_data_path,join=true);
paths,names = vcat(large_scale_paths,default_paths),vcat(large_scale_names,default_names)

thm1_res_pq,thm1_res_pqpv = Dict(),Dict()
thm2_res_pq,thm2_res_pqpv = Dict(),Dict()
thm1_data_pq,thm2_data_pq = Dict(),Dict()
thm1_bound_pq,thm1_bound_pqpv = Dict(),Dict()
thm1_data_pqpv,thm2_data_pqpv = Dict(),Dict()
for (path,name) in zip(paths,names)
    #Make a network model
    net = PM.make_basic_network(PM.parse_file(path))
    #Don't test the large scale networks for radiality
    if name ∉ large_scale_names 
            try
            #Solve the power flow equations
            PM.compute_ac_pf!(net)

            #Compute the jacobians
            J_pq,J_pq_pv = PS.calc_jacobian_matrix(net,[1]),PS.calc_jacobian_matrix(net,[1,2])
            
            #Compute Theorem data
            thm1_res_pq[name],thm1_res_pqpv[name] = test_thm1(net),test_thm1(net,[1,2])
            thm2_res_pq[name],thm2_res_pqpv[name] = test_observability(net),test_observability(net,[1,2])
            thm1_data_pq[name],thm1_data_pqpv[name] = calc_thm1_data(net,[1]),calc_thm1_data(net,[1,2])
            thm2_data_pq[name],thm2_data_pqpv[name] = calc_thm2_data(net,[1]),calc_thm2_data(net,[1,2])
            thm1_bound_pq[name],thm1_bound_pqpv[name] = calc_Δk_bound(net,[1]),calc_Δk_bound(net,[1,2])
            
            #On the first run of the session, plot the network models.
            if test_run_counter == 1
                save_path = "/home/sam/github/PowerSensitivities.jl/figures/spring_22/"
                plot = PMA.plot_network(net,filename=save_path*name*".pdf")
                #push!(plots,plot)
            end
        catch
            continue
        end
    end
end

test_run_counter+=1
