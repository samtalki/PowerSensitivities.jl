#Large scale test(s) of voltage sensitivity estimation conditions.
import PowerModels as PM
import PowerSensitivities as PS
import PowerModelsAnalytics as  PMA
using Plots
using LaTeXStrings


#--- Initalize Parameters
include_industry = true
include_academic = false #Boolean on whether to include the academic radial cases or not
#plots = [] # plots

#- paths and names
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/figure_cases/" #Folder with radial-only systems
rts_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/RTS_GMLC.m"
texas_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/ACTIVSg2000.m"
#tenk_path = "/home/sam/github/PowerSensitivities.jl/data/PowerSystems_data/matpower/case_ACTIVSg10k.m"

large_scale_paths = [rts_path, texas_path] #tenk_path] #east_us_70k_path]
large_scale_names = ["rts-gmlc", "texas"] #"tenk"]#, "70k-east-us"]
default_names = readdir(network_data_path);
default_paths = readdir(network_data_path,join=true);

if include_industry==true && include_academic==true
    paths,names = vcat(large_scale_paths,default_paths),vcat(large_scale_names,default_names)
elseif include_industry==false && include_academic==true
    paths,names = default_paths,default_names
else
    paths,names = large_scale_paths,default_paths
end
