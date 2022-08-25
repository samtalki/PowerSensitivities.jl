#Large scale test(s) of voltage sensitivity estimation conditions.
import PowerModels as PM
import PowerSensitivities as PS
import PowerModelsAnalytics as  PMA
using Plots,ProgressMeter
using LaTeXStrings

include("../test/thm1.jl")
include("../util/matrix.jl")
include("phobs.jl")



"""
Given a PowerModels network data dict, return a 2x2 subplot containing the cummulative and normalized singular values of the 
∂v/∂p and ∂v/∂q inverse power flow Jacobian (voltage sensitivity) submatrices.
"""
function plot_spectral_analysis(net::Dict{String,<:Any})
    name = net["name"] #network name
    spec_dict = PS.calc_spectral_analysis(net)
    svp_spec,svq_spec = spec_dict["svp"],spec_dict["svq"]
    n_singular_vals = length(svp_spec.scuml)
    idx = [i for i in 1:n_singular_vals]
    
    #Check the lengths
    @assert length(svp_spec.scuml) == length(svq_spec.scuml) == length(idx)
    
    #Plot the cummulative singular values of the submatrices
    #dvdp plot
    plt_svp = plot(
        idx,svp_spec.scuml,
        ls=:dashdot,
        ylabel="Value",
        xlabel="Singular Value",
        label="Cummulative Normalized",
        #label=L"\frac{\partial v}{\partial p} "*name
    )
    plot!(
        idx,svp_spec.snorm,
        ls=:dash,
        ylabel="Value",
        xlabel="Singular Value",
        label="Normalized",
        title=L"\frac{\partial v}{\partial q} "*name
    )
    plt_svq = plot(
        idx,svq_spec.snorm,
        ls=:dashdot,
        ylabel="Value",
        xlabel="Singular Value",
        label="Normalized",
        #label=L"\frac{\partial v}{\partial p} \ "*name
        )
    plot!(idx,svq_spec.scuml,
        ls=:dashdot,
        ylabel="Value",
        xlabel="Singular Value",
        label="Cummulative Normalized",
        title=L"\frac{\partial v}{\partial q} \ "*name
    )
    return plot(plt_svp,plt_svq)
end

"""
Given a SpectralAnalysis, plot the normalized and cummulative singular values.
"""
function plot_spectral_analysis(spec::PS.SpectralAnalysis)
    #Make indexes to plot
    n_singular_vals = length(spec.scuml)
    idx = [i for i in 1:n_singular_vals]
    pcuml = plot(
        idx,spec.scuml,
        ylabel="Cummulative Normalized",
        xlabel="Singular Value"
        )
    pnorm = plot(
        idx,spec.snorm,
        ylabel="Normalized",
        xlabel="Singular Value"
        )
    return plot(pcuml,pnorm)
end


"""
Given a list of network model paths and names...
"""
function plot_spectral_analysis(paths::AbstractArray{String},names::AbstractArray{String})
    cuml_plot,norm_plot = plot(),plot()
    @showprogress 1 "Plotting..." for (k,(name,path)) in enumerate(zip(names,paths))
        sleep(0.1)
        #Make a network model
        net = PM.make_basic_network(PM.parse_file(path))
        #Solve the power flow equations
        PM.compute_ac_pf!(net)
        #Compute the spectral analysis
        spec_dict = PS.calc_spectral_analysis(net)
        svp_spec,svq_spec = spec_dict["svp"],spec_dict["svq"]
        n_singular_vals = length(svp_spec.scuml)
        idx = [i for i in 1:n_singular_vals]
        plot!(cuml_plot,idx,svp_spec.scuml,xlabel="Singular Value",ylabel="Cummulative Normalized",label=false,ls=:dash)
        plot!(cuml_plot,idx,svq_spec.scuml,xlabel="Singular Value",ylabel="Cummulative Normalized",label=false,ls=:dash)
        plot!(norm_plot,idx,svp_spec.snorm,xlabel="Singular Value",ylabel="Normalized",label=L"\frac{\partial v}{\partial p} \ "*name,ls=:dashdot)
        plot!(norm_plot,idx,svq_spec.snorm,xlabel="Singular Value",ylabel="Normalized",label=L"\frac{\partial v}{\partial q} \ "*name,ls=:dashdot)
    end
    return plot(cuml_plot,norm_plot)
end



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
    paths,names = large_scale_paths,large_scale_names
    
end

plt = plot_spectral_analysis(paths,names)
savefig(plt,"figures/spring_22/spectral_analysis/large_scale_spectral.pdf")
savefig(plt,"figures/spring_22/spectral_analysis/large_scale_spectral.png")