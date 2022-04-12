## Calculate the minimum power factor
include("../PowerSensitivities.jl")
include("../util/matrix.jl")
import .PowerSensitivities
using PowerModels: parse_file,make_basic_network
using Plots, LaTeXStrings
using LinearAlgebra
theme(:ggplot2)

"""
Given:
    a network data dict under study,
    a chosen maximum power factor pf_max,
    and the SELECTED BUS TYPES under study
Compute:
    the minimum power factor pf_min such that the complex power injections are observable
"""
function calc_pf_min(network::Dict{String,<:Any},sel_bus_types=[1,2],pf_max::Real=1)
    #Compute the indeces that will be considered
    bad_idx = PowerSensitivities.calc_bad_idx(network)
    idx_sel_bus_types = PowerSensitivities.calc_bus_idx_of_type(network,sel_bus_types)
    study_idx = [i for i in 1:length(network["bus"]) if((i ∉ bad_idx) && (i ∈ idx_sel_bus_types))]
    #Compute the matrices of interest
    K = PowerSensitivities.calc_K_matrix(network)[study_idx,study_idx]
    ∂p∂θ = PowerSensitivities.calc_pth_jacobian(network)[study_idx,study_idx]
    ∂q∂θ = PowerSensitivities.calc_qth_jacobian(network)[study_idx,study_idx]
    k_max = maximum(diag(K))
    M = k_max*∂p∂θ - ∂q∂θ
    #M = PowerSensitivities.k(pf_max)*∂p∂θ - ∂q∂θ
    #Check the sizes
    @assert size(∂p∂θ,1) == length(study_idx) && size(∂p∂θ,2) == length(study_idx)
    @assert size(∂q∂θ,1) == length(study_idx) && size(∂q∂θ,2) == length(study_idx)
    @assert size(K,1) == length(study_idx) && size(K,2) == length(study_idx)
    #Compute the bound on Δk
    if pf_max==1
        pf_min = PowerSensitivities.kinv(
            (opnorm(inv(M))^-1)*(opnorm(∂p∂θ)^-1))
    else
        pf_min = PowerSensitivities.kinv(
            PowerSensitivities.k(pf_max) + (opnorm(inv(M))^-1)*(opnorm(∂p∂θ)^-1))
        if pf_min>pf_max
            pf_min = NaN
        end
    end
    return pf_min
end


"""
Given chosen bus types under study, a maximum power factor, and folder of test cases under study,
Compute the minimum power factor implied by Theorem 1 of Talkington and Turizo et al. for the feeder models in network_data_path
"""
function test_pf_min(sel_bus_types=[1,2],pf_max::Real=1,network_data_path=network_data_path)
    names = readdir(network_data_path);
    paths = readdir(network_data_path,join=true);
    pf_min = Dict()
    for (name,path) in zip(names,paths)
        if name ∈ study_network_names
            network =  try
                make_basic_network(parse_file(path));
            catch
                println("PM cannot parse "*name)
                continue
            end
            if PowerSensitivities.is_radial(network) || allow_mesh
                pf_min[name] = try
                    calc_pf_min(network,sel_bus_types,pf_max)
                catch
                    pf_min[name] = NaN
                end
            end
        end
    end
    return pf_min
end

#Test case path and parameters
network_data_path="/home/sam/github/PowerSensitivities.jl/data/matpower/" #Folder with meshed and radial systems
#network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
allow_mesh = true #Whether to allow meshed test cases/require test feeder is radial

#Select networks to plot

#Radial cases
#global study_network_names = ["case4_dist.m" "case16ci.m" "case34sa.m" "case136ma.m"]

#Mesh cases
global study_network_names = ["case5.m" "case9.m" "case14.m" "case24.m"]


#Select maximum power factors to test
global pf_max = 0.75:0.005:1


begin
    #Compute the minimum bus power factors at unity max power factor
    pf_min_unity = test_pf_min([1],1)
    #Compute the minimum bus power factors for a range of power factors - store as a matrix
    pf_min_nonunity = zeros((length(pf_max),length(study_network_names))) #Store as a matrix
    for (i,pf_max_i) in enumerate(pf_max)
        pf_min = test_pf_min([1],pf_max_i)
        for (j,net_name) in enumerate(study_network_names) #Fill up the matrix for potting
            pf_min_nonunity[i,j] = pf_min[net_name]
        end
    end
    #Plot the minimum power factors
    plot(pf_max,pf_min_nonunity,
        title="Feasible Power Factors "*L"\{\alpha : \alpha_{\rm min} \leq \alpha \leq 1\}"*" Given "*L"\alpha_{\rm max}",
        label=study_network_names,
        fillrange=1,
        fillalpha=0.15,
        lw=3,alpha=0.6,
        linestyle=:dot,
        markershape=:square,
        markersize=2,
        markeralpha=0.3,
        xaxis=("Chosen "*L"\alpha_{\rm max}",font(pointsize=12)),
        yaxis=(L"\alpha_{\rm min} = k^{-1}(k(\alpha_{\rm max}) + \Delta k_{\rm max})",font(pointsize=11)),
        legend_position=:bottomleft,
        legend_font_pointsize=12
        #dpi=300,
        #size=(floor(3.5*400),floor((3.5/1.61828)*400))
        )
    savefig("figures/spring_22/mesh_pf_min_current_kmax.png")
    savefig("figures/spring_22/mesh_pf_min_current_kmax.pdf")
end




