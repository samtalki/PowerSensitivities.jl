import PowerModels as PM
using Plots,LaTeXStrings
theme(:dao)
using LinearAlgebra
using DataFrames
using PowerSensitivities

"""
Given a network data dict plot a single 1x2 subplot plot of the q->p and p->q encoding
"""
function plot_network_encoding(network::Dict{String,<:Any},pf_range=(0,1))
    #gr()

    #calculate the parameters to be ploted
    pf = calc_basic_power_factor(network)
    K = calc_K_matrix(network)
    Kinv = inv(K)

    #Plot Q->P matrix values
    active = plot(
        pf,diag(K),
        color_palette = palette(:Greens),
        label=L"$K(\alpha_i)$",
        seriestype = :scatter
        )



    plot!(pf,-1*diag(K),
        color_palette = palette(:Greens),
        label=L"$-K(\alpha_i)$",
        seriestype = :scatter)
    xlabel!(L"Bus Power Factor $\alpha_i$")
    ylabel!(L"$K_{ii} = \pm \frac{1}{\alpha_i} \sqrt{1-\alpha_i^2}$")
    
    #Plot P->Q matrix values
    reactive = plot(pf,diag(Kinv),
        color_palette = palette(:Oranges),
        label=L"$K^{-1}(\alpha_i)$",
        seriestype = :scatter)
    plot!(pf,-1*diag(Kinv),
        color_palette = palette(:Oranges),
        label=L"$-K^{-1}(\alpha_i)$",
        seriestype = :scatter)
    xlabel!(L"Power Factor $\alpha_i$")
    ylabel!(L"$K^{-1}_{ii} = \pm \frac{\alpha_i}{\sqrt{1-\alpha_i^2}}$")
    
    #Plot them together in subplots
    #ylabel!(L"Active Power $p_i(q_i;\alpha_i)")
    plt = plot(active,reactive)
    
    return plt
end

function test_power_encodings(network)
    K,Kinv = compute_power_encodings(network)
    for (i,kii) in enumerate(K) 
        @assert kii>0 "Unity power factor on bus"
    end

end

"""
Given a network data dict, compute the power factor encoding matrix K
"""
function compute_power_encodings(network)
    K = calc_K_matrix(network)
    Kinv = inv(K)
    bad_idx = calc_bad_idx(network)
    good_idx = [i for i in 1:length(network["bus"]) if i ∉ bad_idx]
    return diag(K)[good_idx],diag(Kinv)[good_idx]
end


function plot_network_encoding!(network::Dict{String,<:Any},pf_range=(0,1))
    compute_ac_pf!(network)
    return plot_network_encoding(network,pf_range)
end

#----

function make_network_encoding_dataframe(network::Dict)
    n_bus = length(nethwork["bus"])
    name = network["case_name"]
end



"""
Given a path to network models, plot all the embedding matrices on a polar plot together
"""
function plot_all_case_encodings(network_path="/home/sam/github/PowerSensitivities.jl/data/figure_cases/")
    names = readdir(network_path);
    paths = readdir(network_path,join=true);
    Ks,Kinvs = Dict(),Dict()
    plt_K,plt_inv_K = plot(),plot()
    for (name,path) in zip(names,paths)
        network =  try
            make_basic_network(parse_file(path));
        catch
            println("PM cannot parse "*name)
            continue
        end
        if PowerSensitivities.is_radial(network)
            try  ### Compute the AC power flow solution first!
                compute_ac_pf!(network)  
                Ks[name],Kinvs[name] = compute_power_encodings(network) 
                pf = calc_basic_power_factor(network)
                
                plot!(plt_K,
                    pf,Ks[name],
                    #label=nothing,
                    label=name[5:end-2],
                    #label=L"$K(\alpha)$ "*name,
                    #color_palette = palette(:Greens),
                    seriestype = :scatter
                    )
                xlabel!(L"Power Factor $\alpha_i$")
                ylabel!(L"$K_{ii} = \pm \frac{1}{\alpha_i} \sqrt{1-\alpha_i^2}$")
                plot!(plt_inv_K,
                    pf,Kinvs[name],
                    #label=nothing,
                    label=name[5:end-2],
                    #label=L"$K^{-1}(\alpha)$ "*name,
                    #color_palette = palette(:Oranges),
                    seriestype = :scatter
                    )
                xlabel!(L"Power Factor $\alpha_i$")
                ylabel!(L"$K^{-1}_{ii} = \pm \frac{\alpha_i}{\sqrt{1-\alpha_i^2}}$")

               
            catch
                println("AC Power Flow solution failed for: ",name)
                continue
            end

        end
    end
    plt = plot(plt_K,plt_inv_K)
    # for name in names
    #     #Draw the opnorm
    #     Ki = Diagonal(Ks[name]) #construct K
    #     λmax = opnorm(Ki)
    #     vline(λmax,
    #         #linewidth=3,ls="dashed"
    #         )
    # end
    return plt

end



#----


network_data_path = "/home/sam/github/PowerSensitivities.jl/data/figure_cases/" #Folder with radial-only systems
names = readdir(network_data_path);
paths = readdir(network_data_path,join=true);


plots = [] 

Ks,Kinvs = Dict(),Dict()

for (name,path) in zip(names,paths)
    network =  try
        make_basic_network(parse_file(path));
    catch
        println("PM cannot parse "*name)
        continue
    end
    if PowerSensitivities.is_radial(network)
        try  ### Compute the AC power flow solution first!
            compute_ac_pf!(network)  
        catch
            println("AC Power Flow solution failed for: ",name)
            continue
        end
            Ks[name],Kinvs[name] = compute_power_encodings(network)
            encoding_plot = plot_network_encoding(network)
            #Add title
            title!(name)
            push!(plots,encoding_plot)
            savefig(encoding_plot,"/home/sam/github/PowerSensitivities.jl/figures/spring_22/K_matrices/" *"K_"* name[1:end-2]*".pdf")
            savefig(encoding_plot,"/home/sam/github/PowerSensitivities.jl/figures/spring_22/K_matrices/" *"K_"* name[1:end-2]*".png")		
    
    end


end
