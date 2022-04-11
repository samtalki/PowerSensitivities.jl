#Compute the eigenvalues of case 69 and plot them
include("../PowerSensitivities.jl")
using .PowerSensitivities
using PowerModels:make_basic_network,parse_file
using LinearAlgebra
using DataFrames
using Plots,LaTeXStrings


"""
Given a network data dictionary, a selected bus type ([1] (pq only) or [1,2] (pq and pv))
Compute the jacobian of the network, all of its submatrices, 
and return the eigenvalues up to `num_eigvals` of the Jacobian and its submatrices as  DataFrames.
"""
function calc_jacobian_eigvals_unsorted(network::Dict,sel_bus_types=[1],ϵ=1e-6)
    
    #Compute the indeces that will be considered
    study_idx = calc_bus_idx_of_type(network,sel_bus_types)
    
    #Study_idx
    n_bus,n_study_bus = length(network["bus"]),length(study_idx)
    jac_study_idx = [study_idx; study_idx .+ n_bus]
    
    #Compute the matrices of interest
    J = calc_jacobian_matrix(network)
    ∂pθ,∂qθ = J.pth[study_idx,study_idx],J.qth[study_idx,study_idx] #Power-to angle sensitivities
    ∂pv,∂qv = J.pv[study_idx,study_idx],J.qv[study_idx,study_idx] #Power-to voltage sensitivities
    J = J.matrix[jac_study_idx,jac_study_idx] #full matrix

    #Check the sizes
    @assert size(∂pv,1) == length(study_idx) && size(∂pv,2) == length(study_idx)
    @assert size(∂qv,1) == length(study_idx) && size(∂qv,2) == length(study_idx)
    @assert size(∂pθ,1) == length(study_idx) && size(∂pθ,2) == length(study_idx)
    @assert size(∂qθ,1) == length(study_idx) && size(∂qθ,2) == length(study_idx)
    @assert size(J,1)/2 == size(∂qθ,1) == size(∂qθ,1) == size(∂pv,1) == size(∂qv,1) 
    @assert size(J,2)/2 == size(∂pθ,2) == size(∂qθ,2) == size(∂pv,2) == size(∂qv,2) 
    
    #Compute the eigenvalues of the Jacobian
    eigs_J = Dict("J" => eigvals(Matrix(J)))
    #Compute the eigenvalues of the Jacobian submatrices
    eigs_submatrices = Dict()
    for (M,matrix_name) in zip([∂pv,∂qv,∂pθ,∂qθ],["pv","qv","pth","qth"])
        eigs_submatrices[matrix_name] = eigvals(Matrix(M))
    end 

    return DataFrame(eigs_J),DataFrame(eigs_submatrices)
end

#Plot for case69
begin
    net = make_basic_network(parse_file("data/radial_test/case69.m"))
    eigs_J,eigs_submatrices = calc_jacobian_eigvals_unsorted(net)
    plot(Matrix(eigs_submatrices), labels=permutedims(names(eigs_submatrices)), legend=:topleft,marker=:auto,markeralpha=0.5,markersize=3)
    xlabel!("Bus "*L"i")
    ylabel!("Eigenvalue")
    title!("Case 69 Eigenvalues")
    savefig("figures/spring_22/sub_eigs_case69.pdf")
    plot(real.(eigs_J[!,"J"]), label= L"\Re(\lambda_i)", legend=:topleft)
    plot!(abs.(eigs_J[!,"J"]), label= L"\|\lambda_i \|")
    xlabel!("Bus "*L"i")
    ylabel!("Eigenvalue")
    title!("Case 69 Eigenvalues")
    savefig("figures/spring_22/full_eigs_case69.pdf")
end