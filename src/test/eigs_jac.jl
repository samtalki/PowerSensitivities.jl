#Compute the eigenvalues of all of the Jacobians at the AC power flow solutions
include("../PowerSensitivities.jl")
using .PowerSensitivities
using PowerModels:make_basic_network,parse_file
using LinearAlgebra,CSV,DataFrames
using Tables

"""
Given a network data dictionary, a selected bus type ([1] (pq only) or [1,2] (pq and pv))
Compute the jacobian of the network, all of its submatrices, 
and return the eigenvalues up to `num_eigvals` of the Jacobian and its submatrices as  DataFrames.
"""
function calc_jacobian_eigvals(network::Dict,sel_bus_types=[1],ϵ=1e-6)
    
    #Compute the indeces that will be considered
    study_idx = calc_bus_idx_of_type(network,sel_bus_types)
    
    #Study_idx
    n_total_bus,n_study_bus = length(network["bus"]),length(study_idx)
    jac_study_idx = [study_idx; study_idx .+ n_total_bus]
    
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
        eigs_submatrices[matrix_name] = reverse(eigvals(Matrix(M)))[1:end]
    end 

    return DataFrame(eigs_J),DataFrame(eigs_submatrices)
end

#Compute the leading eigenvalues for every test case
#Folder with Test case systems
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" 
allow_mesh = false #Whether to allow meshed test cases/require test feeder is radial
sel_bus_types = [1] #The bus types to cosnider (PQ by default)
#The maximum number of leading eigenvalues to show (will be less if the case has less buses than this divided by 2)
case_names = readdir(network_data_path);
paths = readdir(network_data_path,join=true);
#Dictionary of minimum eigenvalues for each test case
radial_minimum_eig = Dict()
for (case_name,case_path) in zip(case_names,paths)
    network = try
        make_basic_network(parse_file(case_path)); 
    catch
        println("PM cannot parse "*case_name)
        continue
    end
    if PowerSensitivities.is_radial(network) || allow_mesh
        #Calculate the eigenvalues of the jacobian and submatrices
        eigs_J,eigs_sub = calc_jacobian_eigvals(network,sel_bus_types)
        #Save as CSVs
        file_name = splitext(case_name)[1]
        CSV.write("src/test/results/eigs_jac/"*"J_"*file_name*".csv",eigs_J)
        CSV.write("src/test/results/eigs_jac/"*"sub_"*file_name*".csv",eigs_sub)

        #Calculate full jacobian dataframe and save as csv
        J = calc_jacobian_matrix(network,sel_bus_types)
        J_matrix = Matrix(J.matrix)
        CSV.write("src/test/results/radial/J_full/"*file_name*".csv",Tables.table(J_matrix),writeheader=false)
        
        #Save minimum eigenvalues of the jacobian
        radial_minimum_eig[case_name] = minimum(real.(eigs_J[!,"J"]))

    end
end

#Test for non-radial cases
####
allow_mesh = true #Whether to allow meshed test cases/require test feeder is radial
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/pm_matpower/" 
case_names = readdir(network_data_path);
paths = readdir(network_data_path,join=true);
#Dictionary of minimum eigenvalues for each test case
mesh_minimum_eig = Dict()
for (case_name,case_path) in zip(case_names,paths)
    network = try
        make_basic_network(parse_file(case_path)); 
    catch
        println("PM cannot parse "*case_name)
        continue
    end
    if PowerSensitivities.is_radial(network) || allow_mesh
        #Calculate the eigenvalues of the jacobian and submatrices
        eigs_J,eigs_sub = calc_jacobian_eigvals(network,sel_bus_types)
        #Save as CSVs
        file_name = splitext(case_name)[1]
        CSV.write("src/test/results/eigs_jac/"*"J_"*file_name*".csv",eigs_J)
        CSV.write("src/test/results/eigs_jac/"*"sub_"*file_name*".csv",eigs_sub)

        #Calculate full jacobian dataframe and save as csv
        J = calc_jacobian_matrix(network,sel_bus_types)
        J_matrix = Matrix(J.matrix)
        CSV.write("src/test/results/mesh/J_full/"*file_name*".csv",Tables.table(J_matrix),writeheader=false)
        
        #Save minimum eigenvalues of the jacobian
        mesh_minimum_eig[case_name] = minimum(real.(eigs_J[!,"J"]))

    end
end

# """
# Given a network data dict calculate the jacobain matrix as a DataFrame.
# """
# function calc_jacobian_df(network::Dict,sel_bus_types=[1])
#     #Calculate full Jacobian
#     J = calc_jacobian_matrix(network,sel_bus_types)
#     J_matrix = Matrix(J.matrix)
#     #cols = [J_j for J_j in eachcol(J_matrix)]
#     #bus_names = vcat(["bus_"*string(idx) for idx in 1:size(J_matrix,1)/2], ["bus_"*string(idx) for idx in 1:size(J_matrix,1)/2])
#     #J_df = DataFrame(columns=cols,names=bus_names)
#     return J_df
# end
