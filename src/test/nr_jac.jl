##### Test the jacobian solutions and the voltage sensitivities
using ForwardDiff
using LinearAlgebra
import PowerSensitivities as PS
include("/home/sam/github/PowerSensitivities.jl/src/sens/nr_sens.jl")




"""
Given a Jacobian J and vector of sel_bus_idx, return a sliced Jacobian where all 4 blocks are sliced according to sel_bus_idx.
"""
function slice_jacobian(J::AbstractMatrix,sel_bus_idx::AbstractArray)
    num_bus = Int(size(J,1)/2) #number of buses
    J_sel_idx = [sel_bus_idx; sel_bus_idx .+ num_bus] #Shift up matrix indeces to cover all blocks
    return J[J_sel_idx, J_sel_idx]
end

"""
Given a network data dict test the accuracy of the autodiff jacobian
"""
function test_AD_jacobian(data::Dict{String,<:Any})

    J_base = PM.calc_basic_jacobian_matrix(data)  #Calculate the analytical Jacobian using the power flow equations 
    s_base = PM.calc_basic_bus_injection(data)
    v_base = PM.calc_basic_bus_injection(data) #Rectangular bus voltages v_r + j v_i at basecase
    v_ph_base = calc_phasor_bus_voltage(data)  #Compute the phasor voltages  [θ ; vmag] at basecase

    #Solve the AC Power flow solution and calculate solution data
    PM.compute_ac_pf!(data)
    J_sol = PM.calc_basic_jacobian_matrix(data) #Calculate the analytical Jacobian using the power flow equations 
    s_sol = PM.calc_basic_bus_injection(data)
    v_sol = PM.calc_basic_bus_voltage(data) #Rectangular bus voltages v_r + j v_i at the solution
    v_ph_sol = calc_phasor_bus_voltage(data) #Compute the phasor voltages  [θ ; vmag] at the solution

    #Get the PQ and PQ+PV bus indeces
    pq_idx = PS.calc_bus_idx_of_type(data,[1])
    pq_pv_idx = PS.calc_bus_idx_of_type(data,[1,2])

    #Calculate the PFJ with SensitivityModel
    model = SensitivityModel(data) #Create a sensitivity model
    J_auto = model.J(v_ph_sol)

    # Check they are the same on the specified bus types (PQ only)
    return norm(slice_jacobian(J_sol,pq_idx)-slice_jacobian(J_auto,pq_idx)) < 1e-3
end




#----------------------------------
#Vmag sens



#----
#Test on all Radial cases
#Test case path and parameters
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
names = readdir(network_data_path);
paths = readdir(network_data_path,join=true);

begin
        
    #Results dicts for automatic indexing.
    #PQ bus only, unsuitable indexes removed.
    results = Dict()

    #Save PowerModels networks
    radial_network_dicts = Dict() 

    # """
    # Radial cases
    # Test case path and parameters
    # """
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

            #Check autodiff accuracy
            results[name] = test_AD_jacobian(network)
            
            #Save network
            radial_network_dicts[name] = network
        end
    end

end