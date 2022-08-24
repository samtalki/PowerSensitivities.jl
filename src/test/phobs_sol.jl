include("phobs.jl")

#########
#Radial cases
#Test case path and parameters
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
names = readdir(network_data_path);
paths = readdir(network_data_path,join=true);

begin
		
	#Results dicts for automatic indexing.
	#PQ bus only, unsuitable indexes removed.
	radial_results = Dict()
	meshed_results = Dict()

	#Eigenvalue and other radial_plots
	radial_plots = Dict()
	meshed_plots = Dict()

	#Save PowerModels networks
	radial_network_dicts = Dict() 
	meshed_network_dicts = Dict()

	#Larger minimum eigenvalue
	radial_larger_min_eigval = Dict()
	meshed_larger_min_eigval = Dict()


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
	    if PowerSensitivities.is_radial(network) #Radial network tests
	        try  ### Compute the AC power flow solution first!
	            compute_ac_pf!(network)  
	        catch
	            println("AC Power Flow solution failed for: ",name)
	            continue
	        end

	        #Check observablity
			radial_results[name] = test_observability(network)
            
            #Plot eigenvalues, save additional experiment data
            radial_plots[name] = plot_eigvals(network)
			try
				#savefig(radial_plots[name],"/home/sam/github/PowerSensitivities.jl/figures/spring_22/ph_observability/" * name[1:end-2]*"_eigplot.pdf")
				#savefig(radial_plots[name],"/home/sam/github/PowerSensitivities.jl/figures/spring_22/ph_observability/" * name[1:end-2]*"_eigplot.png")		
			
				#Save network
				radial_network_dicts[name] = network
		
				#Save the largest minimum eigenvalue
				if radial_results[name] !== nothing
					eigmin_A,eigmin_B = radial_results[name]["eigmin_A"],radial_results[name]["eigmin_B"]
					if eigmin_A>eigmin_B
						radial_larger_min_eigval[name] = Dict(
							"matrix" =>"A",
							"eigmin" => eigmin_A
						)
					elseif eigmin_B > eigmin_A
						radial_larger_min_eigval[name] = Dict(
							"matrix" => "B",
							"eigmin" => eigmin_B
						)
					else
						radial_larger_min_eigval[name] = nothing
					end
				end
			catch
				continue
			end
		else #Meshed network tests
			try  ### Compute the AC power flow solution first!
	            compute_ac_pf!(network)  
	        catch
	            println("AC Power Flow solution failed for: ",name)
	            continue
	        end
			try
				#Check observablity
				meshed_results[name] = test_observability(network)
				
				#Plot eigenvalues, save additional experiment data
				meshed_plots[name] = plot_eigvals(network)
				#savefig(meshed_plots[name],"/home/sam/github/PowerSensitivities.jl/figures/spring_22/ph_observability/" * name[1:end-2]*"_eigplot.pdf")
				#savefig(meshed_plots[name],"/home/sam/github/PowerSensitivities.jl/figures/spring_22/ph_observability/" * name[1:end-2]*"_eigplot.png")		
			
				#Save network
				meshed_network_dicts[name] = network
		
				#Save the largest minimum eigenvalue
				if meshed_results[name] !== nothing
					eigmin_A,eigmin_B = meshed_results[name]["eigmin_A"],meshed_results[name]["eigmin_B"]
					if eigmin_A>eigmin_B
						meshed_larger_min_eigval[name] = Dict(
							"matrix" =>"A",
							"eigmin" => eigmin_A
						)
					elseif eigmin_B > eigmin_A
						meshed_larger_min_eigval[name] = Dict(
							"matrix" => "B",
							"eigmin" => eigmin_B
						)
					else
						meshed_larger_min_eigval[name] = nothing
					end
				end
			catch
				continue
			end
		end
	end

end


# #----- Plots

# begin
	
# 	for (name,path) in zip(names,paths)
# 		plot_eigvals()

# 	end
# end