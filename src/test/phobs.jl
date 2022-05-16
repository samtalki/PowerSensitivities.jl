###
#Test Th2 at the AC power flow solution.
###
include("../PowerSensitivities.jl")
include("../test/thm1.jl")
include("../util/matrix.jl")
using PowerModels
using LinearAlgebra
using Plots
using LaTeXStrings
import .PowerSensitivities

# function is_pos_def(A::AbstractMatrix)
# 	#Hermitian part test for positive definitineness of real/complex matrices
# 	Aherm = (1/2)*(A + conj(transpose(A))
# 		)
# 	return eigmin(Aherm) > 0
# end


struct ObservabilityData
	eigmin_A::Union{Real,Complex}
	eigmin_B::Union{Real,Complex}
	A::AbstractMatrix
	B::AbstractMatrix
	K::AbstractMatrix
	vp::AbstractMatrix
	vq::AbstractMatrix
	observable::Bool
end

function test_observability(data::Dict{String,<:Any})
	K = PowerSensitivities.calc_K_matrix(data)
	S = PowerSensitivities.calc_voltage_sensitivity_matrix(data)
	@assert size(K) == size(S.vp) == size(S.vq) "Mismatch in sensitivity dimensions"
	
	#Test matrices
	A = S.vp + S.vq*K
	B = S.vp*inv(K) + S.vq
	is_observable = try
		eigmin(A)>0 || eigmin(B)>0#Test
	catch
		is_observable = false
	end
	#is_observable = is_pos_def(A) || is_pos_def(B)
	if is_observable
		return ObservabilityData(
			eigmin(A),
			eigmin(B),A,B,K,
			S.vp,S.vq,is_observable
			)
	else
		return false
	end
end

function calc_ac_sol_data!(data::Dict{String,<:Any})
	compute_ac_pf!(data)
	Y = calc_basic_admittance_matrix(data)
	s = calc_basic_bus_injection(data)
	v = calc_basic_bus_voltage(data)
	p,q = real.(s),imag.(s)
	J = PowerSensitivities.calc_jacobian_matrix(data)
	S = PowerSensitivities.calc_voltage_sensitivity_matrix(data)
	return Dict(
		"Y" => Y,
		"p" => p,
		"q" => q,
		"J" => J,
		"S" => S
		)
end

"""
Given a network data dict, plot the eigenvalues of the phobs matrices
"""
function plot_eigvals(data::Dict{String,<:Any})
	K = PowerSensitivities.calc_K_matrix(data)
	S = PowerSensitivities.calc_voltage_sensitivity_matrix(data)
	ξ = PowerSensitivities.calc_basic_power_factor(data)
	
	#Phaseless observability matrices
	A = S.vp + S.vq*K
	B = S.vp*inv(K) + S.vq

	pA = plot_eigvals(A)
	# title!(L"$\text{eigs}(\frac{
	# 	\partial \boldsymbol{v}
	# 	}{
	# 	\partial \boldsymbol{p}
	# 	} 
	# 	+ \frac{
	# 	\partial \boldsymbol{v}
	# 	}{
	# 	\partial \boldsymbol{q}
	# 	} 
	# 	\boldsymbol{K}  
	# 	)$")
	pB = plot_eigvals(B)
	# title!(L"$\text{eigs}(
	# 	\frac{
	# 	\partial \boldsymbol{v}
	# 	}{
	# 	\partial \boldsymbol{p}
	# 	}
	# 	\boldsymbol{K}^{-1} 
	# 	+ \frac{
	# 	\partial \boldsymbol{v}
	# 	}{
	# 	\partial \boldsymbol{q}
	# 	}   
	# 	)$")
	p = plot(pA,pB)
	return p
	#title!("Phaseless Observabilities")
	#case_name = data["name"]
	#savefig(p,"/home/sam/github/PowerSensitivities.jl/figures/spring_22/" * case_name*"_eigplot.pdf")
	#savefig(p,"/home/sam/github/PowerSensitivities.jl/figures/spring_22/" * case_name*"_eigplot.png")
end

function plot_eigvals!(data::Dict{String,<:Any})
	compute_ac_pf!(data)
	return plot_eigvals(data)
end

function plot_eigvals(A::AbstractMatrix;sorted=true)
	if sorted
		λ = sort(
			real.(
				eigvals(
					Matrix(A)
					)
				)
			)
	else 
		λ = real.(eigvals(Matrix(A)))
	end
	return plot(
		λ,
		xlabel="Bus Index",
		ylabel="Eigenvalue"
		)
end

#########
#Radial cases
#Test case path and parameters
network_data_path = "/home/sam/github/PowerSensitivities.jl/data/radial_test/" #Folder with radial-only systems
names = readdir(network_data_path);
paths = readdir(network_data_path,join=true);

begin
		
	#Results dicts for automatic indexing.
	#PQ bus only, unsuitable indexes removed.
	results = Dict()

	#Eigenvalue and other plots
	plots = Dict()

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

	        #Check observable
            results[name] = test_observability(network)
            
            #Plot eigenvalues
            plots[name] = plot_eigvals(network)
            savefig(plots[name],"/home/sam/github/PowerSensitivities.jl/figures/spring_22/"*name*".pdf")
            savefig(plots[name],"/home/sam/github/PowerSensitivities.jl/figures/spring_22/"*name*".png")
            

            #Save network
            radial_network_dicts[name] = network
	    end
	end

end


# #----- Plots

# begin
	
# 	for (name,path) in zip(names,paths)
# 		plot_eigvals()

# 	end
# end