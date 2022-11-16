using Plots,LaTeXStrings
"""
Given a matrix A and integer i, return the elements A[i,k] for all k not equal to i
"""
function col_offdiag(A::AbstractMatrix,i::Integer)
	C = []
	for (k,Aᵢₖ) in enumerate(A[i,:])
		if i != k 
			push!(C,Aᵢₖ)
		end
	end
	return C
end

"""
Given a matrix A and integer i, return the elements A[k,i] for all k not equal to i
"""
function row_offdiag(A::AbstractMatrix,i::Integer)
	R = []
	for (k,Aₖᵢ) in enumerate(A[:,i])
		if i != k 
			push!(R,Aₖᵢ)
		end
	end
	return R
end
    
"""
Method 1: Given a matrix A and a diagonal index i construct the discs corresponding to the row/column of a square matrix i
"""
function gershdisc(A::AbstractMatrix,i::Integer)
    Aᵢᵢ = A[i,i]
    #Radii corresponding to the sum of the offdiagonal rowwise/colwise elements of A
    r_col = sum(abs.(col_offdiag(A,i)))
    r_row = sum(abs.(row_offdiag(A,i)))
    return Dict(
        "center" => Aᵢᵢ,
        "r_col" => r_col,
        "r_row" => r_row
    )
end

"""
Method 2: Given a matrix A return all of the gershgorin discs of A
"""
function gershdisc(A::AbstractMatrix)
    #Check if the matrix is square
    @assert size(A)[1] == size(A)[2]
    disc_params = []
    for i in 1:size(A)[1]
        push!(disc_params,gershdisc(A,i))
    end
    return disc_params
end



#----- Plotting functions for the gershdic

"""
Given a basic network data dict, plot the gershdiscs
"""
function plot_gershdisc(network::Dict;row_disc=true,col_disc=true,npts=500)
    J = calc_jacobian_matrix(network,[1,2])
    return plot_gershdisc(J,row_disc=row_disc,col_disc=col_disc,npts=npts)
end

"""
Given a basic network data dict, solve the AC power flow equations and plot the gershdiscs.
"""
function plot_gershdisc!(network::Dict;row_disc=true,col_disc=true,npts=500)
    compute_ac_pf!(network)
    return plot_gershdisc(network,row_disc=row_disc,col_disc=col_disc,npts=npts)
end

"""
Method 1: Given a matrix A and an index i, plot the gershgorin row/col gershdiscs centered at A[i,i].
"""
function plot_gershdisc(A,i;row_disc=true,col_disc=true,npts=500)
    #Construct and unpack disc parameters
    disc = gershdisc(A,i)
    center,r_col,r_row = disc["center"],disc["r_col"],disc["r_row"]
    #Set of all angles from 0-2π
    θ = LinRange(0,2*π,npts)
    #Center of the circle
    c_re,c_im = real(center),imag(center)
    if row_disc #The row disc ∈ R^2 or "∈ C"
        D_row = (c_re .+ r_row * sin.(θ), c_im .+ r_row * cos.(θ))
        p = plot(D_row, seriestype=[:shape],
            fillalpha=0.2,aspect_ratio=1,lw=0.5,
            #c=:blue,linecolor=:black,
            label=L"D_{r}("*string(i)*")")
    end
    if col_disc #The col disc ∈ R^2 or "∈ C"
        D_col = (c_re .+ r_col * sin.(θ), c_im .+ r_col * cos.(θ))
        plot!(p,D_row, seriestype=[:shape],
            fillalpha=0.2,aspect_ratio=1,lw=0.5,
            #c=:blue,linecolor=:black,
            label=L"D_{c}("*string(i)*")")
    end
    return p
end
"""
Method 2: Given a matrix A plot all of the Gershgorin discs of A
"""
function plot_gershdisc(A;row_disc=true,col_disc=true,npts=500)
    @assert size(A)[1] == size(A)[2]
    
    #Set of all disc params
    disc_params = gershdisc(A) 
    
    #Set of all angles from 0-2π
    θ = LinRange(0,2*π,npts)

    p_row = plot()
    p_col = plot()
    for (i,disc_i) in enumerate(disc_params) #For each disc parameters
        center,r_col,r_row = disc_i["center"],disc_i["r_col"],disc_i["r_row"]
        #Center of the circle
        c_re,c_im = real(center),imag(center)
        #For i: the row disc
        D_row = (c_re .+ r_row * sin.(θ), c_im .+ r_row * cos.(θ))
        #For i: the col disc ∈ R^2 or "∈ C"
        D_col = (c_re .+ r_col * sin.(θ), c_im .+ r_col * cos.(θ))
        if row_disc #For i: row disc
            plot!(p_row,D_row,seriestype=[:shape],
            fillalpha=0.2,aspect_ratio=1,lw=0.5,
            #c=:blue,linecolor=:black,
            label=L"D_{r}["*string(i)*"]",
            title="Row Discs",
            xlabel=L"\textsf{Re}(S_{ii})",ylabel=L"\textsf{Im}(S_{ii})")
        end
        if col_disc
            plot!(p_col,D_col,seriestype=[:shape],
            fillalpha=0.2,aspect_ratio=1,lw=0.5,
            #c=:blue,linecolor=:black,
            label=L"D_{c}["*string(i)*"]",
            title="Column Discs",
            xlabel=L"\textsf{Re}(S_{ii})",ylabel=L"\textsf{Im}(S_{ii})")
        end
    end
    
    return plot(p_row,p_col)
end

