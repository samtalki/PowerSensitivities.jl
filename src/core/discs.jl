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

