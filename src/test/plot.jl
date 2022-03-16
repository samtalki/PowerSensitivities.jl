using Gadfly

function plot_jacobian(J::AbstractMatrix)
    p1 = spy(J.sqth,Guide.xlabel("Bus Index"),Guide.ylabel("Bus Index"))
	p2 = spy(J.spth,Guide.xlabel("Bus Index"),Guide.ylabel("Bus Index"))
	title(hstack(p1,p2),"Estimating Angle Submatrices")
end

