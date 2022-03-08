using Flux
using LinearAlgebra


function armijo(f,x_0,a,β,ϵ=1e-3,max_iter=1e2)
    """
    Steepest descent algorithm.
    Find the minimizer of a function f via gradient descent with the Armijo Stepsize.
        x_0: Initial condition
        a,β∈(0,1): Armijo parameters
        ϵ: Stopping/convergence criterion to bound ||∇f||
        max_iter: Maximum iterations before deemed failed to converge
    """
    ∇f(x) = gradient(f,x)[1]
    x = x_0
    i = 0 #Descent iterations
    k = 0 #Line search iterations
    f_ = []
    x_ = []
    push!(f_,f(x))
    push!(x_,x)
    while norm(∇f(x))>ϵ
        Δf = f(x - β^k*∇f(x)) - f(x)
        if Δf <= -a*β^k*norm(∇f(x))^2 # Stepsize selection criterion
            λ = -β^k 
            x = x + λ*∇f(x) #Gradient update
            i+=1; 
            k=0 #Restart line search iterations
            push!(f_,f(x))
            push!(x_,x)
        #elseif i>max_iter
        #    return x,f_,x_
        else
            k+=1
        end
    end
    return x,f_,x_
end