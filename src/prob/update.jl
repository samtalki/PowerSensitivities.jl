function cholupdate(L, x)
    n = length(x);
    for k = 1:n
        r = sqrt(L(k, k)^2 + x(k)^2);
        c = r / L(k, k);
        s = x(k) / L(k, k);
        L[k, k] = r;
        if k < n
            L[(k+1):n, k] = (L[(k+1):n, k] + s * x((k+1):n)) / c;
            x[(k+1):n] = c * x[(k+1):n] - s * L((k+1):n, k);
        end
    end
end
