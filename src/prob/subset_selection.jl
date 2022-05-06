# Find the subset of bus indeces that satify Theorm 1
using JuMP
using Ipopt


"""
Given a network data dict, calculate the largest subset of buses that satify Theorem 1.
"""
function calc_bus_subset(network::Dict)

    function bus_subset_model(network)
        n_bus = length(network["bus"])
        model = Model(Ipopt.Optimizer);
        @variable(model, 0 <= n_sel <= n_bus);
        @objective(model, Max, n_sel);
        @constraint(model, pf_max>=pf_min);
        @NLconstraint(model, (sqrt(1-pf_min^2)/pf_min) - (sqrt(1-pf_max^2)/pf_max) <= Δk_max  );
    end

    
    optimize!(model);
    pf_min,pf_max = value(pf_min),value(pf_max);
    Δpf_max = pf_max - pf_min;
    return Δpf_max
end