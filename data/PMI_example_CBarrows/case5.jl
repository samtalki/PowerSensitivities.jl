
using PowerSystems
using PowerModelsInterface
using Ipopt
using Dates
using DataFrames

# Make a System
sys = System("case5.m")

# add time series for active power demand (note: this adds reactive power profiles from
# same csv, but PowerModelsInterface doesn't use them, it only uses the active power time series)
add_time_series!(sys, "timeseries_pointers_da.json")

pm_map = PowerModelsInterface.get_pm_map(sys)

buses = get_components(Bus, sys)
T = get_forecast_horizon(sys)
results = Dict(
    zip(
        get_name.(buses),
        [
            Dict("va" => zeros(T), "vm" => zeros(T), "pg" => zeros(T), "qg" => zeros(T)) for b = 1:length(buses)
        ],
    ),
)
for t = 1:get_forecast_horizon(sys)
    #run the PF
    soln = run_pf(
        sys,
        ACPPowerModel,
        Ipopt.Optimizer,
        start_time = get_forecast_initial_timestamp(sys),
        period = t,
    )
    #extract the results
    for (id, b) in results
        b["va"][t] += soln["solution"]["bus"][id]["va"]
        b["vm"][t] += soln["solution"]["bus"][id]["vm"]
        for (gid, g) in pm_map["gen"]
            if get_name(get_bus(g)) == id
                b["pg"][t] += soln["solution"]["gen"][gid]["pg"]
                b["qg"][t] += soln["solution"]["gen"][gid]["qg"]
            end
        end
    end
end

#convert results into a dataframe
results_dfs = Dict()
for (b, r) in results
    results_dfs[b] = DataFrame(r)
end

results_dfs["1"]

# for the OPF approach, we can just run one multi network OPF
mn_soln = run_mn_opf(
    sys,
    ACPPowerModel,
    Ipopt.Optimizer,
    start_time = get_forecast_initial_timestamp(sys),
    time_periods = 1:T,
)

#TODO: extract the data from mn_solution
