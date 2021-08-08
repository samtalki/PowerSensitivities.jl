using PowerSystems
using TimeSeries
using Dates
DATA_DIR = "../../data" #hide
system = System(joinpath(DATA_DIR, "matpower/case5.m"))

function qsts(system::System)
    resolution = Dates.Hour(1)
    dates = range(DateTime("2020-01-01T00:00:00"),step=resolution,length=24)
    data = TimeArray(dates,ones(24))
    time_series = SingleTimeSeries("max_active_power",data)
    return time_series    
end


new_renewable = RenewableDispatch(
        name = "WindPowerNew",
        available = true,
        bus = get_component(Bus, system, "3"),
        active_power = 2.0,
        reactive_power = 1.0,
        rating = 1.2,
        prime_mover = PrimeMovers.WT,
        reactive_power_limits = (min = 0.0, max = 0.0),
        base_power = 100.0,
        operation_cost = TwoPartCost(22.0, 0.0),
        power_factor = 1.0
    )

add_component!(system, new_renewable)

ts_data = [0.98, 0.99, 0.99, 1.0, 0.99, 0.99, 0.99, 0.98, 0.95, 0.92, 0.90, 0.88, 0.84, 0.76,
           0.65, 0.52, 0.39, 0.28, 0.19, 0.15, 0.13, 0.11, 0.09, 0.06,]
time_stamps = range(DateTime("2020-01-01"); step = Hour(1), length = 24)
time_series_data_raw = TimeArray(time_stamps, ts_data)
time_series = SingleTimeSeries(name = "active_power", data = time_series_data_raw)

#Add the forecast to the system and component
add_time_series!(system, new_renewable, time_series)
