
mutable struct PerturbSensitivityModel
    model::LinearSensitivityModel
    sp::AbstractArray #Real power voltage magnitude sensitivity matrix
    sq::AbstractArray #Reactive power " "
    spq::AbstractArray #Interleaved matrix
    ΔV::AbstractArray #Estimated Voltage Deviations
    function PerturbSensitivityModel(system_model::System)
        model = LinearSensitivityModel(system_model)
        sp = calc_dv_dp_matrix(system_model)
        sq = cal_dv_dq_matrix(system_model)
        spq = calc_sens_mat_perturb(system_model)
        return new(model,sp,sq,spq,ΔV) 
    end
end

function calc_sens_mat_perturb(sensors,injections,system_model::System)
    """Makes a perturb-and-observe sensitivity matrix for a model"""
    M = size(sensors)
    L = size(injections)
    S = zeros(M,L)
    V_0 = get_voltages
    for (sensor,injection) in zip(sensors,injections)
        1+2
        #voltages_inj = injection_voltages()
    end
end
