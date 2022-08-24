module PowerSensitivities
using PowerModels
using JuMP, Ipopt 
using TimeSeries
using Statistics
using LinearAlgebra
using SparseArrays
using Flux
PowerModels.silence()


#AMI Functionality
include("data/ami.jl")


#Sensitivity functions
include("sens/voltage.jl")
include("sens/angles.jl")

#Core functions
include("core/power_factor.jl")
include("core/jacobian_matrix.jl")
include("core/spectral_analysis.jl")

#Data utilities
include("util/network.jl")
include("util/matrix.jl")
include("util/bus_index.jl")

#Problem statements
include("prob/matrix_completion.jl")


#--- Jacobian-like matrix utilities
#- Classical AC Power flow Jacobian utilities
export PowerFlowJacobian, calc_jacobian_matrix
export calc_pth_jacobian, calc_qth_jacobian
export calc_pv_jacobian, calc_qv_jacobian
#- Voltage Sensitivity Matrix (underdetermined inverse AC Power Flow Jacobian) utilities
export VoltageSensitivityMatrix
export calc_voltage_sensitivity_matrix
#- Spectral Analysis utilities
export SpectralAnalysis
export calc_spectral_analysis

#--- Network data utilities
export make_ami_dataset,make_timeseries_dataset,calc_finite_differences
export set_network_load,is_radial

#--- Functions for testing the Theorems of Talkington and Turizo, et al.
export calc_basic_power_factor
export calc_K_matrix, calc_M_matrix
export calc_delta_k, calc_delta_k_max, calc_max_pf_distance
export calc_vmag_condition

#--- Matrix utilities
export ispd,isnd,isnsd,isinvertible,symmetricdiff,symmetric_part_nd,symmetric_part_nsd,symmetric_part_pd

#--- Special bus indexing utilities
export calc_bus_idx_of_type,calc_bad_idx,calc_study_idx
export calc_net_capacitive_idx,calc_net_inductive_idx

end
