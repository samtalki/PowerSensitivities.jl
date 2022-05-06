# PowerSensitivities.jl
Tools for measurement-based modeling of electric power networks. 

The key goals of this package are:
1. Provide an extensible and useful set of fundamental tools for *approximating the power flow equations from measurement data*. 
2. Provide a platform for applications with reduced reliance on the electric power network model.
3. Provide interfaces for multiple different power system measurements, including advanced metering infrastructure (AMI) data and phasor measurement unit (PMU) data.
4. Pursue interoperability with existing model-based frameworks, including [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) and [PowerSystems.jl]https://github.com/NREL-SIIP/PowerSystems.jl).

This software is under heavy development. It is currently compatible with *basic network models* in [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl).

## Setup
In Julia, run the following
```
import Pkg
Pkg.add("https://github.com/samtalki/PowerSensitivities.jl.git")
```