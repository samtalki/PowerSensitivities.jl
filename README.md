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

## Reproducible Research
The following scripts are included for reproducing the results from "Conditions for Estimation of Sensitivities of Voltage Magnitudes to Complex Power Injections" by Talkington, Turizo, Grijalva, Fernandez, and Molzahn:

- Analytical Results (Julia, PowerModels.jl):
    - Assumption 1 and 2 validity:
    ```
    include("src/test/assum1.jl")
    ```
    - Theorem results at default operating point:
    ```
    include("src/test/thm1_default.jl")
    ```
    - Theorem results at AC power flow solution:
    ```
    include("src/test/thm1_sol.jl")
    ```
    - Power Factor Bounds:
    ```
    include("src/test/pf_min.jl")
    ```
    - Jacobian eigenvalues:
    ```
    include("src/test/eigs_jac.jl")
    ```
- Numerical Results (Python, CVXPY):
    - Matrix Recovery
    ```
    py/SensitivityModelingCVXPY.ipynb
    ```
    - Estimation of complex power from voltage magnitudes
    ```
    py/ImplSMatrix.ipynb
    ```


## Todo

- Port matrix completion/recovery algorithms to native Julia.
- Complete automatic differentiation interfaces.




