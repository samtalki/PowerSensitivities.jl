# PowerSensitivities.jl
Tools for measurement-based modeling of electric power networks.

## Setup
In Julia,

## Reproducible Research
The following scripts are included for reproducing the results from "Conditions for Estimation of Sensitivities of Voltage Magnitudes to Complex Power Injections":

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




