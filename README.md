The goal of this package is to help with differentiating electric power system data and studying the results.

This package is under development.

## Installation

### Clone the repository
```
git clone --recurse-submodules
```

### Activate the Julia environment
```
pkg> activate .
```

## Reproducible Research
We invite you to reproduce the results of our research through the following steps. 

### Analytical Results (Julia):

After activating the Julia environment, reproduce the analytical results with the corresponding commands in the Julia REPL:
- Assumption 1 and 2 validity:
```
include("src/test/assum1.jl")
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

### Numerical Results (Python):
Use your favorite Python environment with the following packages:
- Numpy
- Numba
- CVXPY
- Matplotlib
- Seaborn
- Pandas
- OpenDSSDirect

Reproduce the numerical results with the following notebooks:
- Matrix Recovery
```
py/SensitivityModelingCVXPY.ipynb
```
- Complex deviation estimation
```
py/ImplSMatrix.ipynb
```


