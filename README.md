# PowerSensitivities.jl
Tools for modeling and analyzing the various Jacobian matrices of electric power networks. 


This software is under development. It is currently compatible with *basic network models* in [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl).

## Setup
In Julia, run the following
```
import Pkg
Pkg.add("https://github.com/samtalki/PowerSensitivities.jl.git")
```

## Paper and Reproducibility

This package corresponds to paper "Conditions for Estimation of Sensitivities of Voltage Magnitudes to Complex Power Injections", S. Talkington, D. Turizo, S. Grijalva, J. Fernandez, and D. K. Molzahn, arXiv:2212.01471 [eess.SY] 2022.

[arxiv preprint](https://arxiv.org/abs/2212.01471)

If you find the ideas in this package or this paper useful, please consider adding the following citation
```
@misc{https://doi.org/10.48550/arxiv.2212.01471,
  doi = {10.48550/ARXIV.2212.01471},
  
  url = {https://arxiv.org/abs/2212.01471},
  
  author = {Talkington, Samuel and Turizo, Daniel and Grijalva, Santiago and Fernandez, Jorge and Molzahn, Daniel K.},
  
  keywords = {Systems and Control (eess.SY), Optimization and Control (math.OC), FOS: Electrical engineering, electronic engineering, information engineering, FOS: Electrical engineering, electronic engineering, information engineering, FOS: Mathematics, FOS: Mathematics},
  
  title = {Conditions for Estimation of Sensitivities of Voltage Magnitudes to Complex Power Injections},
  
  publisher = {arXiv},
  
  year = {2022},
  
  copyright = {Creative Commons Attribution 4.0 International}
}
```

The results in this research can be reproduced through the following steps. 

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


