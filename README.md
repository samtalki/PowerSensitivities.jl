This package is under development.

## Reproducible Research
- Analytical Results (Julia, PowerModels.jl):
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
- Numerical Results (Python, CVXPY):
    - Matrix Recovery
    ```
    py/SensitivityModelingCVXPY.ipynb
    ```
    - Estimation of complex power from voltage magnitudes
    ```
    py/ImplSMatrix.ipynb
    ```


