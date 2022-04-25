module dssSensitivity

using OpenDSSDirect

struct PerturbSensitivityModel
    ∂vp::AbstractMatrix
    ∂vq::AbstractMatrix
end

filename = Pkg.dir("OpenDSSDirect", "examples", "8500-Node", "Master.dss")
dss("""
    clear
    compile $filename
    solve
""")

loadnumber = Loads.First()
kWsum = 0.0
kvarsum = 0.0
while loadnumber > 0
    kWsum += Loads.kW()
    kvarsum += Loads.kvar()
    loadnumber = Loads.Next()
end

function PerturbSensitivityModel(path_to_file::String)
    dss(
    """
    clear
    compile $path_to_file
    """
    )
end

end