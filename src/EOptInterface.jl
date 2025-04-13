module EOptInterface

# ---- imports ----
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

using JuMP

# ---- includes ----
include("basefuncs.jl")
include("userfuncs.jl")

# ---- exports ----
export decision_vars, full_solutions, register_nlsystem, register_odesystem

end # module
