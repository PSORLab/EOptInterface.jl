#
# Smallest DMC example
#
# This script shows the public DMC path:
# 1. create a numeric step response,
# 2. register the DMC block in JuMP,
# 3. inspect the created variables,
# 4. update the fixed measurement and history values.
#
# This is the first DMC script a new user should read.
# It does not use a full MTK ODE or DAE model.
# Instead, it starts from the simpler DMC idea:
# - the user already has a numeric step response;
# - the package turns that response into JuMP prediction equations.

import Pkg

Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using EOptInterface
using JuMP

step_response = [0.35, 0.62, 0.8, 0.9, 0.95]
P = 4
M = 2

# Plain JuMP model for the DMC block.
# DMC registration only needs an ordinary JuMP model.
# It does not need the tracking-MPC setup used for mechanistic ODE/DAE models.
model = Model()

# Register the DMC prediction block.
# After this call:
# - `ctx.u` holds the future manipulated-input variables;
# - `ctx.y_pred` holds the predicted output variables;
# - `ctx.y_meas` and `ctx.u_hist` store the fixed current data.
ctx = register_dmcsystem(
    model,
    step_response,
    P;
    M = M,
    y_meas0 = [0.4],
    u_hist0 = fill(0.3, length(step_response) - 1),
    u_bounds = (0.0, 1.0),
    base_name = "demo_dmc",
)

println("u count: ", length(ctx.u))
println("prediction size: ", size(ctx.y_pred))
println("move count: ", length(ctx.du))
println("stored in model.ext[:dmc]: ", haskey(model.ext, :dmc))

# Update the current measured output and input history.
# This is the normal online DMC update step.
# The prediction structure stays the same.
# Only the newest measured values are replaced.
update_dmc_state!(ctx; y_meas = [0.55], u_hist = fill(0.45, length(step_response) - 1))
println("updated y_meas fixed values: ", [JuMP.fix_value(v) for v in ctx.y_meas])
println("updated u_hist fixed values: ", [JuMP.fix_value(v) for v in ctx.u_hist])
