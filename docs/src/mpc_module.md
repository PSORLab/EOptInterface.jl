# Tracking MPC With EOptInterface

This page explains the ModelingToolkit-to-JuMP MPC workflow used by the
examples in this repository. It is written as a reproducible research note:
what model is used, what the controller changes, what output is tracked, and
what gets updated at each sampling time.

## Basic Workflow

For a mechanistic plant model, the normal workflow is:

1. write the plant in `ModelingToolkit`;
2. choose manipulated variables, such as aeration or flow rate;
3. choose tracked outputs and setpoints;
4. build one JuMP MPC problem;
5. during simulation, update the current state and solve again.

The main functions are:

- `build_tracking_mpc(...)`: build the JuMP prediction and tracking problem;
- `solve_tracking_mpc!(...)`: update the current state, solve, and read back the first control move;
- `update_stage_parameter!(...)`: update a disturbance or feed preview;
- `update_tracking_targets!(...)`: change setpoints or soft bounds during a run.

The DMC step-response path is separate and starts from `register_dmcsystem(...)`.

## Minimal Example

```julia
using EOptInterface
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using JuMP, Ipopt

@parameters u = 0.0 d = 1.0
@variables x(t)
@named sys = ODESystem([D(x) ~ -0.8 * x + u + d], t, [x], [u, d])

model = Model(Ipopt.Optimizer)
set_silent(model)

cfg = TrackingMPCConfig(
    PH = 5,
    CH = 2,
    dt = 1.0,
    integrator = "IE",
    system_kind = :ode,
    state_lower = -5.0,
    state_upper = 5.0,
)

ctrl = build_tracking_mpc(
    model,
    sys;
    control_specs = [
        MPCControlSpec(
            sym = sys.u,
            lower = 0.0,
            upper = 2.0,
            delta_max = 0.4,
            move_weight = 0.1,
            first_move_weight = 1.0,
        ),
    ],
    output_specs = [
        MPCOutputSpec(
            sym = sys.x,
            setpoint = 1.0,
            track_weight = 1.0,
            terminal_weight = 2.0,
            lower_soft = 0.9,
            upper_soft = 1.1,
            slack_weight = 100.0,
        ),
    ],
    stage_param_defaults = Dict(sys.d => fill(1.0, cfg.PH + 1)),
    config = cfg,
)

result = solve_tracking_mpc!(ctrl, Dict(sys.x => 0.0), Dict(sys.u => 0.0))
u_apply = result.controls[sys.u][1]
```

`result.controls` contains the optimized input trajectory.
`result.predictions` contains the predicted states.
`result.metrics` contains objective terms that are useful for reports.

## Meaning of the Main Settings

`MPCControlSpec` describes one manipulated variable.

- `sym`: the ModelingToolkit symbol used as the control input;
- `lower`, `upper`: hard bounds;
- `delta_max`: maximum move between samples;
- `move_weight`: penalty on input movement over the horizon;
- `first_move_weight`: penalty on changing from the previously applied input.

`MPCOutputSpec` describes one tracked output.

- `sym`: state or algebraic quantity to track;
- `setpoint`: target value;
- `track_weight`: tracking-error weight;
- `terminal_weight`: extra weight on the final predicted point;
- `lower_soft`, `upper_soft`, `slack_weight`: optional soft operating zone.

`TrackingMPCConfig` holds the numerical settings.

- `PH`: prediction horizon length;
- `CH`: number of free control moves;
- `dt`: sample time;
- `integrator`: transcription rule, such as `"IE"` or `"RK4"`;
- `system_kind`: `:ode` or `:dae`;
- `state_lower`, `state_upper`: default state bounds.

## Closed-Loop Use

Inside a simulation callback or a hand-written loop, the repeated step is:

```julia
state_values = current_state_map(integ, sys)
update_stage_parameter!(ctrl, sys.d, disturbance_forecast)
update_tracking_targets!(ctrl; setpoints = Dict(sys.x => 1.2))
result = solve_tracking_mpc!(ctrl, state_values, Dict(sys.u => integ.ps[sys.u]))
u_apply = result.controls[sys.u][1]
```

The JuMP model is not rebuilt at each sample. Only the current state,
previously applied input, setpoints, and previews are updated.

For a simple live line while the controller runs:

```julia
result = solve_tracking_mpc!(
    ctrl,
    state_values,
    previous_controls;
    show_status = true,
)
```

The NDMC notebook adds its own status line with the variables used in that
case. That extra printing is part of the example, not a different controller.

## DAE and DMC Notes

Use `TrackingMPCConfig(system_kind = :dae)` when the prediction model includes
algebraic equations that should remain in the optimization problem. The small
DAE example is `examples/dae_registration_demo.jl`.

Use `register_dmcsystem(...)` when the plant has already been reduced to step
responses and a classic DMC formulation is enough. The small DMC example is
`examples/dmc_registration_demo.jl`.

## Suggested Reading Order

For the tracking-MPC path:

1. `examples/tracking_mpc_demo.jl`;
2. this page;
3. `src/trackingmpc.jl`;
4. `examples/ndmc_conductivity_mpc_demo.jl` or the NDMC notebook.

The NDMC case now has one repository implementation:
`examples/NDMCExample.jl` loads the model and simulation code from
`examples/ndmc_case.jl` and the plotting code from `examples/ndmc_plots.jl`.
The repository notebook,
`notebooks/ndmc_conductivity_mpc_simple.ipynb`, runs the same case.

For debugging a failed solve, the helper routines in `src/mpcutils.jl` can
print constraints, check initial-condition synchronization, and record MPC
predictions. They are supporting tools, not the main controller definition.
