#
# Small DAE registration example
#
# This script shows direct DAE registration:
# 1. build a tiny DAE that stays a DAE after `structural_simplify`,
# 2. register it with different integrators,
# 3. solve the resulting JuMP models,
# 4. save the summaries to CSV files.
#
# This example is for users who want to call `register_daesystem(...)` directly,
# without building a tracking-MPC experiment around it.
# The point is not process realism.
# The point is to show:
# - what kind of DAE the package expects,
# - how `register_daesystem(...)` is called,
# - and what kind of solved trajectories come back.

import Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using EOptInterface
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using JuMP, Ipopt
using DataFrames, CSV

const GENERATED_DIR = joinpath(@__DIR__, "generated")
const DEMO_TARGET = 0.8
const DEMO_DT = 0.2
const DEMO_HORIZON = 11

function classify_system(sys)
    eqs = ModelingToolkit.full_equations(sys)
    diff_eq_count = 0
    alg_eq_count = 0
    for eq in eqs
        lhs = eq.lhs
        is_diff = try
            Symbolics.operation(lhs) isa ModelingToolkit.Differential
        catch
            false
        end
        if is_diff
            diff_eq_count += 1
        else
            alg_eq_count += 1
        end
    end
    return (
        diff_eq_count = diff_eq_count,
        alg_eq_count = alg_eq_count,
        unknown_count = length(ModelingToolkit.unknowns(sys)),
        equation_count = length(eqs),
    )
end

function find_state_key(x_vars::AbstractDict, label::AbstractString)
    for sym in keys(x_vars)
        if endswith(string(sym), label)
            return sym
        end
    end
    error("Could not find state key ending with $(label). Available keys: $(collect(keys(x_vars)))")
end

function find_parameter_key(sys, label::AbstractString)
    for sym in ModelingToolkit.parameters(sys)
        if endswith(string(sym), label)
            return Num(sym)
        end
    end
    error("Could not find parameter key ending with $(label). Available keys: $(ModelingToolkit.parameters(sys))")
end

function build_demo_dae()
    # This DAE keeps one differential state and one algebraic state.
    # The algebraic equation is implicit, so `structural_simplify(...)` cannot
    # remove it and turn the whole model into a pure ODE.
    @parameters u = 0.0
    ModelingToolkit.@variables x(t) = 0.0 z(t) = 0.0
    @named dae_keep = System([
        D(x) ~ -x + z + u,
        0 ~ z + sin(z) - x,
    ], t)
    return dae_keep
end

function build_integrator_case(sys, integrator::String; target::Float64 = DEMO_TARGET, dt::Float64 = DEMO_DT, horizon::Int = DEMO_HORIZON)
    # Build one JuMP model for one integrator choice.
    # Each integrator gets its own fresh model so the comparison is clean.
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "tol", 1e-8)
    set_optimizer_attribute(model, "max_cpu_time", 30.0)
    set_optimizer_attribute(model, "print_level", 0)

    ctx = decision_vars(sys, Num[]; model = model, horizon = horizon, build_state_trajs = true, lb = -4.0, ub = 4.0)
    # `decision_vars(...; build_state_trajs=true)` allocates one JuMP
    # trajectory per MTK state and returns those trajectories in `ctx.x_vars`.
    x_key = find_state_key(ctx.x_vars, "x(t)")
    z_key = find_state_key(ctx.x_vars, "z(t)")
    u_key = find_parameter_key(sys, "u")
    x_traj = ctx.x_vars[x_key]
    z_traj = ctx.x_vars[z_key]

    for var in x_traj
        set_start_value(var, 0.0)
    end
    for var in z_traj
        set_start_value(var, 0.0)
    end

    fix(x_traj[1], 0.0; force = true)
    fix(z_traj[1], 0.0; force = true)

    u_traj = @variable(model, [1:horizon], lower_bound = -2.0, upper_bound = 2.0, start = 0.0, base_name = "u_stage")
    for var in u_traj
        set_start_value(var, 0.5)
    end

    # Register the DAE over the time grid.
    # This is the key call in the example.
    # It writes the time-stepping equations and algebraic residual equations
    # into the JuMP model.
    register_daesystem(
        model,
        sys,
        (0.0, (horizon - 1) * dt),
        dt,
        integrator;
        p_disc = [u_key],
        p_disc_vars = Dict(u_key => u_traj),
        x_vars = ctx.x_vars,
    )

    @constraint(model, x_traj[end] == target)

    @objective(
        model,
        Min,
        5.0 * sum((x_traj[k] - target)^2 for k in 2:horizon) +
        1e-2 * sum(u_traj[k]^2 for k in 1:horizon)
    )

    optimize!(model)

    summary = (
        integrator = integrator,
        termination_status = string(JuMP.termination_status(model)),
        primal_status = string(JuMP.primal_status(model)),
        objective = objective_value_or_nan(model),
        variable_count = JuMP.num_variables(model),
        constraint_count = JuMP.num_constraints(model; count_variable_in_set_constraints = false),
        x_final = JuMP.value(x_traj[end]),
        z_final = JuMP.value(z_traj[end]),
        u_first = JuMP.value(u_traj[1]),
        u_last = JuMP.value(u_traj[end]),
    )

    times = collect(0.0:dt:((horizon - 1) * dt))
    traj = DataFrame(
        integrator = fill(integrator, horizon),
        stage = collect(1:horizon),
        time = times,
        x = JuMP.value.(x_traj),
        z = JuMP.value.(z_traj),
        u = JuMP.value.(u_traj),
    )

    return traj, summary
end

function main()
    sys = build_demo_dae()
    sys_s = structural_simplify(sys)
    # We save both the original and simplified equation structures so the user
    # can verify that the system is still a DAE after simplification.

    structure_rows = DataFrame([
        (; system = "original", classify_system(sys)...),
        (; system = "simplified", classify_system(sys_s)...),
    ])

    equation_rows = DataFrame(
        system = String[],
        equation_index = Int[],
        lhs_type = String[],
        equation = String[],
    )
    for (system_name, system_obj) in [("original", sys), ("simplified", sys_s)]
        for (idx, eq) in enumerate(ModelingToolkit.full_equations(system_obj))
            lhs = eq.lhs
            lhs_type = try
                Symbolics.operation(lhs) isa ModelingToolkit.Differential ? "differential" : "algebraic"
            catch
                "algebraic"
            end
            push!(equation_rows, (system_name, idx, lhs_type, string(eq)))
        end
    end

    integrator_rows = DataFrame(
        integrator = String[],
        termination_status = String[],
        primal_status = String[],
        objective = Float64[],
        variable_count = Int[],
        constraint_count = Int[],
        x_final = Float64[],
        z_final = Float64[],
        u_first = Float64[],
        u_last = Float64[],
    )
    trajectory_rows = DataFrame(
        integrator = String[],
        stage = Int[],
        time = Float64[],
        x = Float64[],
        z = Float64[],
        u = Float64[],
    )

    for integrator in ("IE", "RK4", "IRK4")
        # Run the same DAE with each supported integrator.
        traj, summary = build_integrator_case(sys, integrator)
        push!(integrator_rows, summary)
        append!(trajectory_rows, traj)
    end

    mkpath(GENERATED_DIR)
    structure_path = joinpath(GENERATED_DIR, "dae_registration_structure_summary.csv")
    equations_path = joinpath(GENERATED_DIR, "dae_registration_equations.csv")
    integrator_path = joinpath(GENERATED_DIR, "dae_registration_integrator_summary.csv")
    trajectory_path = joinpath(GENERATED_DIR, "dae_registration_trajectories.csv")

    CSV.write(structure_path, structure_rows)
    CSV.write(equations_path, equation_rows)
    CSV.write(integrator_path, integrator_rows)
    CSV.write(trajectory_path, trajectory_rows)

    println("Saved DAE demo outputs to:")
    println("  ", structure_path)
    println("  ", equations_path)
    println("  ", integrator_path)
    println("  ", trajectory_path)
end

main()
