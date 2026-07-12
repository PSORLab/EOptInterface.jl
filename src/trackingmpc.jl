# Tracking MPC for ModelingToolkit plants.
#
# The code below follows the same order used in a closed-loop experiment:
# 1. choose the manipulated inputs and measured outputs;
# 2. build one JuMP prediction problem from the ModelingToolkit model;
# 3. at each sample time, write the current plant state into that problem;
# 4. solve the MPC problem and apply the first control move.
#
# The public functions keep this workflow explicit:
# - `build_tracking_mpc(...)` creates the JuMP MPC problem once;
# - `update_tracking_targets!(...)` and `update_stage_parameter!(...)` change
#   setpoints or disturbance previews between solves;
# - `solve_tracking_mpc!(...)` performs one online MPC update.

"""
    MPCControlSpec

One manipulated input in the MPC problem.

For example, this could be an aeration rate, recycle flow, or any model
parameter that the controller is allowed to change. The fields give its bounds,
move limit, and move penalty.
"""
Base.@kwdef struct MPCControlSpec
    sym::Num
    lower::Float64 = -Inf
    upper::Float64 = Inf
    delta_max::Union{Nothing, Float64} = nothing
    move_weight::Float64 = 0.0
    first_move_weight::Float64 = 0.0
    base_name::Union{Nothing, String} = nothing
end

"""
    MPCOutputSpec

One output that the MPC tries to keep near a target.

The output can have a setpoint and optional soft lower/upper limits. Soft
limits are penalized through `slack_weight` instead of making the solve fail
immediately.
"""
Base.@kwdef struct MPCOutputSpec
    sym::Num
    setpoint::Float64 = 0.0
    track_weight::Float64 = 1.0
    terminal_weight::Float64 = 0.0
    lower_soft::Union{Nothing, Float64} = nothing
    upper_soft::Union{Nothing, Float64} = nothing
    slack_weight::Float64 = 0.0
    base_name::Union{Nothing, String} = nothing
end

"""
    TrackingMPCConfig

Basic numerical settings for one tracking MPC experiment.

`PH` is the prediction horizon, `CH` is the number of free control moves, and
`dt` is the sample time used in the prediction model.
"""
Base.@kwdef struct TrackingMPCConfig
    PH::Int = 20
    CH::Int = 10
    dt::Float64 = 1.0
    integrator::String = "IE"
    system_kind::Symbol = :ode
    state_lower = 0.0
    state_upper = 1e6
    rhs0 = 0.0
    track_start_index::Int = 2
    lower_clip = nothing
    upper_clip = nothing
    time_lower::Float64 = 0.0
    time_upper::Union{Nothing, Float64} = nothing
    enforce_alg_at_t1::Bool = true
end

"""
    TrackingMPCController

Stored MPC problem used during closed-loop simulation.

Users normally create this with `build_tracking_mpc(...)` and then pass it to
`solve_tracking_mpc!(...)` at each sampling time.
"""
mutable struct TrackingMPCController
    cfg::TrackingMPCConfig
    sys
    model::JuMP.Model
    N::Int
    Nc::Int
    x_vars::Dict{Num, Vector{JuMP.VariableRef}}
    c_ic::Dict{Num, JuMP.ConstraintRef}
    control_specs::Vector{MPCControlSpec}
    output_specs::Vector{MPCOutputSpec}
    control_vars::Dict{Num, Vector{JuMP.VariableRef}}
    stage_param_vars::Dict{Num, Vector{JuMP.VariableRef}}
    control_names::Dict{Num, String}
    output_names::Dict{Num, String}
    t_stage::JuMP.VariableRef
    first_move_anchors::Dict{Num, JuMP.ConstraintRef}
    delta_vars::Dict{Num, Vector{JuMP.VariableRef}}
    target_vars::Dict{Num, JuMP.VariableRef}
    lower_ref_vars::Dict{Num, JuMP.VariableRef}
    upper_ref_vars::Dict{Num, JuMP.VariableRef}
    lower_slack_vars::Dict{Num, Vector{JuMP.VariableRef}}
    upper_slack_vars::Dict{Num, Vector{JuMP.VariableRef}}
    track_terms::Dict{Num, Any}
    soft_terms::Dict{Num, Any}
    move_terms::Dict{Num, Any}
    first_move_terms::Dict{Num, Any}
    objective_expr::Any
    last_controls::Dict{Num, Vector{Float64}}
end

"""
    _resolved_base_name(sym, custom, suffix)

Build the JuMP base name for one MPC variable family.

The code uses the same naming rule for controls, outputs, and stage parameters.
"""
function _resolved_base_name(sym, custom::Union{Nothing, String}, suffix::AbstractString)::String
    return resolve_mpc_base_name(sym; custom=custom, suffix=suffix)
end

"""
    _metric_key(prefix, label)

Build a metric key for the `metrics` dictionary returned by
`solve_tracking_mpc!(...)`.

Examples are keys such as `:track_C1_y` or `:move_Q_air_u`.
"""
function _metric_key(prefix::AbstractString, label::AbstractString)::Symbol
    return Symbol(sanitize_mpc_name(string(prefix, "_", label)))
end

"""
    _canonical_bound_spec(sys, spec, pool)

Normalize bound dictionaries onto the standard MTK keys.

If `spec` is not a dictionary, it is returned unchanged.
This is used when `build_tracking_mpc(...)` receives state bounds or initial
values keyed by user-facing symbols.
"""
function _canonical_bound_spec(sys, spec, pool::Symbol)
    if spec isa AbstractDict
        out = Dict{Num, Any}()
        for (sym, value) in pairs(spec)
            # Convert user input names to the standard MTK keys.
            out[canonical_system_symbol(sys, sym; pool=pool)] = value
        end
        return out
    end
    return spec
end

"""
    _canonical_value_map(sys, values, pool)

Convert an input dictionary to a `Num => Float64` map with standard MTK keys.

This is used before writing current plant values or previous controls into the
MPC model.
"""
function _canonical_value_map(sys, values::AbstractDict, pool::Symbol)::Dict{Num, Float64}
    out = Dict{Num, Float64}()
    for (sym, value) in pairs(values)
        # Store online updates under standard MTK keys.
        out[canonical_system_symbol(sys, sym; pool=pool)] = float(value)
    end
    return out
end

"""
    _store_result_alias!(dest, original, canonical, value)

Store one result under both the user-facing key and the standard key.

This makes the returned `controls` and `predictions` dictionaries easier to use.
"""
function _store_result_alias!(dest::AbstractDict, original, canonical, value)
    # Keep both the user-facing key and the standard key.
    dest[original] = value
    if string(original) != string(canonical)
        dest[canonical] = value
    end
    return dest
end

"""
    _stage_values(values, N, label)

Turn one stage-parameter input into a length-`N` vector.

It accepts scalars, tuples, and vectors.
This is used by `build_tracking_mpc(...)` and `update_stage_parameter!(...)`.
"""
function _stage_values(values, N::Int, label::AbstractString)::Vector{Float64}
    if values isa AbstractVector
        vals = Float64.(collect(values))
        length(vals) == N || error("Stage parameter $(label) must have length $(N), got $(length(vals)).")
        return vals
    elseif values isa Tuple
        vals = Float64.(collect(values))
        length(vals) == N || error("Stage parameter $(label) must have length $(N), got $(length(vals)).")
        return vals
    else
        return fill(float(values), N)
    end
end

"""
    _validate_unique_syms(sys, specs, pool, label)

Check whether a list of control or output specs contains duplicate MTK symbols.

`build_tracking_mpc(...)` uses this before any JuMP variables are created.
"""
function _validate_unique_syms(sys, specs, pool::Symbol, label::AbstractString)
    seen = Set{String}()
    for spec in specs
        key = string(canonical_system_symbol(sys, spec.sym; pool=pool))
        key in seen && error("Duplicate $(label) symbol detected for $(spec.sym).")
        push!(seen, key)
    end
    return nothing
end

"""
    _value_or_float(x)

Return a plain `Float64`.

If `x` is already numeric, it is converted directly.
If `x` is a JuMP expression, the helper reads its current value.
"""
function _value_or_float(x)::Float64
    return x isa Real ? float(x) : float(JuMP.value(x))
end

"""
    _register_tracking_system!(cfg, model, sys, p_disc, p_disc_vars, x_vars, t_stage)

Add the ODE or DAE prediction equations to the JuMP model.

It calls:
- `register_odesystem(...)` when `system_kind = :ode`
- `register_daesystem(...)` when `system_kind = :dae`

`build_tracking_mpc(...)` calls this after the state, control, and stage
trajectories are in place.
"""
function _register_tracking_system!(ctrl::TrackingMPCConfig,
                                    model::JuMP.Model,
                                    sys,
                                    p_disc::Vector{Num},
                                    p_disc_vars::Dict{Num, Vector{JuMP.VariableRef}},
                                    x_vars::AbstractDict,
                                    t_stage::JuMP.VariableRef)
    iv_mtk = ModelingToolkit.get_iv(sys)
    iv_map = Dict(iv_mtk => t_stage)
    kind = lowercase(string(ctrl.system_kind))
    if kind == "ode"
        # ODE path.
        register_odesystem(
            model,
            sys,
            (0.0, ctrl.PH * ctrl.dt),
            ctrl.dt,
            ctrl.integrator;
            p_disc = p_disc,
            p_disc_vars = p_disc_vars,
            x_vars = x_vars,
            t_map = iv_map,
        )
    elseif kind == "dae"
        # DAE path.
        register_daesystem(
            model,
            sys,
            (0.0, ctrl.PH * ctrl.dt),
            ctrl.dt,
            ctrl.integrator;
            p_disc = p_disc,
            p_disc_vars = p_disc_vars,
            x_vars = x_vars,
            t_map = iv_map,
            enforce_alg_at_t1 = ctrl.enforce_alg_at_t1,
        )
    else
        error("Unsupported system_kind=$(ctrl.system_kind). Use :ode or :dae.")
    end
    return nothing
end

"""
    build_tracking_mpc(model, sys; control_specs, output_specs, config=..., stage_param_defaults=...)

Build one tracking MPC problem around a ModelingToolkit ODE or DAE system.

This function builds the JuMP model once.
It adds state trajectories, control trajectories, optional preview parameters,
tracking terms, and move penalties.

This is the main setup function for the tracking MPC examples.

What this function creates:
- one JuMP trajectory for each model state;
- one JuMP trajectory for each control input;
- optional stage-wise preview trajectories, such as disturbances or prices;
- one objective made of tracking penalties, soft-bound penalties, and move
  penalties;
- one `TrackingMPCController` value that stores the JuMP references needed
  during simulation.

Why this matters during simulation:
- the expensive model construction happens once here;
- later online MPC steps only update values and re-solve;
- this avoids rebuilding the full JuMP model at every sample.

Algorithm:
1. Check horizon settings, control specifications, and output specifications.
2. Normalize all control, output, and stage-parameter names to one MTK key.
3. Build or reuse the state trajectories and IC constraints.
4. Create one JuMP trajectory for each control and each stage parameter.
5. Call `_register_tracking_system!(...)` to add the ODE or DAE prediction
   equations.
6. Build the tracking terms, soft-bound terms, and move penalties.
7. Store the resulting JuMP references and naming maps in a
   `TrackingMPCController` so later online solves only update values and
   re-optimize.
"""
function build_tracking_mpc(model::JuMP.Model,
                            sys;
                            control_specs::AbstractVector{MPCControlSpec},
                            output_specs::AbstractVector{MPCOutputSpec}=MPCOutputSpec[],
                            config::TrackingMPCConfig=TrackingMPCConfig(),
                            stage_param_defaults::AbstractDict=Dict{Num, Any}(),
                            x_vars::Union{Nothing, AbstractDict}=nothing,
                            c_ic::Union{Nothing, AbstractDict}=nothing,
                            store_ext::Bool=true)
    config.PH >= 1 || error("TrackingMPCConfig.PH must be at least 1.")
    config.CH >= 1 || error("TrackingMPCConfig.CH must be at least 1.")

    N = config.PH + 1
    Nc = config.CH
    Nc <= N || error("TrackingMPCConfig.CH=$(config.CH) exceeds the horizon size N=$(N).")
    1 <= config.track_start_index <= N || error("track_start_index must be in 1:N.")

    control_vec = collect(control_specs)
    output_vec = collect(output_specs)
    isempty(control_vec) && error("At least one control spec is required.")

    _validate_unique_syms(sys, control_vec, :parameters, "control")
    _validate_unique_syms(sys, output_vec, :unknowns, "output")

    canonical_controls = Dict{Num, Num}()
    for spec in control_vec
        # Save one standard MTK key for each control.
        canonical_controls[_as_num(spec.sym)] = canonical_system_parameter(sys, spec.sym)
    end

    canonical_outputs = Dict{Num, Num}()
    for spec in output_vec
        canonical_outputs[_as_num(spec.sym)] = canonical_system_unknown(sys, spec.sym)
    end

    stage_specs = NamedTuple[]
    control_keys = Set(string(sym) for sym in values(canonical_controls))
    seen_stage_keys = Set{String}()
    for (sym, values) in pairs(stage_param_defaults)
        original = _as_num(sym)
        canonical = canonical_system_parameter(sys, sym)
        key = string(canonical)
        key in control_keys && error("Stage parameter $(sym) is already defined as a control variable.")
        key in seen_stage_keys && error("Duplicate stage parameter symbol detected for $(sym).")
        push!(seen_stage_keys, key)
        push!(stage_specs, (original = original, canonical = canonical, values = values))
    end

    if (x_vars === nothing) != (c_ic === nothing)
        error("Pass both x_vars and c_ic together, or pass neither.")
    end

    p_disc_syms = Num[canonical_controls[_as_num(spec.sym)] for spec in control_vec]
    append!(p_disc_syms, Num[entry.canonical for entry in stage_specs])

    state_trajs = if x_vars === nothing
        # Default path: build one trajectory and one IC constraint per state.
        decision_ctx = decision_vars(
            sys,
            p_disc_syms;
            model = model,
            horizon = N,
            build_state_trajs = true,
            lb = _canonical_bound_spec(sys, config.state_lower, :unknowns),
            ub = _canonical_bound_spec(sys, config.state_upper, :unknowns),
            rhs0 = _canonical_bound_spec(sys, config.rhs0, :unknowns),
            store_ext = store_ext,
        )
        decision_ctx.x_vars
    else
        Dict{Num, Vector{JuMP.VariableRef}}(
            canonical_system_unknown(sys, var) => traj for (var, traj) in pairs(x_vars)
        )
    end

    ic_map = if c_ic === nothing
        model.ext[:c_ic_by_sym]
    else
        Dict{Num, JuMP.ConstraintRef}(
            canonical_system_unknown(sys, var) => cref for (var, cref) in pairs(c_ic)
        )
    end

    control_vars = Dict{Num, Vector{JuMP.VariableRef}}()
    control_names = Dict{Num, String}()
    for spec in control_vec
        canonical_sym = canonical_controls[_as_num(spec.sym)]
        base = _resolved_base_name(spec.sym, spec.base_name, "u")
        # Build one JuMP control vector for each control variable.
        vec = @variable(model, [1:N], lower_bound=spec.lower, upper_bound=spec.upper, base_name=base)
        control_vars[canonical_sym] = vec
        _store_result_alias!(control_names, spec.sym, canonical_sym, base)
    end

    stage_param_vars = Dict{Num, Vector{JuMP.VariableRef}}()
    for entry in stage_specs
        base = _resolved_base_name(entry.original, nothing, "stage")
        vec = @variable(model, [1:N], base_name=base)
        vals = _stage_values(entry.values, N, string(entry.original))
        for k in 1:N
            # Stage parameters are fixed values, not optimization variables.
            JuMP.fix(vec[k], vals[k]; force=true)
        end
        stage_param_vars[entry.canonical] = vec
    end

    p_disc = copy(p_disc_syms)
    p_disc_vars = Dict{Num, Vector{JuMP.VariableRef}}()
    for spec in control_vec
        canonical_sym = canonical_controls[_as_num(spec.sym)]
        p_disc_vars[canonical_sym] = control_vars[canonical_sym]
    end
    for entry in stage_specs
        p_disc_vars[entry.canonical] = stage_param_vars[entry.canonical]
    end

    t_upper = isnothing(config.time_upper) ? config.PH * config.dt : config.time_upper
    t_stage = @variable(model, lower_bound=config.time_lower, upper_bound=t_upper, base_name="t_stage")
    JuMP.fix(t_stage, config.time_lower; force=true)

    # Add the plant model after all trajectories have been created.
    _register_tracking_system!(config, model, sys, p_disc, p_disc_vars, state_trajs, t_stage)

    first_move_anchors = Dict{Num, JuMP.ConstraintRef}()
    delta_vars = Dict{Num, Vector{JuMP.VariableRef}}()
    move_terms = Dict{Num, Any}()
    first_move_terms = Dict{Num, Any}()

    for spec in control_vec
        canonical_sym = canonical_controls[_as_num(spec.sym)]
        u = control_vars[canonical_sym]
        if Nc < N
            # After the control horizon, hold the last move constant.
            for k in (Nc + 1):N
                @constraint(model, u[k] == u[Nc])
            end
        end

        d1 = @variable(model, base_name="$(control_names[spec.sym])_d1")
        # Tie the first move to an anchor that is updated online.
        first_move_anchors[canonical_sym] = @constraint(model, u[1] - d1 == 0.0)
        first_move_terms[canonical_sym] = @expression(model, spec.first_move_weight * d1^2)

        if N >= 2
            du = @variable(model, [1:(N - 1)], base_name="$(control_names[spec.sym])_du")
            # Store `du` explicitly so move changes can be inspected or penalized.
            @constraint(model, [k = 2:N], du[k - 1] == u[k] - u[k - 1])
            if spec.delta_max !== nothing
                step = float(spec.delta_max)
                @constraint(model, [k = 1:(N - 1)], -step <= du[k] <= step)
            end
            delta_vars[canonical_sym] = du
            move_terms[canonical_sym] = @expression(model, spec.move_weight * sum(du[k]^2 for k in 1:(N - 1)))
        else
            delta_vars[canonical_sym] = JuMP.VariableRef[]
            move_terms[canonical_sym] = 0.0
        end
    end

    output_names = Dict{Num, String}()
    target_vars = Dict{Num, JuMP.VariableRef}()
    lower_ref_vars = Dict{Num, JuMP.VariableRef}()
    upper_ref_vars = Dict{Num, JuMP.VariableRef}()
    lower_slack_vars = Dict{Num, Vector{JuMP.VariableRef}}()
    upper_slack_vars = Dict{Num, Vector{JuMP.VariableRef}}()
    track_terms = Dict{Num, Any}()
    soft_terms = Dict{Num, Any}()
    objective_terms = Any[]

    for spec in output_vec
        canonical_sym = canonical_outputs[_as_num(spec.sym)]
        haskey(state_trajs, canonical_sym) || error("Output symbol $(spec.sym) is not present in x_vars.")
        base = _resolved_base_name(spec.sym, spec.base_name, "y")
        _store_result_alias!(output_names, spec.sym, canonical_sym, base)

        sp = @variable(model, base_name="$(base)_sp")
        JuMP.fix(sp, spec.setpoint; force=true)
        target_vars[canonical_sym] = sp

        y = state_trajs[canonical_sym]
        # The tracking term is written on the state trajectory.
        track_terms[canonical_sym] = @expression(
            model,
            spec.track_weight * sum((y[k] - sp)^2 for k in config.track_start_index:N) +
            spec.terminal_weight * (y[N] - sp)^2
        )
        push!(objective_terms, track_terms[canonical_sym])

        if (spec.lower_soft !== nothing || spec.upper_soft !== nothing) && spec.slack_weight <= 0.0
            error("Output $(spec.sym) defines soft bounds but has slack_weight <= 0.")
        end

        soft_parts = Any[]
        if spec.lower_soft !== nothing
            lo = @variable(model, base_name="$(base)_lo")
            JuMP.fix(lo, float(spec.lower_soft); force=true)
            s_lo = @variable(model, [1:N], lower_bound=0.0, base_name="$(base)_s_lo")
            # Slack variables soften the bound and add a penalty.
            @constraint(model, [k = config.track_start_index:N], y[k] + s_lo[k] >= lo)
            lower_ref_vars[canonical_sym] = lo
            lower_slack_vars[canonical_sym] = s_lo
            push!(soft_parts, sum(s_lo[k]^2 for k in config.track_start_index:N))
        end

        if spec.upper_soft !== nothing
            hi = @variable(model, base_name="$(base)_hi")
            JuMP.fix(hi, float(spec.upper_soft); force=true)
            s_hi = @variable(model, [1:N], lower_bound=0.0, base_name="$(base)_s_hi")
            @constraint(model, [k = config.track_start_index:N], y[k] - s_hi[k] <= hi)
            upper_ref_vars[canonical_sym] = hi
            upper_slack_vars[canonical_sym] = s_hi
            push!(soft_parts, sum(s_hi[k]^2 for k in config.track_start_index:N))
        end

        soft_terms[canonical_sym] = isempty(soft_parts) ? 0.0 : @expression(model, spec.slack_weight * reduce(+, soft_parts))
        push!(objective_terms, soft_terms[canonical_sym])
    end

    append!(objective_terms, values(move_terms))
    append!(objective_terms, values(first_move_terms))
    objective_expr = isempty(objective_terms) ? 0.0 : reduce(+, objective_terms)
    # Build the full objective from all pieces.
    @objective(model, Min, objective_expr)

    ctrl = TrackingMPCController(
        config,
        sys,
        model,
        N,
        Nc,
        state_trajs,
        ic_map,
        control_vec,
        output_vec,
        control_vars,
        stage_param_vars,
        control_names,
        output_names,
        t_stage,
        first_move_anchors,
        delta_vars,
        target_vars,
        lower_ref_vars,
        upper_ref_vars,
        lower_slack_vars,
        upper_slack_vars,
        track_terms,
        soft_terms,
        move_terms,
        first_move_terms,
        objective_expr,
        Dict{Num, Vector{Float64}}(),
    )

    if store_ext
        model.ext[:tracking_mpc] = ctrl
    end

    return ctrl
end

"""
    _trajectory_base_name(traj)

Return the base JuMP name for one trajectory vector.

If the first variable is named `x[1]`, this returns `x`.
This is used by the trajectory query helpers.
"""
function _trajectory_base_name(traj::AbstractVector{<:JuMP.VariableRef})::Union{Nothing, String}
    isempty(traj) && return nothing
    name = JuMP.name(first(traj))
    isempty(name) && return nothing
    return replace(name, r"\[\d+\]$" => "")
end

"""
    _tracking_aliases(ctrl, sym, traj, kind)

Build the list of names that can refer to one trajectory.

It includes MTK names, display names, sanitized names, and the stored public
base names.
"""
function _tracking_aliases(ctrl::TrackingMPCController,
                           sym::Num,
                           traj::AbstractVector{<:JuMP.VariableRef},
                           kind::Symbol)::Vector{String}
    aliases = String[string(sym), display_mpc_name(sym), sanitize_mpc_name(sym)]
    base = _trajectory_base_name(traj)
    if base !== nothing
        push!(aliases, base)
    end
    if kind === :control && haskey(ctrl.control_names, sym)
        push!(aliases, ctrl.control_names[sym])
    elseif kind === :stage && haskey(ctrl.stage_param_vars, sym)
        push!(aliases, resolve_mpc_base_name(sym; suffix="stage"))
    end
    return unique(filter(!isempty, aliases))
end

"""
    _resolve_tracking_query(ctrl, query, kind)

Convert a state, control, or stage-parameter query into one standard MTK symbol.

It turns a user query into one standard MTK symbol.
The query can be:
- an MTK symbol,
- a display name,
- or a JuMP base name.
"""
function _resolve_tracking_query(ctrl::TrackingMPCController, query, kind::Symbol)::Num
    dict = if kind === :state
        ctrl.x_vars
    elseif kind === :control
        ctrl.control_vars
    elseif kind === :stage
        ctrl.stage_param_vars
    else
        error("Unsupported tracking lookup kind $(kind). Use :state, :control, or :stage.")
    end

    if !(query isa AbstractString)
        # Non-string queries are treated as MTK symbols.
        return kind === :state ?
            canonical_system_unknown(ctrl.sys, query) :
            canonical_system_parameter(ctrl.sys, query)
    end

    matches = Num[]
    for (sym, traj) in pairs(dict)
        query in _tracking_aliases(ctrl, sym, traj, kind) || continue
        push!(matches, sym)
    end

    isempty(matches) && error("No $(kind) trajectory matches \"$(query)\".")
    length(matches) == 1 && return only(matches)

    labels = join(display_mpc_name.(matches), ", ")
    error("Ambiguous $(kind) trajectory query \"$(query)\". Matches: $(labels)")
end

"""
    state_traj(ctrl, query)

Return the JuMP state trajectory that matches `query`.

You can query by MTK symbol, display name, or JuMP base name.
"""
function state_traj(ctrl::TrackingMPCController, query)
    sym = _resolve_tracking_query(ctrl, query, :state)
    return ctrl.x_vars[sym]
end

"""
    control_traj(ctrl, query)

Return the JuMP control trajectory that matches `query`.

You can query by MTK symbol, display name, or JuMP base name.
"""
function control_traj(ctrl::TrackingMPCController, query)
    sym = _resolve_tracking_query(ctrl, query, :control)
    return ctrl.control_vars[sym]
end

"""
    stage_param_traj(ctrl, query)

Return the fixed stage-parameter trajectory that matches `query`.
"""
function stage_param_traj(ctrl::TrackingMPCController, query)
    sym = _resolve_tracking_query(ctrl, query, :stage)
    return ctrl.stage_param_vars[sym]
end

"""
    _tracking_selected_syms(ctrl, queries, kind)

Choose which output or control symbols appear in `print_tracking_status(...)`.

It decides which output or control symbols should appear in the printed status
line.
"""
function _tracking_selected_syms(ctrl::TrackingMPCController, queries, kind::Symbol)::Vector{Num}
    if queries === nothing
        if kind === :state
            return [canonical_system_unknown(ctrl.sys, spec.sym) for spec in ctrl.output_specs]
        elseif kind === :control
            return [canonical_system_parameter(ctrl.sys, spec.sym) for spec in ctrl.control_specs]
        else
            error("Unsupported tracking selection kind $(kind).")
        end
    end

    out = Num[]
    seen = Set{String}()
    lookup_kind = kind === :state ? :state : :control
    for query in queries
        sym = _resolve_tracking_query(ctrl, query, lookup_kind)
        key = string(sym)
        key in seen && continue
        push!(seen, key)
        push!(out, sym)
    end
    return out
end

"""
    _format_tracking_scalar(value; digits=3)

Format one number for status printing.

It prints `NaN` cleanly and rounds ordinary numbers to a small number of digits.
"""
function _format_tracking_scalar(value::Real; digits::Int=3)::String
    return isnan(value) ? "NaN" : string(round(float(value); digits=digits))
end

"""
    _tracking_state_snapshot(ctrl, state_values)

Normalize one `state => value` dictionary onto the standard MTK state keys.

`print_tracking_status(...)` uses this so it can show the current plant state
with the same symbols used inside the controller.
"""
function _tracking_state_snapshot(ctrl::TrackingMPCController,
                                  state_values::Union{Nothing, AbstractDict})
    state_values === nothing && return Dict{Num, Float64}()
    return _canonical_value_map(ctrl.sys, state_values, :unknowns)
end

"""
    print_tracking_status(io, ctrl, result; ...)

Print one short line that summarizes the current MPC solve.

The line can include time, solve status, objective value, tracked outputs, and
the first control move.

This function exists for users who want live feedback while the simulation is
running.
It is not part of the optimization itself.
It is only a reporting helper.

What it prints:
- the current time, if the caller supplies it;
- the JuMP termination status;
- whether that status is treated as acceptable for MPC;
- the total objective value;
- selected tracked outputs together with their setpoints;
- selected control values, usually the first move to be applied.

How it connects to other functions:
- `solve_tracking_mpc!(...)` can call this automatically with `show_status=true`;
- `state_traj(...)`, `control_traj(...)`, and the naming helpers make it easier
  to choose which variables should appear in the printed line.
"""
function print_tracking_status(io::IO,
                               ctrl::TrackingMPCController,
                               result;
                               time=nothing,
                               state_values::Union{Nothing, AbstractDict}=nothing,
                               output_syms=nothing,
                               control_syms=nothing,
                               digits::Int=3,
                               prefix::AbstractString="[MPC]",
                               accepted_statuses=default_mpc_accepted_statuses())
    summary = String[prefix]
    if time !== nothing
        push!(summary, "t=" * _format_tracking_scalar(time; digits=digits))
    end

    push!(summary, "status=" * string(result.status))
    push!(summary, "accepted=" * string(is_accepted_mpc_status(result.status; accepted_statuses=accepted_statuses)))
    objective = get(result.metrics, :objective, NaN)
    push!(summary, "obj=" * _format_tracking_scalar(objective; digits=digits))

    snapshot = _tracking_state_snapshot(ctrl, state_values)
    selected_outputs = _tracking_selected_syms(ctrl, output_syms, :state)
    if !isempty(selected_outputs)
        parts = String[]
        for sym in selected_outputs
            label = display_mpc_name(sym)
            current = if haskey(snapshot, sym)
                snapshot[sym]
            elseif haskey(result.predictions, sym) && !isempty(result.predictions[sym])
                result.predictions[sym][1]
            else
                NaN
            end
            target = haskey(ctrl.target_vars, sym) ? JuMP.fix_value(ctrl.target_vars[sym]) : NaN
            push!(parts, string(label, "=", _format_tracking_scalar(current; digits=digits),
                                " (sp=", _format_tracking_scalar(target; digits=digits), ")"))
        end
        push!(summary, "y: " * join(parts, ", "))
    end

    selected_controls = _tracking_selected_syms(ctrl, control_syms, :control)
    if !isempty(selected_controls)
        parts = String[]
        for sym in selected_controls
            label = haskey(ctrl.control_names, sym) ? ctrl.control_names[sym] : display_mpc_name(sym)
            seq = get(result.controls, sym, Float64[])
            applied = isempty(seq) ? NaN : seq[1]
            push!(parts, string(label, "=", _format_tracking_scalar(applied; digits=digits)))
        end
        push!(summary, "u: " * join(parts, ", "))
    end

    head, tail... = summary
    if isempty(tail)
        println(io, head)
    else
        println(io, string(head, " ", join(tail, " | ")))
    end
    return nothing
end

print_tracking_status(ctrl::TrackingMPCController, result; kwargs...) =
    print_tracking_status(stdout, ctrl, result; kwargs...)

"""
    update_tracking_targets!(ctrl; setpoints=..., lower_soft=..., upper_soft=...)

Update the fixed reference values in an existing tracking MPC controller.

This function changes the optimization targets without rebuilding the whole
controller.
That is important because rebuilding the JuMP model is usually much more
expensive than simply updating fixed reference values.

Typical use cases:
- change the setpoint during a setpoint-tracking study;
- change a lower soft limit or upper soft limit during a scenario switch;
- keep the same controller structure but study different operating goals.

How it connects to other functions:
- `build_tracking_mpc(...)` creates the fixed JuMP variables that store
  setpoints and soft limits;
- this function updates those fixed values in-place;
- the next call to `solve_tracking_mpc!(...)` then uses the new targets.
"""
function update_tracking_targets!(ctrl::TrackingMPCController;
                                  setpoints::AbstractDict=Dict{Any, Any}(),
                                  lower_soft::AbstractDict=Dict{Any, Any}(),
                                  upper_soft::AbstractDict=Dict{Any, Any}())
    for (sym, value) in pairs(setpoints)
        canonical_sym = canonical_system_unknown(ctrl.sys, sym)
        haskey(ctrl.target_vars, canonical_sym) || error("No setpoint variable is registered for $(sym).")
        # The reference is a fixed scalar JuMP variable.
        JuMP.fix(ctrl.target_vars[canonical_sym], float(value); force=true)
    end

    for (sym, value) in pairs(lower_soft)
        canonical_sym = canonical_system_unknown(ctrl.sys, sym)
        haskey(ctrl.lower_ref_vars, canonical_sym) || error("No lower soft bound is registered for $(sym).")
        JuMP.fix(ctrl.lower_ref_vars[canonical_sym], float(value); force=true)
    end

    for (sym, value) in pairs(upper_soft)
        canonical_sym = canonical_system_unknown(ctrl.sys, sym)
        haskey(ctrl.upper_ref_vars, canonical_sym) || error("No upper soft bound is registered for $(sym).")
        JuMP.fix(ctrl.upper_ref_vars[canonical_sym], float(value); force=true)
    end

    return ctrl
end

"""
    update_stage_parameter!(ctrl, sym, values)

Update one fixed stage parameter trajectory, such as a disturbance preview.

This function is how the controller receives a preview that changes from stage
to stage across the prediction horizon.
Common examples are disturbance previews, feed forecasts, or time-varying
prices.

Why this helper exists:
- the JuMP variables for stage parameters are created once in
  `build_tracking_mpc(...)`;
- each new MPC solve may need new preview values;
- this function updates those values without rebuilding the controller.

How it connects to other functions:
- `build_tracking_mpc(...)` creates the fixed stage-parameter trajectories;
- this function overwrites those fixed values;
- `solve_tracking_mpc!(...)` then solves with the new preview.
"""
function update_stage_parameter!(ctrl::TrackingMPCController, sym, values)
    canonical_sym = canonical_system_parameter(ctrl.sys, sym)
    haskey(ctrl.stage_param_vars, canonical_sym) || error("No stage parameter trajectory is registered for $(sym).")
    vals = _stage_values(values, ctrl.N, string(sym))
    for k in 1:ctrl.N
        # Preview values are fixed JuMP variables.
        JuMP.fix(ctrl.stage_param_vars[canonical_sym][k], vals[k]; force=true)
    end
    return ctrl
end

"""
    current_state_map(integ, i_state)

Read the current DifferentialEquations state into a `var => value` dictionary.

This method uses a precomputed `state => index` map.

Use this when you already know the mapping from each MTK state to the matching
position in `integ.u`.
This is slightly more efficient than rebuilding the mapping every time.

How it connects to other functions:
- `make_state_index_map(...)` builds the `state => index` map;
- `prepare_tracking_mpc_step!(...)` and other online helpers often need the
  dictionary returned by this function.
"""
function current_state_map(integ, i_state::AbstractDict)
    values = Dict{Any, Float64}()
    for (var, idx) in pairs(i_state)
        values[var] = float(integ.u[idx])
    end
    return values
end

"""
    current_state_map(integ, sys)

Read the current DifferentialEquations state into a `var => value` dictionary.

This method builds the lookup from `unknowns(sys)`.

This is the more convenient version of `current_state_map(...)`.
It is easier to use in quick scripts because you only need the integrator and
the MTK system.
The tradeoff is that it rebuilds the lookup map each time you call it.
"""
function current_state_map(integ, sys)
    vars = collect(ModelingToolkit.unknowns(sys))
    i_state = make_state_index_map(vars)
    return current_state_map(integ, i_state)
end

"""
    seed_mpc_log!(logctx, state_values, t0; control_values=...)

Store one initial log entry before the closed-loop run starts.

Use this right before the first simulated control step.
It records the initial plant condition so plots and later analysis start from
the true initial point instead of from the first callback time.

How it connects to other functions:
- `make_mpc_log(...)` creates the empty log;
- `seed_mpc_log!(...)` writes the first row;
- `log_mpc_state!(...)` then appends all later rows.
"""
function seed_mpc_log!(logctx::MPCLog,
                       state_values::AbstractDict,
                       t0::Real;
                       control_values::AbstractDict=Dict{Any, Any}())
    # Add an initial point before the first callback fires.
    push!(logctx.ts, float(t0))
    for (var, hist) in logctx.Xhist
        push!(hist, float(state_values[var]))
    end
    for (key, hist) in logctx.Controlhist
        push!(hist, float(control_values[key]))
    end
    return nothing
end

"""
    prepare_tracking_mpc_step!(ctrl, state_values, previous_controls; ...)

Update the ICs, warm starts, and first-move bounds before solving.

This is the last setup step right before one online MPC solve.
It does not solve the optimization problem by itself.
Its only job is to make sure the already-built JuMP model matches the newest
plant information.

Inputs:
- `state_values`: the current plant state estimate or measurement;
- `previous_controls`: the control values that are currently active at the
  plant.

Why both are needed:
- the plant state tells the optimizer where the prediction starts;
- the previous control tells the optimizer what "next move" is allowed under
  move-rate limits.

Algorithm:
1. Normalize the incoming plant state and previous control names.
2. Clean very small values so the MPC starts from a numerically tidy point.
3. Write the current plant state into the first point of each state trajectory.
4. Fill the rest of each state trajectory with a warm start.
5. Update the first-move anchor to the last applied control value.
6. Shift the previous optimal control sequence when possible.
7. Reapply control bounds and first-move delta limits before the next solve.
"""
function prepare_tracking_mpc_step!(ctrl::TrackingMPCController,
                                    state_values::AbstractDict,
                                    previous_controls::AbstractDict;
                                    lower_clip=ctrl.cfg.lower_clip,
                                    upper_clip=ctrl.cfg.upper_clip,
                                    atol::Real=1e-12,
                                    missing::Symbol=:error)
    # Normalize the incoming plant values before writing them into JuMP.
    prepared = _canonical_value_map(ctrl.sys, state_values, :unknowns)
    previous = _canonical_value_map(ctrl.sys, previous_controls, :parameters)
    zero_small_values!(prepared; atol=atol)

    sync_state_trajectories!(
        ctrl.c_ic,
        ctrl.x_vars,
        prepared;
        lower_clip = lower_clip,
        upper_clip = upper_clip,
        set_start = true,
        fill_horizon = true,
        missing = missing,
    )

    for spec in ctrl.control_specs
        canonical_sym = canonical_system_parameter(ctrl.sys, spec.sym)
        haskey(previous, canonical_sym) || error("Missing previous control value for $(spec.sym).")
        u_prev = previous[canonical_sym]
        # Update the anchor tied to the last applied control.
        JuMP.set_normalized_rhs(ctrl.first_move_anchors[canonical_sym], u_prev)

        if haskey(ctrl.last_controls, canonical_sym)
            # Best warm start: shift the last optimal sequence by one step.
            shift_warm_start!(ctrl.control_vars[canonical_sym], ctrl.last_controls[canonical_sym])
        else
            # First solve: use a constant warm start.
            warm_start_with_constant!(ctrl.control_vars[canonical_sym], u_prev)
        end

        set_first_move_bounds!(
            ctrl.control_vars[canonical_sym],
            u_prev;
            lower = spec.lower,
            upper = spec.upper,
            delta_max = spec.delta_max,
        )
    end

    return ctrl
end

_elapsed_wall_s(start_ns::Integer)::Float64 = (time_ns() - start_ns) / 1.0e9

"""
    solve_tracking_mpc!(ctrl, state_values, previous_controls; ...)

Prepare and solve one MPC step.

This is the main online call used inside a closed-loop simulation or callback.
For most users, this is the function that "runs MPC once."

What it returns:
- the optimized control trajectories;
- the predicted state or output trajectories;
- a small dictionary of scalar metrics;
- the solver termination status.

Why it is separated from `build_tracking_mpc(...)`:
- `build_tracking_mpc(...)` is the expensive one-time construction step;
- this function is the repeated online step that uses the already-built model.

Algorithm:
1. Call `prepare_tracking_mpc_step!(...)` to refresh ICs, warm starts, and move
   bounds.
2. Solve the JuMP model.
3. Check the termination status against the accepted MPC statuses.
4. Read back the optimized control trajectories.
5. Read back the predicted state or output trajectories.
6. Compute metrics such as the total objective and its main pieces.
7. Optionally print a compact live status line for the user.
8. Save the last control sequence so the next MPC step can shift it as a warm
   start.

If `return_timing=true`, the returned named tuple also includes a `timing`
field with a wall-clock breakdown of the online MPC step.
"""
function solve_tracking_mpc!(ctrl::TrackingMPCController,
                             state_values::AbstractDict,
                             previous_controls::AbstractDict;
                             lower_clip=ctrl.cfg.lower_clip,
                             upper_clip=ctrl.cfg.upper_clip,
                             atol::Real=1e-12,
                             include_all_states::Bool=false,
                             show_status::Bool=false,
                             status_io::IO=stdout,
                             status_time=nothing,
                             status_output_syms=nothing,
                             status_control_syms=nothing,
                             status_digits::Int=3,
                             status_prefix::AbstractString="[MPC]",
                             return_timing::Bool=false,
                             accepted_statuses=(MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED))
    total_start_ns = time_ns()
    # Refresh ICs, warm starts, and first-move bounds before each solve.
    prepare_start_ns = time_ns()
    prepare_tracking_mpc_step!(
        ctrl,
        state_values,
        previous_controls;
        lower_clip = lower_clip,
        upper_clip = upper_clip,
        atol = atol,
    )
    prepare_tracking_mpc_step_s = _elapsed_wall_s(prepare_start_ns)

    optimize_start_ns = time_ns()
    JuMP.optimize!(ctrl.model)
    jump_optimize_s = _elapsed_wall_s(optimize_start_ns)

    postprocess_start_ns = time_ns()
    status = JuMP.termination_status(ctrl.model)
    status in accepted_statuses || error("Tracking MPC solve failed with status $(status).")

    controls = Dict{Any, Vector{Float64}}()
    for spec in ctrl.control_specs
        canonical_sym = canonical_system_parameter(ctrl.sys, spec.sym)
        # Store both the user key and the standard key.
        seq = collect(float.(JuMP.value.(ctrl.control_vars[canonical_sym])))
        _store_result_alias!(controls, spec.sym, canonical_sym, seq)
        ctrl.last_controls[canonical_sym] = copy(seq)
    end

    predictions = Dict{Any, Vector{Float64}}()
    for spec in ctrl.output_specs
        canonical_sym = canonical_system_unknown(ctrl.sys, spec.sym)
        _store_result_alias!(predictions, spec.sym, canonical_sym, collect(float.(JuMP.value.(ctrl.x_vars[canonical_sym]))))
    end
    if include_all_states
        for (sym, traj) in ctrl.x_vars
            predictions[sym] = collect(float.(JuMP.value.(traj)))
        end
    end
    for spec in ctrl.control_specs
        canonical_sym = canonical_system_parameter(ctrl.sys, spec.sym)
        _store_result_alias!(predictions, spec.sym, canonical_sym, controls[spec.sym])
    end

    metrics = Dict{Symbol, Float64}(:objective => float(JuMP.objective_value(ctrl.model)))
    for spec in ctrl.output_specs
        canonical_sym = canonical_system_unknown(ctrl.sys, spec.sym)
        label = ctrl.output_names[spec.sym]
        # Save objective pieces separately for later inspection.
        metrics[_metric_key("track", label)] = _value_or_float(ctrl.track_terms[canonical_sym])
        metrics[_metric_key("soft", label)] = _value_or_float(ctrl.soft_terms[canonical_sym])
    end
    for spec in ctrl.control_specs
        canonical_sym = canonical_system_parameter(ctrl.sys, spec.sym)
        label = ctrl.control_names[spec.sym]
        metrics[_metric_key("move", label)] = _value_or_float(ctrl.move_terms[canonical_sym])
        metrics[_metric_key("move1", label)] = _value_or_float(ctrl.first_move_terms[canonical_sym])
    end
    solve_postprocess_s = _elapsed_wall_s(postprocess_start_ns)

    status_print_s = 0.0
    if show_status
        status_print_start_ns = time_ns()
        print_tracking_status(
            status_io,
            ctrl,
            (controls = controls, predictions = predictions, metrics = metrics, status = status);
            time = status_time,
            state_values = state_values,
            output_syms = status_output_syms,
            control_syms = status_control_syms,
            digits = status_digits,
            prefix = status_prefix,
            accepted_statuses = accepted_statuses,
        )
        status_print_s = _elapsed_wall_s(status_print_start_ns)
    end

    if return_timing
        result = (
            controls = controls,
            predictions = predictions,
            metrics = metrics,
            status = status,
            timing = (
                prepare_tracking_mpc_step_s = prepare_tracking_mpc_step_s,
                jump_optimize_s = jump_optimize_s,
                solve_postprocess_s = solve_postprocess_s,
                status_print_s = status_print_s,
                total_s = _elapsed_wall_s(total_start_ns),
            ),
        )
        return result
    end

    result = (controls = controls, predictions = predictions, metrics = metrics, status = status)
    return result
end
