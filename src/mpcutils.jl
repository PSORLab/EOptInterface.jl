# Helper routines used by the MPC examples.
#
# Most users do not need to read this file first. The main scientific workflow
# is in `trackingmpc.jl` and the examples. The functions here handle practical
# tasks that come up during a closed-loop run: naming JuMP variables, copying
# measured states into initial-condition constraints, warm-starting controls,
# logging trajectories, and inspecting a JuMP model when a solve looks wrong.

"""
    dump_all_constraints(model)

Print all JuMP constraints in a model.

This is a debugging aid for small models, not part of the normal MPC workflow.
"""
function dump_all_constraints(model::JuMP.Model)
    println("==== Constraint types & counts ====")
    for (func_type, set_type) in JuMP.list_of_constraint_types(model)
        println("[$(JuMP.num_constraints(model, func_type, set_type))]  F=$(func_type), S=$(set_type)")
    end

    println("==== All constraints (one by one) ====")
    for (func_type, set_type) in JuMP.list_of_constraint_types(model)
        for cref in JuMP.all_constraints(model, func_type, set_type)
            co = JuMP.constraint_object(cref)
            println("* ", cref, " :: ", typeof(co.func), " in ", typeof(co.set))
        end
    end

    return nothing
end

"""
    _constraint_uses_any_var(func, target_idxs, target_names)

Check whether a JuMP/MOI expression contains one of the variables of interest.
"""
function _constraint_uses_any_var(func, target_idxs::Set, target_names::Set{String})
    if func isa JuMP.VariableRef
        return JuMP.index(func) in target_idxs
    elseif func isa JuMP.AffExpr
        return any(JuMP.index(var) in target_idxs for (var, _) in _affexpr_var_coeff_pairs(func))
    elseif func isa MOI.ScalarAffineFunction{Float64}
        return any(term.variable_index in target_idxs for term in func.terms)
    elseif func isa MOI.ScalarQuadraticFunction{Float64}
        return any(term.variable_index in target_idxs for term in func.affine_terms) ||
               any(term.variable_1 in target_idxs || term.variable_2 in target_idxs for term in func.quadratic_terms)
    end

    rendered = sprint(show, func)
    return any(occursin(name, rendered) for name in target_names)
end

"""
    find_constraints_with_vars(model, vars; constraint_sets=(MOI.EqualTo{Float64},), verbose=true)

Return constraints that contain one of the variables in `vars`.

This is useful when a state or control appears to be fixed by an unexpected
constraint.
"""
function find_constraints_with_vars(model::JuMP.Model,
                                    vars::AbstractVector{<:JuMP.VariableRef};
                                    constraint_sets::Tuple=(MOI.EqualTo{Float64},),
                                    verbose::Bool=true)
    target_idxs = Set(JuMP.index.(vars))
    target_names = Set(filter(!isempty, String[JuMP.name(var) for var in vars]))
    hits = JuMP.ConstraintRef[]

    for (func_type, set_type) in JuMP.list_of_constraint_types(model)
        any(set_type <: allowed for allowed in constraint_sets) || continue
        for cref in JuMP.all_constraints(model, func_type, set_type)
            co = JuMP.constraint_object(cref)
            if _constraint_uses_any_var(co.func, target_idxs, target_names)
                push!(hits, cref)
                verbose && println("HIT -> ", cref, " :: ", co)
            end
        end
    end

    return hits
end

"""
    split_mtk_state_path(var)

Split a ModelingToolkit-style state name such as `reactor1₊S_O(t)` into:
- a subsystem path vector, e.g. `[:reactor1]`
- a terminal state symbol, e.g. `:S_O`

Top-level states such as `x(t)` return `(Symbol[], :x)`.
"""
function split_mtk_state_path(var)
    text = string(var)
    indexed_match = match(r"^\((.*)\)\[(\d+)\]$", text)
    index_suffix = ""
    if indexed_match !== nothing
        text = indexed_match.captures[1]
        index_suffix = "_" * indexed_match.captures[2]
    end
    text = replace(text, r"\(t\)" => "")
    parts = split(text, '₊')
    if length(parts) < 2
        return Symbol[], Symbol(text * index_suffix)
    end
    return Symbol.(parts[1:end-1]), Symbol(parts[end] * index_suffix)
end

"""
    sanitize_mpc_name(sym)

Convert a ModelingToolkit symbol or label into a safe JuMP name.

Examples:
- `reactor1₊S_O(t)` -> `"reactor1_S_O"`
- `clarifier.outlet_stream` -> `"clarifier_outlet_stream"`
"""
function sanitize_mpc_name(sym)::String
    raw = replace(string(sym), r"\(t\)" => "")
    # Replace MTK separators and punctuation with underscores.
    raw = replace(raw, '₊' => '_')
    raw = replace(raw, '.' => '_')
    raw = replace(raw, r"[^A-Za-z0-9_]+" => "_")
    raw = strip(raw, '_')
    return isempty(raw) ? "var" : raw
end

"""
    display_mpc_name(sym)

Convert an MTK symbol into a dotted display name.

Examples:
- `clarifier₊outlet_stream₊S_N2O(t)` -> `"clarifier.outlet_stream.S_N2O"`
- `x(t)` -> `"x"`
"""
function display_mpc_name(sym)::String
    raw = replace(string(sym), r"\(t\)" => "")
    return replace(raw, '₊' => '.')
end

"""
    resolve_mpc_base_name(sym; custom=nothing, suffix=nothing)

Return the JuMP base name for an MPC variable group.

If `custom` is given, use it.
Otherwise build the name from `sym`.
If `suffix` is given, append it.
"""
function resolve_mpc_base_name(sym;
                               custom::Union{Nothing, AbstractString}=nothing,
                               suffix::Union{Nothing, AbstractString}=nothing)::String
    if custom !== nothing && !isempty(custom)
        return String(custom)
    end
    base = sanitize_mpc_name(sym)
    if suffix === nothing || isempty(suffix)
        return base
    end
    return string(base, "_", suffix)
end

"""
    state_trajectory_base_name(sys, var)

Build the JuMP base name for a state trajectory.

Nested MTK variables use the subsystem path as a prefix.
Top-level variables use the root system name when possible.
"""
function state_trajectory_base_name(sys, var)::String
    path_syms, state_sym = split_mtk_state_path(var)
    if isempty(path_syms)
        # Top-level states use the root system name when possible.
        prefix = hasfield(typeof(sys), :name) ? string(getfield(sys, :name)) : ""
    else
        prefix = join(string.(path_syms), "_")
    end

    if isempty(prefix)
        return sanitize_mpc_name(state_sym)
    end
    return sanitize_mpc_name(string(prefix, "_", state_sym))
end

_as_num(sym)::Num = sym isa Num ? sym : Num(sym)

"""
    _root_symbol_prefix(sys)

Return the MTK root-name prefix for a system.

For example, if the system name is `sys`, this function returns `"sys₊"`.
It is used when the package tries to match `sys.x` with the unprefixed symbol
`x(t)`.
"""
function _root_symbol_prefix(sys)::Union{Nothing, String}
    hasfield(typeof(sys), :name) || return nothing
    name = string(getfield(sys, :name))
    isempty(name) && return nothing
    return string(name, '₊')
end

"""
    _system_symbol_aliases(sys, sym)

Build a list of valid string aliases for one MTK symbol.

This lets the code accept both prefixed and unprefixed MTK names.
"""
function _system_symbol_aliases(sys, sym)::Vector{String}
    aliases = String[string(sym)]
    prefix = _root_symbol_prefix(sys)
    if prefix !== nothing
        text = aliases[1]
        if startswith(text, prefix)
            stripped = replace(text, prefix => ""; count=1)
            stripped in aliases || push!(aliases, stripped)
        else
            prefixed = string(prefix, text)
            prefixed in aliases || push!(aliases, prefixed)
        end
    end
    return aliases
end

"""
    _system_symbol_lookup(sys, pool)

Build the string-to-symbol lookup table used by `canonical_system_symbol(...)`.

`pool` decides whether the lookup is built from unknowns, parameters, or both.
"""
function _system_symbol_lookup(sys, pool::Symbol)::Dict{String, Num}
    candidates = if pool === :unknowns
        collect(ModelingToolkit.unknowns(sys))
    elseif pool === :parameters
        collect(ModelingToolkit.parameters(sys))
    elseif pool === :all
        vcat(collect(ModelingToolkit.unknowns(sys)), collect(ModelingToolkit.parameters(sys)))
    else
        error("Unsupported symbol pool $(pool). Use :unknowns, :parameters, or :all.")
    end

    lookup = Dict{String, Num}()
    for candidate in candidates
        canon = _as_num(candidate)
        for alias in _system_symbol_aliases(sys, candidate)
            existing = get(lookup, alias, nothing)
            if existing === nothing || string(existing) == string(canon)
                lookup[alias] = canon
            end
        end
    end
    return lookup
end

"""
    canonical_system_symbol(sys, sym; pool=:all, allow_missing=false)

Resolve a user-facing MTK symbol such as `sys.x` or `sys.u` to the standard MTK
symbol used inside the package.
"""
function canonical_system_symbol(sys, sym; pool::Symbol=:all, allow_missing::Bool=false)::Num
    lookup = _system_symbol_lookup(sys, pool)
    key = string(sym)
    if haskey(lookup, key)
        # Fast path: the symbol already matches a known alias.
        return lookup[key]
    end
    if allow_missing
        return _as_num(sym)
    end
    pool_label = pool === :all ? "unknowns/parameters" : String(Symbol(pool))
    sys_label = hasfield(typeof(sys), :name) ? string(getfield(sys, :name)) : string(typeof(sys))
    error("Symbol $(sym) is not present in $(pool_label) for system $(sys_label).")
end

"""
    canonical_system_unknown(sys, sym; allow_missing=false)

Resolve `sym` to the standard state symbol in `unknowns(sys)`.
"""
function canonical_system_unknown(sys, sym; allow_missing::Bool=false)::Num
    return canonical_system_symbol(sys, sym; pool=:unknowns, allow_missing=allow_missing)
end

"""
    canonical_system_parameter(sys, sym; allow_missing=false)

Resolve `sym` to the standard parameter symbol in `parameters(sys)`.
"""
function canonical_system_parameter(sys, sym; allow_missing::Bool=false)::Num
    return canonical_system_symbol(sys, sym; pool=:parameters, allow_missing=allow_missing)
end

"""
    canonicalize_system_symbols(sys, syms; pool=:all, allow_missing=false)

Normalize a list of system symbols.
Keep the first-seen order.
Drop duplicates after normalization.
"""
function canonicalize_system_symbols(sys,
                                     syms;
                                     pool::Symbol=:all,
                                     allow_missing::Bool=false)::Vector{Num}
    out = Num[]
    seen = Set{String}()
    for sym in syms
        canon = canonical_system_symbol(sys, sym; pool=pool, allow_missing=allow_missing)
        key = string(canon)
        key in seen && continue
        push!(seen, key)
        push!(out, canon)
    end
    return out
end

"""
    _split_path_and_state(var)

Return the subsystem path and terminal state name for one variable.

This keeps the path-handling code short in the functions below.
"""
function _split_path_and_state(var)
    return split_mtk_state_path(var)
end

"""
    _subsystem_path_syms(var)

Return only the subsystem-path part of an MTK variable name.

For `clarifier₊outlet_stream₊S_N2O(t)`, this returns the path symbols and drops
the final state name.
"""
function _subsystem_path_syms(var)::Vector{Symbol}
    path_syms, _ = _split_path_and_state(var)
    return path_syms
end

"""
    _resolve_path(sys, path_syms)

Walk through nested subsystem names and return the final object.

This is an internal helper for flowsheet-style MTK models with nested
subsystems.
"""
function _resolve_path(sys, path_syms::Vector{Symbol})
    obj = sys
    start_idx = 1
    if !isempty(path_syms) && hasfield(typeof(sys), :name) && string(path_syms[1]) == string(getfield(sys, :name))
        start_idx = 2
    end
    for name in path_syms[start_idx:end]
        obj = getproperty(obj, name)
    end
    return obj
end

"""
    collect_subsystems(sys, vars; skip_invalid=true)

Resolve the subsystem path referenced by each variable in `vars` and return the
deduplicated subsystem objects in order of first appearance.
"""
function collect_subsystems(sys, vars; skip_invalid::Bool=true)
    seen = Set{String}()
    subsystems = Any[]

    for var in vars
        path = _subsystem_path_syms(var)
        isempty(path) && continue
        key = join(string.(path), ".")
        key in seen && continue
        try
            push!(subsystems, _resolve_path(sys, path))
            push!(seen, key)
        catch err
            skip_invalid || rethrow(err)
        end
    end

    return subsystems
end

"""
    pretty_subsystems(sys, subs)

Format subsystem objects as human-readable `sys.a.b` strings when the subsystem
exposes a `name` field. Falls back to `string(obj)` otherwise.
"""
function pretty_subsystems(sys, subs)::Vector{String}
    _ = sys
    names = String[]
    for obj in subs
        name_str = if hasproperty(obj, :name)
            string(getproperty(obj, :name))
        else
            string(obj)
        end
        push!(names, "sys." * replace(name_str, '₊' => '.'))
    end
    return names
end

"""
    group_state_symbols(vars)

Group variables by their leading subsystem name. For example,
`reactor2₊S_NO3(t)` is grouped under `:reactor2` with state `:S_NO3`.
"""
function group_state_symbols(vars::Vector)
    groups = Dict{Symbol, Vector{Symbol}}()
    for var in vars
        path_syms, state_sym = _split_path_and_state(var)
        isempty(path_syms) && continue
        states = get!(groups, first(path_syms), Symbol[])
        state_sym in states || push!(states, state_sym)
    end
    return Dict(key => Tuple(values) for (key, values) in groups)
end

"""
    build_state_trajs_from_vars!(model, sys, vars, N; lb=0.0, ub=1e6, rhs0=0.0, store_ext=true)

Create one JuMP state trajectory for each variable in `vars`.
Returns `(x_vars, c_ic)` where:
- `x_vars[var]` is the trajectory vector for the MTK variable or resolved leaf object.
- `c_ic[var]` is the initial-condition equality on the first time point.

`lb`, `ub`, and `rhs0` may be scalars or dictionaries keyed by the MTK state
variables. Top-level state variables are supported in addition to nested
subsystem variables.

Algorithm:
1. Read each MTK state name and split it into a subsystem path and a leaf state
   name.
2. Group states by subsystem so naming and lookup stay consistent.
3. For each state, choose the lower bound, upper bound, and IC right-hand side.
4. Create one JuMP vector of length `N` for that state.
5. Add one equality constraint to pin the first time point to the chosen IC.
6. Store the trajectories and IC constraints in dictionaries keyed by the MTK
   state symbol.
7. Optionally save those dictionaries in `model.ext` for later MPC helpers.
"""
function _bound_for_var(spec, var, default)
    # `spec` can be a scalar, a dictionary, or `nothing`.
    # This helper turns all of those cases into one numeric value.
    if spec isa AbstractDict
        return get(spec, var, default)
    elseif spec === nothing
        return default
    else
        return spec
    end
end

function build_state_trajs_from_vars!(model::JuMP.Model,
                                      sys,
                                      vars::Vector,
                                      N::Int;
                                      lb=0.0,
                                      ub=1e6,
                                      rhs0=0.0,
                                      store_ext::Bool=true)
    grouped = Dict{String, Tuple{Any, Set{Symbol}}}()
    original_var_lookup = Dict{Tuple{String, Symbol}, Any}()
    for var in vars
        path_syms, state_sym = _split_path_and_state(var)
        state_sym == Symbol("") && continue
        # Group states by subsystem path.
        key = isempty(path_syms) ? "__root__" : join(string.(path_syms), ".")
        if !haskey(grouped, key)
            subsys = isempty(path_syms) ? sys : _resolve_path(sys, path_syms)
            grouped[key] = (subsys, Set{Symbol}())
        end
        push!(grouped[key][2], state_sym)
        original_var_lookup[(key, state_sym)] = var
    end

    x_vars = Dict{Any, Vector{JuMP.VariableRef}}()
    c_ic = Dict{Any, JuMP.ConstraintRef}()
    ic_map = Dict{JuMP.VariableRef, JuMP.ConstraintRef}()

    for (path_key, (subsys, state_syms)) in grouped
        for state_sym in sort!(collect(state_syms); by=string)
            name_source = original_var_lookup[(path_key, state_sym)]
            var_key = name_source isa Union{AbstractString, Symbol} ? getproperty(subsys, state_sym) : name_source
            lb_var = _bound_for_var(lb, var_key, 0.0)
            ub_var = _bound_for_var(ub, var_key, 1e6)
            rhs0_var = _bound_for_var(rhs0, var_key, 0.0)
            # Build one JuMP vector for this MTK state.
            traj = @variable(model, [1:N], lower_bound=lb_var, upper_bound=ub_var, base_name=state_trajectory_base_name(sys, name_source))
            x_vars[var_key] = traj
            # Store the first-point equality separately.
            ic = @constraint(model, traj[1] == rhs0_var)
            c_ic[var_key] = ic
            ic_map[traj[1]] = ic
        end
    end

    if store_ext
        model.ext[:x_vars] = x_vars
        model.ext[:c_ic_by_sym] = c_ic
        model.ext[:ic_map] = ic_map
        model.ext[:ic_built] = true
    end

    return x_vars, c_ic
end

"""
    zero_small_values!(data; atol=1e-6)

Set very small numbers to zero in-place.
Works on dictionaries and arrays.
"""
function zero_small_values!(data::AbstractDict; atol::Real=1e-6)
    for (key, value) in data
        if value isa Real && abs(float(value)) <= atol
            data[key] = 0.0
        end
    end
    return data
end

function zero_small_values!(data::AbstractArray{<:Real}; atol::Real=1e-6)
    for idx in eachindex(data)
        if abs(float(data[idx])) <= atol
            data[idx] = zero(eltype(data))
        end
    end
    return data
end

"""
    _apply_value_clips(value; lower_clip=nothing, upper_clip=nothing)

Clip one numeric value to optional lower and upper limits.

This is used when the current plant state is copied into MPC initial conditions.
"""
function _apply_value_clips(value::Real; lower_clip=nothing, upper_clip=nothing)
    clipped = float(value)
    if lower_clip !== nothing
        clipped = max(clipped, float(lower_clip))
    end
    if upper_clip !== nothing
        clipped = min(clipped, float(upper_clip))
    end
    return clipped
end

"""
    sync_state_trajectories!(c_ic, x_vars, values; set_start=true, fill_horizon=true,
                             lower_clip=nothing, upper_clip=nothing, missing=:error)

Update the initial-condition constraints from `values`.
You can also copy the same value across the whole trajectory as a start value.

This function is one of the most important online-MPC helpers.
The controller already has a JuMP model and already has state trajectories.
What changes from one solve to the next is the current plant state.
This function copies that newest plant state into the model.

Inputs:
- `c_ic`: the dictionary of initial-condition constraints, usually one per
  state trajectory.
- `x_vars`: the dictionary of JuMP state trajectories.
- `values`: the newest measured or estimated plant state values.

What this function does:
- it updates the right-hand side of each IC equality;
- it optionally writes warm-start values into the JuMP state trajectory;
- it can also fill the whole horizon with the current state value, which is a
  simple and robust default warm start.

How it connects to other functions:
- `build_state_trajs_from_vars!(...)` creates the original trajectories and IC
  constraints;
- `prepare_tracking_mpc_step!(...)` calls this before every MPC solve;
- `solve_tracking_mpc!(...)` depends on this update so the optimization problem
  starts from the current plant state, not from stale data.
"""
function sync_state_trajectories!(c_ic::AbstractDict,
                                  x_vars::Union{Nothing, AbstractDict},
                                  values::AbstractDict;
                                  set_start::Bool=true,
                                  fill_horizon::Bool=true,
                                  lower_clip=nothing,
                                  upper_clip=nothing,
                                  missing::Symbol=:error)
    missing in (:error, :skip) || error("Unsupported missing policy $(missing). Use :error or :skip.")

    for (var, cref) in c_ic
        if !haskey(values, var)
            missing === :skip && continue
            error("Missing value for state $(var).")
        end

        x0 = _apply_value_clips(values[var]; lower_clip=lower_clip, upper_clip=upper_clip)
        # Update the equality on the first point.
        JuMP.set_normalized_rhs(cref, x0)

        if set_start && x_vars !== nothing
            haskey(x_vars, var) || (missing === :skip && continue)
            haskey(x_vars, var) || error("Missing trajectory for state $(var).")
            traj = x_vars[var]
            isempty(traj) && continue
            JuMP.set_start_value(traj[1], x0)
            if fill_horizon
                # Fill the whole horizon with the current value.
                for k in 2:length(traj)
                    JuMP.set_start_value(traj[k], x0)
                end
            end
        end
    end

    return nothing
end

function sync_state_trajectories!(c_ic::AbstractDict, values::AbstractDict; kwargs...)
    return sync_state_trajectories!(c_ic, nothing, values; set_start=false, kwargs...)
end

"""
    warm_start_with_constant!(controls, value)

Fill every control start value with the same scalar.

This is the simplest warm-start strategy.
It is useful on the very first MPC solve, when no previous optimal sequence
exists yet.

How it connects to other functions:
- `prepare_tracking_mpc_step!(...)` uses this on the first solve;
- later solves usually switch to `shift_warm_start!(...)`, which reuses the
  last optimized control sequence.
"""
function warm_start_with_constant!(controls::AbstractVector{<:JuMP.VariableRef}, value::Real)
    for control in controls
        JuMP.set_start_value(control, float(value))
    end
    return nothing
end

"""
    shift_warm_start!(controls, previous_values; repeat_last=true)

Shift a previous optimal sequence by one step and use it as the new warm start.

This is the standard receding-horizon warm-start strategy.
Suppose the previous optimal sequence was `[u1, u2, u3, u4]`.
On the next solve, the first move `u1` has already been applied to the plant.
So the new warm start becomes something like `[u2, u3, u4, u4]`.

Why this helps:
- the new optimization problem is usually similar to the last one;
- the last solution is therefore a good starting guess;
- good starting guesses often reduce solver work.

How it connects to other functions:
- `solve_tracking_mpc!(...)` stores the last optimal control sequence;
- `prepare_tracking_mpc_step!(...)` reads that stored sequence and calls this
  helper before the next solve.
"""
function shift_warm_start!(controls::AbstractVector{<:JuMP.VariableRef},
                           previous_values::AbstractVector{<:Real};
                           repeat_last::Bool=true)
    length(controls) == length(previous_values) ||
        error("controls and previous_values must have the same length.")

    isempty(controls) && return nothing

    # Shift the old sequence one step to the left.
    for k in 1:(length(controls) - 1)
        JuMP.set_start_value(controls[k], float(previous_values[k + 1]))
    end

    # Hold the tail value constant by default.
    tail_value = repeat_last ? previous_values[end] : previous_values[max(end - 1, 1)]
    JuMP.set_start_value(controls[end], float(tail_value))
    return nothing
end

"""
    set_first_move_bounds!(controls, previous_value; lower, upper, delta_max)

Apply a bound to the first move.
Return the `(lower_bound, upper_bound)` pair that was applied.

This helper is how the package enforces the idea that the next control move
cannot jump too far away from the control value currently applied to the plant.

Example:
- if the plant is currently at `u = 100`;
- and `delta_max = 20`;
- then the first optimized move can only lie in `[80, 120]`,
  after the ordinary lower and upper bounds are also respected.

How it connects to other functions:
- `MPCControlSpec` stores `lower`, `upper`, and `delta_max`;
- `prepare_tracking_mpc_step!(...)` calls this before every solve;
- `solve_tracking_mpc!(...)` therefore always solves with a first move that
  respects both hard bounds and move-rate limits.
"""
function set_first_move_bounds!(controls::AbstractVector{<:JuMP.VariableRef},
                                previous_value::Real;
                                lower::Real=-Inf,
                                upper::Real=Inf,
                                delta_max::Union{Nothing, Real}=nothing)
    isempty(controls) && error("controls must contain at least one JuMP variable.")

    applied_lb = float(lower)
    applied_ub = float(upper)
    if delta_max !== nothing
        # Shrink the first-step bounds to enforce the move-rate limit.
        step = float(delta_max)
        applied_lb = max(applied_lb, float(previous_value) - step)
        applied_ub = min(applied_ub, float(previous_value) + step)
    end

    JuMP.set_lower_bound(controls[1], applied_lb)
    JuMP.set_upper_bound(controls[1], applied_ub)
    return applied_lb, applied_ub
end

"""
    make_state_index_map(state_vars)

Build a `state -> index` lookup from `state_vars`.
"""
function make_state_index_map(state_vars)
    vars = collect(state_vars)
    return Dict{Any, Int}(var => i for (i, var) in pairs(vars))
end

"""
    MPCLog

Simple container for closed-loop MPC logging.

It stores:
- state indices,
- logged times,
- state histories,
- control histories,
- prediction times,
- prediction histories,
- and scalar metric histories.

Think of this as the package's "black box recorder" for one MPC run.
It does not solve anything by itself.
It simply keeps all the time histories in one place so later code can inspect,
plot, or compare them.
"""
mutable struct MPCLog
    i_state::Dict{Any, Int}
    ts::Vector{Float64}
    Xhist::Dict{Any, Vector{Float64}}
    Controlhist::Dict{Any, Vector{Float64}}
    pred_times::Vector{Float64}
    Predhist::Dict{Any, Vector{Vector{Float64}}}
    Metrichist::Dict{Symbol, Vector{Float64}}
end

"""
    make_mpc_log(state_vars; control_keys=Any[], predicted_keys=Any[], metric_keys=Symbol[])

Create a generic closed-loop MPC log.

Use this when you already know which states, controls, predicted channels, and
metrics you want to record.

How it connects to other functions:
- `log_mpc_state!(...)` appends measured plant data into this object;
- `record_mpc_prediction!(...)` appends the prediction from each MPC solve;
- `record_mpc_metrics!(...)` appends scalar objective terms and diagnostics;
- `compute_step_prediction_errors(...)` later reads this log to compare plant
  reality against earlier MPC predictions.
"""
function make_mpc_log(state_vars::AbstractVector;
                      control_keys=Any[],
                      predicted_keys=Any[],
                      metric_keys=Symbol[])
    vars = collect(state_vars)
    # Pre-create the history containers.
    i_state = make_state_index_map(vars)
    xhist = Dict{Any, Vector{Float64}}(var => Float64[] for var in vars)
    controlhist = Dict{Any, Vector{Float64}}(key => Float64[] for key in control_keys)
    predhist = Dict{Any, Vector{Vector{Float64}}}(key => Vector{Float64}[] for key in predicted_keys)
    metrichist = Dict{Symbol, Vector{Float64}}(key => Float64[] for key in metric_keys)
    return MPCLog(i_state, Float64[], xhist, controlhist, Float64[], predhist, metrichist)
end

"""
    make_mpc_log(sys; kwargs...)

Convenience method that builds `MPCLog` directly from `unknowns(sys)`.

Use this when you already have a ModelingToolkit system and do not want to
manually pass the state vector.
"""
function make_mpc_log(sys::ModelingToolkit.System; kwargs...)
    return make_mpc_log(collect(unknowns(sys)); kwargs...)
end

"""
    log_mpc_state!(logctx, integ; control_values=Dict(), missing_control=:error)

Record the current time, states, and control values.

This is the main "plant-side" logging function.
It reads the current DifferentialEquations integrator state and appends one new
record to the log.

Typical use in a closed-loop script:
1. the plant integrator advances in time;
2. the MPC callback computes or updates the control signal;
3. this function stores the new plant state and the applied control.

How it connects to other functions:
- `current_state_map(...)` can build a `state => value` dictionary from the same
  integrator;
- `seed_mpc_log!(...)` writes the initial entry before the loop starts;
- this function then appends one more entry at each later time step.
"""
function log_mpc_state!(logctx::MPCLog,
                        integ;
                        control_values::AbstractDict=Dict{Any, Any}(),
                        missing_control::Symbol=:error)
    missing_control in (:error, :nan, :hold, :skip) ||
        error("Unsupported missing_control policy $(missing_control).")

    push!(logctx.ts, float(integ.t))
    for (var, idx) in logctx.i_state
        # Use the saved index map to read the integrator state.
        push!(logctx.Xhist[var], float(integ.u[idx]))
    end

    for (key, hist) in logctx.Controlhist
        if haskey(control_values, key)
            push!(hist, float(control_values[key]))
        elseif missing_control === :nan
            push!(hist, NaN)
        elseif missing_control === :hold
            push!(hist, isempty(hist) ? NaN : hist[end])
        elseif missing_control === :skip
            continue
        else
            error("Missing control value for key $(key).")
        end
    end

    return nothing
end

"""
    record_mpc_prediction!(logctx, t, predictions; missing=:error)

Store the prediction trajectories from one MPC solve at time `t`.

This is the "optimizer-side" logging function.
Unlike `log_mpc_state!(...)`, which stores what the plant actually did, this
function stores what the MPC predicted at one solve time.

Why this matters:
- after the run, you can ask whether the MPC predictions were accurate;
- you can compare the predicted next step with the later measured plant value;
- that makes debugging much easier.

How it connects to other functions:
- `solve_tracking_mpc!(...)` returns a `predictions` dictionary;
- this function stores that dictionary in `MPCLog`;
- `compute_step_prediction_errors(...)` later compares those stored predictions
  against the logged plant state history.
"""
function record_mpc_prediction!(logctx::MPCLog,
                                t::Real,
                                predictions::AbstractDict;
                                missing::Symbol=:error)
    missing in (:error, :skip) || error("Unsupported missing policy $(missing). Use :error or :skip.")

    previous_len = length(logctx.pred_times)
    for key in keys(predictions)
        if !haskey(logctx.Predhist, key)
            # Backfill new prediction channels with empty records.
            logctx.Predhist[key] = [Float64[] for _ in 1:previous_len]
        end
    end

    push!(logctx.pred_times, float(t))
    for (key, hist) in logctx.Predhist
        if haskey(predictions, key)
            push!(hist, collect(float.(predictions[key])))
        elseif missing === :skip
            push!(hist, Float64[])
        else
            error("Missing prediction for key $(key).")
        end
    end

    return nothing
end

"""
    record_mpc_metrics!(logctx, metrics; missing=:error)

Append scalar MPC diagnostics such as tracking cost or move cost.

Use this when you want more than just states and controls.
For example, you may want to record:
- the full objective value;
- the tracking penalty;
- the soft-constraint penalty;
- the control-move penalty.

How it connects to other functions:
- `solve_tracking_mpc!(...)` already builds a `metrics` dictionary;
- this function appends those scalar values into the log;
- later analysis scripts can plot how each objective component changed over
  time.
"""
function record_mpc_metrics!(logctx::MPCLog,
                             metrics::AbstractDict{<:Symbol, <:Real};
                             missing::Symbol=:error)
    missing in (:error, :skip) || error("Unsupported missing policy $(missing). Use :error or :skip.")

    for key in keys(metrics)
        get!(logctx.Metrichist, key, Float64[])
    end

    for (key, hist) in logctx.Metrichist
        if haskey(metrics, key)
            push!(hist, float(metrics[key]))
        elseif missing === :skip
            continue
        else
            error("Missing metric for key $(key).")
        end
    end

    return nothing
end

"""
    default_mpc_accepted_statuses()

Return the default JuMP/MOI solve statuses treated as acceptable for MPC.
"""
function default_mpc_accepted_statuses()
    return (
        MOI.OPTIMAL,
        MOI.LOCALLY_SOLVED,
        MOI.ALMOST_OPTIMAL,
        MOI.ALMOST_LOCALLY_SOLVED,
    )
end

"""
    is_accepted_mpc_status(status; accepted_statuses=default_mpc_accepted_statuses())

Return `true` when `status` is in the accepted MPC status set.
"""
function is_accepted_mpc_status(status; accepted_statuses=default_mpc_accepted_statuses())
    return status in accepted_statuses
end

"""
    objective_value_or_nan(model)

Return the JuMP objective value as `Float64`.
Return `NaN` if the objective is not available.
"""
function objective_value_or_nan(model::JuMP.Model)::Float64
    JuMP.has_values(model) || return NaN
    return try
        float(JuMP.objective_value(model))
    catch
        NaN
    end
end

"""
    summarize_mpc_solve(model; accepted_statuses=default_mpc_accepted_statuses())

Return a short solve summary `(status, accepted, objective)`.
"""
function summarize_mpc_solve(model::JuMP.Model; accepted_statuses=default_mpc_accepted_statuses())
    status = JuMP.termination_status(model)
    return (
        status = status,
        accepted = is_accepted_mpc_status(status; accepted_statuses=accepted_statuses),
        objective = objective_value_or_nan(model),
    )
end

"""
    _time_index(times, target; match=:nearest, atol=1e-8)

Find the index of a target time in a logged time vector.

Use `match=:nearest` for the closest entry or `match=:exact` when you want the
time to match within a tolerance.
"""
function _time_index(times::AbstractVector{<:Real}, target::Real; match::Symbol=:nearest, atol::Real=1e-8)
    isempty(times) && error("times must contain at least one entry.")
    match in (:nearest, :exact) || error("Unsupported match policy $(match). Use :nearest or :exact.")

    diffs = abs.(float.(times) .- float(target))
    idx = argmin(diffs)
    if match === :exact && diffs[idx] > atol
        error("No logged time matches target=$(target) within atol=$(atol).")
    end
    return idx
end

"""
    compute_step_prediction_errors(logctx, key; actual_key=key, step_index=2, match=:nearest, atol=1e-8)

Compare one predicted step against the later logged actual value.
Return times, predicted values, actual values, and errors.

This helper is useful after a closed-loop run is finished.
It answers a simple question:
"When the MPC predicted the plant state several steps ahead, how wrong was that
prediction when the future actually arrived?"

Example:
- at solve time `t_k`, the MPC predicts `x[k+1]`, `x[k+2]`, and so on;
- later, the plant actually reaches `t_{k+1}`, `t_{k+2}`, and so on;
- this helper lines up those two records and computes the difference.

How it connects to other functions:
- `record_mpc_prediction!(...)` must already have stored the prediction
  histories;
- `log_mpc_state!(...)` must already have stored the actual plant histories;
- this function reads both logs and returns aligned arrays for later plotting or
  statistics.
"""
function compute_step_prediction_errors(logctx::MPCLog,
                                        key;
                                        actual_key=key,
                                        step_index::Int=2,
                                        match::Symbol=:nearest,
                                        atol::Real=1e-8)
    step_index >= 1 || error("step_index must be at least 1.")
    haskey(logctx.Predhist, key) || error("No prediction history stored for key $(key).")
    haskey(logctx.Xhist, actual_key) || error("No state history stored for key $(actual_key).")

    pred_hist = logctx.Predhist[key]
    max_records = min(length(pred_hist), length(logctx.pred_times)) - (step_index - 1)

    times = Float64[]
    predicted = Float64[]
    actual = Float64[]
    errors = Float64[]

    for k in 1:max(0, max_records)
        traj = pred_hist[k]
        length(traj) >= step_index || continue
        t_target = logctx.pred_times[k + step_index - 1]
        # Match the predicted step to the logged plant value at that time.
        idx = _time_index(logctx.ts, t_target; match=match, atol=atol)
        pred_val = float(traj[step_index])
        act_val = float(logctx.Xhist[actual_key][idx])
        push!(times, t_target)
        push!(predicted, pred_val)
        push!(actual, act_val)
        push!(errors, act_val - pred_val)
    end

    return (times=times, predicted=predicted, actual=actual, errors=errors)
end

"""
    _affexpr_var_coeff_pairs(func)

Extract `(variable, coefficient)` pairs from a `JuMP.AffExpr`.

This is used by the constraint-audit functions later in this file.
"""
function _affexpr_var_coeff_pairs(func::JuMP.AffExpr)
    pairs = Tuple{JuMP.VariableRef, Float64}[]

    for term in JuMP.linear_terms(func)
        if term isa Tuple && length(term) == 2
            first_term, second_term = term
            if second_term isa JuMP.VariableRef
                push!(pairs, (second_term, float(first_term)))
            elseif first_term isa JuMP.VariableRef
                push!(pairs, (first_term, float(second_term)))
            end
        elseif term isa Pair
            if term.second isa JuMP.VariableRef
                push!(pairs, (term.second, float(term.first)))
            elseif term.first isa JuMP.VariableRef
                push!(pairs, (term.first, float(term.second)))
            end
        elseif hasproperty(term, :variable) && hasproperty(term, :coefficient)
            push!(pairs, (getproperty(term, :variable), float(getproperty(term, :coefficient))))
        elseif hasproperty(term, :first) && hasproperty(term, :second)
            if getproperty(term, :second) isa JuMP.VariableRef
                push!(pairs, (getproperty(term, :second), float(getproperty(term, :first))))
            elseif getproperty(term, :first) isa JuMP.VariableRef
                push!(pairs, (getproperty(term, :first), float(getproperty(term, :second))))
            end
        end
    end

    return pairs
end

"""
    _single_var_from_affexpr(func)

Return the one variable term from an affine expression when it has exactly one.

If an affine expression contains exactly one variable term, this returns that
variable and its coefficient. Otherwise it returns `(nothing, 0.0)`.
"""
function _single_var_from_affexpr(func::JuMP.AffExpr)
    pairs = _affexpr_var_coeff_pairs(func)
    return length(pairs) == 1 ? first(pairs) : (nothing, 0.0)
end

"""
    find_conflicts(model; atol=1e-8)

Check a JuMP model for contradictory single-variable affine equalities and bounds.
Returns a vector of human-readable problem strings.
"""
function find_conflicts(model::JuMP.Model; atol::Float64=1e-8)
    problems = String[]
    eq_rhs_by_var = Dict{JuMP.VariableRef, Vector{Float64}}()
    lbs = Dict{JuMP.VariableRef, Float64}()
    ubs = Dict{JuMP.VariableRef, Float64}()

    for cref in JuMP.all_constraints(model, JuMP.AffExpr, MOI.EqualTo{Float64})
        func = JuMP.constraint_object(cref).func
        var, coef = _single_var_from_affexpr(func)
        isnothing(var) && continue
        isapprox(coef, 1.0; atol=atol) || continue
        rhs = JuMP.constraint_object(cref).set.value
        push!(get!(eq_rhs_by_var, var, Float64[]), rhs)
    end

    for cref in JuMP.all_constraints(model, JuMP.VariableRef, MOI.GreaterThan{Float64})
        var = JuMP.constraint_object(cref).func
        lb = JuMP.constraint_object(cref).set.lower
        lbs[var] = haskey(lbs, var) ? max(lbs[var], lb) : lb
    end

    for cref in JuMP.all_constraints(model, JuMP.VariableRef, MOI.LessThan{Float64})
        var = JuMP.constraint_object(cref).func
        ub = JuMP.constraint_object(cref).set.upper
        ubs[var] = haskey(ubs, var) ? min(ubs[var], ub) : ub
    end

    for cref in JuMP.all_constraints(model, JuMP.AffExpr, MOI.GreaterThan{Float64})
        func = JuMP.constraint_object(cref).func
        var, coef = _single_var_from_affexpr(func)
        isnothing(var) && continue
        isapprox(coef, 1.0; atol=atol) || continue
        lb = JuMP.constraint_object(cref).set.lower
        lbs[var] = haskey(lbs, var) ? max(lbs[var], lb) : lb
    end

    for cref in JuMP.all_constraints(model, JuMP.AffExpr, MOI.LessThan{Float64})
        func = JuMP.constraint_object(cref).func
        var, coef = _single_var_from_affexpr(func)
        isnothing(var) && continue
        isapprox(coef, 1.0; atol=atol) || continue
        ub = JuMP.constraint_object(cref).set.upper
        ubs[var] = haskey(ubs, var) ? min(ubs[var], ub) : ub
    end

    for cref in JuMP.all_constraints(model, JuMP.AffExpr, MOI.Interval{Float64})
        func = JuMP.constraint_object(cref).func
        var, coef = _single_var_from_affexpr(func)
        isnothing(var) && continue
        isapprox(coef, 1.0; atol=atol) || continue
        set = JuMP.constraint_object(cref).set
        lbs[var] = haskey(lbs, var) ? max(lbs[var], set.lower) : set.lower
        ubs[var] = haskey(ubs, var) ? min(ubs[var], set.upper) : set.upper
    end

    for (var, rhs_vals) in eq_rhs_by_var
        uniq_rhs = unique(round.(rhs_vals; digits=10))
        if length(uniq_rhs) >= 2
            push!(problems, "Conflict: $(JuMP.name(var)) has multiple equalities: " * join(string.(uniq_rhs), ", "))
        end
    end

    for var in union(keys(lbs), keys(ubs))
        if haskey(lbs, var) && haskey(ubs, var) && lbs[var] > ubs[var] + atol
            push!(problems, "Conflict: infeasible bounds on $(JuMP.name(var)): lb=$(lbs[var]) > ub=$(ubs[var]).")
        end
    end

    for (var, rhs_vals) in eq_rhs_by_var
        lb = get(lbs, var, -Inf)
        ub = get(ubs, var, Inf)
        for rhs in rhs_vals
            if rhs < lb - atol
                push!(problems, "Conflict: $(JuMP.name(var)) == $(rhs) violates lb=$(lb).")
            elseif rhs > ub + atol
                push!(problems, "Conflict: $(JuMP.name(var)) == $(rhs) violates ub=$(ub).")
            end
        end
    end

    return problems
end

"""
    print_conflicts(model; atol=1e-8)

Pretty-print the output of `find_conflicts`.
"""
function print_conflicts(model::JuMP.Model; atol::Float64=1e-8)
    problems = find_conflicts(model; atol=atol)
    if isempty(problems)
        println("Constraint audit: clean")
    else
        println("Constraint audit: found $(length(problems)) issue(s):")
        for problem in problems
            println(" - ", problem)
        end
    end
    return nothing
end

"""
    dump_single_var_affines(model; atol=1e-8)

List single-variable affine equalities and bounds grouped by variable.
"""
function dump_single_var_affines(model::JuMP.Model; atol::Float64=1e-8)
    by_var = Dict{JuMP.VariableRef, Vector{String}}()

    for cref in JuMP.all_constraints(model, JuMP.AffExpr, MOI.EqualTo{Float64})
        func = JuMP.constraint_object(cref).func
        var, coef = _single_var_from_affexpr(func)
        isnothing(var) && continue
        isapprox(coef, 1.0; atol=atol) || continue
        rhs = JuMP.constraint_object(cref).set.value
        push!(get!(by_var, var, String[]), "$(JuMP.name(var)) == $(rhs)")
    end

    for cref in JuMP.all_constraints(model, JuMP.AffExpr, MOI.GreaterThan{Float64})
        func = JuMP.constraint_object(cref).func
        var, coef = _single_var_from_affexpr(func)
        isnothing(var) && continue
        isapprox(coef, 1.0; atol=atol) || continue
        lb = JuMP.constraint_object(cref).set.lower
        push!(get!(by_var, var, String[]), "$(JuMP.name(var)) >= $(lb)")
    end

    for cref in JuMP.all_constraints(model, JuMP.AffExpr, MOI.LessThan{Float64})
        func = JuMP.constraint_object(cref).func
        var, coef = _single_var_from_affexpr(func)
        isnothing(var) && continue
        isapprox(coef, 1.0; atol=atol) || continue
        ub = JuMP.constraint_object(cref).set.upper
        push!(get!(by_var, var, String[]), "$(JuMP.name(var)) <= $(ub)")
    end

    for cref in JuMP.all_constraints(model, JuMP.AffExpr, MOI.Interval{Float64})
        func = JuMP.constraint_object(cref).func
        var, coef = _single_var_from_affexpr(func)
        isnothing(var) && continue
        isapprox(coef, 1.0; atol=atol) || continue
        set = JuMP.constraint_object(cref).set
        push!(get!(by_var, var, String[]), "$(set.lower) <= $(JuMP.name(var)) <= $(set.upper)")
    end

    for (var, lines) in by_var
        println("* ", JuMP.name(var))
        for line in lines
            println("  - ", line)
        end
    end

    return nothing
end

"""
    _value_affine(model, func)

Internal fallback evaluator for `MOI.ScalarAffineFunction`.

This is used when `JuMP.value(func)` is not available for a particular affine
function object.
"""
function _value_affine(model::JuMP.Model, func::MOI.ScalarAffineFunction{Float64})
    value = func.constant
    for term in func.terms
        var = JuMP.VariableRef(model, term.variable_index)
        value += term.coefficient * JuMP.value(var)
    end
    return value
end

"""
    worst_residuals(model; topk=20)

Return the largest affine equality residuals and bound violations after a solve.
"""
function worst_residuals(model::JuMP.Model; topk::Int=20)
    residuals = String[]
    pairs = Tuple{Float64, JuMP.ConstraintRef}[]

    for cref in JuMP.all_constraints(model, MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64})
        func = JuMP.constraint_object(cref).func
        rhs = JuMP.constraint_object(cref).set.value
        value = try
            JuMP.value(func)
        catch
            _value_affine(model, func)
        end
        push!(pairs, (abs(value - rhs), cref))
    end

    sort!(pairs; by=first, rev=true)
    for (residual, cref) in Iterators.take(pairs, topk)
        push!(residuals, "equality residual $(residual) : $(cref)")
    end

    for var in JuMP.all_variables(model)
        value = JuMP.value(var)
        if JuMP.has_lower_bound(var) && value < JuMP.lower_bound(var) - 1e-8
            push!(residuals, "lower-bound violation $(JuMP.lower_bound(var) - value) : $(var)")
        end
        if JuMP.has_upper_bound(var) && value > JuMP.upper_bound(var) + 1e-8
            push!(residuals, "upper-bound violation $(value - JuMP.upper_bound(var)) : $(var)")
        end
    end

    return residuals
end

"""
    summarize_ic_uniqueness(model)

Count how many first-index trajectory variables have exactly one affine equality
constraint and report any duplicates.
"""
function summarize_ic_uniqueness(model::JuMP.Model)
    counts = Dict{String, Int}()

    for cref in JuMP.all_constraints(model, JuMP.AffExpr, MOI.EqualTo{Float64})
        func = JuMP.constraint_object(cref).func
        var, coef = _single_var_from_affexpr(func)
        isnothing(var) && continue
        isapprox(coef, 1.0; atol=1e-8) || continue
        name = JuMP.name(var)
        occursin(r"\[\s*1\s*\]", name) || continue
        base = replace(name, r"\[\s*1\s*\]" => "[1]")
        counts[base] = get(counts, base, 0) + 1
    end

    duplicates = filter(pair -> pair[2] != 1, collect(counts))
    return (total_ic_vars=length(counts), non_unique=duplicates)
end

"""
    check_ic_sync(current_values, ic_constraints; atol=1e-8)

Compare a dictionary of current values against the RHS stored in initial-condition
constraints. Returns a named tuple with the maximum absolute error and all per-key
differences.
"""
function check_ic_sync(current_values::AbstractDict, ic_constraints::AbstractDict; atol::Real=1e-8)
    diffs = Dict{Any, Float64}()

    for key in intersect(Set(keys(current_values)), Set(keys(ic_constraints)))
        rhs = JuMP.normalized_rhs(ic_constraints[key])
        diffs[key] = abs(float(current_values[key]) - rhs)
    end

    max_err = isempty(diffs) ? 0.0 : maximum(values(diffs))
    return (max_err=max_err, diffs=diffs, ok=max_err <= atol)
end
