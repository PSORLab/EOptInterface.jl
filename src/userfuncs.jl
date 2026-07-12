# ------------------------------------------------------------------------------
# ModelingToolkit-to-JuMP registration routines
#
# This file turns a ModelingToolkit system into a JuMP model.
# It does not choose controls or run the online loop. It takes one MTK model
# and writes the matching JuMP variables and constraints.
#
# Main flow:
# 1. `decision_vars(...)` picks the MTK symbols used in the problem.
# 2. `register_odesystem(...)` or `register_daesystem(...)` adds the dynamics.
# 3. `build_tracking_mpc(...)` builds the full MPC problem on top of that.
# 4. `full_solutions(...)` maps solved values back to MTK symbols.
#
# Good reading order for new users:
# 1. Read `decision_vars(...)` first.
#    It answers the simple question: "Which MTK symbols become optimization
#    variables?"
# 2. Then read `register_odesystem(...)` or `register_daesystem(...)`.
#    These functions answer the next question: "How do those variables satisfy
#    the dynamic model over time?"
# 3. Read `full_solutions(...)` last.
#    It answers the reporting question: "How do I translate solved JuMP values
#    back into MTK notation?"
# ------------------------------------------------------------------------------

"""
    _canonical_parameter_var_dict(sys, data)

Rewrite a dictionary so all keys are the standard parameter symbols from
`parameters(sys)`.

This matters because users may pass different but equivalent MTK names.
The rest of the package wants one consistent key form before it writes JuMP
variables into the MPC model.
"""
function _canonical_parameter_var_dict(sys, data::AbstractDict)
    out = Dict{Num, Vector{JuMP.VariableRef}}()
    for (sym, vars) in pairs(data)
        # Convert each input name to the standard MTK parameter key.
        out[canonical_system_parameter(sys, sym)] = vars
    end
    return out
end

"""
    _canonical_state_var_dict(sys, data)

Rewrite a dictionary so all keys are the standard state symbols from
`unknowns(sys)`.

This is the state-side companion to `_canonical_parameter_var_dict(...)`.
It keeps later IC and trajectory updates keyed by one canonical MTK symbol set.
"""
function _canonical_state_var_dict(sys, data::AbstractDict)
    out = Dict{Num, Vector{JuMP.VariableRef}}()
    for (sym, vars) in pairs(data)
        # Do the same thing for states.
        out[canonical_system_unknown(sys, sym)] = vars
    end
    return out
end

# original decision_vars function that just calls the more general one with defaults
"""
    decision_vars(sys)
    decision_vars(sys, ps; model=nothing, horizon=nothing, build_state_trajs=false, ...)

Identify the ModelingToolkit symbols that should participate in an MPC or other
dynamic optimization problem.

What this function returns depends on `build_state_trajs`:
- `build_state_trajs=false`:
  returns a flat vector of symbols, consisting of
  1. all MTK states (`unknowns(sys)`), and
  2. the parameters in `ps`, plus any parameters without defaults.
- `build_state_trajs=true`:
  also creates one JuMP trajectory for every state and returns a named tuple with
  `all`, `states`, `params`, `x_vars`, and `c_ic`.

Why this function exists:
- it gives the package one place to decide which MTK symbols go into the NLP;
- it normalizes symbols, so callers can use either raw MTK symbols or names
  like `sys.x`;
- it can also build the state trajectories used by the ODE and DAE functions.

How it connects to other functions:
- `build_tracking_mpc(...)` calls this to get `x_vars` and `c_ic`;
- `register_odesystem(...)` and `register_daesystem(...)` use those
  trajectories;
- `full_solutions(...)` uses the same symbol order when it reads values back.

Algorithm:
1. Read all MTK states from `unknowns(sys)`.
2. Normalize the user-supplied parameter names in `ps`.
3. Add any parameters that have no default value, because the NLP still needs
   them.
4. If `build_state_trajs=false`, return the final symbol list.
5. If `build_state_trajs=true`, call `build_state_trajs_from_vars!(...)` to
   allocate one JuMP trajectory and one IC constraint per state.
6. Return the symbol list together with the trajectory containers.
"""
function decision_vars(sys::ModelingToolkit.System)
    return decision_vars(sys, Num[];
                         model=nothing,
                         horizon=nothing,
                         build_state_trajs=false,
                         lb=0.0,
                         ub=1e6,
                         rhs0=0.0,
                         store_ext=true)
end

function decision_vars(sys::ModelingToolkit.System, ps::Vector{Num};
                       model::Union{Nothing, JuMP.Model}=nothing,
                       horizon::Union{Nothing, Int}=nothing,
                       build_state_trajs::Bool=false,
                       lb=0.0,
                       ub=1e6,
                       rhs0=0.0,
                       store_ext::Bool=true)
    # All MTK unknowns become state trajectories.
    state_vars = collect(ModelingToolkit.unknowns(sys))
    # Normalize parameter names first.
    ps_canonical = canonicalize_system_symbols(sys, ps; pool=:parameters)
    # Parameters without defaults must also be optimized.
    params_to_opt = collect(union(ps_canonical, setdiff(ModelingToolkit.parameters(sys),
                                                        keys(ModelingToolkit.defaults(sys)))))
    all_vars = vcat(state_vars, params_to_opt)

    build_state_trajs || return all_vars
    isnothing(model) && error("decision_vars(...; build_state_trajs=true) requires `model`.")
    isnothing(horizon) && error("decision_vars(...; build_state_trajs=true) requires `horizon`.")

    # Build the JuMP state trajectories in the shared helper.
    x_vars, c_ic = build_state_trajs_from_vars!(
        model,
        sys,
        state_vars,
        horizon;
        lb = lb,
        ub = ub,
        rhs0 = rhs0,
        store_ext = store_ext,
    )

    return (
        all = all_vars,
        states = state_vars,
        params = params_to_opt,
        x_vars = x_vars,
        c_ic = c_ic,
    )
end
"""
    register_nlsystem(model, sys, obj, ineqs)

Register a purely algebraic ModelingToolkit system inside a JuMP model.

Use this when the model has no time dynamics.
It adds the MTK equations, the inequality constraints, and the objective.

How it connects to other functions:
- unlike the ODE and DAE functions, this one does not build trajectories;
- it uses the same MTK-to-JuMP conversion tools;
- it is separate from `build_tracking_mpc(...)`, which is for time-based MPC.
"""
function register_nlsystem(model::JuMP.Model, sys::ModelingToolkit.System, obj::Symbolics.Num, ineqs::Vector{Symbolics.Num})
    # Turn MTK expressions into callable JuMP expressions.
    h = EOptInterface.mtk_generate_model_equations(sys)
    f = EOptInterface.mtk_generate_reduced_expression(obj, sys)
    g = []
    for i in eachindex(ineqs)
        gi = EOptInterface.mtk_generate_reduced_expression(ineqs[i], sys)
        push!(g, gi)
    end
    JuMP.@constraint(model, [i in eachindex(h)], h[i](JuMP.all_variables(model)...) == 0)
    JuMP.@constraint(model, [i in eachindex(g)], g[i](JuMP.all_variables(model)...) >= 0)
    JuMP.@objective(model, Min, f(JuMP.all_variables(model)...))
end

"""
    ensure_ic!(model, ic_map, x_vars, u, val)

Ensure that the first point of a state trajectory has the chosen value.

This function either updates an existing initial-condition constraint or
creates it.

How it connects to other functions:
- `decision_vars(...; build_state_trajs=true)` and
  `build_state_trajs_from_vars!(...)` create the initial-condition constraints
  stored in `c_ic`;
- `prepare_tracking_mpc_step!(...)` updates those constraints before each solve;
- `get_ic_constraint!(...)` checks that the IC constraint is present.
"""
function ensure_ic!(model::JuMP.Model,
                    ic_map::Dict{Any, JuMP.ConstraintRef},
                    x_vars::Dict, u, val)
    if haskey(ic_map, u)
        # Normal path: update the existing IC constraint.
        JuMP.set_normalized_rhs(ic_map[u], val)
    else
        # Fallback path: create the IC constraint if it is missing.
        con = @constraint(model, x_vars[u][1] == val)
        ic_map[u] = con
    end
    return nothing
end

"""
    to_num_key_map(u0_dict)

Return a copy of `u0_dict` with all keys converted to `Num`.

This is a small compatibility helper for older code paths that may still store
initial-condition dictionaries under other symbolic key types.
"""
function to_num_key_map(u0_dict)
    # u0_dict :: Dict{BasicSymbolic{Real}, Float64}
    # return :: Dict{Num, Float64}
    Dict( Num(k) => v for (k,v) in u0_dict )
end

"""
    get_ic_constraint!(model, x; idx=1, rhs_if_missing=nothing)

Return the equality constraint for the first point of a trajectory.

This helper looks for an equality on `x[idx]`.
It returns the constraint if it exists.
It can create the constraint if it is missing.
It can also delete duplicates.

How it connects to other functions:
- `build_state_trajs_from_vars!(...)` creates the original IC constraints;
- `ensure_ic!(...)` updates or adds them during online MPC use;
- this helper is mainly for checking and repair.
"""
function get_ic_constraint!(model::JuMP.Model,
                            x::AbstractVector{<:JuMP.VariableRef};
                            idx::Int=1,
                            rhs_if_missing::Union{Nothing,Float64}=nothing)
    @assert 1 <= idx <= length(x) "idx out of bounds"
    target = x[idx]
    target_idx = JuMP.index(target)

    hits = JuMP.ConstraintRef[]
    for c in JuMP.all_constraints(model, JuMP.AffExpr, MOI.EqualTo{Float64})
        co = JuMP.constraint_object(c)
        terms = co.func.terms
        # Only keep simple one-variable equalities.
        # Support both dict-like and vector-like term storage (JuMP version differences).
        if length(terms) == 1
            v = begin
                if terms isa AbstractDict
                    first(keys(terms))
                else
                    first(terms).first
                end
            end
            if JuMP.index(v) == target_idx
                push!(hits, c)
            end
        end
    end

    if isempty(hits)
        isnothing(rhs_if_missing) && error("No scalar IC equality found for $(JuMP.name(target))")
        # Create the missing IC if the caller asks for it.
        return @constraint(model, target == rhs_if_missing)
    end

    # Keep one IC constraint and delete duplicates.
    for c in hits[2:end]
        JuMP.delete(model, c)
    end
    return first(hits)
end

"""
    register_odesystem(model, sys, tspan, tstep, integrator; ...)

Register a ModelingToolkit ODE system as JuMP constraints over time.

This function turns a continuous ODE into the step equations used by MPC.

Supported integrators:
- `"EE"`: Explicit Euler
- `"IE"` / `"BDF1"`: Implicit Euler
- `"RK4"`: Classical fourth-order Runge-Kutta
- `"IRK4"`: two-stage Gauss-Legendre implicit Runge-Kutta

How it connects to other functions:
- `decision_vars(...; build_state_trajs=true)` or `build_state_trajs_from_vars!(...)`
  creates the `x_vars` trajectories used here;
- `_register_tracking_system!(...)` in `trackingmpc.jl` calls this when the user
  selects `system_kind=:ode`;
- `prepare_tracking_mpc_step!(...)` and `solve_tracking_mpc!(...)` do not
  rebuild these constraints. They only update values and solve again.

Algorithm:
1. Build the time grid from `tspan` and `tstep`.
2. Normalize state and parameter names so all later lookups use one MTK key.
3. Turn each MTK right-hand side into a callable function of states, parameters,
   and time.
4. Read or infer the JuMP state trajectories `x_vars`.
5. For each time step, add one set of discretization equations.
6. The exact equations depend on `integrator`, for example Euler, RK4, or IRK4.
7. Leave the built constraints in the model for repeated MPC solves.
"""

function register_odesystem(model::JuMP.Model,
                            odesys::ModelingToolkit.System,
                            tspan::Tuple{Real,Real},
                            tstep::Real,
                            integrator::String;
                            p_disc::Vector{Num}=Num[],
                            p_disc_vars::Dict{Num,Vector{JuMP.VariableRef}}=Dict{Num,Vector{JuMP.VariableRef}}(),
                            x_vars::AbstractDict=Dict(),
                            t_map::Dict{SymbolicUtils.BasicSymbolic{Real},JuMP.VariableRef}=Dict{SymbolicUtils.BasicSymbolic{Real},JuMP.VariableRef}())

    # Time grid and dimensions
    N = Int(floor((tspan[2] - tspan[1]) / tstep)) + 1
    # State variable dimensions
    V = length(ModelingToolkit.unknowns(odesys))
    # Normalize all user input names once at the start.
    p_disc = canonicalize_system_symbols(odesys, p_disc; pool=:parameters)
    p_disc_vars = _canonical_parameter_var_dict(odesys, p_disc_vars)
    x_vars = _canonical_state_var_dict(odesys, x_vars)

    # Normalize integrator key
    intg = uppercase(integrator)

    # 1) Build callable right-hand-side functions.
    param_dict = copy(ModelingToolkit.defaults(odesys))
    for var in ModelingToolkit.unknowns(odesys)
        pop!(param_dict, var, nothing)
    end
    for p in p_disc
        pop!(param_dict, p, nothing)
    end

    t_MTK = ModelingToolkit.get_iv(odesys)
    # Preserve the legacy no-keyword use case: if no JuMP time variable is
    # passed in, evaluate the MTK RHS at numeric times on the discretization grid.
    t_base = haskey(t_map, t_MTK) ? t_map[t_MTK] : tspan[1]
    t_at(i::Int; c::Real=0.0) = t_base + (i-1 + c) * tstep

    params_to_opt = setdiff(decision_vars(odesys, p_disc), ModelingToolkit.unknowns(odesys))
    static_param_syms = setdiff(params_to_opt, p_disc)

    dx = Vector{Function}(undef, V)
    dx_exprs = Vector{Symbolics.Num}(undef, V)
    full_eqs = ModelingToolkit.full_equations(odesys)
    for j in 1:V
        # Turn each MTK right-hand side into a callable function.
        dxj_expr = full_eqs[j].rhs
        while !isempty(intersect(Symbolics.get_variables(dxj_expr), keys(param_dict)))
            dxj_expr = SymbolicUtils.substitute(dxj_expr, param_dict)
        end
        dx_exprs[j] = dxj_expr
        dx[j] = build_function(dxj_expr,
                            decision_vars(odesys, p_disc)..., t_MTK;
                            expression = Val{false})
    end

    # Jacobian function for Rosenbrock methods.
    Jfun = nothing
    if intg in ("ROS2", "ROS", "ROSENBROCK")
        J_expr = Symbolics.jacobian(dx_exprs, ModelingToolkit.unknowns(odesys))
        Jfun = build_function(J_expr, decision_vars(odesys, p_disc)..., t_MTK; expression = Val{false})
    end

    # 2) Check the discretized parameter arrays.
    for p in p_disc
        @assert haskey(p_disc_vars, p) "Missing p_disc_vars[$p] for $p (should have length N)"
        @assert length(p_disc_vars[p]) == N "p_disc_vars[$p] must have length N=$N"
    end
    flat_pvars = reduce(vcat, values(p_disc_vars); init=JuMP.VariableRef[])

    # 3) Build the state matrix `xs`.
    unknowns = ModelingToolkit.unknowns(odesys)
    xs = Array{JuMP.VariableRef}(undef, V, N)
    static_param_vars = JuMP.VariableRef[]
    if !isempty(x_vars)
        for (j, u) in enumerate(unknowns)
            @assert haskey(x_vars, u) "Missing x_vars[$u]"
            @assert length(x_vars[u]) == N "x_vars[$u] must have length N=$N"
            xs[j, :] = x_vars[u]
        end
        if !isempty(static_param_syms)
            state_var_refs = vec(xs)
            extra_vars = setdiff(JuMP.all_variables(model), vcat(state_var_refs, flat_pvars))
            @assert length(extra_vars) >= length(static_param_syms) "register_odesystem could not infer static parameter variables from the model. Pass p_disc_vars=... or build the static parameter variables explicitly before calling register_odesystem."
            static_param_vars = extra_vars[1:length(static_param_syms)]
        end
    else
        # Old fallback path. Passing `x_vars` is safer.
        # Preserve the original examples that:
        # 1. add the state trajectories first, and
        # 2. append free design parameters afterwards.
        all_vars = JuMP.all_variables(model)
        state_var_count = V * N
        fallback_vars = setdiff(all_vars, flat_pvars)
        if length(fallback_vars) == state_var_count
            xs = reshape(fallback_vars, V, N)
        else
            legacy_param_count = length(setdiff(decision_vars(odesys), unknowns))
            if length(all_vars) >= state_var_count + legacy_param_count
                xs = reshape(all_vars[1:state_var_count], V, N)
                if !isempty(static_param_syms)
                    @assert length(all_vars) >= state_var_count + length(static_param_syms) "register_odesystem legacy fallback could not infer static parameter variables."
                    static_param_vars = all_vars[(state_var_count + 1):(state_var_count + length(static_param_syms))]
                end
            else
                error("register_odesystem fallback could not infer a $(V)x$(N) state matrix from $(length(all_vars)) JuMP variables. Pass x_vars=... explicitly.")
            end
        end
        if isempty(static_param_vars) && !isempty(static_param_syms)
            extra_vars = setdiff(all_vars, vcat(vec(xs), flat_pvars))
            @assert length(extra_vars) >= length(static_param_syms) "register_odesystem legacy fallback could not infer static parameter variables."
            static_param_vars = extra_vars[1:length(static_param_syms)]
        end
        @warn "register_odesystem: inferring state ordering from all_variables; pass x_vars=... to avoid mismatches."
        # Preserve the legacy API path used by the older examples: if the caller
        # did not pass explicit state trajectories, also create/repair the
        # initial-condition equalities from MTK defaults.
        defaults = ModelingToolkit.defaults(odesys)
        for (j, u) in enumerate(unknowns)
            haskey(defaults, u) || continue
            get_ic_constraint!(model, xs[j, :]; idx=1, rhs_if_missing=float(defaults[u]))
        end
    end


    # 4) Add the ODE step constraints.
    for i in 1:(N-1)
        p_args_i   = [p_disc_vars[p][i] for p in p_disc]
        all_p_args_i = vcat(p_args_i, static_param_vars)

        if intg == "EE"
            for j in 1:V
                # Explicit Euler uses the current slope.
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + tstep * dx[j](xs[:, i]..., all_p_args_i..., t_at(i)))
            end

        elseif intg == "IE"
            for j in 1:V
                # Implicit Euler uses the next slope.
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + tstep * dx[j](xs[:, i+1]..., all_p_args_i..., t_at(i, c=1.0)))
            end

        elseif intg == "RK4"
            # RK4 uses four slope estimates.
            k1 = [dx[j](xs[:, i]..., all_p_args_i..., t_at(i)) for j in 1:V]
            xk1_half = [xs[j, i] + 0.5 * tstep * k1[j] for j in 1:V]

            k2 = [dx[j](xk1_half..., all_p_args_i..., t_at(i, c=0.5)) for j in 1:V]
            xk2_half = [xs[j, i] + 0.5 * tstep * k2[j] for j in 1:V]

            k3 = [dx[j](xk2_half..., all_p_args_i..., t_at(i, c=0.5)) for j in 1:V]
            xk3_full = [xs[j, i] + tstep * k3[j] for j in 1:V]

            k4 = [dx[j](xk3_full..., all_p_args_i..., t_at(i, c=1.0)) for j in 1:V]

            for j in 1:V
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + (tstep/6) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]))
            end

        elseif intg == "IRK4"
            # Implicit RK order 4 (2-stage Gauss-Legendre)
            c1 = 0.5 - sqrt(3)/6
            c2 = 0.5 + sqrt(3)/6
            a11 = 1/4
            a12 = 1/4 - sqrt(3)/6
            a21 = 1/4 + sqrt(3)/6
            a22 = 1/4
            b1 = 1/2
            b2 = 1/2

            # In IRK4 the stage slopes are also decision variables.
            ks = JuMP.@variable(model, [1:V, 1:2])

            y1 = [xs[j, i] + tstep*(a11*ks[j,1] + a12*ks[j,2]) for j in 1:V]
            for j in 1:V
                JuMP.@constraint(model, ks[j,1] == dx[j](y1..., all_p_args_i..., t_at(i, c=c1)))
            end

            y2 = [xs[j, i] + tstep*(a21*ks[j,1] + a22*ks[j,2]) for j in 1:V]
            for j in 1:V
                JuMP.@constraint(model, ks[j,2] == dx[j](y2..., all_p_args_i..., t_at(i, c=c2)))
            end

            for j in 1:V
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + tstep*(b1*ks[j,1] + b2*ks[j,2]))
            end
        else
            error("Available integrators: EE, IE, RK4, IRK4 (Gauss-Legendre 2-stage)")
        end
    end
end

"""
    full_solutions(model, sys)
    full_solutions(model, sys, ps)

Read solved JuMP values back into MTK symbols.

This is a reporting helper. It substitutes solved values into the MTK model and
returns the observed MTK expressions.

How it connects to other functions:
- `decision_vars(...)` determines the symbol ordering used here;
- it is useful after `register_nlsystem(...)`, `register_odesystem(...)`, or
  `register_daesystem(...)`;
- it lets case-study scripts stay close to the original MTK notation.
"""
function full_solutions(model::JuMP.Model, sys::ModelingToolkit.System)
    return full_solutions(model, sys, Num[])
end

function full_solutions(model::JuMP.Model, sys::ModelingToolkit.System, ps::Vector{Num})
    # Rebuild the same symbol order used during registration.
    vars = decision_vars(sys, ps)
    sub_dict = ModelingToolkit.defaults(sys)
    for i in eachindex(vars)
        sub_dict[vars[i]] = JuMP.value.(JuMP.all_variables(model)[i])
    end
    for eqn in ModelingToolkit.observed(sys)
        sub_dict[eqn.lhs] = eqn.rhs
    end
    soln_dict = Dict()
    for i in eachindex(ModelingToolkit.observed(sys))
        rhsExpr = ModelingToolkit.observed(sys)[i].rhs
        # Keep substituting until the observed expression is fully reduced.
        while !isempty(intersect(Symbolics.get_variables(rhsExpr), keys(sub_dict)))
            rhsExpr = SymbolicUtils.substitute(rhsExpr, sub_dict)
        end
        soln_dict[ModelingToolkit.observed(sys)[i].lhs] = rhsExpr
    end
    return soln_dict
end

"""
    register_daesystem(model, sys, tspan, tstep, integrator; ...)

Register a ModelingToolkit DAE system as JuMP constraints over time.

Compared with `register_odesystem(...)`, this function also keeps the algebraic
equations in the model.

Assumed equation forms:
- `D(x) ~ f(...)`
- `0 ~ g(...)`
- `z ~ h(...)` (converted internally to a residual equation)

Supported integrators for the differential part:
- `"EE"`: Explicit Euler
- `"IE"`: Implicit Euler
- `"RK4"`: Classical explicit Runge-Kutta order 4
- `"IRK4"`: two-stage Gauss-Legendre implicit Runge-Kutta order 4

How it connects to other functions:
- `decision_vars(...; build_state_trajs=true)` provides the state trajectories
  consumed here;
- `_register_tracking_system!(...)` in `trackingmpc.jl` chooses this function
  when the user sets `system_kind=:dae`;
- use this for models like ASM2d where the algebraic equations must stay in the
  prediction model.

Algorithm:
1. Build the time grid and normalize all user-supplied names.
2. Split the MTK equations into two groups:
   differential equations for stepped states, and algebraic residual equations
   that must stay equal to zero.
3. Build callable functions for the differential right-hand sides and algebraic
   residuals.
4. Read or infer the JuMP state trajectories `x_vars`.
5. At each time step, add the chosen time-stepping equations for the
   differential states.
6. Also enforce the algebraic residuals, either at the nodes only or at the
   internal RK stages when the method needs them.
7. Keep both parts in the model so the prediction stays on the DAE manifold.
"""
function register_daesystem(model::JuMP.Model,
                            sys::ModelingToolkit.System,
                            tspan::Tuple{Real,Real},
                            tstep::Real,
                            integrator::String;
                            p_disc::Vector{Num}=Num[],
                            p_disc_vars::Dict{Num,Vector{JuMP.VariableRef}}=Dict{Num,Vector{JuMP.VariableRef}}(),
                            x_vars::AbstractDict=Dict(),
                            t_map::Dict{SymbolicUtils.BasicSymbolic{Real},JuMP.VariableRef}=Dict{SymbolicUtils.BasicSymbolic{Real},JuMP.VariableRef}(),
                            enforce_alg_at_t1::Bool=true)

    # -----------------------------
    # Time grid and dimensions
    # -----------------------------
    N = Int(floor((tspan[2] - tspan[1]) / tstep)) + 1
    # Normalize all user input names once at the start.
    p_disc = canonicalize_system_symbols(sys, p_disc; pool=:parameters)
    p_disc_vars = _canonical_parameter_var_dict(sys, p_disc_vars)
    x_vars = _canonical_state_var_dict(sys, x_vars)
    unknowns = ModelingToolkit.unknowns(sys)
    V = length(unknowns)

    intg = uppercase(integrator)
    @assert intg in ("EE", "IE", "RK4", "IRK4") "register_daesystem: integrator must be one of EE, IE, RK4, IRK4"

    # If no time variable is passed in, use numeric time.
    t_MTK = ModelingToolkit.get_iv(sys)
    t_base = haskey(t_map, t_MTK) ? t_map[t_MTK] : tspan[1]
    t_at(i::Int; c::Real=0.0) = t_base + (i-1 + c) * tstep

    # -----------------------------
    # Substitute default parameters, except states and discretized parameters.
    # -----------------------------
    param_dict = copy(ModelingToolkit.defaults(sys))
    for u in unknowns
        pop!(param_dict, u, nothing)
    end
    for p in p_disc
        pop!(param_dict, p, nothing)
    end
    sub_defaults(expr) = begin
        ex = expr
        while !isempty(intersect(Symbolics.get_variables(ex), keys(param_dict)))
            ex = SymbolicUtils.substitute(ex, param_dict)
        end
        ex
    end

    # -----------------------------
    # Split equations into differential and algebraic parts.
    # -----------------------------
    idx_of = Dict(u => i for (i,u) in enumerate(unknowns))
    full_eqs = ModelingToolkit.full_equations(sys)

    # differential: state index -> right-hand side
    f_expr = Dict{Int,Symbolics.Num}()

    # algebraic: residual expressions that must equal zero
    alg_residuals = Symbolics.Num[]

    is_diff_lhs(lhs) = begin
        try
            op = Symbolics.operation(lhs)
            op isa ModelingToolkit.Differential
        catch
            false
        end
    end

    diff_var(lhs) = Symbolics.arguments(lhs)[1]

    for eq in full_eqs
        lhs = eq.lhs
        rhs = eq.rhs

        if is_diff_lhs(lhs)
            # Differential equation: this is stepped forward in time.
            v = diff_var(lhs)
            @assert haskey(idx_of, v) "DAE parse: differential var $v not in unknowns(sys)"
            j = idx_of[v]
            f_expr[j] = sub_defaults(rhs)
        else
            # Algebraic equation: turn it into a residual equal to zero.
            push!(alg_residuals, sub_defaults(lhs - rhs))
        end
    end

    diff_idx = sort(collect(keys(f_expr)))
    diff_set = Set(diff_idx)
    alg_idx = [j for j in 1:V if !(j in diff_set)]
    @assert !isempty(diff_idx) "register_daesystem: no differential equations detected"
    # It is okay if no algebraic residuals remain after simplification.

    make_stage_state(diff_vals::AbstractVector, alg_vals::AbstractVector) = begin
        @assert length(diff_vals) == length(diff_idx)
        @assert length(alg_vals) == length(alg_idx)
        stage = Vector{Any}(undef, V)
        for (pos, j) in enumerate(diff_idx)
            stage[j] = diff_vals[pos]
        end
        for (pos, j) in enumerate(alg_idx)
            stage[j] = alg_vals[pos]
        end
        stage
    end

    # -----------------------------
    # Build callable functions for the differential and algebraic parts.
    # -----------------------------
    # Use the same input order as `register_odesystem(...)`.
    dv = decision_vars(sys, p_disc)

    f_fun = Dict{Int,Function}()
    for j in diff_idx
        f_fun[j] = build_function(f_expr[j], dv..., t_MTK; expression=Val{false})
    end

    g_fun = Function[]
    for r in alg_residuals
        # Build the algebraic residual functions once and reuse them.
        push!(g_fun, build_function(r, dv..., t_MTK; expression=Val{false}))
    end

    # -----------------------------
    # Check the discretized parameter arrays.
    # -----------------------------
    for p in p_disc
        @assert haskey(p_disc_vars, p) "Missing p_disc_vars[$p] for $p (should have length N)"
        @assert length(p_disc_vars[p]) == N "p_disc_vars[$p] must have length N=$N"
    end

    flat_pvars = reduce(vcat, values(p_disc_vars); init=JuMP.VariableRef[])

    # -----------------------------
    # Build the state matrix `xs` with size V x N.
    # -----------------------------
    xs = Array{JuMP.VariableRef}(undef, V, N)
    if !isempty(x_vars)
        for (j, u) in enumerate(unknowns)
            @assert haskey(x_vars, u) "Missing x_vars[$u]"
            @assert length(x_vars[u]) == N "x_vars[$u] must have length N=$N"
            xs[j, :] = x_vars[u]
        end
    else
        # Old fallback path. Passing `x_vars` is safer.
        xs = reshape(setdiff(JuMP.all_variables(model), flat_pvars), V, N)
        @warn "register_daesystem: inferring state ordering from all_variables; pass x_vars=... to avoid mismatches."
    end

    # -----------------------------
    # Step only the differential states.
    # -----------------------------
    for i in 1:(N-1)
        p_args_i = [p_disc_vars[p][i] for p in p_disc]

        if intg == "EE"
            for j in diff_idx
                # Only differential states are advanced here.
                JuMP.@constraint(model,
                    xs[j, i+1] == xs[j, i] + tstep * f_fun[j](xs[:, i]..., p_args_i..., t_at(i))
                )
            end

        elseif intg == "IE"
            for j in diff_idx
                JuMP.@constraint(model,
                    xs[j, i+1] == xs[j, i] + tstep * f_fun[j](xs[:, i+1]..., p_args_i..., t_at(i, c=1.0))
                )
            end

        elseif intg == "RK4"
            k1 = [f_fun[j](xs[:, i]..., p_args_i..., t_at(i)) for j in diff_idx]

            # Add temporary algebraic stage variables for the RK points.
            y2_alg = isempty(alg_idx) ? Any[] : Any[@variable(model) for _ in alg_idx]
            y2_diff = [xs[j, i] + 0.5 * tstep * k1[pos] for (pos, j) in enumerate(diff_idx)]
            y2 = make_stage_state(y2_diff, y2_alg)
            for gg in g_fun
                JuMP.@constraint(model, gg(y2..., p_args_i..., t_at(i, c=0.5)) == 0)
            end
            k2 = [f_fun[j](y2..., p_args_i..., t_at(i, c=0.5)) for j in diff_idx]

            y3_alg = isempty(alg_idx) ? Any[] : Any[@variable(model) for _ in alg_idx]
            y3_diff = [xs[j, i] + 0.5 * tstep * k2[pos] for (pos, j) in enumerate(diff_idx)]
            y3 = make_stage_state(y3_diff, y3_alg)
            for gg in g_fun
                JuMP.@constraint(model, gg(y3..., p_args_i..., t_at(i, c=0.5)) == 0)
            end
            k3 = [f_fun[j](y3..., p_args_i..., t_at(i, c=0.5)) for j in diff_idx]

            y4_alg = isempty(alg_idx) ? Any[] : Any[@variable(model) for _ in alg_idx]
            y4_diff = [xs[j, i] + tstep * k3[pos] for (pos, j) in enumerate(diff_idx)]
            y4 = make_stage_state(y4_diff, y4_alg)
            for gg in g_fun
                JuMP.@constraint(model, gg(y4..., p_args_i..., t_at(i, c=1.0)) == 0)
            end
            k4 = [f_fun[j](y4..., p_args_i..., t_at(i, c=1.0)) for j in diff_idx]

            for (pos, j) in enumerate(diff_idx)
                JuMP.@constraint(model,
                    xs[j, i+1] == xs[j, i] + (tstep / 6.0) * (k1[pos] + 2 * k2[pos] + 2 * k3[pos] + k4[pos])
                )
            end

        elseif intg == "IRK4"
            c1 = 0.5 - sqrt(3) / 6
            c2 = 0.5 + sqrt(3) / 6
            a11 = 1 / 4
            a12 = 1 / 4 - sqrt(3) / 6
            a21 = 1 / 4 + sqrt(3) / 6
            a22 = 1 / 4
            b1 = 1 / 2
            b2 = 1 / 2

            # The implicit stage slopes and algebraic stage variables are solved
            # together.
            ks = JuMP.@variable(model, [1:length(diff_idx), 1:2])
            y1_alg = isempty(alg_idx) ? Any[] : Any[@variable(model) for _ in alg_idx]
            y2_alg = isempty(alg_idx) ? Any[] : Any[@variable(model) for _ in alg_idx]

            y1_diff = [xs[j, i] + tstep * (a11 * ks[pos, 1] + a12 * ks[pos, 2]) for (pos, j) in enumerate(diff_idx)]
            y2_diff = [xs[j, i] + tstep * (a21 * ks[pos, 1] + a22 * ks[pos, 2]) for (pos, j) in enumerate(diff_idx)]
            y1 = make_stage_state(y1_diff, y1_alg)
            y2 = make_stage_state(y2_diff, y2_alg)

            for (pos, j) in enumerate(diff_idx)
                JuMP.@constraint(model, ks[pos, 1] == f_fun[j](y1..., p_args_i..., t_at(i, c=c1)))
                JuMP.@constraint(model, ks[pos, 2] == f_fun[j](y2..., p_args_i..., t_at(i, c=c2)))
            end

            for gg in g_fun
                JuMP.@constraint(model, gg(y1..., p_args_i..., t_at(i, c=c1)) == 0)
                JuMP.@constraint(model, gg(y2..., p_args_i..., t_at(i, c=c2)) == 0)
            end

            for (pos, j) in enumerate(diff_idx)
                JuMP.@constraint(model,
                    xs[j, i+1] == xs[j, i] + tstep * (b1 * ks[pos, 1] + b2 * ks[pos, 2])
                )
            end
        end
    end

    # -----------------------------
    # Enforce the algebraic equations at the time nodes.
    # -----------------------------
    if !isempty(g_fun)
        k_start = enforce_alg_at_t1 ? 1 : 2
        for k in k_start:N
            p_args_k = [p_disc_vars[p][k] for p in p_disc]
            for gg in g_fun
                # These constraints keep the prediction on the algebraic manifold.
                JuMP.@constraint(model, gg(xs[:, k]..., p_args_k..., t_at(k)) == 0)
            end
        end
    end

    return nothing
end
