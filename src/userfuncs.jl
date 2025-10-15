"""
    decision_vars(sys, ps)

Displays the decision variables for optimization problem of a ModelingToolkit model.
`ps` specifies which MTK parameters (e.g., u_in_valve(t)) will be treated as decision variables (and can be discretized over time).
"""
function decision_vars(sys::ModelingToolkit.System, ps::Vector{Num}=Num[])
    # Parameters to optimize = user-specified ps ∪ parameters without default values
    params_to_opt = union(ps, setdiff(ModelingToolkit.parameters(sys),
                                      keys(ModelingToolkit.defaults(sys))))
    return vcat(ModelingToolkit.unknowns(sys), collect(params_to_opt))
end

"""
    register_nlsystem(model, sys, obj, ineqs)

Registers a ModelingToolkit algebraic model, objective function, and inequality constraints as algebraic constraints in JuMP.

# Arguments
- `model::Model`: the JuMP model
- `sys::System`: the ModelingToolkit model
- `obj::Num`: a symbolic expression of the objective function using the ModelingToolkit model variables
- `ineqs::Vector{Num}`: a vector of symbolic expressions of inequality constraints using the ModelingToolkit model variables
"""
function register_nlsystem(model::JuMP.Model, sys::ModelingToolkit.System, obj::Symbolics.Num, ineqs::Vector{Symbolics.Num})
    h = EOptInterface.mtk_generate_model_equations(sys)
    f = EOptInterface.mtk_generate_reduced_expression(obj, sys)
    g = []
    for i in eachindex(ineqs)
        gi = EOptInterface.mtk_generate_reduced_expression(ineqs[i], sys)
        push!(g, gi)
    end
    JuMP.@constraint(model, [i in eachindex(h)], h[i](JuMP.all_variables(model)...) == 0)
    JuMP.@constraint(model, [i in eachindex(g)], g[i](JuMP.all_variables(model)...) ≤ 0)
    JuMP.@objective(model, Min, f(JuMP.all_variables(model)...))
end

"""
    register_odesystem(model, sys, tspan, tstep, integrator)

Registers a ModelingToolkit dynamic model as algebraic constraints in JuMP by discretizing ODEs as a system of algebraic equations.

- Supported integrators:
  - "EE"       : Explicit Euler
  - "IE"|"BDF1": Implicit Euler
  - "RK4"      : Classical 4th-order Runge–Kutta
  - "IRK4"     : 2-stage Gauss–Legendre implicit RK (order 4)
  - "RADAU"|"RADAU5"|"RADAUIIA": 3-stage Radau IIA implicit RK (order 5)
  - "BDF"|"BDF2": BDF2 with IE startup
  - "RK5"|"DP5"|"DOPRI5": Dormand–Prince 5th-order explicit RK
  - "ROS2"|"ROS"|"Rosenbrock": 2-stage Rosenbrock–W method (linearly implicit, order 2)
  - "TSIT5"    : Tsitouras 5(4) explicit RK (7-stage, FSAL; step uses 6 stages)
  - "TRBDF2"   : 2-stage SDIRK approximation of TR-BDF2 (L-stable, order 2)
# Arguments
- `model::Model`: the JuMP model
- `sys::System`: the ModelingToolkit model
- `tspan::Type{Number,Number}`: the time span over which the dynamic model is simulated
- `tstep::Number`: the time step used in the integration scheme
- `integrator::String`: integration scheme used in discretization (see list above)
"""
function register_odesystem(model::JuMP.Model,
                            odesys::ModelingToolkit.System,
                            tspan::Tuple{Real,Real},
                            tstep::Real,
                            integrator::String;
                            p_disc::Vector{Num}=Num[],
                            p_disc_vars::Dict{Num,Vector{JuMP.VariableRef}}=Dict(),
                            x_vars::Dict{Num,Vector{JuMP.VariableRef}}=Dict(),
                            fallback_map::Dict{Num,Real}=Dict(),
                            pure_time_inputs::Dict{Num,Function}=Dict(),  # ← 新增
                            t0::Real = nothing)

    # Time grid and dimensions
    N = Int(floor((tspan[2] - tspan[1]) / tstep)) + 1
    times = (t0 .+ (0:N-1) .* tstep)               # ← 新增
    pure_syms = collect(keys(pure_time_inputs))    # ← 新增
    V = length(ModelingToolkit.unknowns(odesys))

    # Normalize integrator key
    intg = uppercase(integrator)

    # 1) Generate RHS functions f_j(x, p) and (for Rosenbrock) Jacobian J(x, p)
    param_dict = copy(ModelingToolkit.defaults(odesys))
    for var in ModelingToolkit.unknowns(odesys)
        pop!(param_dict, var, nothing)
    end
    for p in p_disc
        pop!(param_dict, p, nothing)
    end
    for s in pure_syms
        pop!(param_dict, s, nothing)
    end

    dx = Vector{Function}(undef, V)
    dx_exprs = Vector{Symbolics.Num}(undef, V)
    full_eqs = ModelingToolkit.full_equations(odesys)
    for j in 1:V
        dxj_expr = full_eqs[j].rhs
        while !isempty(intersect(Symbolics.get_variables(dxj_expr), keys(param_dict)))
            dxj_expr = SymbolicUtils.substitute(dxj_expr, param_dict)
        end
        dx_exprs[j] = dxj_expr
        dx[j] = build_function(dxj_expr,
                            decision_vars(odesys, vcat(p_disc, pure_syms))...;
                            expression = Val{false})
    end

    # Jacobian function for Rosenbrock methods (evaluated at x_i)
    Jfun = nothing
    if intg in ("ROS2", "ROS", "ROSENBROCK")
        J_expr = Symbolics.jacobian(dx_exprs, ModelingToolkit.unknowns(odesys))
        Jfun = build_function(J_expr, decision_vars(odesys, p_disc)...; expression = Val{false})
    end

    # 2) Validate discretized parameter arrays
    for p in p_disc
        @assert haskey(p_disc_vars, p) "Missing p_disc_vars[$p] for $p (should have length N)"
        @assert length(p_disc_vars[p]) == N "p_disc_vars[$p] must have length N=$N"
    end
    flat_pvars = reduce(vcat, values(p_disc_vars); init=JuMP.VariableRef[])

    # 3) Build the state matrix xs
    unknowns = ModelingToolkit.unknowns(odesys)
    xs = Array{JuMP.VariableRef}(undef, V, N)
    if !isempty(x_vars)
        for (j, u) in enumerate(unknowns)
            @assert haskey(x_vars, u) "Missing x_vars[$u]"
            @assert length(x_vars[u]) == N "x_vars[$u] must have length N=$N"
            xs[j, :] = x_vars[u]
        end
    else
        xs = reshape(setdiff(JuMP.all_variables(model), flat_pvars), V, N)
        @warn "register_odesystem: inferring state ordering from all_variables; pass x_vars=... to avoid mismatches."
    end

    # 4) Initial-condition constraints
    dflt = ModelingToolkit.defaults(odesys)

    ivals = [ get(dflt, u) do
                get(fallback_map, u) do
                    error("No initial value for $u")
                end
            end
            for u in unknowns ]
    for j in 1:V
        # println("Initial condition for $(unknowns[j]): $(ivals[j])")
        JuMP.@constraint(model, xs[j, 1] == ivals[j])
    end
    println(V)
    println(xs[102, :])
    # 5) ODE discretization constraints
    for i in 1:(N-1)
        p_args_i   = [p_disc_vars[p][i] for p in p_disc]
        p_pure_i   = [pure_time_inputs[s](times[i]) for s in pure_syms]   # ← 新增

        if intg == "EE"
            for j in 1:V
                println(j)
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + tstep * dx[j](xs[:, i]..., p_args_i..., p_pure_i...))
            end

        elseif intg == "IE" || intg == "BDF1"
            for j in 1:V
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + tstep * dx[j](xs[:, i+1]..., p_args_i..., p_pure_i...))
            end

        elseif intg == "RK4"
            k1 = [dx[j](xs[:, i]..., p_args_i..., p_pure_i...) for j in 1:V]
            xk1_half = [xs[j, i] + 0.5 * tstep * k1[j] for j in 1:V]

            k2 = [dx[j](xk1_half..., p_args_i..., p_pure_i...) for j in 1:V]
            xk2_half = [xs[j, i] + 0.5 * tstep * k2[j] for j in 1:V]

            k3 = [dx[j](xk2_half..., p_args_i..., p_pure_i...) for j in 1:V]
            xk3_full = [xs[j, i] + tstep * k3[j] for j in 1:V]

            k4 = [dx[j](xk3_full..., p_args_i..., p_pure_i...) for j in 1:V]

            for j in 1:V
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + (tstep/6) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]))
            end

        elseif intg == "RK5" || intg == "DP5" || intg == "DOPRI5"
            # Dormand–Prince 5(4) explicit RK (7 stages), params constant during step
            # Butcher tableau (c, A, b) for the 5th-order solution
            # c
            c2 = 1//5
            c3 = 3//10
            c4 = 4//5
            c5 = 8//9
            c6 = 1//1
            c7 = 1//1
            # A (lower triangular)
            a21 = 1//5

            a31 = 3//40;   a32  = 9//40

            a41 = 44//45;  a42  = -56//15;  a43  = 32//9

            a51 = 19372//6561; a52 = -25360//2187; a53 = 64448//6561; a54 = -212//729

            a61 = 9017//3168;  a62 = -355//33;    a63 = 46732//5247; a64 = 49//176; a65 = -5103//18656

            a71 = 35//384;     a72 = 0//1;        a73 = 500//1113;   a74 = 125//192; a75 = -2187//6784; a76 = 11//84
            # b (5th order)
            b1 = 35//384; b2 = 0//1; b3 = 500//1113; b4 = 125//192; b5 = -2187//6784; b6 = 11//84; b7 = 0//1

            k1 = [dx[j](xs[:, i]..., p_args_i..., p_pure_i...) for j in 1:V]

            y2 = [xs[j, i] + tstep*(a21*k1[j]) for j in 1:V]
            k2 = [dx[j](y2..., p_args_i..., p_pure_i...) for j in 1:V]

            y3 = [xs[j, i] + tstep*(a31*k1[j] + a32*k2[j]) for j in 1:V]
            k3 = [dx[j](y3..., p_args_i..., p_pure_i...) for j in 1:V]

            y4 = [xs[j, i] + tstep*(a41*k1[j] + a42*k2[j] + a43*k3[j]) for j in 1:V]
            k4 = [dx[j](y4..., p_args_i..., p_pure_i...) for j in 1:V]

            y5 = [xs[j, i] + tstep*(a51*k1[j] + a52*k2[j] + a53*k3[j] + a54*k4[j]) for j in 1:V]
            k5 = [dx[j](y5..., p_args_i..., p_pure_i...) for j in 1:V]

            y6 = [xs[j, i] + tstep*(a61*k1[j] + a62*k2[j] + a63*k3[j] + a64*k4[j] + a65*k5[j]) for j in 1:V]
            k6 = [dx[j](y6..., p_args_i..., p_pure_i...) for j in 1:V]

            y7 = [xs[j, i] + tstep*(a71*k1[j] + a72*k2[j] + a73*k3[j] + a74*k4[j] + a75*k5[j] + a76*k6[j]) for j in 1:V]
            k7 = [dx[j](y7..., p_args_i..., p_pure_i...) for j in 1:V]

            for j in 1:V
                JuMP.@constraint(model,
                    xs[j, i+1] == xs[j, i] + tstep*(b1*k1[j] + b2*k2[j] + b3*k3[j] + b4*k4[j] + b5*k5[j] + b6*k6[j] + b7*k7[j])
                )
            end

        elseif intg == "TSIT5"
            # Tsitouras 5(4) explicit RK (7 stages, FSAL). Step uses 6 stages.
            k1 = [dx[j](xs[:, i]..., p_args_i..., p_pure_i...) for j in 1:V]

            y2 = [xs[j, i] + tstep*(0.1610000000000000*k1[j]) for j in 1:V]
            k2 = [dx[j](y2..., p_args_i..., p_pure_i...) for j in 1:V]

            y3 = [xs[j, i] + tstep*(-0.008480655492356989*k1[j] + 0.3354806554923570*k2[j]) for j in 1:V]
            k3 = [dx[j](y3..., p_args_i..., p_pure_i...) for j in 1:V]

            y4 = [xs[j, i] + tstep*(2.8971530571054935*k1[j] - 6.3594484899750750*k2[j] + 4.3622954328695815*k3[j]) for j in 1:V]
            k4 = [dx[j](y4..., p_args_i..., p_pure_i...) for j in 1:V]

            y5 = [xs[j, i] + tstep*(5.3258648284392570*k1[j] -11.7488835640628280*k2[j] + 7.4955393428898365*k3[j] -0.09249506636175525*k4[j]) for j in 1:V]
            k5 = [dx[j](y5..., p_args_i..., p_pure_i...) for j in 1:V]

            y6 = [xs[j, i] + tstep*(5.8614554429464200*k1[j] -12.9209693178471100*k2[j] + 8.1593678985761590*k3[j] -0.07158497328140100*k4[j] -0.02826905039406838*k5[j]) for j in 1:V]
            k6 = [dx[j](y6..., p_args_i..., p_pure_i...) for j in 1:V]

            # Update uses a7 row (FSAL): y_{i+1} = y_i + h*sum(a7j * kj, j=1..6)
            for j in 1:V
                JuMP.@constraint(model,
                    xs[j, i+1] == xs[j, i] + tstep*(0.09646076681806523*k1[j] + 0.01000000000000000*k2[j] +
                                                    0.47988965041449960*k3[j] + 1.3790085741037420*k4[j] +
                                                   -3.2900695154360810*k5[j] + 2.3247105240997740*k6[j])
                )
            end

        elseif intg == "IRK4"
            # Implicit RK order 4 (2-stage Gauss–Legendre)
            c1 = 0.5 - sqrt(3)/6
            c2 = 0.5 + sqrt(3)/6
            a11 = 1/4
            a12 = 1/4 - sqrt(3)/6
            a21 = 1/4 + sqrt(3)/6
            a22 = 1/4
            b1 = 1/2
            b2 = 1/2

            ks = JuMP.@variable(model, [1:V, 1:2])

            y1 = [xs[j, i] + tstep*(a11*ks[j,1] + a12*ks[j,2]) for j in 1:V]
            for j in 1:V
                JuMP.@constraint(model, ks[j,1] == dx[j](y1..., p_args_i..., p_pure_i...))
            end

            y2 = [xs[j, i] + tstep*(a21*ks[j,1] + a22*ks[j,2]) for j in 1:V]
            for j in 1:V
                JuMP.@constraint(model, ks[j,2] == dx[j](y2..., p_args_i..., p_pure_i...))
            end

            for j in 1:V
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + tstep*(b1*ks[j,1] + b2*ks[j,2]))
            end

        elseif intg == "RADAU" || intg == "RADAU5" || intg == "RADAUIIA"
            # 3-stage Radau IIA (order 5) implicit RK
            s6 = sqrt(6.0)
            # Nodes (not explicitly used here)
            _c1 = (4 - s6)/10
            _c2 = (4 + s6)/10
            _c3 = 1.0
            # A matrix
            a11 = 11/45 - 7*s6/360
            a12 = 37/225 - 169*s6/1800
            a13 = -2/225 + s6/75

            a21 = 37/225 + 169*s6/1800
            a22 = 11/45 + 7*s6/360
            a23 = -2/225 - s6/75

            a31 = 4/9 - s6/36
            a32 = 4/9 + s6/36
            a33 = 1/9
            # b vector
            b1 = a31
            b2 = a32
            b3 = a33

            ks = JuMP.@variable(model, [1:V, 1:3])

            y1 = [xs[j, i] + tstep*(a11*ks[j,1] + a12*ks[j,2] + a13*ks[j,3]) for j in 1:V]
            for j in 1:V
                JuMP.@constraint(model, ks[j,1] == dx[j](y1..., p_args_i..., p_pure_i...))
            end

            y2 = [xs[j, i] + tstep*(a21*ks[j,1] + a22*ks[j,2] + a23*ks[j,3]) for j in 1:V]
            for j in 1:V
                JuMP.@constraint(model, ks[j,2] == dx[j](y2..., p_args_i..., p_pure_i...))
            end

            y3 = [xs[j, i] + tstep*(a31*ks[j,1] + a32*ks[j,2] + a33*ks[j,3]) for j in 1:V]
            for j in 1:V
                JuMP.@constraint(model, ks[j,3] == dx[j](y3..., p_args_i..., p_pure_i...))
            end

            for j in 1:V
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + tstep*(b1*ks[j,1] + b2*ks[j,2] + b3*ks[j,3]))
            end

        elseif intg == "TRBDF2"
            # 2-stage SDIRK approximation (L-stable, order 2)
            γ = 2 - sqrt(2)  # ~0.585786
            a11 = γ
            a21 = 1 - γ
            a22 = γ
            b1 = 1 - γ
            b2 = γ

            ks = JuMP.@variable(model, [1:V, 1:2])

            y1 = [xs[j, i] + tstep*(a11*ks[j,1]) for j in 1:V]
            for j in 1:V
                JuMP.@constraint(model, ks[j,1] == dx[j](y1..., p_args_i..., p_pure_i...))
            end

            y2 = [xs[j, i] + tstep*(a21*ks[j,1] + a22*ks[j,2]) for j in 1:V]
            for j in 1:V
                JuMP.@constraint(model, ks[j,2] == dx[j](y2..., p_args_i..., p_pure_i...))
            end

            for j in 1:V
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + tstep*(b1*ks[j,1] + b2*ks[j,2]))
            end

        elseif intg == "BDF" || intg == "BDF2"
            if i == 1
                for j in 1:V
                    JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + tstep * dx[j](xs[:, i+1]..., p_args_i..., p_pure_i...))
                end
            else
                for j in 1:V
                    JuMP.@constraint(model,
                        xs[j, i+1] == (4/3) * xs[j, i] - (1/3) * xs[j, i-1] + (2/3) * tstep * dx[j](xs[:, i+1]..., p_args_i..., p_pure_i...)
                    )
                end
            end

        elseif intg == "ROS2" || intg == "ROS" || intg == "ROSENBROCK"
            # Rosenbrock–W method (order 2, 2 stages), Jacobian frozen at x_i
            @assert Jfun !== nothing "Jacobian function not built"
            Ji = Jfun(xs[:, i]..., p_args_i...)  # V×V matrix (expressions)

            γ = 1 - 1/sqrt(2)
            a21 = 1.0
            c21 = -2.0  # W-method (Hairer/Wanner), J frozen at x_i
            m1 = 0.5
            m2 = 0.5

            ks = JuMP.@variable(model, [1:V, 1:2])

            # Stage 1: (I - γ h J) k1 = h f(x_i)
            f1 = [dx[j](xs[:, i]..., p_args_i..., p_pure_i...) for j in 1:V]
            for j in 1:V
                JuMP.@constraint(model,
                    ks[j,1] - γ*tstep*sum(Ji[j,m] * ks[m,1] for m in 1:V) == tstep * f1[j]
                )
            end

            # Stage 2: (I - γ h J) k2 = h f(x_i + a21*k1) + h J (c21*k1)
            y2 = [xs[j, i] + a21*ks[j,1] for j in 1:V]
            f2 = [dx[j](y2..., p_args_i..., p_pure_i...) for j in 1:V]
            for j in 1:V
                JuMP.@constraint(model,
                    ks[j,2] - γ*tstep*sum(Ji[j,m] * ks[m,2] for m in 1:V) ==
                        tstep * f2[j] + tstep * sum(Ji[j,m] * (c21*ks[m,1]) for m in 1:V)
                )
            end

            # Update
            for j in 1:V
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + m1*ks[j,1] + m2*ks[j,2])
            end
        else
            error("Available integrators: EE, IE, RK4, IRK4 (Gauss–Legendre 2-stage), RADAU (Radau IIA 3-stage), BDF (BDF2 with IE startup), RK5 (Dormand–Prince), TSIT5 (Tsitouras 5), TRBDF2 (SDIRK2), ROS2/ROS/ROSENBROCK (Rosenbrock–W order 2)")
        end
    end
end

"""
    full_solutions(model, sys)

Returns a dictionary of optimal solution values for the observed variables for a ModelingToolkit algebraic model.
"""
function full_solutions(model::JuMP.Model, sys::ModelingToolkit.System; ps::Vector{Num}=Num[])
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
        while !isempty(intersect(Symbolics.get_variables(rhsExpr), keys(sub_dict)))
            rhsExpr = SymbolicUtils.substitute(rhsExpr, sub_dict)
        end
        soln_dict[ModelingToolkit.observed(sys)[i].lhs] = rhsExpr
    end
    return soln_dict
end
