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

# Arguments
- `model::Model`: the JuMP model
- `sys::System`: the ModelingToolkit model
- `tspan::Type{Number,Number}`: the time span over which the dynamic model is simulated
- `tstep::Number`: the time step used in the integration scheme
- `integrator::String`: integration scheme used in discretization, `"EE"` for Explicit Euler or `"IE"` for Implicit Euler
"""
function register_odesystem(model::JuMP.Model,
                            odesys::ModelingToolkit.System,
                            tspan::Tuple{Real,Real},
                            tstep::Real,
                            integrator::String;
                            p_disc::Vector{Num}=Num[],
                            p_disc_vars::Dict{Num,Vector{JuMP.VariableRef}}=Dict())

    # Time grid and dimensions
    N = Int(floor((tspan[2] - tspan[1]) / tstep)) + 1
    V = length(ModelingToolkit.unknowns(odesys))

    # 1) Generate the RHS function dx[j](unknowns..., p_disc...) for each unknown
    #    - Use defaults to replace parameters not selected for discretization with constants
    #    - unknowns and p_disc remain symbolic, as arguments to dx
    param_dict = copy(ModelingToolkit.defaults(odesys))
    for var in ModelingToolkit.unknowns(odesys)
        pop!(param_dict, var, nothing)
    end
    for p in p_disc
        pop!(param_dict, p, nothing)
    end

    dx = Vector{Function}(undef, V)
    full_eqs = ModelingToolkit.full_equations(odesys)
    for j in 1:V
        dxj_expr = full_eqs[j].rhs
        while !isempty(intersect(Symbolics.get_variables(dxj_expr), keys(param_dict)))
            dxj_expr = SymbolicUtils.substitute(dxj_expr, param_dict)
        end
        dx[j] = build_function(dxj_expr,
                               decision_vars(odesys, p_disc)...;
                               expression = Val{false})
    end

    # 2) Extract the JuMP variable trajectories for discretized parameters
    #    The caller should have already created these: e.g., @variable(model, u_in_valve_k[1:N] ∈ [0,1])
    for p in p_disc
        @assert haskey(p_disc_vars, p) "Missing p_disc_vars[$p] for $p (should have length N)"
        @assert length(p_disc_vars[p]) == N "p_disc_vars[$p] must have length N=$N"
    end
    flat_pvars = reduce(vcat, values(p_disc_vars); init=JuMP.VariableRef[])

    # 3) Extract the state variable matrix xs over time from the model
    #    Method: exclude discretized parameter variables from all_variables, leaving V×N unknowns trajectories
    xs = reshape(setdiff(JuMP.all_variables(model), flat_pvars), V, N)

    # 4) Initial value constraints: xs[:,1] equals the default values of unknowns in MTK defaults
    JuMP.@constraint(model, xs[:, 1] .== [
        ModelingToolkit.defaults(odesys)[ModelingToolkit.unknowns(odesys)[i]]
        for i in eachindex(ModelingToolkit.unknowns(odesys))
    ])

    # 5) ODE discretization constraints (EE/IE), at each time step assemble parameter arguments:
    #    For step i, parameter arguments = [ p_disc_vars[p][i] for p in p_disc ]
    for i in 1:(N-1)
        p_args_i = [p_disc_vars[p][i] for p in p_disc]
        if integrator == "EE"
            for j in 1:V
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + tstep * dx[j](xs[:, i]..., p_args_i...))
            end
        elseif integrator == "IE"
            for j in 1:V
                JuMP.@constraint(model, xs[j, i+1] == xs[j, i] + tstep * dx[j](xs[:, i+1]..., p_args_i...))
            end
        else
            error("Available integrators: EE, IE")
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
