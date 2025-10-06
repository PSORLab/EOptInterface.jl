"""
    decision_vars(sys)

Displays the decision variables for optimization problem of a ModelingToolkit model.
"""
function decision_vars(sys::ModelingToolkit.System)
    return [
        ModelingToolkit.unknowns(sys); 
        setdiff(ModelingToolkit.parameters(sys),keys(ModelingToolkit.defaults(sys)))
        ]
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
function register_odesystem(model::JuMP.Model, odesys::ModelingToolkit.System, tspan::Tuple{Real,Real}, tstep::Real, integrator::String)
    N = Int(floor((tspan[2] - tspan[1])/tstep))+1 # number of discrete time nodes
    V = length(ModelingToolkit.unknowns(odesys)) # number of ode variables
    param_dict = copy(ModelingToolkit.defaults(odesys))
    for var in ModelingToolkit.unknowns(odesys)
        pop!(param_dict, var)
    end
    dx = []
    for j in 1:V
        dxj_expr = ModelingToolkit.full_equations(odesys)[j].rhs
        # Fully subsitute parameters with default values
        while ~isempty(intersect(Symbolics.get_variables(dxj_expr),keys(param_dict)))
            dxj_expr = SymbolicUtils.substitute(dxj_expr, param_dict)
        end
        dxj = build_function(
            dxj_expr,
            EOptInterface.decision_vars(odesys)..., 
            expression = Val{false}
            )
        push!(dx, dxj)
    end
    ps = JuMP.all_variables(model)[end-length(setdiff(EOptInterface.decision_vars(odesys),ModelingToolkit.unknowns(odesys)))+1:end]
    xs = reshape(setdiff(JuMP.all_variables(model),ps), V, N)
    # extracting initial conditions from MTK ODESystem -> algebraic JuMP constraint for x[1:V,1]
    JuMP.fix.(xs[:,1], [ModelingToolkit.defaults(odesys)[ModelingToolkit.unknowns(odesys)[i]] for i in eachindex(ModelingToolkit.unknowns(odesys))], force=true)
    # formulating JuMP constraints of ode discretizations
    if integrator == "EE"
        JuMP.@constraint(model, [j in 1:V,i in 1:(N-1)], xs[j,i+1] == xs[j,i] + tstep*dx[j](xs[:,i]...,ps...))
    elseif integrator == "IE"
        JuMP.@constraint(model, [j in 1:V,i in 1:(N-1)], xs[j,i+1] == xs[j,i] + tstep*dx[j](xs[:,i+1]...,ps...))
    else
        print("Available integrators: EE, IE")
    end
end

"""
    full_solutions(model, sys)

Returns a dictionary of optimal solution values for the observed variables for a ModelingToolkit algebraic model.
"""
function full_solutions(model::JuMP.Model, sys::ModelingToolkit.System)
    vars = decision_vars(sys)
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
        while ~isempty(intersect(Symbolics.get_variables(rhsExpr),keys(sub_dict)))
            rhsExpr = SymbolicUtils.substitute(rhsExpr, sub_dict)
        end
        soln_dict[ModelingToolkit.observed(sys)[i].lhs] = rhsExpr
    end
    return soln_dict
end