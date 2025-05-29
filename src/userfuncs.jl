"""
    decision_vars(sys)

Displays the decision variables for optimization problem of a ModelingToolkit model.
"""
function decision_vars(sys::System)
    x = [unknowns(sys); setdiff(ModelingToolkit.parameters(sys),keys(ModelingToolkit.defaults(sys)))]
    return x
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
function register_nlsystem(model::Model, sys::System, obj::Num, ineqs::Vector{Num})
    h = mtkns_modeleqs(sys)
    f = mtkns_usereqs(obj, sys)
    g = []
    for i in eachindex(ineqs)
        gi = mtkns_usereqs(ineqs[i], sys)
        push!(g, gi)
    end
    @constraint(model, [i in eachindex(h)], h[i](JuMP.all_variables(model)...) == 0)
    @constraint(model, [i in eachindex(g)], g[i](JuMP.all_variables(model)...) ≤ 0)
    @objective(model, Min, f(JuMP.all_variables(model)...))
end

"""
    register_nlsystem(model, sys, tspan, tstep, integrator)

Registers a ModelingToolkit dynamic model as algebraic constraints in JuMP by discretizing ODEs as a system of algebraic equations.

# Arguments
- `model::Model`: the JuMP model
- `sys::System`: the ModelingToolkit model
- `tspan::Type{Number,Number}`: the time span over which the dynamic model is simulated
- `tstep::Number`: the time step used in the integration scheme
- `integrator::String`: integration scheme used in discretization, `"EE"` for Explicit Euler or `"IE"` for Implicit Euler
"""
function register_odesystem(model::Model, odesys::System, tspan::Tuple{Number,Number}, tstep::Number, integrator::String)
    N = Int(floor((tspan[2] - tspan[1])/tstep))+1 # number of discrete time nodes
    V = length(unknowns(odesys)) # number of ode variables
    param_dict = copy(ModelingToolkit.defaults(odesys))
    for var in unknowns(odesys)
        pop!(param_dict, var)
    end
    dx = []
    for j in 1:V
        dxj_expr = full_equations(odesys)[j].rhs
        # Fully subsitute parameters with default values
        while ~isempty(intersect(get_variables(dxj_expr),keys(param_dict)))
            dxj_expr = substitute(dxj_expr, param_dict)
        end
        dxj = build_function(
            dxj_expr,
            unknowns(odesys)...,
            setdiff(ModelingToolkit.parameters(odesys),keys(ModelingToolkit.defaults(odesys)))..., 
            expression = Val{false}
        )
        push!(dx, dxj)
    end
    # extracting initial conditions from MTK ODESystem -> algebraic JuMP constraint for x[1:V,1]
    ps = JuMP.all_variables(model)[1:length(setdiff(decision_vars(odesys),unknowns(odesys)))]
    xs = reshape(setdiff(JuMP.all_variables(model),JuMP.all_variables(model)[1:length(setdiff(decision_vars(odesys),unknowns(odesys)))]), V, N)
    @constraint(model, xs[:,1] == [ModelingToolkit.defaults(odesys)[unknowns(odesys)[i]] for i in eachindex(unknowns(odesys))])
    # formulating JuMP constraints of ode discretizations
    for i in 1:(N-1)
        for j in 1:V
            if integrator == "EE"
                @constraint(model, xs[j,i+1] == xs[j,i] + tstep*dx[j](xs[:,i]...,ps...))
            elseif integrator == "IE"
                @constraint(model, xs[j,i+1] == xs[j,i] + tstep*dx[j](xs[:,i+1]...,ps...))
            else
                print("Available integrators: EE, IE")
                break
            end
        end
    end
end

"""
    full_solutions(model, sys)

Returns a dictionary of optimal solution values for the observed variables for a ModelingToolkit algebraic model.
"""
function full_solutions(model::Model, sys::System)
    vars = [unknowns(sys); setdiff(ModelingToolkit.parameters(sys),keys(ModelingToolkit.defaults(sys)))]
    sub_dict = ModelingToolkit.defaults(sys)
    for i in eachindex(vars)
        sub_dict[vars[i]] = JuMP.value.(JuMP.all_variables(model)[i])
    end
    for eqn in observed(sys)
        sub_dict[eqn.lhs] = eqn.rhs
    end
    soln_dict = Dict()
    for i in eachindex(observed(sys))
        rhsExpr = observed(sys)[i].rhs
        while ~isempty(intersect(get_variables(rhsExpr),keys(sub_dict)))
            rhsExpr = substitute(rhsExpr, sub_dict)
        end
        soln_dict[observed(sys)[i].lhs] = rhsExpr
    end
    return soln_dict
end