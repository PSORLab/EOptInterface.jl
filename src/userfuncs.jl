# MTK: NonlinearSystem or ODESystem decision vars <<<
function decision_vars(sys::ODESystem)
    x = [unknowns(sys); setdiff(ModelingToolkit.parameters(sys),keys(ModelingToolkit.defaults(sys)))]
    return x
end

# full_solutions(sys::ODESystem) <<<
function full_solutions(sys::ODESystem)
    vars = [unknowns(sys); setdiff(ModelingToolkit.parameters(sys),keys(ModelingToolkit.defaults(sys)))]
    sub_dict = ModelingToolkit.defaults(sys)
    for i in eachindex(vars)
        sub_dict[vars[i]] = JuMP.value.(x[i])
    end
    for eqn in observed(sys)
        sub_dict[eqn.lhs] = eqn.rhs
    end
    soln_dict = Dict()
    for i in eachindex(observed(sys))
        rhsExpr = observed(s)[i].rhs
        while ~isempty(intersect(get_variables(rhsExpr),keys(sub_dict)))
            rhsExpr = substitute(rhsExpr, sub_dict)
        end
        soln_dict[observed(s)[i].lhs] = rhsExpr
    end
    return soln_dict
end

# MTK: NLSystem -> JuMP <<<
function register_nlsystem(model::Model, sys::ODESystem, obj::Num, ineqs::Vector{Num})
    h = mtkns_modeleqs(sys)
    f = mtkns_usereqs(obj, sys)
    g = []
    for i in eachindex(ineqs)
        gi = mtkns_usereqs(ineqs[i], sys)
        push!(g, gi)
    end
    @constraint(model, [i in eachindex(h)], h[i](x...) == 0)
    @constraint(model, [i in eachindex(g)], g[i](x...) ≤ 0)
    @objective(model, Min, f(x...))
end

# MTK: ODESystem -> JuMP <<<
function register_odesystem(model::Model, odesys::ODESystem, tspan::Tuple{Number,Number}, tstep::Number, solver::String)
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
    @constraint(model, x[1:V,1] == [ModelingToolkit.defaults(o)[unknowns(o)[i]] for i in eachindex(unknowns(o))])
    # formulating JuMP constraints of ode discretizations
    for i in 1:(N-1)
        for j in 1:V
            if solver == "EE"
                @constraint(model, x[j,i+1] == x[j,i] + tstep*dx[j](x[:,i]...,p...))
            elseif solver == "IE"
                @constraint(model, x[j,i+1] == x[j,i] + tstep*dx[j](x[:,i+1]...,p...))
            elseif solver == "RK4"
                k1 = dx[j](x[:,i]...,p...)
                k2 = dx[j](x[:,i].+tstep/2*k1...,p...)
                k3 = dx[j](x[:,i].+tstep/2*k2...,p...)
                k4 = dx[j](x[:,i].+tstep*k3...,p...)
                @constraint(model, x[j,i+1] == x[j,i] + tstep*(k1 + 2*k2 + 2*k3 + k4)/6)
            elseif solver == "IRK4"
                k1 = dx[j](x[:,i+1]...,p...)
                k2 = dx[j](x[:,i+1].+tstep/2*k1...,p...)
                k3 = dx[j](x[:,i+1].+tstep/2*k2...,p...)
                k4 = dx[j](x[:,i+1].+tstep*k3...,p...)
                @constraint(model, x[j,i+1] == x[j,i] + tstep*(k1 + 2*k2 + 2*k3 + k4)/6)
            else
                print("Available integrators: EE, IE, RK4, IRK4")
                break
            end
        end
    end
end