# Copyright (c) 2025 Joseph Choi, Dimitri Alston, Pengfei Xu, Matthew Stuber,
# and the University of Connecticut (UConn)
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EOptInterface
# An abstraction layer for optimizing equation-oriented/acausal models
# https://github.com/PSORLab/EOptInterface.jl
################################################################################
# src/userfuncs.jl
# Defines user functions for retrieving decision variables from ModelingToolkit
# systems, registering ModelingToolkit equations as constraints in JuMP models,
# and calculating full-space solutions from ModelingToolkit systems.
################################################################################

"""
    decision_vars(sys)

Returns the decision variables for an optimization problem formulated from a ModelingToolkit system.
"""
function decision_vars(sys::ModelingToolkit.System)
    return [
        ModelingToolkit.unknowns(sys); 
        setdiff(ModelingToolkit.parameters(sys),keys(ModelingToolkit.defaults(sys)))
        ]
end

"""
    register_nlsystem(model, sys, obj, ineqs)

Automatically formulates algebraic JuMP constraints and objective function from
an algebraic ModelingToolkit system and user-provided constraints and objective symbolic expressions.

# Arguments
- `model::JuMP.Model`: the JuMP model
- `sys::ModelingToolkit.System`: the ModelingToolkit system
- `obj::Symbolics.Num`: a symbolic expression of the objective function using the ModelingToolkit system variables
- `ineqs::Vector{Symbolics.Num}`: a vector of symbolic expressions of inequality constraints using the ModelingToolkit system variables
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

Automatically applies forward transcription and registers the discretized ODE ModelingToolkit system as algebraic JuMP constraints.

# Arguments
- `model::JuMP.Model`: the JuMP model
- `sys::ModelingToolkit.System`: the ModelingToolkit model
- `tspan::Tuple{Real,Real}`: the time span over which the dynamic model is simulated
- `tstep::Real`: the time step used in the integration scheme
- `integrator::String`: integration scheme used in discretization, `"EE"` for explicit Euler or `"IE"` for implicit Euler
"""
function register_odesystem(model::JuMP.Model, odesys::ModelingToolkit.System, tspan::Tuple{Real,Real}, tstep::Real, integrator::String)
    param_dict = copy(ModelingToolkit.defaults(odesys))
    N = Int(floor((tspan[2] - tspan[1])/tstep)) + 1  # Number of discrete time nodes
    V = length(ModelingToolkit.unknowns(odesys))  # Number of ode variables
    for var in ModelingToolkit.unknowns(odesys)
        pop!(param_dict, var)
    end
    dx = []
    for j in 1:V
        dxj_expr = ModelingToolkit.full_equations(odesys)[j].rhs
        # Fully substitute parameters with default values
        while ~isempty(intersect(Symbolics.get_variables(dxj_expr), keys(param_dict)))
            dxj_expr = SymbolicUtils.substitute(dxj_expr, param_dict)
        end
        dxj = Symbolics.build_function(
            dxj_expr,
            EOptInterface.decision_vars(odesys)..., 
            expression = Val{false}
            )
        push!(dx, dxj)
    end
    ps = JuMP.all_variables(model)[end-length(setdiff(EOptInterface.decision_vars(odesys),ModelingToolkit.unknowns(odesys)))+1:end]
    xs = reshape(setdiff(JuMP.all_variables(model),ps), V, N)
    JuMP.fix.(xs[:,1], [ModelingToolkit.defaults(odesys)[ModelingToolkit.unknowns(odesys)[i]] for i in eachindex(ModelingToolkit.unknowns(odesys))], force=true)
    # Extract initial conditions from the ModelingToolkit system and fix them in the JuMP model for x[1:V,1]
    # Formulate JuMP constraints based on chosen ODE discretization method
    if integrator == "EE"
        JuMP.@constraint(model, [j in 1:V, i in 1:(N-1)], xs[j,i+1] == xs[j,i] + tstep*dx[j](xs[:,i]..., ps...))
    elseif integrator == "IE"
        JuMP.@constraint(model, [j in 1:V, i in 1:(N-1)], xs[j,i+1] == xs[j,i] + tstep*dx[j](xs[:,i+1]..., ps...))
    else
        error("Available integrators: EE, IE")
    end
end

"""
    full_solution(model, sys)

Returns a dictionary of optimal solution values for the observed variables of an algebraic ModelingToolkit system if the JuMP model is solved.
"""
function full_solution(model::JuMP.Model, sys::ModelingToolkit.System)
    vars = EOptInterface.decision_vars(sys)
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
        while ~isempty(intersect(Symbolics.get_variables(rhsExpr), keys(sub_dict)))
            rhsExpr = SymbolicUtils.substitute(rhsExpr, sub_dict)
        end
        soln_dict[ModelingToolkit.observed(sys)[i].lhs] = rhsExpr
    end
    return soln_dict
end
