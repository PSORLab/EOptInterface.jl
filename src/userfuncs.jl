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
    $(DocStringExtensions.TYPEDSIGNATURES)

Returns the decision variables for an optimization problem from a 
`ModelingToolkit.System`.
"""
function decision_vars(sys::ModelingToolkit.System)
    return [
        ModelingToolkit.unknowns(sys); 
        setdiff(ModelingToolkit.parameters(sys), keys(ModelingToolkit.initial_conditions(sys)))
        ]
end

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Automatically formulates and adds user-provided `Symbolics.Num` objective 
function and `Vector{Symbolics.Num}` constraints from a `ModelingToolkit.System` 
to a `JuMP.Model`.
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
    return
end

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Automatically applies specified direct transcription method and registers the 
discretized ODE `ModelingToolkit.System` as algebraic JuMP constraints.
Current supports Explicit Euler "EE" and Implicit Euler "IE".
"""
function register_odesystem(model::JuMP.Model, sys::ModelingToolkit.System, tspan::Tuple{Real,Real}, tstep::Real, integrator::String)
    if integrator != "EE" && integrator != "IE"
        error("Available integrators: EE, IE")
    end
    # Number of discrete time nodes
    N = Int(floor((tspan[2] - tspan[1])/tstep)) + 1
    # Number of ODE variables
    V = length(ModelingToolkit.unknowns(sys))
    param_dict = copy(ModelingToolkit.initial_conditions(sys).dict)
    for var in ModelingToolkit.unknowns(sys)
        pop!(param_dict, var)
    end
    dx = []
    for j in 1:V
        dxj_expr = ModelingToolkit.full_equations(sys)[j].rhs
        dxj_expr = SymbolicUtils.substitute(dxj_expr, ModelingToolkit.bindings(sys))
        dxj_expr = SymbolicUtils.substitute(dxj_expr, ModelingToolkit.bindings(sys))
        # Fully substitute parameters with default values
        while ~isempty(intersect(Symbolics.get_variables(dxj_expr), keys(param_dict)))
            dxj_expr = SymbolicUtils.substitute(dxj_expr, param_dict)
        end
        dxj = Symbolics.build_function(
            dxj_expr,
            EOptInterface.decision_vars(sys)..., 
            expression = Val{false}
            )
        push!(dx, dxj)
    end
    ps = JuMP.all_variables(model)[end-length(setdiff(EOptInterface.decision_vars(sys), ModelingToolkit.unknowns(sys)))+1:end]
    xs = reshape(setdiff(JuMP.all_variables(model), ps), V, N)
    # Extract initial conditions from the ModelingToolkit system and fix them in the JuMP model for x[1:V,1]
    JuMP.fix.(xs[:,1], [ModelingToolkit.initial_conditions(sys)[ModelingToolkit.unknowns(sys)[i]].val for i in eachindex(ModelingToolkit.unknowns(sys))], force=true)
    # Formulate JuMP constraints based on chosen ODE discretization method
    if integrator == "EE"
        JuMP.@constraint(model, [j in 1:V, i in 1:(N-1)], xs[j,i+1] == xs[j,i] + tstep*dx[j](xs[:,i]..., ps...))
    elseif integrator == "IE"
        JuMP.@constraint(model, [j in 1:V, i in 1:(N-1)], xs[j,i+1] == xs[j,i] + tstep*dx[j](xs[:,i+1]..., ps...))
    end
    return
end

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Returns a dictionary of optimal solution values for the observed variables of a 
`ModelingToolkit.System` if the `JuMP.Model` is solved.
"""
function full_solution(model::JuMP.Model, sys::ModelingToolkit.System)
    vars = EOptInterface.decision_vars(sys)
    sub_dict = ModelingToolkit.initial_conditions(sys)
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
