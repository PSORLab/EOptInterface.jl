# Copyright (c) 2025 Joseph Choi, Dimitri Alston, Pengfei Xu, Matthew Stuber,
# and the University of Connecticut (UConn)
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EOptInterface
# An abstraction layer for optimizing equation-oriented/acausal models
# https://github.com/PSORLab/EOptInterface.jl
################################################################################
# src/basefuncs.jl
# Defines base functions for generating symbolic equations from ModelingToolkit
# systems.
################################################################################

function mtk_generate_model_equations(sys::ModelingToolkit.System)
    param_dict = copy(ModelingToolkit.initial_conditions(sys))
    h = []
    for i in eachindex(ModelingToolkit.unknowns(sys))
        expr = (ModelingToolkit.full_equations(ModelingToolkit.expand_connections(sys))[i].rhs 
            - ModelingToolkit.full_equations(ModelingToolkit.expand_connections(sys))[i].lhs)
        while ~isempty(intersect(Symbolics.get_variables(expr), keys(param_dict)))
            expr = SymbolicUtils.substitute(expr, param_dict)
        end
        hi = Symbolics.build_function(
            expr, 
            EOptInterface.decision_vars(sys)..., 
            expression = Val{false}
            )
        push!(h, hi)
    end
    return h
end

function mtk_generate_reduced_expression(expr::Symbolics.Num, sys::ModelingToolkit.System)
    sub_dict = ModelingToolkit.initial_conditions(sys)
    for eqn in ModelingToolkit.observed(sys)
        sub_dict[eqn.lhs] = eqn.rhs
    end
    while ~isempty(intersect(Symbolics.get_variables(expr), keys(sub_dict)))
        expr = SymbolicUtils.substitute(expr, sub_dict)
    end
    return Symbolics.build_function(expr, EOptInterface.decision_vars(sys)..., expression = Val{false})

end
