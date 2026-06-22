# Copyright (c) 2025 Joseph Choi, Dimitri Alston, Pengfei Xu, Matthew Stuber,
# and the University of Connecticut (UConn)
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EOptInterface
# An abstraction layer for optimizing equation-oriented/acausal models
# https://github.com/PSORLab/EOptInterface.jl
################################################################################
# ext/EOIInfiniteOptExt.jl
# An extension for using InfiniteOpt.
################################################################################

module EOIInfiniteOptExt

    # Imports
    import DocStringExtensions
    import EOptInterface
    import InfiniteOpt
    import JuMP
    import ModelingToolkit
    import SymbolicUtils
    import Symbolics

    # Functions
    """
        $(DocStringExtensions.TYPEDSIGNATURES)

    Automatically registers ODEs from a `ModelingToolkit.System` as 
    InfiniteOpt constraints.
    """
    function EOptInterface.register_odesystem(model::InfiniteOpt.InfiniteModel, sys::ModelingToolkit.System)
        t = InfiniteOpt.all_parameters(model, InfiniteOpt.InfiniteParameter)[1]
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
        xs = reshape(setdiff(JuMP.all_variables(model), ps), V)
        # Extract initial conditions from the ModelingToolkit system
        ic_dict = [ModelingToolkit.initial_conditions(sys)[i].val for i in intersect(ModelingToolkit.unknowns(sys), keys(ModelingToolkit.initial_conditions(sys)))]
        # Formulate InfiniteOpt constraints
        InfiniteOpt.@constraint(model, [i in 1:V], xs[i](0) == ic_dict[i])
        InfiniteOpt.@constraint(model, [i in 1:V], InfiniteOpt.∂(xs[i], t) == dx[i](xs..., ps...))
        return
    end

end