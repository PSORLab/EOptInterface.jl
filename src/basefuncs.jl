function mtk_generate_model_equations(sys::ModelingToolkit.System)
    # Create subsitution dictionary for parameters with assigned default values
    param_dict = copy(ModelingToolkit.defaults(sys))
    # Function holder
    funcs = []
    # Create functions for each model equation
    for i in eachindex(ModelingToolkit.unknowns(sys))
        # Full model equation 0 = f(x)
        expr = ModelingToolkit.full_equations(ModelingToolkit.expand_connections(sys))[i].rhs 
            - ModelingToolkit.full_equations(ModelingToolkit.expand_connections(sys))[i].lhs
        # Subsitute all parameter values
        while ~isempty(intersect(Symbolics.get_variables(expr),keys(param_dict)))
            expr = SymbolicUtils.substitute(expr, param_dict)
        end
        # Create function that has inputs: 
        # x = [state variables (unknowns); desgin variables (parameters with no assigned default values)]
        add_func = Symbolics.build_function(
            expr, 
            EOptInterface.decision_vars(sys)..., 
            expression = Val{false}
            )
        # Add function to function holder
        push!(funcs, add_func)
    end
    return funcs
end

function mtk_generate_reduced_expression(expr::Symbolics.Num, sys::ModelingToolkit.System)
    # Creates a dictionary of parameter and observed variable substitutions
    sub_dict = ModelingToolkit.defaults(sys)
    for eqn in ModelingToolkit.observed(sys)
        sub_dict[eqn.lhs] = eqn.rhs
    end
    # Makes substitutions until no substitutions can be made
    while ~isempty(intersect(Symbolics.get_variables(expr), keys(sub_dict)))
        expr = SymbolicUtils.substitute(expr, sub_dict)
    end
    return Symbolics.build_function(expr, EOptInterface.decision_vars(sys)..., expression = Val{false})
end