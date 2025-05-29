function mtkns_modeleqs(sys::System)
    # Create subsitution dictionary for parameters with assigned default values
    param_dict = copy(ModelingToolkit.defaults(sys))
    # Function holder
    func = []
    # Create functions for each model equation
    for i in eachindex(unknowns(sys))
        # Full model equation 0 = f(x)
        expr = full_equations(expand_connections(sys))[i].rhs - full_equations(expand_connections(sys))[i].lhs
        # Subsitute all parameter values
        while ~isempty(intersect(get_variables(expr),keys(param_dict)))
            expr = substitute(expr, param_dict)
        end
        # Create function that has inputs: 
        # x = [state variables (unknowns); desgin variables (parameters with no assigned default values)]
        add_func = build_function(
            expr, 
            unknowns(sys)..., 
            setdiff(ModelingToolkit.parameters(sys),keys(ModelingToolkit.defaults(sys)))..., 
            expression = Val{false}
        )
        # Add function to function holder
        push!(func, add_func)
    end
    return func
end

function mtkns_usereqs(expr::Num, sys::System)
    # Creates a dictionary of parameter and observed variable substitutions
    sub_dict = ModelingToolkit.defaults(sys)
    for eqn in observed(sys)
        sub_dict[eqn.lhs] = eqn.rhs
    end
    # Makes substitutions until no substitutions can be made
    while ~isempty(intersect(string.(get_variables(expr)), string.(keys(sub_dict))))
        expr = substitute(expr, sub_dict)
    end
    # Function inputs (design and state variables)
    x = [unknowns(sys); setdiff(ModelingToolkit.parameters(sys),keys(ModelingToolkit.defaults(sys)))]
    return build_function(expr, x..., expression = Val{false})
end