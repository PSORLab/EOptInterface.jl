function mtk_generate_model_equations(sys::ModelingToolkit.System)
    param_dict = copy(ModelingToolkit.defaults(sys))
    h = []
    for i in eachindex(ModelingToolkit.unknowns(sys))
        expr = (ModelingToolkit.full_equations(ModelingToolkit.expand_connections(sys))[i].rhs 
            - ModelingToolkit.full_equations(ModelingToolkit.expand_connections(sys))[i].lhs)
        while ~isempty(intersect(Symbolics.get_variables(expr),keys(param_dict)))
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
    sub_dict = ModelingToolkit.defaults(sys)
    for eqn in ModelingToolkit.observed(sys)
        sub_dict[eqn.lhs] = eqn.rhs
    end
    while ~isempty(intersect(Symbolics.get_variables(expr), keys(sub_dict)))
        expr = SymbolicUtils.substitute(expr, sub_dict)
    end
    return Symbolics.build_function(expr, EOptInterface.decision_vars(sys)..., expression = Val{false})
end