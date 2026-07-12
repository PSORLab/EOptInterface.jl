"""
    mtk_generate_model_equations(sys)

Build callable Julia functions for each dynamic equation in a
`ModelingToolkit.System`.

Most users will not call this directly. It is useful when the raw residual
functions are needed for inspection or a custom experiment.
Its job is:

1. expand the MTK equations,
2. move each equation into residual form `rhs - lhs`,
3. substitute default parameter values,
4. build one executable function per equation.

The returned vector is useful when you want to inspect or reuse the raw model
equations outside the standard tracking MPC workflow.
"""
function mtk_generate_model_equations(sys::ModelingToolkit.System)
    param_dict = copy(ModelingToolkit.defaults(sys))
    h = []
    for i in eachindex(ModelingToolkit.unknowns(sys))
        # Convert the MTK equation into one residual expression.
        expr = (ModelingToolkit.full_equations(ModelingToolkit.expand_connections(sys))[i].rhs 
            - ModelingToolkit.full_equations(ModelingToolkit.expand_connections(sys))[i].lhs)
        # Replace symbolic defaults with numeric values so the generated function
        # is easier to evaluate.
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

"""
    mtk_generate_reduced_expression(expr, sys)

Substitute defaults and observed-variable definitions into one symbolic
expression, then build an executable Julia function.

Use this helper when you want to evaluate one derived quantity from a
`ModelingToolkit` system without manually expanding all observed equations
yourself.
"""
function mtk_generate_reduced_expression(expr::Symbolics.Num, sys::ModelingToolkit.System)
    sub_dict = ModelingToolkit.defaults(sys)
    for eqn in ModelingToolkit.observed(sys)
        # Observed equations act like named aliases. Add them to the substitution
        # dictionary so the final expression is fully reduced.
        sub_dict[eqn.lhs] = eqn.rhs
    end
    while ~isempty(intersect(Symbolics.get_variables(expr), keys(sub_dict)))
        expr = SymbolicUtils.substitute(expr, sub_dict)
    end
    return Symbolics.build_function(expr, EOptInterface.decision_vars(sys)..., expression = Val{false})

end
