using EOptInterface
using Ipopt
using JuMP
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SciCompDSL

@mtkmodel TankValve begin
    @parameters begin
        A = 1.0     # tank cross-section area  [m²]
        q_in = 2.0     # constant inlet flow      [m³/s]
        k_v                # valve coefficient (FREE – no default)
    end
    @variables begin
        h(t) = 1.0    # liquid level [m]  (ODE state, IC = 1)
        q_out(t), [irreducible = true]           # outlet flow  [m³/s] (algebraic, no IC)
    end
    @equations begin
        # ODE: tank mass balance
        D(h) ~ (q_in - q_out) / A
        # Algebraic constraint: valve equation  (written as 0 ~ rhs)
        q_out ~ k_v * sqrt(h)
    end
end

@mtkcompile sys = TankValve()

tspan = (0.0, 10.0)
tstep = 0.01

N = Int(floor((tspan[2] - tspan[1]) / tstep)) + 1
V = length(ModelingToolkit.unknowns(sys))
n_dvars = length(decision_vars(sys))
P = n_dvars - V                      # number of free parameters

jump_model = JuMP.Model(Ipopt.Optimizer)

JuMP.@variable(jump_model, 0.0 <= xs[1:V, 1:N] <= 100.0)
JuMP.@variable(jump_model, 0.01 <= ps[1:P] <= 10.0)

register_daesystem(jump_model, sys, tspan, tstep, "EE")

h_target = 4.0
JuMP.@objective(jump_model, Min, (xs[1, N] - h_target)^2)

JuMP.optimize!(jump_model)

println("═══ Results ═══")
println("Termination : ", JuMP.termination_status(jump_model))
println("Primal       : ", JuMP.primal_status(jump_model))
println("Solve time   : ", round(JuMP.solve_time(jump_model), digits=4), " s")
println("Objective    : ", round(JuMP.objective_value(jump_model), digits=6))

# Extract optimal k_v
kv_opt = JuMP.value(ps[1])
println("k_v*         : ", round(kv_opt, digits=4))
println("Expected k_v : 1.0  (analytical: q_in / sqrt(h_target) = 2/√4 = 1)")
