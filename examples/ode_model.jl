using ModelingToolkit, JuMP, Ipopt, EOptInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

# Create ModelingToolkit ODE system
@mtkmodel KineticParameterEstimation begin
    @parameters begin
        # Known parameters
        T = 273.0
        K_2 = 46.0*exp(6500.0/T - 18.0)
        K_3 = 2.0*K_2
        k_1 = 53.0
        k_1s = k_1*1e-6
        k_5 = 1.2e-3
        c_O2 = 2e-3

        # Unknown parameters (free design variables)
        k_2f
        k_3f
        k_4
    end
    @variables begin
        # Initial conditions given for differential variables
        x_A(t) = 0.0
        x_B(t) = 0.0
        x_D(t) = 0.0
        x_Y(t) = 0.4
        x_Z(t) = 140.0

        I(t)
    end
    @equations begin
        D(x_A) ~ k_1*x_Z*x_Y - c_O2*(k_2f + k_3f)*x_A + k_2f/K_2*x_D + k_3f/K_3*x_B - k_5*x_A^2
        D(x_B) ~ c_O2*k_3f*x_A - (k_3f/K_3 + k_4)*x_B
        D(x_D) ~ c_O2*k_2f*x_A - k_2f/K_2*x_D
        D(x_Y) ~ -k_1s*x_Z*x_Y
        D(x_Z) ~ -k_1*x_Z*x_Y
        I ~ x_A + 2/21*x_B + 2/21*x_D
    end
end

# Compile system
@mtkcompile system = KineticParameterEstimation()

# Define integration time span
tspan = (0.0, 2.0)
tstep = 0.01
N = Int(floor((tspan[2] - tspan[1])/tstep)) + 1

# Include experimental intensity data
include("kinetic_intensity_data.jl")
# Define intensity function
intensity(x_A, x_B, x_D) = x_A + 2/21*x_B + 2/21*x_D

# Create JuMP model
model = Model(Ipopt.Optimizer)

# Retrieve decision variables from ModelingToolkit system
# Returns [x_Z(t), x_Y(t), x_D(t), x_B(t), x_A(t), k_2f, k_3f, k_4]
decision_vars(system)

# Create discretized state decision variables
# z = [x_Z(t), x_Y(t), x_D(t), x_B(t), x_A(t)]
V = length(unknowns(system))
@variable(model, -75.0 <= z[1:V,1:N] <= 150.0)

# Create free design decision variables
# p = [k_2f, k_3f, k_4]
pL = [10.0, 10.0, 0.001]
pU = [1200.0, 1200.0, 40.0]
@variable(model, pL[i] <= p[i=1:3] <= pU[i])

# Register ModelingToolkit ODE system as constraints
register_odesystem(model, system, tspan, tstep, "IE")

# Define objective function
@objective(model, Min, sum((intensity(z[5,i], z[4,i], z[3,i]) - data[i-1])^2 for i in 2:N))

# Optimize model and retrieve results
JuMP.optimize!(model)
println("Termination Status: $(JuMP.termination_status(model))")
println("Primal Status: $(JuMP.primal_status(model))")
println("Solve Time: $(round.(JuMP.solve_time(model), digits=5))")
println("f^* = $(round(JuMP.objective_value(model), digits=5))")
println("p* = $(round.(JuMP.value.(p), digits=3))")
