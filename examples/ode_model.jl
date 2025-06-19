using ModelingToolkit, JuMP, EOptInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel KineticParameterEstimation begin
    @parameters begin
        T = 273
        K_2 = 46*exp(6500/T-18)
        K_3 = 2*K_2
        k_1 = 53
        k_1s = k_1*1e-6
        k_5 = 1.2e-3
        c_O2 = 2e-3

        k_2f
        k_3f
        k_4
    end
    @variables begin
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

@mtkcompile o = KineticParameterEstimation()

tspan = (0.0,2.0)
tstep = 0.01
include("kinetic_intensity_data.jl")
intensity(x_A,x_B,x_D) = x_A + 2/21*x_B + 2/21*x_D

using EAGO
model = Model(EAGO.Optimizer)
decision_vars(o)
pL = [0.001, 10, 10]
pU = [40, 1200, 1200]
@variable(model, pL[i] <= p[i=1:3] <= pU[i])
N = Int(floor((tspan[2] - tspan[1])/tstep))+1
V = length(unknowns(o))
@variable(model, -75 <= x[1:V,1:N] <= 150.0 )
register_odesystem(model, o, tspan, tstep, "EE")
@objective(model, Min, sum((intensity(x[1,i+1],x[2,i+1],x[3,i+1]) - data[i])^2 for i in 1:(N-1)))
JuMP.optimize!(model)
JuMP.value.(p)
JuMP.value.(x)