# EOptInterface.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://joseph03choi.github.io/EOptInterface.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://joseph03choi.github.io/EOptInterface.jl/dev/)
[![Build Status](https://github.com/joseph03choi/EOptInterface.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/joseph03choi/EOptInterface.jl/actions/workflows/CI.yml?query=branch%3Amaster)

EOptInterface.jl is an abstraction layer for automatically formulating JuMP mathematical programming models from ModelingToolkit equation-oriented/acausal models.

## Feature Summary

```julia
decision_vars(::System)
```
Displays the optimization problem decision variables.

```julia
register_nlsystem(::Model, ::System, obj::Num, ineqs::Vector{Num})
```
Registers algebraic JuMP constraints and objective from ModelingToolkit algebraic models built using `@mtkbuild`.

```julia
full_solutions(::Model, ::System)
```
Returns a dictionary of optimal solution values for all eliminated variables from ModelingToolkit's structural simplification step.

```julia
register_odesystem(::Model, ::System, tspan::Tuple{Number,Number}, tstep::Number, solver::String)
```
Registers algebraic JuMP constraints from ModelingToolkit differential equation models built using `@mtkbuild`. Available integration schemes: `"EE", "IE"`




## Example Usage
Optimizing algebraic `@mtkbuild` models using EAGO solver
```julia
using ModelingToolkit, JuMP, EOptInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

@connector Stream begin
    @variables begin
        F(t),   [input=true]
        y_A(t), [input=true]
        y_B(t), [input=true]
        y_C(t), [input=true]
    end
    @parameters begin
        V_A = 8.937e-2
        V_B = 1.018e-1
        V_C = 1.13e-1
    end
end
@mtkmodel Influent begin
    @components begin
        out = Stream()
    end
    @parameters begin
        F
        y_A = 1
        y_B = 0
        y_C = 0
    end
    @equations begin
        out.F ~ F
        out.y_A ~ y_A
        out.y_B ~ y_B
        out.y_C ~ y_C
    end
end
@mtkmodel Mixer begin
    @components begin
        in1 = Stream()
        in2 = Stream()
        out = Stream()
    end
    @equations begin
        out.F ~ in1.F + in2.F
        out.y_A ~ (in1.y_A*in1.F + in2.y_A*in2.F)/(in1.F + in2.F)
        out.y_B ~ (in1.y_B*in1.F + in2.y_B*in2.F)/(in1.F + in2.F)
        out.y_C ~ (in1.y_C*in1.F + in2.y_C*in2.F)/(in1.F + in2.F)
    end
end
@mtkmodel CSTR begin
    @components begin
        in = Stream()
        out = Stream()
    end
    @parameters begin
        V
        k_1 = 0.4
        k_2 = 0.055
    end
    begin
        r_1 = k_1*out.y_A/(out.y_A*in.V_A + out.y_B*in.V_B + out.y_C*in.V_C)
        r_2 = k_2*out.y_B/(out.y_A*in.V_A + out.y_B*in.V_B + out.y_C*in.V_C)
    end
    @equations begin
        out.F ~ in.F
        out.y_A + out.y_B + out.y_C ~ 1
        out.y_B*out.F ~ in.y_B*in.F + (r_1 - r_2)*V
        out.y_C*out.F ~ in.y_C*in.F + r_2*V
    end
end
@mtkmodel Separator1 begin
    @components begin
        in = Stream()
        outV = Stream()
        outL = Stream()
    end
    @equations begin
        in.F ~ outV.F + outL.F
        in.y_B*in.F ~ outL.y_B*outL.F
        in.y_C*in.F ~ outL.y_C*outL.F
        
        outV.y_A + outV.y_B + outV.y_C ~ 1
        outV.y_C ~ 0
        outV.y_B ~ 0

        outL.y_A + outL.y_B + outL.y_C ~ 1
        outL.y_A ~ 0
    end
end
@mtkmodel Separator2 begin
    @components begin
        in = Stream()
        outV = Stream()
        outL = Stream()
    end
    @equations begin
        in.F ~ outV.F + outL.F
        in.y_B*in.F ~ outV.F

        outV.y_A + outV.y_B + outV.y_C ~ 1
        outV.y_A ~ 0
        outV.y_C ~ 0

        outL.y_A + outL.y_B + outL.y_C ~ 1
        outL.y_A ~ 0
        outL.y_B ~ 0
    end
end
@mtkmodel ReactorSeparatorRecycle begin
    @components begin
        influent = Influent()
        mixer = Mixer()
        cstr = CSTR()
        sep1 = Separator1()
        sep2 = Separator2()
    end
    @equations begin
        connect(influent.out, mixer.in1)
        connect(mixer.out, cstr.in)
        connect(cstr.out, sep1.in)
        connect(sep1.outV, mixer.in2)
        connect(sep1.outL, sep2.in)
    end
end

@mtkcompile s = ReactorSeparatorRecycle()

exprF5 = s.sep2.outV.F
exprTau = s.cstr.V/(s.cstr.out.F*(s.cstr.out.y_A*s.cstr.in.V_A + s.cstr.out.y_B*s.cstr.in.V_B + s.cstr.out.y_C*s.cstr.in.V_C))
f_CSTR = (25764 + 8178*s.cstr.V)/2.5
s1cap = 132718 + s.cstr.out.F*(369*s.cstr.out.y_A - 1113.9*s.cstr.out.y_B)
s2cap = 25000 + s.sep1.outL.F*(6984.5*s.sep1.outL.y_B - 3869.53*s.sep1.outL.y_C^2)
s1op = s.cstr.out.F*(3+36.11*s.cstr.out.y_A + 7.71*s.cstr.out.y_B)*26.32e-3
s2op = s.sep1.outL.F*(26.21 + 29.45*s.sep1.outL.y_B)*26.32e-3;
f_Sep = (s1cap+s2cap)/2.5 + 0.52*(s1op+s2op)
g1 = 25 - exprF5
g2 = 475/3600 - exprTau
obj = f_CSTR + f_Sep

using EAGO
model = Model(EAGO.Optimizer)
decision_vars(s)
xL = zeros(6)
xU = [100, 1, 1, 1, 100, 10]
@variable(model, xL[i] <= x[i=1:6] <= xU[i])
register_nlsystem(model, s, obj, [g1, g2])
JuMP.optimize!(model)
full_solutions(model, s)
```
Optimizing ODE `@mtkbuild` models using EAGO solver
```julia
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
include("kinetic_intensity_data.jl") # see \examples\kinetic_intensity_data.jl
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
```

## References
1. Y. Ma, S. Gowda, R. Anantharaman, C. Laughman, V. Shah, and C. Rackauckas, **ModelingToolkit: A composable graph transformation system for equation-based modeling**, 2021.
2. M. Lubin, O. Dowson, J. Dias Garcia, J. Huchette, B. Legat, and J. P. Vielma, **JuMP 1.0: Recent improvements to a modeling language for mathematical optimization**, *Mathematical Programming Computation*, vol. 15, p. 581-589, 2023.
3. M. Wilhelm and M. Stuber, **EAGO.jl Easy Advanced Global Optimization in Julia**, *Optimization Methods and Software*, vol. 37, no. 2, pp. 425-450, 2022.
