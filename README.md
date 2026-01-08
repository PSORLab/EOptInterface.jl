# EOptInterface.jl

`EOptInterface.jl`, or <ins>**E**</ins>quation-oriented <ins>**Opt**</ins>imization <ins>**Interface**</ins>, is an abstraction layer for automatically formulating mathematical programming `Model`s in `JuMP` directly from componentized equation-oriented models built using `ModelingToolkit`'s high-level `@mtkmodel` interface.

| **PSOR Lab** | **Documentation**                                                 | **Build Status**                                                                                |
|:------------:|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/Developed_by-PSOR_Lab-342674)](https://psor.uconn.edu/) | [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PSORLab.github.io/EOptInterface.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PSORLab.github.io/EOptInterface.jl/dev/) | [![Build Status](https://github.com/PSORLab/EOptInterface.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/PSORLab/EOptInterface.jl/actions/workflows/CI.yml?query=branch%3Amaster) |



## Feature Summary

```julia
decision_vars(::ModelingToolkit.System)
```
Returns the decision variables for an optimization problem formulated from a `ModelingToolkit` system.

```julia
register_nlsystem(::JuMP.Model, ::ModelingToolkit.System, obj::Symbolics.Num, ineqs::Vector{Symbolics.Num})
```
Automatically formulates algebraic `JuMP` constraints and objective function from
an algebraic `ModelingToolkit` system and user-provided constraints and objective symbolic expressions.

```julia
full_solutions(::JuMP.Model, ::ModelingToolkit.System)
```
Returns a dictionary of optimal solution values for the observed variables of an algebraic `ModelingToolkit` system if the `JuMP` model is solved.

```julia
register_odesystem(::JuMP.Model, ::ModelingToolkit.System, tspan::Tuple{Number,Number}, tstep::Number, solver::String)
```
Automatically applies forward transcription and registers the discretized ODE `ModelingToolkit` system as algebraic `JuMP` constraints. Available integration methods: `"EE"` (explicit Euler), `"IE"` (implicit Euler).




## Example Usage: Algebraic Models

An optimal reactor-separator-recycle process design problem originally presented by [Kokosis and Floudas⁵](https://doi.org/10.1016/0009-2509(91)85063-4) is used to demonstrate the use of `register_nlsystem` to formulate and solve a reduced-space `Model` using the deterministic global optimizer `EAGO`.

<p align="center">
  <img width="1111" height="579" alt="image" src="https://github.com/user-attachments/assets/310a47fc-504e-4864-9d77-8f8165e4feb0" />
  Figure 1: Illustration of the reactor-separator-recycle system.
</p>

```julia
using ModelingToolkit, JuMP, EOptInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

# Formulate MTK algebraic System
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
        F       # Free design variable
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
        V       # Free design variable
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

# Define symbolic expressions of constraints and objective
# Use syntax System.Component.Connector.Variable, System.Component.Component.Parameter, or System.Component.Parameter
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

# Solve using EAGO
using EAGO
model = Model(EAGO.Optimizer)
decision_vars(s) # Displays: sep1₊in₊F(t), sep1₊in₊y_B(t), sep1₊in₊y_C(t), sep1₊outL₊y_C(t), influent₊F, cstr₊V
xL = zeros(6) # lower bound on x
xU = [100, 1, 1, 1, 100, 10] # upper bound on x
@variable(model, xL[i] <= x[i=1:6] <= xU[i]) # ̂x = (̂z,p), ̂z = (z(t),...), p = (influent₊F, cstr₊V)
register_nlsystem(model, s, obj, [g1, g2])
JuMP.optimize!(model)
JuMP.value.(x)

# Obtain observed variable solutions
full_solutions(model, s)
```

## Example Usage: ODE Models

A nonlinear kinetic parameter estimation problem originally described by [Taylor⁶](http://hdl.handle.net/1721.1/33716) is used to demonstrate the use of `register_odesystem` to formulate a `Model` from an ODE `System`, solved using `Ipopt`.

```julia
using ModelingToolkit, JuMP, EOptInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

# Formulate MTK ODE System
@mtkmodel KineticParameterEstimation begin
    @parameters begin
        T = 273
        K_2 = 46*exp(6500/T-18)
        K_3 = 2*K_2
        k_1 = 53
        k_1s = k_1*1e-6
        k_5 = 1.2e-3
        c_O2 = 2e-3

        k_2f    # Free design variable
        k_3f    # Free design variable
        k_4     # Free design variable
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

# Solve using Ipopt
using Ipopt
model = Model(Ipopt.Optimizer)
decision_vars(o) # Displays: x_Z(t), x_Y(t), x_D(t), x_B(t), x_A(t), k_2f, k_3f, k_4
# FIRST, create discretized differntial state decision variables "z"
N = Int(floor((tspan[2] - tspan[1])/tstep))+1
V = length(unknowns(o))
zL = zeros(V) # lower bound on z
zU = [140.0, 0.4, 140.0, 140.0, 140.0] # upper bound on z
@variable(model, zL[i] <= z[i in 1:V,1:N] <= zU[i]) # ̇z = (x_Z(t), x_Y(t), x_D(t), x_B(t), x_A(t))
# SECOND, create free design decision variables "p"
pL = [10, 10, 0.001] # lower bound on p
pU = [1200, 1200, 40] # upper bound on p
@variable(model, pL[i] <= p[i=1:3] <= pU[i]) # p = (k_2f, k_3f, k_4)
register_odesystem(model, o, tspan, tstep, "EE") # Try changing between "EE" and "IE"!
@objective(model, Min, sum((intensity(z[5,i],z[4,i],z[3,i]) - data[i-1])^2 for i in 2:N))
JuMP.optimize!(model)
```

## References
1. Y. Ma, S. Gowda, R. Anantharaman, C. Laughman, V. Shah, and C. Rackauckas, **ModelingToolkit: A composable graph transformation system for equation-based modeling**, 2021.
2. M. Lubin, O. Dowson, J. Dias Garcia, J. Huchette, B. Legat, and J. P. Vielma, **JuMP 1.0: Recent improvements to a modeling language for mathematical optimization**, *Mathematical Programming Computation*, vol. 15, p. 581-589, 2023.
3. M. Wilhelm and M. Stuber, **EAGO.jl Easy Advanced Global Optimization in Julia**, *Optimization Methods and Software*, vol. 37, no. 2, p. 425-450, 2022.
4. A. Wächter, & L. T. Biegler, **On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming**, *Mathematical programming*, vol. 106, no. 1, p. 25-57, 2006.
5. A. C. Kokossis, C. A. Floudas, **Synthesis of isothermal reactor-separator-recycle systems**, *Chemical Engineering Science*, vol. 46, p. 1361-1383, 1991.
6. J. W. Taylor, **Direct measurement and analysis of cyclohexadienyl oxidation**, Ph.D. thesis, Massachusetts Institute of Technology, 2005.

