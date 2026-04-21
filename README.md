# EOptInterface.jl

`EOptInterface.jl`, or <ins>**E**</ins>quation-oriented <ins>**Opt**</ins>imization <ins>**Interface**</ins>, is an abstraction layer for optimizing equation-oriented/acausal models built in [`ModelingToolkit.jl`](https://github.com/SciML/ModelingToolkit.jl) [[1](#references)] using [`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) [[2](#references)].

| **PSOR Lab** | **Build Status**                                                                                |
|:------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/Developed_by-PSOR_Lab-342674)](https://psor.uconn.edu/) | [![Build Status](https://github.com/PSORLab/EOptInterface.jl/workflows/CI/badge.svg?branch=main)](https://github.com/PSORLab/EOptInterface.jl/actions?query=workflow%3ACI) [![codecov](https://codecov.io/gh/PSORLab/EOptInterface.jl/graph/badge.svg?token=1x9pOV439N)](https://codecov.io/gh/PSORLab/EOptInterface.jl) |

<!-- 
 **Documentation**                                                 |
 :-----------------------------------------------------------------:|
  [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PSORLab.github.io/EOptInterface.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PSORLab.github.io/EOptInterface.jl/dev/) |
-->

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
full_solution(::JuMP.Model, ::ModelingToolkit.System)
```
Returns a dictionary of optimal solution values for the observed variables of an algebraic `ModelingToolkit` system if the `JuMP` model is solved.

```julia
register_odesystem(::JuMP.Model, ::ModelingToolkit.System, tspan::Tuple{Real,Real}, tstep::Real, solver::String)
```
Automatically applies forward transcription and registers the discretized ODE `ModelingToolkit` system as algebraic `JuMP` constraints. Available integration methods: `"EE"` (explicit Euler), `"IE"` (implicit Euler).

## Examples

The code for these examples can be found in [`src/examples/`](https://github.com/PSORLab/EOptInterface.jl/tree/main/examples).

### [Algebraic System](https://github.com/PSORLab/EOptInterface.jl/blob/main/examples/algebraic_model.jl)

An optimal reactor-separator-recycle process design problem originally presented by [[3](#references)] is used to demonstrate the use of `register_nlsystem` to formulate and solve a reduced-space model using the deterministic global optimizer [`EAGO.jl`](https://github.com/PSORLab/EAGO.jl) [[4](#references)].

### [ODE System](https://github.com/PSORLab/EOptInterface.jl/blob/main/examples/ode_model.jl)

A nonlinear kinetic parameter estimation problem originally described by [[5](#references)] is used to demonstrate the use of `register_odesystem` to formulate and solve an ODE system using `Ipopt` [[6](#references)].

## References
1. Ma, Y., Gowda, S., Anantharaman, R., Laughman, C., Shah, V., and Rackauckas, C. ModelingToolkit: A Composable Graph Transformation System For Equation-Based Modeling. (2022). DOI: [10.48550/arXiv.2103.05244)](https://doi.org/10.48550/arXiv.2103.05244)
2. Lubin, M., Dowson, O., Dias Garcia, J., Huchette, J., Legat, B., and Vielma, J.P. JuMP 1.0: recent improvements to a modeling language for mathematical optimization. *Mathematical Programming Computation.* 15, 581-589 (2023). DOI: [10.1007/s12532-023-00239-3](https://doi.org/10.1007/s12532-023-00239-3)
3. Kokossis, A.C. and Floudas, C.A. Synthesis of isothermal reactor-separator-recycle systems. *Chemical Engineering Science.* 46, 1361-1383 (1991). DOI: [10.1016/0009-2509(91)85063-4](https://doi.org/10.1016/0009-2509(91)85063-4)
4. Wilhelm, M. E. and Stuber, M.D. EAGO.jl: easy advanced global optimization in Julia. *Optimization Methods and Software.* 37(2), 425-450 (2022). DOI: [10.1080/10556788.2020.1786566](https://doi.org/10.1080/10556788.2020.1786566)
5. Taylor, J.W. Direct measurement and analysis of cyclohexadienyl oxidation. Ph.D. thesis, Massachusetts Institute of Technology. (2005). URL: http://hdl.handle.net/1721.1/33716
6. Wächter, A. and Biegler, L.T. On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming. *Mathematical Programming.* 106, 25-57 (2006). DOI: [10.1007/s10107-004-0559-y](https://doi.org/10.1007/s10107-004-0559-y)
