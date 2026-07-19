# EOptInterface.jl

`EOptInterface.jl` connects equation-oriented models written in
[`ModelingToolkit.jl`](https://github.com/SciML/ModelingToolkit.jl) with
optimization problems written in [`JuMP.jl`](https://github.com/jump-dev/JuMP.jl).

This `MPC` branch is a research branch. It contains the tracking MPC, DMC, NDMC,
DAE, and wastewater examples developed on the ModelingToolkit 10 code line. The
released package and newer ModelingToolkit 11 development are maintained on the
repository's `main` branch. Results from this branch should therefore be cited
with the branch name and commit hash.

[![Repository](https://img.shields.io/badge/repository-PSORLab%2FEOptInterface.jl-342674)](https://github.com/PSORLab/EOptInterface.jl)
[![Main documentation](https://img.shields.io/badge/docs-main-blue.svg)](https://PSORLab.github.io/EOptInterface.jl/stable/)

## Research workflow

The main operations in this branch are:

1. register a ModelingToolkit algebraic, ODE, or DAE model in JuMP;
2. build one tracking MPC problem from the dynamic model;
3. update measured states, setpoints, and disturbance previews online;
4. solve the existing JuMP problem and apply the first control move;
5. record trajectories and objective terms for analysis.

The JuMP MPC model is built once. Closed-loop updates change numerical data and
solve the same model again rather than rebuilding it at every sample.

## Installation

Clone the repository and select the research branch:

```bash
git clone https://github.com/PSORLab/EOptInterface.jl.git
cd EOptInterface.jl
git switch MPC
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the package tests:

```bash
julia --project=. test/runtests.jl
```

The examples use a separate environment:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## Main functions

- `register_nlsystem(...)` adds algebraic model equations to JuMP.
- `register_odesystem(...)` discretizes an ODE model with `EE`, `IE`/`BDF1`,
  `RK4`, or `IRK4` equations.
- `register_daesystem(...)` retains differential and algebraic equations in the
  prediction problem.
- `build_tracking_mpc(...)` builds the tracking objective, state trajectories,
  control trajectories, bounds, and move penalties.
- `solve_tracking_mpc!(...)` updates one online MPC step and returns the control
  and predicted-state trajectories.
- `register_dmcsystem(...)` builds a step-response DMC prediction block.

See [`docs/src/mpc_module.md`](docs/src/mpc_module.md) for a compact tracking MPC
example.

## Examples

Small tracking MPC example:

```bash
julia --project=examples examples/tracking_mpc_demo.jl
```

ODE and DAE registration examples:

```bash
julia --project=examples examples/ode_model.jl
julia --project=examples examples/dae_registration_demo.jl
```

The ODE parameter-estimation example uses EAGO for global optimization and can
take several minutes. It is not part of the quick regression test.

NDMC conductivity closed-loop experiment:

```bash
NDMC_SHOW_DETAILED_STATUS=1 \
  julia --project=examples examples/ndmc_conductivity_mpc_demo.jl
```

The NDMC implementation is organized as:

- `examples/ndmc_conductivity_mpc_demo.jl`: command-line entry point;
- `examples/NDMCExample.jl`: names used by the script and notebook;
- `examples/ndmc_case.jl`: plant, controller, and closed-loop simulation;
- `examples/ndmc_profile.jl`: optional timing and benchmarking reports;
- `examples/ndmc_plots.jl`: NDMC figures;
- `notebooks/ndmc_conductivity_mpc_simple.ipynb`: interactive experiment.

## Canonical NDMC case

| Setting | Value |
|---|---:|
| Simulation interval | 0-4000 s |
| MPC sample time | 20 s |
| Saved trajectory interval | 10 s |
| Prediction span | 400 s |
| Move span | 60 s |
| Conductivity setpoint | 280 |
| Aeration bounds | 0-800 |
| Influent shock | 2100-2250 s |
| Shock value in zone 3 | 320 |

The canonical result files in `examples/generated/` use this case. Timing and
machine-information files are intentionally not versioned because they depend
on the computer and Julia session used for the run.

## Notebook use

Start the notebook from the repository root so its first cell can locate
`notebooks/bootstrap_examples_notebook_env.jl`. The first cell activates and
instantiates the shared `examples/` environment. The main experiment cells are
separate from optional timing cells so scientific results can be reproduced
without running the benchmarking workflow.

## Documentation

Build the local documentation with:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Reproducibility note

When reporting numerical results, record the branch, commit hash, Julia version,
solver version, controller settings, and whether the run used the canonical
10 s saved-output interval. The legacy DMC comparison uses a deterministic Ipopt
initialization setting.

## References

1. Y. Ma et al., "ModelingToolkit: A composable graph transformation system for
   equation-based modeling," 2021.
2. M. Lubin et al., "JuMP 1.0: Recent improvements to a modeling language for
   mathematical optimization," *Mathematical Programming Computation*, 2023.
