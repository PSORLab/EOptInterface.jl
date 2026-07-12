#
# NDMC conductivity MPC example
#
# Run the repository NDMC conductivity MPC experiment.

import Pkg

Pkg.activate(@__DIR__; io=devnull)
Pkg.instantiate(; io=devnull)

using EOptInterface

include(joinpath(@__DIR__, "NDMCExample.jl"))

NDMCExample.main(; output_dir=joinpath(@__DIR__, "generated"))
