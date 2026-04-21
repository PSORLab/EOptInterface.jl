# Copyright (c) 2025 Joseph Choi, Dimitri Alston, Pengfei Xu, Matthew Stuber,
# and the University of Connecticut (UConn)
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EOptInterface
# An abstraction layer for optimizing equation-oriented/acausal models
# https://github.com/PSORLab/EOptInterface.jl
################################################################################
# src/EOptInterface.jl
# The main file for EOptInterface.
################################################################################

module EOptInterface

# ---- imports ----
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

using JuMP

# ---- includes ----
include("basefuncs.jl")
include("userfuncs.jl")

# ---- exports ----
export decision_vars, full_solutions, register_nlsystem, register_odesystem

end # module
