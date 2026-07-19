module NDMCExample

using EOptInterface
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using JuMP, Ipopt
using OrdinaryDiffEq, DiffEqCallbacks
using DataFrames, CSV
using Dates
using Printf
using Sockets
using Plots

include("ndmc_case.jl")
include("ndmc_profile.jl")
include("ndmc_plots.jl")

# Main experiment settings and output filenames.
export NDMC_APPLIED_CONTROL_FILENAME,
       NDMC_CLOSED_LOOP_FILENAME,
       NDMC_DEFAULT_OUTPUT_DIR,
       NDMC_NAME_MAP_FILENAME,
       NDMCCaseResult,
       NDMCMPCConfig,
       build_ndmc_controller,
       build_ndmc_system,
       load_config_from_env,
       main,
       ndmc_case_output_paths,
       run_ndmc_case,
       save_ndmc_case_outputs!,

       # Plotting used by the notebook.
       NDMC_NOTEBOOK_PLOT_CONFIG,
       NDMC_PLOT_CONFIG,
       apply_eoi_publication_style!,
       ndmc_plot_closed_loop_input,
       ndmc_plot_closed_loop_outputs,
       ndmc_plot_difference_pair,
       ndmc_plot_full_horizon_pair,
       ndmc_plot_legacy_comparison_pair,
       ndmc_plot_penalty_terms,
       ndmc_plot_shock_window_c3,
       ndmc_plot_shock_window_pair,
       ndmc_plot_shock_window_qair,
       ndmc_plot_supporting_states

# Timing functions remain available as `NDMCExample.<function>` but are not
# imported by default into the scientific notebook namespace.

end # module
