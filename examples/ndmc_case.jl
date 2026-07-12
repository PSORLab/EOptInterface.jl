# NDMC conductivity MPC case.
#
# Reading path for the scientific example:
# 1. constants and rate expressions define the plant;
# 2. `NDMCPlant` gives the ModelingToolkit equations;
# 3. `build_ndmc_controller(...)` defines the MPC problem;
# 4. `run_ndmc_case(...)` simulates the closed loop.
#
# Timing and BenchmarkTools helpers are kept in this same file because the
# notebook reports them, but they are optional reporting code rather than part
# of the controller definition.

using BenchmarkTools
using Statistics
using XLSX

# Values carried over from the legacy NDMC script so the comparison uses the
# same sample time and horizons.
const NDMC_LEGACY_DT = 20.0
const NDMC_LEGACY_PSPAN = 400.0
const NDMC_LEGACY_MSPAN = 60.0
const NDMC_LEGACY_PH = Int(round(NDMC_LEGACY_PSPAN / NDMC_LEGACY_DT))
const NDMC_LEGACY_CH = Int(round(NDMC_LEGACY_MSPAN / NDMC_LEGACY_DT))

# Physical and empirical parameters used by the NDMC plant model.
const NDMC_M_N = 14.0067
const NDMC_R_AOMAX = 0.67
const NDMC_K_OAO = 0.3
const NDMC_X_AO = 0.505
const NDMC_PHI_OAO = 2.5
const NDMC_SOTE = 0.1
const NDMC_COS = 9.1
const NDMC_VOLUME_TOTAL = 1000.0
const NDMC_ZONE_VOLUME = NDMC_VOLUME_TOTAL / 4.0
const NDMC_FLOW = NDMC_VOLUME_TOTAL / 240.0 / 4.0
const NDMC_LA0 = 149.6
const NDMC_A = 60.2
const NDMC_B = 0.229
const NDMC_SMOOTH_POS_EPS = 1e-8

const NDMC_DEFAULT_OUTPUT_DIR = normpath(joinpath(@__DIR__, "generated"))
const NDMC_CLOSED_LOOP_FILENAME = "ndmc_conductivity_eopt_closed_loop.csv"
const NDMC_APPLIED_CONTROL_FILENAME = "ndmc_conductivity_eopt_applied_control.csv"
const NDMC_NAME_MAP_FILENAME = "ndmc_conductivity_eopt_name_map.csv"
const NDMC_PROFILE_SUMMARY_FILENAME = "ndmc_conductivity_eopt_profile_summary.csv"
const NDMC_PROFILE_STEPS_FILENAME = "ndmc_conductivity_eopt_profile_steps.csv"
const NDMC_PROFILE_ENV_FILENAME = "ndmc_conductivity_eopt_profile_env.csv"
const NDMC_TIMING_WORKBOOK_FILENAME = "ndmc_conductivity_eopt_timing_tables.xlsx"

_elapsed_s(start_ns::Integer)::Float64 = (time_ns() - start_ns) / 1.0e9

# ---------------------------------------------------------------------------
# Optional timing/reporting helpers
# ---------------------------------------------------------------------------

"""
    NDMCTimingProfile

Optional wall-clock timing table for the NDMC closed-loop example.

The controller can run without this object. It is only used when the notebook
or script is collecting timing tables for reporting.
"""
Base.@kwdef mutable struct NDMCTimingProfile
    enabled::Bool = false
    run_label::String = "single"
    offline_prewarm_s::Float64 = NaN
    build_system_s::Float64 = NaN
    build_tracking_mpc_s::Float64 = NaN
    build_controller_total_s::Float64 = NaN
    initial_mpc_step_total_s::Float64 = NaN
    closed_loop_simulation_s::Float64 = NaN
    closed_loop_mpc_callback_s::Float64 = 0.0
    closed_loop_log_callback_s::Float64 = 0.0
    closed_loop_disturbance_callback_s::Float64 = 0.0
    closed_loop_integrator_only_s::Float64 = NaN
    total_script_s::Float64 = NaN
    step_index::Vector{Int} = Int[]
    time::Vector{Float64} = Float64[]
    state_capture_s::Vector{Float64} = Float64[]
    forecast_update_s::Vector{Float64} = Float64[]
    prepare_tracking_mpc_step_s::Vector{Float64} = Float64[]
    jump_optimize_s::Vector{Float64} = Float64[]
    solve_postprocess_s::Vector{Float64} = Float64[]
    status_print_s::Vector{Float64} = Float64[]
    solve_tracking_mpc_s::Vector{Float64} = Float64[]
    step_total_s::Vector{Float64} = Float64[]
    status::Vector{String} = String[]
    objective::Vector{Float64} = Float64[]
    q_prev::Vector{Float64} = Float64[]
    q_apply::Vector{Float64} = Float64[]
end

"""
    NDMCCaseResult

Results from one NDMC closed-loop run.

The fields store the model, controller, solution, and data tables used by the
notebook figures.
"""
struct NDMCCaseResult
    cfg
    sys
    ctrl::TrackingMPCController
    logctx::MPCLog
    solution
    closed_loop_df::DataFrame
    applied_df::DataFrame
    name_map_df::DataFrame
    profile::NDMCTimingProfile
end

_env_flag(name::AbstractString; default::AbstractString="false") =
    lowercase(strip(get(ENV, name, default))) in ("1", "true", "yes", "on")

_float_or_nan(x) = x isa Real ? Float64(x) : NaN

ndmc_profiling_enabled() = _env_flag("NDMC_ENABLE_PROFILING")
ndmc_offline_prewarm_enabled() = _env_flag("NDMC_OFFLINE_PREWARM")

function _ndmc_profile_mode(profile_mode)::Symbol
    mode = profile_mode isa Symbol ? profile_mode : Symbol(lowercase(strip(String(profile_mode))))
    mode in (:single, :cold_and_warm) ||
        error("NDMC profile mode must be `:single` or `:cold_and_warm`.")
    return mode
end

function ndmc_profile_mode()::Symbol
    return _ndmc_profile_mode(get(ENV, "NDMC_PROFILE_MODE", "single"))
end

function ndmc_case_output_paths(output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR)
    return (
        closed_loop = joinpath(output_dir, NDMC_CLOSED_LOOP_FILENAME),
        applied = joinpath(output_dir, NDMC_APPLIED_CONTROL_FILENAME),
        name_map = joinpath(output_dir, NDMC_NAME_MAP_FILENAME),
    )
end

function ndmc_profile_output_paths(output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR)
    return (
        summary = joinpath(output_dir, NDMC_PROFILE_SUMMARY_FILENAME),
        steps = joinpath(output_dir, NDMC_PROFILE_STEPS_FILENAME),
        env = joinpath(output_dir, NDMC_PROFILE_ENV_FILENAME),
    )
end

function ndmc_timing_workbook_path(output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR)
    return joinpath(output_dir, NDMC_TIMING_WORKBOOK_FILENAME)
end

"""
    save_ndmc_case_outputs!(result; output_dir=NDMC_DEFAULT_OUTPUT_DIR, announce=true)

Write the NDMC case outputs to CSV.
"""
function save_ndmc_case_outputs!(result::NDMCCaseResult;
                                 output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR,
                                 announce::Bool=true)
    mkpath(output_dir)
    paths = ndmc_case_output_paths(output_dir)
    CSV.write(paths.closed_loop, result.closed_loop_df)
    CSV.write(paths.applied, result.applied_df)
    CSV.write(paths.name_map, result.name_map_df)
    if announce
        println("Saved NDMC example outputs to:")
        println("  ", paths.closed_loop)
        println("  ", paths.applied)
        println("  ", paths.name_map)
    end
    return paths
end

"""
    save_ndmc_profile_outputs!(summary_df, steps_df, env_df; output_dir=..., announce=true)

Write the NDMC profiling tables to CSV.
"""
function save_ndmc_profile_outputs!(summary_df::DataFrame,
                                    steps_df::DataFrame,
                                    env_df::DataFrame;
                                    output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR,
                                    announce::Bool=true)
    mkpath(output_dir)
    paths = ndmc_profile_output_paths(output_dir)
    CSV.write(paths.summary, summary_df)
    CSV.write(paths.steps, steps_df)
    CSV.write(paths.env, env_df)
    if announce
        println("Saved NDMC profiling outputs to:")
        println("  ", paths.summary)
        println("  ", paths.steps)
        println("  ", paths.env)
    end
    return paths
end

"""
    save_ndmc_timing_workbook!(sheet_tables...; output_dir=..., filename=..., announce=true)

Write multiple NDMC timing tables to one Excel workbook.

Each argument must be a `"sheet_name" => dataframe` pair. This is intended for
the notebook workflow, where several readable timing tables and isolated
BenchmarkTools tables are written into one workbook for reporting.
"""
function save_ndmc_timing_workbook!(sheet_tables::Pair...;
                                    output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR,
                                    filename::AbstractString=NDMC_TIMING_WORKBOOK_FILENAME,
                                    announce::Bool=true)
    isempty(sheet_tables) && error("No NDMC timing tables were provided for workbook export.")

    normalized_tables = Pair{String,DataFrame}[]
    for (sheet_name, table) in sheet_tables
        sheet_name isa AbstractString ||
            error("Workbook sheet names must be strings. Got $(typeof(sheet_name)).")
        table isa AbstractDataFrame ||
            error("Workbook tables must be DataFrames or AbstractDataFrames. Got $(typeof(table)) for sheet \"$(sheet_name)\".")

        ncodeunits(sheet_name) <= 31 ||
            error("Excel sheet name \"$(sheet_name)\" is longer than Excel's 31-character limit.")
        push!(normalized_tables, String(sheet_name) => DataFrame(table))
    end

    mkpath(output_dir)
    workbook_path = joinpath(output_dir, filename)
    XLSX.writetable(workbook_path, normalized_tables...; overwrite=true)

    if announce
        println("Saved NDMC timing workbook to:")
        println("  ", workbook_path)
    end
    return workbook_path
end

"""
    record_ndmc_step_timing!(profile; kwargs...)

Append one MPC-step timing row to the NDMC profiling collector.
"""
function record_ndmc_step_timing!(profile::NDMCTimingProfile;
                                  step_index::Integer,
                                  time::Real,
                                  state_capture_s::Real,
                                  forecast_update_s::Real,
                                  prepare_tracking_mpc_step_s::Real,
                                  jump_optimize_s::Real,
                                  solve_postprocess_s::Real,
                                  status_print_s::Real,
                                  solve_tracking_mpc_s::Real,
                                  step_total_s::Real,
                                  status,
                                  objective,
                                  q_prev::Real,
                                  q_apply::Real)
    profile.enabled || return profile
    push!(profile.step_index, Int(step_index))
    push!(profile.time, Float64(time))
    push!(profile.state_capture_s, Float64(state_capture_s))
    push!(profile.forecast_update_s, Float64(forecast_update_s))
    push!(profile.prepare_tracking_mpc_step_s, Float64(prepare_tracking_mpc_step_s))
    push!(profile.jump_optimize_s, Float64(jump_optimize_s))
    push!(profile.solve_postprocess_s, Float64(solve_postprocess_s))
    push!(profile.status_print_s, Float64(status_print_s))
    push!(profile.solve_tracking_mpc_s, Float64(solve_tracking_mpc_s))
    push!(profile.step_total_s, Float64(step_total_s))
    push!(profile.status, string(status))
    push!(profile.objective, _float_or_nan(objective))
    push!(profile.q_prev, Float64(q_prev))
    push!(profile.q_apply, Float64(q_apply))
    return profile
end

function ndmc_step_timing_dataframe(profile::NDMCTimingProfile)
    return DataFrame(
        run_label = fill(profile.run_label, length(profile.step_index)),
        step_index = profile.step_index,
        time = profile.time,
        state_capture_s = profile.state_capture_s,
        forecast_update_s = profile.forecast_update_s,
        prepare_tracking_mpc_step_s = profile.prepare_tracking_mpc_step_s,
        jump_optimize_s = profile.jump_optimize_s,
        solve_postprocess_s = profile.solve_postprocess_s,
        status_print_s = profile.status_print_s,
        solve_tracking_mpc_s = profile.solve_tracking_mpc_s,
        step_total_s = profile.step_total_s,
        status = profile.status,
        objective = profile.objective,
        q_prev = profile.q_prev,
        q_apply = profile.q_apply,
    )
end

_sum_selected(v::AbstractVector{<:Real}, idxs::AbstractVector{<:Integer}) =
    isempty(idxs) ? 0.0 : sum(@view v[idxs])

function ndmc_step_rollups(profile::NDMCTimingProfile)
    initial_idxs = findall(==(0), profile.step_index)
    online_idxs = findall(>(0), profile.step_index)
    return (
        initial_state_capture_s = _sum_selected(profile.state_capture_s, initial_idxs),
        initial_forecast_update_s = _sum_selected(profile.forecast_update_s, initial_idxs),
        initial_prepare_tracking_mpc_step_s = _sum_selected(profile.prepare_tracking_mpc_step_s, initial_idxs),
        initial_jump_solve_s = _sum_selected(profile.jump_optimize_s, initial_idxs),
        initial_solve_postprocess_s = _sum_selected(profile.solve_postprocess_s, initial_idxs),
        initial_status_print_s = _sum_selected(profile.status_print_s, initial_idxs),
        initial_solve_tracking_mpc_total_s = _sum_selected(profile.solve_tracking_mpc_s, initial_idxs),
        initial_step_total_recorded_s = _sum_selected(profile.step_total_s, initial_idxs),
        online_state_capture_s = _sum_selected(profile.state_capture_s, online_idxs),
        online_forecast_update_s = _sum_selected(profile.forecast_update_s, online_idxs),
        online_prepare_tracking_mpc_step_s = _sum_selected(profile.prepare_tracking_mpc_step_s, online_idxs),
        online_jump_solve_s = _sum_selected(profile.jump_optimize_s, online_idxs),
        online_solve_postprocess_s = _sum_selected(profile.solve_postprocess_s, online_idxs),
        online_status_print_s = _sum_selected(profile.status_print_s, online_idxs),
        online_solve_tracking_mpc_total_s = _sum_selected(profile.solve_tracking_mpc_s, online_idxs),
        online_step_total_recorded_s = _sum_selected(profile.step_total_s, online_idxs),
        state_capture_total_s = sum(profile.state_capture_s),
        forecast_update_total_s = sum(profile.forecast_update_s),
        prepare_tracking_mpc_step_total_s = sum(profile.prepare_tracking_mpc_step_s),
        jump_solve_total_s = sum(profile.jump_optimize_s),
        solve_postprocess_total_s = sum(profile.solve_postprocess_s),
        status_print_total_s = sum(profile.status_print_s),
        solve_tracking_mpc_total_s = sum(profile.solve_tracking_mpc_s),
        step_total_recorded_s = sum(profile.step_total_s),
    )
end

function ndmc_summary_dataframe(profile::NDMCTimingProfile,
                                cfg;
                                expected_total_solves::Integer)
    step_rollups = ndmc_step_rollups(profile)
    callback_other_s = profile.closed_loop_log_callback_s + profile.closed_loop_disturbance_callback_s
    online_mpc_preparation_s =
        step_rollups.online_state_capture_s +
        step_rollups.online_forecast_update_s +
        step_rollups.online_prepare_tracking_mpc_step_s
    mpc_callback_residual_s = max(
        0.0,
        profile.closed_loop_mpc_callback_s -
        step_rollups.online_state_capture_s -
        step_rollups.online_forecast_update_s -
        step_rollups.online_prepare_tracking_mpc_step_s -
        step_rollups.online_jump_solve_s -
        step_rollups.online_solve_postprocess_s -
        step_rollups.online_status_print_s,
    )
    offline_prewarm_s = isfinite(profile.offline_prewarm_s) ? profile.offline_prewarm_s : 0.0
    return DataFrame(
        run_label = [profile.run_label],
        offline_prewarm_s = [profile.offline_prewarm_s],
        build_system_s = [profile.build_system_s],
        build_controller_total_s = [profile.build_controller_total_s],
        build_tracking_mpc_s = [profile.build_tracking_mpc_s],
        initial_mpc_step_total_s = [profile.initial_mpc_step_total_s],
        closed_loop_simulation_s = [profile.closed_loop_simulation_s],
        closed_loop_wall_s = [profile.closed_loop_simulation_s],
        closed_loop_mpc_callback_s = [profile.closed_loop_mpc_callback_s],
        closed_loop_log_callback_s = [profile.closed_loop_log_callback_s],
        closed_loop_disturbance_callback_s = [profile.closed_loop_disturbance_callback_s],
        closed_loop_other_callback_s = [callback_other_s],
        closed_loop_integrator_only_s = [profile.closed_loop_integrator_only_s],
        initial_state_capture_s = [step_rollups.initial_state_capture_s],
        initial_forecast_update_s = [step_rollups.initial_forecast_update_s],
        initial_prepare_tracking_mpc_step_s = [step_rollups.initial_prepare_tracking_mpc_step_s],
        initial_jump_solve_s = [step_rollups.initial_jump_solve_s],
        initial_solve_postprocess_s = [step_rollups.initial_solve_postprocess_s],
        initial_status_print_s = [step_rollups.initial_status_print_s],
        initial_solve_tracking_mpc_total_s = [step_rollups.initial_solve_tracking_mpc_total_s],
        initial_step_total_recorded_s = [step_rollups.initial_step_total_recorded_s],
        online_state_capture_s = [step_rollups.online_state_capture_s],
        online_forecast_update_s = [step_rollups.online_forecast_update_s],
        online_prepare_tracking_mpc_step_s = [step_rollups.online_prepare_tracking_mpc_step_s],
        online_jump_solve_s = [step_rollups.online_jump_solve_s],
        online_solve_postprocess_s = [step_rollups.online_solve_postprocess_s],
        online_status_print_s = [step_rollups.online_status_print_s],
        online_solve_tracking_mpc_total_s = [step_rollups.online_solve_tracking_mpc_total_s],
        online_mpc_preparation_s = [online_mpc_preparation_s],
        online_step_total_recorded_s = [step_rollups.online_step_total_recorded_s],
        state_capture_total_s = [step_rollups.state_capture_total_s],
        forecast_update_total_s = [step_rollups.forecast_update_total_s],
        prepare_tracking_mpc_step_total_s = [step_rollups.prepare_tracking_mpc_step_total_s],
        jump_solve_total_s = [step_rollups.jump_solve_total_s],
        solve_postprocess_total_s = [step_rollups.solve_postprocess_total_s],
        status_print_total_s = [step_rollups.status_print_total_s],
        solve_tracking_mpc_total_s = [step_rollups.solve_tracking_mpc_total_s],
        jump_model_assembly_s = [profile.build_tracking_mpc_s],
        jump_total_s = [profile.build_tracking_mpc_s + step_rollups.jump_solve_total_s],
        mpc_callback_residual_s = [mpc_callback_residual_s],
        total_script_s = [profile.total_script_s],
        end_to_end_total_including_offline_prewarm_s = [offline_prewarm_s + profile.total_script_s],
        control_interval_s = [cfg.dt],
        simulation_start_s = [cfg.simulation_span[1]],
        simulation_end_s = [cfg.simulation_span[2]],
        expected_total_solves = [Int(expected_total_solves)],
        recorded_total_solves = [length(profile.step_index)],
    )
end

function _ndmc_profile_env_dataframe(profile_mode::Symbol,
                                     offline_prewarm::Bool,
                                     profiles::AbstractVector{<:NDMCTimingProfile})
    return DataFrame(
        timestamp = [string(now())],
        julia_version = [string(VERSION)],
        active_project = [Base.active_project()],
        nthreads = [Threads.nthreads()],
        cpu_threads = [Sys.CPU_THREADS],
        kernel = [string(Sys.KERNEL)],
        machine = [string(Sys.MACHINE)],
        arch = [string(Sys.ARCH)],
        hostname = [Sockets.gethostname()],
        profile_mode = [string(profile_mode)],
        offline_prewarm_enabled = [offline_prewarm],
        run_count = [length(profiles)],
    )
end

function ndmc_profile_dataframes(profiles::AbstractVector{<:NDMCTimingProfile},
                                 cfg;
                                 expected_total_solves::Integer,
                                 profile_mode::Symbol,
                                 offline_prewarm::Bool=false)
    summary_df = reduce(
        vcat,
        [ndmc_summary_dataframe(profile, cfg; expected_total_solves=expected_total_solves) for profile in profiles],
    )
    steps_df = reduce(vcat, [ndmc_step_timing_dataframe(profile) for profile in profiles], cols=:union)
    env_df = _ndmc_profile_env_dataframe(profile_mode, offline_prewarm, profiles)
    return (summary_df=summary_df, steps_df=steps_df, env_df=env_df)
end

# ---------------------------------------------------------------------------
# NDMC plant model
# ---------------------------------------------------------------------------

# Smooth positive-part helper. This avoids a hard kink at zero.
smooth_nonnegative(x) = 0.5 * (x + sqrt(x * x + NDMC_SMOOTH_POS_EPS^2))

function ndmc_rate_term(cO)
    cO_pos = smooth_nonnegative(cO)
    return NDMC_R_AOMAX * (cO_pos / (NDMC_K_OAO + cO_pos)) * NDMC_X_AO * 0.001 / 60.0 / NDMC_M_N
end

function ndmc_conductivity_sink(cO)
    rate = ndmc_rate_term(cO)
    return (NDMC_LA0 - (NDMC_A + NDMC_B * NDMC_LA0) * sqrt(rate)) * rate * 1e3
end

function ndmc_oxygen_sink(cO)
    cO_pos = smooth_nonnegative(cO)
    return (NDMC_R_AOMAX * (cO_pos / (NDMC_K_OAO + cO_pos))) * NDMC_PHI_OAO * NDMC_X_AO / 60.0
end

function ndmc_k_vector()
    k1_hi = 2.37570423107917e-03
    k2_hi = 1.41065048786101e-03
    k3_hi = 1.50410563455869e-03
    k4_hi = 9.47661458622374e-01

    k1_me = 1.38465080722995e-03
    k2_me = 2.91428968468164e-03
    k3_me = 2.57662674484102e-03
    k4_me = 1.89896205118280e+00

    k1 = 0.5 * (k1_hi + k1_me)
    k2 = 0.5 * (k2_hi + k2_me)
    k3 = 0.5 * (k3_hi + k3_me)
    k4 = 0.5 * (k4_hi + k4_me)
    return (1000.0 / 0.38) .* [k1, k2, k3, k4]
end

"""
    NDMCMPCConfig

Numerical settings for the NDMC closed-loop experiment.

The defaults reproduce the legacy NDMC sample time, horizons, influent shock,
and one shared aeration input.
"""
Base.@kwdef struct NDMCMPCConfig
    simulation_span::Tuple{Float64, Float64} = (0.0, 4000.0)
    dt::Float64 = NDMC_LEGACY_DT
    PH::Int = NDMC_LEGACY_PH
    CH::Int = NDMC_LEGACY_CH
    setpoint::Float64 = 280.0
    move_weight::Float64 = 0.0
    first_move_weight::Float64 = 0.0
    terminal_weight::Float64 = 0.0
    Q_init::Float64 = 168.0
    Q_min::Float64 = 0.0
    Q_max::Float64 = 800.0
    delta_max::Union{Nothing, Float64} = nothing
    use_disturbance_forecast::Bool = false
    disturbance_start::Float64 = 2100.0
    disturbance_stop::Float64 = 2250.0
    Cs::Float64 = 285.0
    Cin_spike::NTuple{3, Float64} = (285.0, 285.0, 320.0)
    save_dt::Float64 = 10.0
    initial_state::NTuple{5, Float64} = (280.0, 280.0, 280.0, 280.0, 0.0)
    state_lower::NTuple{5, Float64} = (0.0, 0.0, 0.0, 0.0, 0.0)
    state_upper::NTuple{5, Float64} = (600.0, 600.0, 600.0, 600.0, 20.0)
    ipopt_tol::Float64 = 1e-6
    ipopt_max_cpu_time::Float64 = 60.0
    show_detailed_status::Bool = true
    show_generic_status::Bool = false
    status_digits::Int = 3
    status_prefix::String = "[NDMC MPC]"
    status_stream::Symbol = :stdout
end

function load_config_from_env()
    dt = parse(Float64, get(ENV, "NDMC_DT", string(NDMC_LEGACY_DT)))
    pspan = parse(Float64, get(ENV, "NDMC_PSPAN", string(NDMC_LEGACY_PSPAN)))
    mspan = parse(Float64, get(ENV, "NDMC_MSPAN", string(NDMC_LEGACY_MSPAN)))
    t_end = parse(Float64, get(ENV, "NDMC_T_END", "4000.0"))
    ch = max(1, Int(round(mspan / dt)))
    ph = max(ch, Int(round(pspan / dt)))
    return NDMCMPCConfig(
        simulation_span = (0.0, t_end),
        dt = dt,
        PH = ph,
        CH = ch,
        setpoint = parse(Float64, get(ENV, "NDMC_SETPOINT", "280.0")),
        move_weight = parse(Float64, get(ENV, "NDMC_MOVE_WEIGHT", "0.0")),
        first_move_weight = parse(Float64, get(ENV, "NDMC_FIRST_MOVE_WEIGHT", "0.0")),
        terminal_weight = parse(Float64, get(ENV, "NDMC_TERMINAL_WEIGHT", "0.0")),
        Q_init = parse(Float64, get(ENV, "NDMC_Q_INIT", "168.0")),
        Q_min = parse(Float64, get(ENV, "NDMC_Q_MIN", "0.0")),
        Q_max = parse(Float64, get(ENV, "NDMC_Q_MAX", "800.0")),
        delta_max = haskey(ENV, "NDMC_DELTA_MAX") ? parse(Float64, ENV["NDMC_DELTA_MAX"]) : nothing,
        use_disturbance_forecast = _env_flag("NDMC_USE_DISTURBANCE_FORECAST"),
        disturbance_start = parse(Float64, get(ENV, "NDMC_DIST_START", "2100.0")),
        disturbance_stop = parse(Float64, get(ENV, "NDMC_DIST_STOP", "2250.0")),
        Cs = parse(Float64, get(ENV, "NDMC_CS", "285.0")),
        Cin_spike = (
            parse(Float64, get(ENV, "NDMC_CIN1", "285.0")),
            parse(Float64, get(ENV, "NDMC_CIN2", "285.0")),
            parse(Float64, get(ENV, "NDMC_CIN3", "320.0")),
        ),
        save_dt = parse(Float64, get(ENV, "NDMC_SAVE_DT", string(dt))),
        initial_state = (
            parse(Float64, get(ENV, "NDMC_C1_INIT", "280.0")),
            parse(Float64, get(ENV, "NDMC_C2_INIT", "280.0")),
            parse(Float64, get(ENV, "NDMC_C3_INIT", "280.0")),
            parse(Float64, get(ENV, "NDMC_CMIX_INIT", "280.0")),
            parse(Float64, get(ENV, "NDMC_CO_INIT", "0.0")),
        ),
        state_lower = (
            parse(Float64, get(ENV, "NDMC_C_LO", "0.0")),
            parse(Float64, get(ENV, "NDMC_C_LO", "0.0")),
            parse(Float64, get(ENV, "NDMC_C_LO", "0.0")),
            parse(Float64, get(ENV, "NDMC_CMIX_LO", "0.0")),
            parse(Float64, get(ENV, "NDMC_CO_LO", "0.0")),
        ),
        state_upper = (
            parse(Float64, get(ENV, "NDMC_C_HI", "600.0")),
            parse(Float64, get(ENV, "NDMC_C_HI", "600.0")),
            parse(Float64, get(ENV, "NDMC_C_HI", "600.0")),
            parse(Float64, get(ENV, "NDMC_CMIX_HI", "600.0")),
            parse(Float64, get(ENV, "NDMC_CO_HI", "20.0")),
        ),
        ipopt_tol = parse(Float64, get(ENV, "NDMC_IPOPT_TOL", "1e-6")),
        ipopt_max_cpu_time = parse(Float64, get(ENV, "NDMC_IPOPT_MAX_CPU_TIME", "60.0")),
        show_detailed_status = _env_flag("NDMC_SHOW_DETAILED_STATUS"; default="true"),
        show_generic_status = _env_flag("NDMC_SHOW_GENERIC_STATUS"),
        status_digits = parse(Int, get(ENV, "NDMC_STATUS_DIGITS", "3")),
        status_prefix = get(ENV, "NDMC_STATUS_PREFIX", "[NDMC MPC]"),
        status_stream = begin
            stream = Symbol(lowercase(strip(get(ENV, "NDMC_STATUS_STREAM", "stdout"))))
            stream in (:stdout, :stderr) || error("NDMC_STATUS_STREAM must be `stdout` or `stderr`.")
            stream
        end,
    )
end

function _ndmc_status_io(cfg::NDMCMPCConfig)
    if cfg.status_stream === :stdout
        return stdout
    elseif cfg.status_stream === :stderr
        return stderr
    end
    error("Unsupported NDMC status stream $(cfg.status_stream). Use :stdout or :stderr.")
end

function _fmt_ndmc_scalar(x; digits::Int=3)
    if x isa Real && isfinite(float(x))
        return @sprintf("%.*f", digits, float(x))
    end
    return "NaN"
end

function print_ndmc_mpc_status(io::IO,
                               cfg::NDMCMPCConfig,
                               sys,
                               result,
                               state_values::AbstractDict,
                               q_prev::Real,
                               q_apply::Real,
                               cin_now,
                               t_now::Real,
                               step::Integer)
    accepted = is_accepted_mpc_status(result.status; accepted_statuses=default_mpc_accepted_statuses())
    pred_c3 = get(result.predictions, sys.C3, Float64[])
    pred_c3_next = length(pred_c3) >= 2 ? pred_c3[2] : NaN
    pred_c3_terminal = isempty(pred_c3) ? NaN : pred_c3[end]

    parts = String[
        cfg.status_prefix,
        "step=" * string(step),
        "t=" * _fmt_ndmc_scalar(t_now; digits=cfg.status_digits),
        "status=" * string(result.status),
        "accepted=" * string(accepted),
        "obj=" * _fmt_ndmc_scalar(get(result.metrics, :objective, NaN); digits=cfg.status_digits),
        "C1=" * _fmt_ndmc_scalar(get(state_values, sys.C1, NaN); digits=cfg.status_digits) * " (sp=" * _fmt_ndmc_scalar(cfg.setpoint; digits=cfg.status_digits) * ")",
        "C2=" * _fmt_ndmc_scalar(get(state_values, sys.C2, NaN); digits=cfg.status_digits) * " (sp=" * _fmt_ndmc_scalar(cfg.setpoint; digits=cfg.status_digits) * ")",
        "C3=" * _fmt_ndmc_scalar(get(state_values, sys.C3, NaN); digits=cfg.status_digits) * " (sp=" * _fmt_ndmc_scalar(cfg.setpoint; digits=cfg.status_digits) * ")",
        "Q_prev=" * _fmt_ndmc_scalar(q_prev; digits=cfg.status_digits),
        "Q_apply=" * _fmt_ndmc_scalar(q_apply; digits=cfg.status_digits),
        "Cin1=" * _fmt_ndmc_scalar(cin_now[1]; digits=cfg.status_digits),
        "Cin2=" * _fmt_ndmc_scalar(cin_now[2]; digits=cfg.status_digits),
        "Cin3=" * _fmt_ndmc_scalar(cin_now[3]; digits=cfg.status_digits),
        "pred_C3_next=" * _fmt_ndmc_scalar(pred_c3_next; digits=cfg.status_digits),
        "pred_C3_terminal=" * _fmt_ndmc_scalar(pred_c3_terminal; digits=cfg.status_digits),
    ]

    println(io, string(first(parts), " ", join(parts[2:end], " | ")))
    flush(io)
    return nothing
end

# Five-state NDMC plant used by both the closed-loop simulation and the MPC
# prediction model. `C1`, `C2`, and `C3` are the tracked conductivities.
@mtkmodel NDMCPlant begin
    @parameters begin
        Cin1 = 285.0
        Cin2 = 285.0
        Cin3 = 320.0
        Q_air = 168.0
        k1 = 1.0
        k2 = 1.0
        k3 = 1.0
        k4 = 1.0
    end
    @variables begin
        C1(t) = 280.0
        C2(t) = 280.0
        C3(t) = 280.0
        Cmix(t) = 280.0
        cO(t) = 0.0
    end
    @equations begin
        D(C1) ~ (k1 * (Cmix - C1) + NDMC_FLOW * Cin1 - NDMC_FLOW * C1) / NDMC_ZONE_VOLUME - ndmc_conductivity_sink(cO)
        D(C2) ~ (k2 * (Cmix - C2) + NDMC_FLOW * Cin2 - NDMC_FLOW * C2) / NDMC_ZONE_VOLUME - ndmc_conductivity_sink(cO)
        D(C3) ~ (k3 * (Cmix - C3) + NDMC_FLOW * Cin3 - NDMC_FLOW * C3) / NDMC_ZONE_VOLUME - ndmc_conductivity_sink(cO)
        D(Cmix) ~ (k4 * (C1 + C2 + C3 - 3.0 * Cmix)) / NDMC_ZONE_VOLUME - ndmc_conductivity_sink(cO)
        D(cO) ~ ((NDMC_SOTE * 0.2967 * Q_air) / NDMC_COS / NDMC_VOLUME_TOTAL) * (NDMC_COS - cO) - ndmc_oxygen_sink(cO)
    end
end

# Build the plant with the nominal inflow and initial aeration.
function build_ndmc_system(cfg::NDMCMPCConfig)
    k = ndmc_k_vector()
    @mtkbuild sys = NDMCPlant(
        Cin1 = cfg.Cs,
        Cin2 = cfg.Cs,
        Cin3 = cfg.Cs,
        Q_air = cfg.Q_init,
        k1 = k[1],
        k2 = k[2],
        k3 = k[3],
        k4 = k[4],
    )
    return sys
end

# Initial measured state supplied to the first MPC solve.
function ndmc_initial_state(sys, cfg::NDMCMPCConfig)
    return Dict(
        sys.C1 => cfg.initial_state[1],
        sys.C2 => cfg.initial_state[2],
        sys.C3 => cfg.initial_state[3],
        sys.Cmix => cfg.initial_state[4],
        sys.cO => cfg.initial_state[5],
    )
end

# Influent disturbance used in the closed-loop benchmark.
function disturbance_triplet(t_now::Real, cfg::NDMCMPCConfig)
    if cfg.disturbance_start <= t_now < cfg.disturbance_stop
        return cfg.Cin_spike
    end
    return (cfg.Cs, cfg.Cs, cfg.Cs)
end

# ---------------------------------------------------------------------------
# MPC setup and one-step update
# ---------------------------------------------------------------------------

# Build one JuMP tracking-MPC problem for the NDMC plant.
#
# The controller tracks C1, C2, and C3 to the same setpoint and manipulates
# only Q_air. The Cin parameters are stored as stage-wise values so the
# disturbance preview can be refreshed before each solve.
function build_ndmc_controller(sys;
                               cfg::NDMCMPCConfig,
                               timing_profile::Union{Nothing, NDMCTimingProfile}=nothing)
    controller_start_ns = time_ns()
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", cfg.ipopt_tol)
    set_optimizer_attribute(model, "max_cpu_time", cfg.ipopt_max_cpu_time)
    set_optimizer_attribute(model, "print_level", 0)
    set_silent(model)

    mpc_cfg = TrackingMPCConfig(
        PH = cfg.PH,
        CH = cfg.CH,
        dt = cfg.dt,
        integrator = "IE",
        system_kind = :ode,
        state_lower = Dict(
            sys.C1 => cfg.state_lower[1],
            sys.C2 => cfg.state_lower[2],
            sys.C3 => cfg.state_lower[3],
            sys.Cmix => cfg.state_lower[4],
            sys.cO => cfg.state_lower[5],
        ),
        state_upper = Dict(
            sys.C1 => cfg.state_upper[1],
            sys.C2 => cfg.state_upper[2],
            sys.C3 => cfg.state_upper[3],
            sys.Cmix => cfg.state_upper[4],
            sys.cO => cfg.state_upper[5],
        ),
        rhs0 = 0.0,
        track_start_index = 2,
        lower_clip = 0.0,
    )

    build_tracking_mpc_start_ns = time_ns()
    ctrl = build_tracking_mpc(
        model,
        sys;
        control_specs = [
            MPCControlSpec(
                sym = sys.Q_air,
                lower = cfg.Q_min,
                upper = cfg.Q_max,
                delta_max = cfg.delta_max,
                move_weight = cfg.move_weight,
                first_move_weight = cfg.first_move_weight,
                base_name = "Q_air",
            ),
        ],
        output_specs = [
            MPCOutputSpec(sym = sys.C1, setpoint = cfg.setpoint, track_weight = 1.0, terminal_weight = cfg.terminal_weight, base_name = "C1"),
            MPCOutputSpec(sym = sys.C2, setpoint = cfg.setpoint, track_weight = 1.0, terminal_weight = cfg.terminal_weight, base_name = "C2"),
            MPCOutputSpec(sym = sys.C3, setpoint = cfg.setpoint, track_weight = 1.0, terminal_weight = cfg.terminal_weight, base_name = "C3"),
        ],
        stage_param_defaults = Dict(
            sys.Cin1 => fill(cfg.Cs, cfg.PH + 1),
            sys.Cin2 => fill(cfg.Cs, cfg.PH + 1),
            sys.Cin3 => fill(cfg.Cs, cfg.PH + 1),
        ),
        config = mpc_cfg,
    )

    if timing_profile !== nothing && timing_profile.enabled
        timing_profile.build_tracking_mpc_s = _elapsed_s(build_tracking_mpc_start_ns)
        timing_profile.build_controller_total_s = _elapsed_s(controller_start_ns)
    end

    return ctrl
end

# Write the current influent values into the ODE integrator.
function set_plant_disturbance!(integ, sys, cfg::NDMCMPCConfig, t_now::Real)
    cin1, cin2, cin3 = disturbance_triplet(t_now, cfg)
    integ.ps[sys.Cin1] = cin1
    integ.ps[sys.Cin2] = cin2
    integ.ps[sys.Cin3] = cin3
    return nothing
end

# Write the current or previewed influent values into the MPC prediction model.
function update_disturbance_forecast!(ctrl::TrackingMPCController, sys, cfg::NDMCMPCConfig, t_now::Real)
    cin_now = disturbance_triplet(t_now, cfg)
    for (sym, idx) in zip((sys.Cin1, sys.Cin2, sys.Cin3), 1:3)
        vals = if cfg.use_disturbance_forecast
            [disturbance_triplet(t_now + (k - 1) * cfg.dt, cfg)[idx] for k in 1:ctrl.N]
        else
            fill(cin_now[idx], ctrl.N)
        end
        update_stage_parameter!(ctrl, sym, vals)
    end
    return ctrl
end

# Run one receding-horizon update and return the aeration value to apply.
function solve_ndmc_step!(ctrl::TrackingMPCController,
                          sys,
                          cfg::NDMCMPCConfig,
                          state_values::AbstractDict,
                          q_prev::Real,
                          t_now::Real;
                          step::Union{Nothing, Integer}=nothing,
                          state_capture_s::Real=0.0,
                          timing_profile::Union{Nothing, NDMCTimingProfile}=nothing)
    step_start_ns = time_ns()
    forecast_update_start_ns = time_ns()
    update_disturbance_forecast!(ctrl, sys, cfg, t_now)
    forecast_update_s = _elapsed_s(forecast_update_start_ns)
    cin_now = disturbance_triplet(t_now, cfg)
    collect_mpc_timing = timing_profile !== nothing && timing_profile.enabled
    solve_tracking_mpc_start_ns = time_ns()
    result = solve_tracking_mpc!(
        ctrl,
        state_values,
        Dict(sys.Q_air => q_prev);
        lower_clip = 0.0,
        show_status = cfg.show_generic_status,
        status_io = _ndmc_status_io(cfg),
        status_time = t_now,
        status_output_syms = [sys.C1, sys.C2, sys.C3],
        status_control_syms = [sys.Q_air],
        status_digits = cfg.status_digits,
        status_prefix = cfg.status_prefix,
        return_timing = collect_mpc_timing,
    )
    solve_tracking_mpc_s = _elapsed_s(solve_tracking_mpc_start_ns)
    prepare_tracking_mpc_step_s = collect_mpc_timing ? result.timing.prepare_tracking_mpc_step_s : 0.0
    jump_optimize_s = collect_mpc_timing ? result.timing.jump_optimize_s : 0.0
    solve_postprocess_s = collect_mpc_timing ? result.timing.solve_postprocess_s : 0.0
    status_print_s = collect_mpc_timing ? result.timing.status_print_s : 0.0
    q_apply = result.controls[sys.Q_air][1]
    if abs(q_apply - q_prev) < 1e-8
        q_apply = q_prev
    end
    step_total_s = _elapsed_s(step_start_ns)

    if timing_profile !== nothing && timing_profile.enabled
        record_ndmc_step_timing!(
            timing_profile;
            step_index = something(step, -1),
            time = t_now,
            state_capture_s = state_capture_s,
            forecast_update_s = forecast_update_s,
            prepare_tracking_mpc_step_s = prepare_tracking_mpc_step_s,
            jump_optimize_s = jump_optimize_s,
            solve_postprocess_s = solve_postprocess_s,
            status_print_s = status_print_s,
            solve_tracking_mpc_s = solve_tracking_mpc_s,
            step_total_s = step_total_s,
            status = result.status,
            objective = get(result.metrics, :objective, NaN),
            q_prev = q_prev,
            q_apply = q_apply,
        )
    end

    if cfg.show_detailed_status
        print_ndmc_mpc_status(
            _ndmc_status_io(cfg),
            cfg,
            sys,
            result,
            state_values,
            q_prev,
            q_apply,
            cin_now,
            t_now,
            something(step, -1),
        )
    end
    return (
        q_apply = q_apply,
        result = result,
        timing = (
            state_capture_s = state_capture_s,
            forecast_update_s = forecast_update_s,
            prepare_tracking_mpc_step_s = prepare_tracking_mpc_step_s,
            jump_optimize_s = jump_optimize_s,
            solve_postprocess_s = solve_postprocess_s,
            status_print_s = status_print_s,
            solve_tracking_mpc_s = solve_tracking_mpc_s,
            step_total_s = step_total_s,
        ),
    )
end

# ---------------------------------------------------------------------------
# Closed-loop simulation
# ---------------------------------------------------------------------------

function metric_column(logctx::MPCLog, key::Symbol)
    return get(logctx.Metrichist, key, fill(NaN, length(logctx.pred_times)))
end

function ndmc_expected_total_solves(cfg::NDMCMPCConfig)::Int
    times_ctrl = collect((cfg.simulation_span[1] + cfg.dt):cfg.dt:cfg.simulation_span[2])
    return 1 + length(times_ctrl)
end

function ndmc_config_with(cfg::NDMCMPCConfig; kwargs...)
    names = fieldnames(NDMCMPCConfig)
    values = NamedTuple{names}(Tuple(getfield(cfg, name) for name in names))
    return NDMCMPCConfig(; values..., kwargs...)
end

function ndmc_offline_prewarm_config(cfg::NDMCMPCConfig)
    t0, t1 = cfg.simulation_span
    warmup_end = min(t1, t0 + cfg.dt)
    return ndmc_config_with(
        cfg;
        simulation_span = (t0, warmup_end),
        save_dt = max(eps(Float64), min(cfg.save_dt, max(cfg.dt, warmup_end - t0))),
        show_detailed_status = false,
        show_generic_status = false,
    )
end

function run_ndmc_offline_prewarm!(cfg::NDMCMPCConfig)::Float64
    warmup_start_ns = time_ns()
    run_ndmc_case(
        ndmc_offline_prewarm_config(cfg);
        timing_profile = nothing,
        save_outputs = false,
        announce_outputs = false,
        perform_offline_prewarm = false,
    )
    return _elapsed_s(warmup_start_ns)
end

"""
    run_ndmc_case(cfg; kwargs...)

Run one NDMC closed-loop simulation.

This is the main experiment function used by both the notebook and the script.
It builds the plant, builds the MPC problem, runs the ODE simulation with MPC
callbacks, and writes the standard CSV outputs when requested.
"""
function run_ndmc_case(cfg::NDMCMPCConfig;
                       timing_profile::Union{Nothing, NDMCTimingProfile}=nothing,
                       save_outputs::Bool=true,
                       announce_outputs::Bool=true,
                       perform_offline_prewarm::Bool=false,
                       output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR)
    profile = isnothing(timing_profile) ? NDMCTimingProfile(enabled=false) : timing_profile
    if perform_offline_prewarm
        announce_outputs && println("Running one offline NDMC warm-up pass before the saved case...")
        prewarm_elapsed_s = run_ndmc_offline_prewarm!(cfg)
        if profile.enabled
            profile.offline_prewarm_s = prewarm_elapsed_s
        end
        announce_outputs && println("Offline NDMC warm-up pass finished in $(round(prewarm_elapsed_s; digits=3)) s.")
    end
    total_start_ns = time_ns()

    build_system_start_ns = time_ns()
    sys = build_ndmc_system(cfg)
    if profile.enabled
        profile.build_system_s = _elapsed_s(build_system_start_ns)
    end
    state0 = ndmc_initial_state(sys, cfg)
    ctrl = build_ndmc_controller(sys; cfg=cfg, timing_profile=profile)
    mpc_step = Ref(0)

    initial = solve_ndmc_step!(
        ctrl,
        sys,
        cfg,
        state0,
        cfg.Q_init,
        cfg.simulation_span[1];
        step = mpc_step[],
        timing_profile = profile,
    )
    if profile.enabled
        profile.initial_mpc_step_total_s = initial.timing.step_total_s
    end
    mpc_step[] += 1
    initial_q = initial.q_apply
    initial_cin = disturbance_triplet(cfg.simulation_span[1], cfg)

    track_keys = Dict(sym => Symbol("track_" * ctrl.output_names[sym]) for sym in (sys.C1, sys.C2, sys.C3))
    move_key = Symbol("move_" * ctrl.control_names[sys.Q_air])
    move1_key = Symbol("move1_" * ctrl.control_names[sys.Q_air])

    logctx = make_mpc_log(
        sys;
        control_keys = [sys.Q_air],
        predicted_keys = [sys.C1, sys.C2, sys.C3, sys.Q_air],
        metric_keys = vcat(collect(values(track_keys)), [move_key, move1_key, :objective]),
    )
    record_mpc_prediction!(logctx, cfg.simulation_span[1], initial.result.predictions)
    record_mpc_metrics!(logctx, initial.result.metrics; missing=:skip)
    seed_mpc_log!(logctx, state0, cfg.simulation_span[1]; control_values=Dict(sys.Q_air => initial_q))

    p0 = Dict(
        sys.Q_air => initial_q,
        sys.Cin1 => initial_cin[1],
        sys.Cin2 => initial_cin[2],
        sys.Cin3 => initial_cin[3],
    )

    prob_cl = ODEProblem(sys, merge(state0, p0), cfg.simulation_span)

    times_ctrl = collect((cfg.simulation_span[1] + cfg.dt):cfg.dt:cfg.simulation_span[2])
    times_log = collect((cfg.simulation_span[1] + cfg.save_dt):cfg.save_dt:cfg.simulation_span[2])
    disturbance_times = unique(filter(t -> cfg.simulation_span[1] <= t <= cfg.simulation_span[2], [cfg.disturbance_start, cfg.disturbance_stop]))

    function mpc_affect!(integ)
        callback_start_ns = time_ns()
        set_plant_disturbance!(integ, sys, cfg, integ.t)
        state_capture_start_ns = time_ns()
        state_values = current_state_map(integ, sys)
        q_prev = integ.ps[sys.Q_air]
        state_capture_s = _elapsed_s(state_capture_start_ns)
        solved = solve_ndmc_step!(
            ctrl,
            sys,
            cfg,
            state_values,
            q_prev,
            integ.t;
            step = mpc_step[],
            state_capture_s = state_capture_s,
            timing_profile = profile,
        )
        mpc_step[] += 1
        q_apply = solved.q_apply
        integ.ps[sys.Q_air] = q_apply
        record_mpc_prediction!(logctx, integ.t, solved.result.predictions)
        record_mpc_metrics!(logctx, solved.result.metrics; missing=:skip)
        if profile.enabled
            profile.closed_loop_mpc_callback_s += _elapsed_s(callback_start_ns)
        end
        return nothing
    end

    function disturbance_affect!(integ)
        callback_start_ns = time_ns()
        set_plant_disturbance!(integ, sys, cfg, integ.t)
        if profile.enabled
            profile.closed_loop_disturbance_callback_s += _elapsed_s(callback_start_ns)
        end
        return nothing
    end

    function log_affect!(integ)
        callback_start_ns = time_ns()
        log_mpc_state!(logctx, integ; control_values=Dict(sys.Q_air => integ.ps[sys.Q_air]), missing_control=:hold)
        if profile.enabled
            profile.closed_loop_log_callback_s += _elapsed_s(callback_start_ns)
        end
        return nothing
    end

    mpc_cb = PresetTimeCallback(times_ctrl, mpc_affect!; save_positions=(false, false))
    dist_cb = isempty(disturbance_times) ? nothing : PresetTimeCallback(disturbance_times, disturbance_affect!; save_positions=(false, false))
    log_cb = PresetTimeCallback(times_log, log_affect!; save_positions=(false, false))
    cbs = dist_cb === nothing ? CallbackSet(mpc_cb, log_cb) : CallbackSet(dist_cb, mpc_cb, log_cb)

    solve_start_ns = time_ns()
    sol = solve(
        prob_cl,
        Rodas5P();
        adaptive = true,
        abstol = 1e-9,
        reltol = 1e-7,
        callback = cbs,
        tstops = sort(unique(vcat(times_ctrl, times_log, disturbance_times))),
        saveat = times_log,
    )
    if profile.enabled
        profile.closed_loop_simulation_s = _elapsed_s(solve_start_ns)
        profile.closed_loop_integrator_only_s = max(
            0.0,
            profile.closed_loop_simulation_s -
            profile.closed_loop_mpc_callback_s -
            profile.closed_loop_log_callback_s -
            profile.closed_loop_disturbance_callback_s,
        )
    end

    closed_loop_df = DataFrame(
        time = logctx.ts,
        C1 = logctx.Xhist[sys.C1],
        C2 = logctx.Xhist[sys.C2],
        C3 = logctx.Xhist[sys.C3],
        Cmix = logctx.Xhist[sys.Cmix],
        cO = logctx.Xhist[sys.cO],
        Q_air = logctx.Controlhist[sys.Q_air],
    )
    applied_df = DataFrame(
        time = logctx.pred_times,
        Q_applied = [traj[1] for traj in logctx.Predhist[sys.Q_air]],
        track_C1 = metric_column(logctx, track_keys[sys.C1]),
        track_C2 = metric_column(logctx, track_keys[sys.C2]),
        track_C3 = metric_column(logctx, track_keys[sys.C3]),
        move_Q_air = metric_column(logctx, move_key),
        move1_Q_air = metric_column(logctx, move1_key),
        objective = metric_column(logctx, :objective),
    )
    name_map_df = DataFrame(
        category = ["control", "output", "output", "output"],
        symbol = [string(sys.Q_air), string(sys.C1), string(sys.C2), string(sys.C3)],
        base_name = [
            ctrl.control_names[sys.Q_air],
            ctrl.output_names[sys.C1],
            ctrl.output_names[sys.C2],
            ctrl.output_names[sys.C3],
        ],
    )

    if profile.enabled
        profile.total_script_s = _elapsed_s(total_start_ns)
    end

    result = NDMCCaseResult(
        cfg,
        sys,
        ctrl,
        logctx,
        sol,
        closed_loop_df,
        applied_df,
        name_map_df,
        profile,
    )

    if save_outputs
        save_ndmc_case_outputs!(result; output_dir=output_dir, announce=announce_outputs)
    end

    return result
end

# ---------------------------------------------------------------------------
# Optional profiling and BenchmarkTools timings
# ---------------------------------------------------------------------------

"""
    run_ndmc_profile_cases(cfg; kwargs...)

Run the NDMC case with timing collection enabled.

This is for reporting wall-clock timing. It does not change the controller
formulation.
"""
function run_ndmc_profile_cases(cfg::NDMCMPCConfig;
                                profile_mode::Union{Symbol, AbstractString}=:single,
                                offline_prewarm::Bool=false,
                                save_outputs::Bool=true,
                                announce_outputs::Bool=true,
                                output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR)
    mode = _ndmc_profile_mode(profile_mode)
    expected_total_solves = ndmc_expected_total_solves(cfg)
    results = NDMCCaseResult[]

    if mode === :cold_and_warm && offline_prewarm
        println("Ignoring offline_prewarm=true because profile_mode=:cold_and_warm already measures the cold and warm paths separately.")
        offline_prewarm = false
    end

    if mode === :cold_and_warm
        push!(
            results,
            run_ndmc_case(
                cfg;
                timing_profile = NDMCTimingProfile(enabled=true, run_label="cold"),
                save_outputs = save_outputs,
                announce_outputs = false,
                perform_offline_prewarm = offline_prewarm,
                output_dir = output_dir,
            ),
        )
        push!(
            results,
            run_ndmc_case(
                cfg;
                timing_profile = NDMCTimingProfile(enabled=true, run_label="warm"),
                save_outputs = save_outputs,
                announce_outputs = announce_outputs,
                perform_offline_prewarm = offline_prewarm,
                output_dir = output_dir,
            ),
        )
    else
        push!(
            results,
            run_ndmc_case(
                cfg;
                timing_profile = NDMCTimingProfile(enabled=true, run_label="single"),
                save_outputs = save_outputs,
                announce_outputs = announce_outputs,
                perform_offline_prewarm = offline_prewarm,
                output_dir = output_dir,
            ),
        )
    end

    profiles = [result.profile for result in results]
    tables = ndmc_profile_dataframes(
        profiles,
        cfg;
        expected_total_solves = expected_total_solves,
        profile_mode = mode,
        offline_prewarm = offline_prewarm,
    )
    paths = save_outputs ?
        save_ndmc_profile_outputs!(tables.summary_df, tables.steps_df, tables.env_df; output_dir=output_dir, announce=announce_outputs) :
        nothing
    return (
        results = results,
        profiles = profiles,
        summary_df = tables.summary_df,
        steps_df = tables.steps_df,
        env_df = tables.env_df,
        paths = paths,
    )
end

ndmc_default_btime_snapshot_time(cfg::NDMCMPCConfig)::Float64 =
    min(cfg.simulation_span[2], cfg.simulation_span[1] + 2.0 * cfg.dt)

function _ndmc_snapshot_probe_config(cfg::NDMCMPCConfig, snapshot_time::Real)
    t0, t1 = cfg.simulation_span
    snap_t = clamp(Float64(snapshot_time), t0, t1)
    snap_dt = max(eps(Float64), snap_t - t0)
    return ndmc_config_with(
        cfg;
        simulation_span = (t0, snap_t),
        save_dt = snap_dt,
        show_detailed_status = false,
        show_generic_status = false,
    ), snap_t
end

function ndmc_btime_snapshot_values(cfg::NDMCMPCConfig;
                                    snapshot_time::Real=ndmc_default_btime_snapshot_time(cfg))
    probe_cfg, snap_t = _ndmc_snapshot_probe_config(cfg, snapshot_time)
    result = run_ndmc_case(
        probe_cfg;
        save_outputs = false,
        announce_outputs = false,
        perform_offline_prewarm = false,
    )
    snapshot_row = result.closed_loop_df[end, :]
    cin_now = disturbance_triplet(snap_t, cfg)
    return (
        snapshot_time = Float64(snapshot_row.time),
        C1 = Float64(snapshot_row.C1),
        C2 = Float64(snapshot_row.C2),
        C3 = Float64(snapshot_row.C3),
        Cmix = Float64(snapshot_row.Cmix),
        cO = Float64(snapshot_row.cO),
        Q_air = Float64(snapshot_row.Q_air),
        Cin1 = Float64(cin_now[1]),
        Cin2 = Float64(cin_now[2]),
        Cin3 = Float64(cin_now[3]),
    )
end

function _ndmc_snapshot_state_map(sys, snapshot)
    return Dict(
        sys.C1 => snapshot.C1,
        sys.C2 => snapshot.C2,
        sys.C3 => snapshot.C3,
        sys.Cmix => snapshot.Cmix,
        sys.cO => snapshot.cO,
    )
end

function _ndmc_snapshot_previous_controls(sys, snapshot)
    return Dict(sys.Q_air => snapshot.Q_air)
end

function _ndmc_trial_summary_row(stage::AbstractString,
                                 includes::AbstractString,
                                 excludes::AbstractString,
                                 closest_profile_quantity::AbstractString,
                                 trial::BenchmarkTools.Trial)
    times_ns = Float64.(trial.times)
    return (
        stage = String(stage),
        includes = String(includes),
        excludes = String(excludes),
        closest_profile_quantity = String(closest_profile_quantity),
        samples = length(times_ns),
        evals = trial.params.evals,
        minimum_s = minimum(times_ns) / 1.0e9,
        median_s = Statistics.median(times_ns) / 1.0e9,
        mean_s = Statistics.mean(times_ns) / 1.0e9,
        p95_s = Statistics.quantile(times_ns, 0.95) / 1.0e9,
        memory_bytes = trial.memory,
        allocs = trial.allocs,
    )
end

function _ndmc_run_benchmark(bench::BenchmarkTools.Benchmark, samples::Integer, evals::Integer, seconds::Real)
    bench.params.samples = Int(samples)
    bench.params.evals = Int(evals)
    bench.params.seconds = Float64(seconds)
    return BenchmarkTools.run(bench)
end

"""
    ndmc_btime_profile_tables(cfg; snapshot_time=..., samples=10, evals=1, seconds=2.0)

Run isolated BenchmarkTools timings for the NDMC example.

These are small hot-stage timings. They complement the closed-loop timing
profile and should not be interpreted as full end-to-end runtime.
"""
function ndmc_btime_profile_tables(cfg::NDMCMPCConfig;
                                   snapshot_time::Real=ndmc_default_btime_snapshot_time(cfg),
                                   samples::Integer=10,
                                   evals::Integer=1,
                                   seconds::Real=2.0)
    bench_cfg = ndmc_config_with(
        cfg;
        show_detailed_status = false,
        show_generic_status = false,
    )
    snapshot = ndmc_btime_snapshot_values(bench_cfg; snapshot_time=snapshot_time)
    snap_t = snapshot.snapshot_time
    snap_C1 = snapshot.C1
    snap_C2 = snapshot.C2
    snap_C3 = snapshot.C3
    snap_Cmix = snapshot.Cmix
    snap_cO = snapshot.cO
    snap_Q_air = snapshot.Q_air
    accepted_statuses = default_mpc_accepted_statuses()

    build_system_bench = BenchmarkTools.@benchmarkable begin
        build_ndmc_system($bench_cfg)
    end
    build_system_trial = _ndmc_run_benchmark(build_system_bench, samples, evals, seconds)

    sys_build = build_ndmc_system(bench_cfg)
    build_controller_bench = BenchmarkTools.@benchmarkable begin
        build_ndmc_controller($sys_build; cfg=$bench_cfg)
    end
    build_controller_trial = _ndmc_run_benchmark(build_controller_bench, samples, evals, seconds)

    forecast_update_bench = BenchmarkTools.@benchmarkable begin
        update_disturbance_forecast!(ctrl, sys, $bench_cfg, $snap_t)
    end setup=(
        sys = build_ndmc_system($bench_cfg);
        ctrl = build_ndmc_controller(sys; cfg=$bench_cfg);
    )
    forecast_update_trial = _ndmc_run_benchmark(forecast_update_bench, samples, evals, seconds)

    prepare_online_step_bench = BenchmarkTools.@benchmarkable begin
        prepare_tracking_mpc_step!(ctrl, state_values, previous_controls)
    end setup=(
        sys = build_ndmc_system($bench_cfg);
        ctrl = build_ndmc_controller(sys; cfg=$bench_cfg);
        state_values = Dict(
            sys.C1 => $snap_C1,
            sys.C2 => $snap_C2,
            sys.C3 => $snap_C3,
            sys.Cmix => $snap_Cmix,
            sys.cO => $snap_cO,
        );
        previous_controls = Dict(sys.Q_air => $snap_Q_air);
    )
    prepare_online_step_trial = _ndmc_run_benchmark(prepare_online_step_bench, samples, evals, seconds)

    optimize_prepared_problem_bench = BenchmarkTools.@benchmarkable begin
        JuMP.optimize!(ctrl.model)
        status = JuMP.termination_status(ctrl.model)
        status in $accepted_statuses || error("Prepared NDMC benchmark solve failed with status $(status).")
        nothing
    end setup=(
        sys = build_ndmc_system($bench_cfg);
        ctrl = build_ndmc_controller(sys; cfg=$bench_cfg);
        state_values = Dict(
            sys.C1 => $snap_C1,
            sys.C2 => $snap_C2,
            sys.C3 => $snap_C3,
            sys.Cmix => $snap_Cmix,
            sys.cO => $snap_cO,
        );
        previous_controls = Dict(sys.Q_air => $snap_Q_air);
        update_disturbance_forecast!(ctrl, sys, $bench_cfg, $snap_t);
        prepare_tracking_mpc_step!(ctrl, state_values, previous_controls);
    )
    optimize_prepared_problem_trial = _ndmc_run_benchmark(optimize_prepared_problem_bench, samples, evals, seconds)

    solve_tracking_step_bench = BenchmarkTools.@benchmarkable begin
        solve_tracking_mpc!(
            ctrl,
            state_values,
            previous_controls;
            lower_clip = 0.0,
            show_status = false,
            return_timing = false,
        )
    end setup=(
        sys = build_ndmc_system($bench_cfg);
        ctrl = build_ndmc_controller(sys; cfg=$bench_cfg);
        state_values = Dict(
            sys.C1 => $snap_C1,
            sys.C2 => $snap_C2,
            sys.C3 => $snap_C3,
            sys.Cmix => $snap_Cmix,
            sys.cO => $snap_cO,
        );
        previous_controls = Dict(sys.Q_air => $snap_Q_air);
        update_disturbance_forecast!(ctrl, sys, $bench_cfg, $snap_t);
    )
    solve_tracking_step_trial = _ndmc_run_benchmark(solve_tracking_step_bench, samples, evals, seconds)

    single_control_update_bench = BenchmarkTools.@benchmarkable begin
        solve_ndmc_step!(
            ctrl,
            sys,
            $bench_cfg,
            state_values,
            $snap_Q_air,
            $snap_t;
            step = 1,
            timing_profile = nothing,
        )
    end setup=(
        sys = build_ndmc_system($bench_cfg);
        ctrl = build_ndmc_controller(sys; cfg=$bench_cfg);
        state_values = Dict(
            sys.C1 => $snap_C1,
            sys.C2 => $snap_C2,
            sys.C3 => $snap_C3,
            sys.Cmix => $snap_Cmix,
            sys.cO => $snap_cO,
        );
    )
    single_control_update_trial = _ndmc_run_benchmark(single_control_update_bench, samples, evals, seconds)

    summary_df = DataFrame([
        _ndmc_trial_summary_row(
            "build dynamic model",
            "Build the MTK ODE model for C1, C2, C3, Cmix, and cO after the package code is already loaded.",
            "Building the MPC optimization problem, solving for Q_air, callbacks, and closed-loop simulation.",
            "MTK model build time only.",
            build_system_trial,
        ),
        _ndmc_trial_summary_row(
            "build controller and optimization problem",
            "Build the MPC optimization problem for future Q_air moves and predicted C1/C2/C3 trajectories on an already-built model.",
            "Dynamic-model build, online solves for Q_air, callbacks, and closed-loop simulation.",
            "Controller setup incl optimization model assembly.",
            build_controller_trial,
        ),
        _ndmc_trial_summary_row(
            "update future inflow values for Cin1/Cin2/Cin3",
            "Fill the prediction horizon with the future Cin1, Cin2, and Cin3 values seen by the MPC at one fixed simulation time.",
            "Reading the current plant state, writing state/warm-start data into the MPC, optimize!, reading back predictions/Q_air, and closed-loop progression.",
            "Online future-inflow update for Cin1/Cin2/Cin3 only.",
            forecast_update_trial,
        ),
        _ndmc_trial_summary_row(
            "write current state into MPC and refresh Q_air warm-start",
            "Write the fixed current C1/C2/C3/Cmix/cO values into the first MPC state point, shift the stored Q_air sequence when available, and refresh first-move bounds.",
            "Updating future Cin1/Cin2/Cin3, optimize!, reading back predictions/Q_air, simulation progression, and callbacks outside this preparation stage.",
            "Online write-current-state + refresh-Q_air-warm-start time only.",
            prepare_online_step_trial,
        ),
        _ndmc_trial_summary_row(
            "solve prepared optimization problem",
            "Run one nonlinear optimize! call after future Cin1/Cin2/Cin3 and the current C1/C2/C3/Cmix/cO state are already written into the MPC problem, then check solver status.",
            "Future-inflow update, writing the current state, Q_air warm-start refresh, reading back predictions/metrics, simulation progression, and callbacks outside optimize! itself.",
            "Online JuMP optimize! time only.",
            optimize_prepared_problem_trial,
        ),
        _ndmc_trial_summary_row(
            "compute next Q_air from one fixed current state",
            "Starting from one fixed C1/C2/C3/Cmix/cO state and previous Q_air value, run the local MPC solve step: write state, refresh Q_air warm-start, optimize!, and read back the next Q_air and predicted trajectories.",
            "Updating future Cin1/Cin2/Cin3 before the solve, reading the plant state before entering this function, and overall closed-loop simulation progression.",
            "Local MPC solve-step = write current state + refresh Q_air warm-start + optimize! + readback/status, excluding future-inflow update.",
            solve_tracking_step_trial,
        ),
        _ndmc_trial_summary_row(
            "full fixed-point MPC update for next Q_air",
            "For one fixed simulation point, update future Cin1/Cin2/Cin3 and then compute the next Q_air from the current C1/C2/C3/Cmix/cO state.",
            "Reading the plant state before entering this timing block, logging callbacks, inflow-shock event callbacks outside this timing block, and ODE integrator progression.",
            "Future-inflow update + local MPC solve-step for one fixed-point control update.",
            single_control_update_trial,
        ),
    ])

    summary_df.memory_mib = summary_df.memory_bytes ./ 2.0^20

    snapshot_df = DataFrame(
        snapshot_time = [snap_t],
        C1 = [snap_C1],
        C2 = [snap_C2],
        C3 = [snap_C3],
        Cmix = [snap_Cmix],
        cO = [snap_cO],
        Q_air = [snap_Q_air],
        Cin1 = [snapshot.Cin1],
        Cin2 = [snapshot.Cin2],
        Cin3 = [snapshot.Cin3],
        samples = [Int(samples)],
        evals = [Int(evals)],
        seconds_budget = [Float64(seconds)],
    )

    return (summary_df=summary_df, snapshot_df=snapshot_df)
end

_ndmc_pretty_run_label(label::AbstractString) = get(
    Dict(
        "cold" => "cold / first call",
        "warm" => "warm / second run",
        "prewarmed_main" => "main run after offline prewarm",
    ),
    label,
    label,
)

function _ndmc_sort_run_rows(df::DataFrame)
    run_order = Dict("cold" => 1, "warm" => 2, "prewarmed_main" => 3)
    :run_label in names(df) || return df
    ordered = copy(df)
    insertcols!(ordered, 1, :_run_order => [get(run_order, string(label), 99) for label in ordered.run_label])
    sort!(ordered, :_run_order)
    select!(ordered, Not(:_run_order))
    return ordered
end

function _ndmc_subset_run(df::DataFrame, label::AbstractString)
    view = filter(:run_label => ==(label), df)
    nrow(view) == 1 || error("Expected exactly one row for $(label), found $(nrow(view)).")
    return view
end

function _ndmc_stage_values(row)
    return [
        isfinite(row.offline_prewarm_s) ? row.offline_prewarm_s : 0.0,
        row.build_system_s,
        row.build_controller_total_s,
        row.build_tracking_mpc_s,
        row.initial_mpc_step_total_s,
        row.closed_loop_wall_s,
        row.total_script_s,
    ]
end

function _ndmc_profile_summary_view(df::DataFrame)
    view = _ndmc_sort_run_rows(select(
        df,
        :run_label,
        :offline_prewarm_s,
        :build_system_s,
        :build_controller_total_s,
        :build_tracking_mpc_s,
        :initial_mpc_step_total_s,
        :closed_loop_wall_s,
        :total_script_s,
        :expected_total_solves,
        :recorded_total_solves,
    ))
    transform!(view, :run_label => ByRow(_ndmc_pretty_run_label) => :run_label)
    rename!(view, Dict(
        :run_label => :mode,
        :offline_prewarm_s => Symbol("offline warm-up time before main run only (s)"),
        :build_system_s => Symbol("MTK model build time only (s)"),
        :build_controller_total_s => Symbol("controller setup time incl building the MPC problem for future Q_air and predicted C1/C2/C3 trajectories (s)"),
        :build_tracking_mpc_s => Symbol("optimization model assembly time only = create JuMP variables/constraints/objective for future Q_air and predicted C1/C2/C3 (subset of controller setup) (s)"),
        :initial_mpc_step_total_s => Symbol("initial t = 0 control update incl future Cin1/Cin2/Cin3 fill + write current C1/C2/C3/Cmix/cO + refresh Q_air warm-start + optimize! + read next Q_air/predictions (s)"),
        :closed_loop_wall_s => Symbol("closed-loop solve wall time incl ODE integration + MPC callback + logging callback + inflow-shock callback (s)"),
        :total_script_s => Symbol("total measured runtime incl model build + controller build + initial control update + closed-loop solve, excl offline warm-up (s)"),
        :expected_total_solves => Symbol("expected total NMPC solves"),
        :recorded_total_solves => Symbol("recorded total NMPC solves"),
    ))
    return view
end

function _ndmc_runtime_split_view(df::DataFrame)
    view = _ndmc_sort_run_rows(select(
        df,
        :run_label,
        :offline_prewarm_s,
        :build_system_s,
        :build_tracking_mpc_s,
        :initial_jump_solve_s,
        :online_jump_solve_s,
        :jump_solve_total_s,
        :solve_tracking_mpc_total_s,
        :closed_loop_integrator_only_s,
        :closed_loop_mpc_callback_s,
        :closed_loop_other_callback_s,
        :total_script_s,
        :end_to_end_total_including_offline_prewarm_s,
    ))
    transform!(view, :run_label => ByRow(_ndmc_pretty_run_label) => :run_label)
    rename!(view, Dict(
        :run_label => :mode,
        :offline_prewarm_s => Symbol("offline warm-up time before main run only (s)"),
        :build_system_s => Symbol("MTK model build time only (s)"),
        :build_tracking_mpc_s => Symbol("JuMP model assembly time only = create future-Q_air / predicted-C1-C2-C3 variables, constraints, and objective (subset of controller setup) (s)"),
        :initial_jump_solve_s => Symbol("initial JuMP optimize! time only at t = 0 after future Cin1/Cin2/Cin3 and current C1/C2/C3/Cmix/cO are already written (subset of initial control update) (s)"),
        :online_jump_solve_s => Symbol("online JuMP optimize! time total only across online steps after each local problem is already filled (s)"),
        :jump_solve_total_s => Symbol("total JuMP optimize! time only = initial optimize! + online optimize! calls (s)"),
        :solve_tracking_mpc_total_s => Symbol("total MPC solve-step time = write current C1/C2/C3/Cmix/cO + refresh Q_air warm-start + optimize! + read next Q_air/predictions, excl future Cin1/Cin2/Cin3 update (s)"),
        :closed_loop_integrator_only_s => Symbol("ODE integration time only = closed-loop wall - MPC callback - logging callback - inflow-shock callback (s)"),
        :closed_loop_mpc_callback_s => Symbol("online MPC callback wall time incl read current state + fill future Cin1/Cin2/Cin3 + write state + refresh Q_air warm-start + optimize! + readback + leftover callback work (s)"),
        :closed_loop_other_callback_s => Symbol("other callback wall time = plant-data logging callback + inflow-shock callback (s)"),
        :total_script_s => Symbol("measured runtime excl offline warm-up, incl model/controller build + initial control update + closed-loop solve (s)"),
        :end_to_end_total_including_offline_prewarm_s => Symbol("end-to-end runtime incl offline warm-up pass + measured runtime (s)"),
    ))
    return view
end

function _ndmc_callback_breakdown_view(df::DataFrame)
    view = _ndmc_sort_run_rows(select(
        df,
        :run_label,
        :online_state_capture_s,
        :online_forecast_update_s,
        :online_prepare_tracking_mpc_step_s,
        :online_mpc_preparation_s,
        :online_jump_solve_s,
        :online_solve_postprocess_s,
        :online_status_print_s,
        :closed_loop_mpc_callback_s,
        :mpc_callback_residual_s,
        :closed_loop_integrator_only_s,
    ))
    transform!(view, :run_label => ByRow(_ndmc_pretty_run_label) => :run_label)
    rename!(view, Dict(
        :run_label => :mode,
        :online_state_capture_s => Symbol("online read current plant state C1/C2/C3/Cmix/cO and current Q_air only (subset of MPC callback) (s)"),
        :online_forecast_update_s => Symbol("online fill future Cin1/Cin2/Cin3 values along the prediction horizon only (subset of MPC callback) (s)"),
        :online_prepare_tracking_mpc_step_s => Symbol("online write current C1/C2/C3/Cmix/cO into MPC and refresh Q_air warm-start / first-move bounds only (subset of MPC callback) (s)"),
        :online_mpc_preparation_s => Symbol("online MPC preparation total = read current state + fill future Cin1/Cin2/Cin3 + write state into MPC + refresh Q_air warm-start (s)"),
        :online_jump_solve_s => Symbol("online JuMP optimize! time only after future Cin1/Cin2/Cin3 and current state are already written (subset of MPC callback) (s)"),
        :online_solve_postprocess_s => Symbol("online readback time only = read next Q_air, predicted trajectories, and objective pieces after optimize! (subset of MPC callback) (s)"),
        :online_status_print_s => Symbol("online status-print time only = print one MPC log line after solve (subset of MPC callback) (s)"),
        :closed_loop_mpc_callback_s => Symbol("online MPC callback wall time incl read current state + fill future Cin1/Cin2/Cin3 + write state + refresh Q_air warm-start + optimize! + readback + leftover callback work (s)"),
        :mpc_callback_residual_s => Symbol("remaining MPC callback time not explained by the named pieces = callback wall - state read - future inflow fill - state/warm-start write - optimize! - readback - status print (s)"),
        :closed_loop_integrator_only_s => Symbol("ODE integration time only = closed-loop wall - MPC callback - logging callback - inflow-shock callback (s)"),
    ))
    return view
end

function _ndmc_online_step_stats(df::DataFrame)
    online_steps = filter(:step_index => >(0), df)
    stats = combine(
        groupby(online_steps, :run_label),
        nrow => :online_solve_count,
        :solve_tracking_mpc_s => mean => :mean_online_solve_s,
        :solve_tracking_mpc_s => median => :median_online_solve_s,
        :solve_tracking_mpc_s => (x -> quantile(x, 0.95)) => :p95_online_solve_s,
        :solve_tracking_mpc_s => maximum => :max_online_solve_s,
        :step_total_s => mean => :mean_online_step_total_s,
        :step_total_s => median => :median_online_step_total_s,
        :step_total_s => (x -> quantile(x, 0.95)) => :p95_online_step_total_s,
        :step_total_s => maximum => :max_online_step_total_s,
    )
    stats = _ndmc_sort_run_rows(stats)
    transform!(stats, :run_label => ByRow(_ndmc_pretty_run_label) => :run_label)
    rename!(stats, Dict(
        :run_label => :mode,
        :online_solve_count => Symbol("online NMPC solves"),
        :mean_online_solve_s => Symbol("mean online MPC solve-step time = write current state + refresh Q_air warm-start + optimize! + read next Q_air/predictions, excl future Cin1/Cin2/Cin3 update (s)"),
        :median_online_solve_s => Symbol("median online MPC solve-step time = write current state + refresh Q_air warm-start + optimize! + read next Q_air/predictions, excl future Cin1/Cin2/Cin3 update (s)"),
        :p95_online_solve_s => Symbol("p95 online MPC solve-step time = write current state + refresh Q_air warm-start + optimize! + read next Q_air/predictions, excl future Cin1/Cin2/Cin3 update (s)"),
        :max_online_solve_s => Symbol("max online MPC solve-step time = write current state + refresh Q_air warm-start + optimize! + read next Q_air/predictions, excl future Cin1/Cin2/Cin3 update (s)"),
        :mean_online_step_total_s => Symbol("mean online full MPC update time = fill future Cin1/Cin2/Cin3 + solve for next Q_air, excl reading plant state before this step (s)"),
        :median_online_step_total_s => Symbol("median online full MPC update time = fill future Cin1/Cin2/Cin3 + solve for next Q_air, excl reading plant state before this step (s)"),
        :p95_online_step_total_s => Symbol("p95 online full MPC update time = fill future Cin1/Cin2/Cin3 + solve for next Q_air, excl reading plant state before this step (s)"),
        :max_online_step_total_s => Symbol("max online full MPC update time = fill future Cin1/Cin2/Cin3 + solve for next Q_air, excl reading plant state before this step (s)"),
    ))
    return stats
end

function _ndmc_profile_comparison_view(profile_summary::DataFrame)
    cold_summary = _ndmc_subset_run(profile_summary, "cold")
    warm_summary = _ndmc_subset_run(profile_summary, "warm")
    prewarmed_summary = _ndmc_subset_run(profile_summary, "prewarmed_main")

    comparison_view = DataFrame(
        timing_item = [
            "offline warm-up pass",
            "MTK model build only",
            "controller setup incl building the MPC problem for future Q_air and predicted C1/C2/C3 trajectories",
            "optimization model assembly only = create future-Q_air / predicted-C1-C2-C3 variables, constraints, and objective (subset of controller setup)",
            "initial t = 0 control update incl future Cin1/Cin2/Cin3 fill + write current C1/C2/C3/Cmix/cO + refresh Q_air warm-start + optimize! + read next Q_air/predictions",
            "closed-loop solve wall time incl ODE integration + MPC callback + logging callback + inflow-shock callback",
            "total measured runtime incl model/controller build + initial control update + closed-loop solve, excl offline warm-up",
        ],
        cold_first_call_s = [
            0.0,
            cold_summary.build_system_s[1],
            cold_summary.build_controller_total_s[1],
            cold_summary.build_tracking_mpc_s[1],
            cold_summary.initial_mpc_step_total_s[1],
            cold_summary.closed_loop_wall_s[1],
            cold_summary.total_script_s[1],
        ],
        warm_second_run_s = [
            0.0,
            warm_summary.build_system_s[1],
            warm_summary.build_controller_total_s[1],
            warm_summary.build_tracking_mpc_s[1],
            warm_summary.initial_mpc_step_total_s[1],
            warm_summary.closed_loop_wall_s[1],
            warm_summary.total_script_s[1],
        ],
        main_run_after_offline_prewarm_s = [
            prewarmed_summary.offline_prewarm_s[1],
            prewarmed_summary.build_system_s[1],
            prewarmed_summary.build_controller_total_s[1],
            prewarmed_summary.build_tracking_mpc_s[1],
            prewarmed_summary.initial_mpc_step_total_s[1],
            prewarmed_summary.closed_loop_wall_s[1],
            prewarmed_summary.total_script_s[1],
        ],
    )
    comparison_view.estimated_first_call_overhead_s =
        comparison_view.cold_first_call_s .- comparison_view.warm_second_run_s
    rename!(comparison_view, Dict(
        :timing_item => Symbol("timing item"),
        :cold_first_call_s => Symbol("cold / first call (s)"),
        :warm_second_run_s => Symbol("warm / second run (s)"),
        :main_run_after_offline_prewarm_s => Symbol("main run after offline prewarm (s)"),
        :estimated_first_call_overhead_s => Symbol("estimated first-call overhead (cold - warm) (s)"),
    ))
    return comparison_view
end

function _ndmc_profile_budget_view(profile_summary::DataFrame)
    cold_summary = _ndmc_subset_run(profile_summary, "cold")
    warm_summary = _ndmc_subset_run(profile_summary, "warm")
    prewarmed_summary = _ndmc_subset_run(profile_summary, "prewarmed_main")
    budget_view = DataFrame(
        scenario = [
            "cold run total",
            "warm run total",
            "main run after offline prewarm",
            "offline warm-up + main run",
        ],
        seconds = [
            cold_summary.total_script_s[1],
            warm_summary.total_script_s[1],
            prewarmed_summary.total_script_s[1],
            prewarmed_summary.offline_prewarm_s[1] + prewarmed_summary.total_script_s[1],
        ],
    )
    rename!(budget_view, Dict(
        :scenario => Symbol("timing budget"),
        :seconds => Symbol("seconds"),
    ))
    return budget_view
end

function _ndmc_profile_env_view(coldwarm_env::DataFrame, prewarm_env::DataFrame)
    env_view = DataFrame(
        run_group = ["cold and warm run", "single run with offline prewarm"],
        timestamp = [coldwarm_env.timestamp[1], prewarm_env.timestamp[1]],
        julia_version = [coldwarm_env.julia_version[1], prewarm_env.julia_version[1]],
        julia_threads = [coldwarm_env.nthreads[1], prewarm_env.nthreads[1]],
        cpu_threads = [coldwarm_env.cpu_threads[1], prewarm_env.cpu_threads[1]],
        kernel = [coldwarm_env.kernel[1], prewarm_env.kernel[1]],
        machine = [coldwarm_env.machine[1], prewarm_env.machine[1]],
        architecture = [coldwarm_env.arch[1], prewarm_env.arch[1]],
        active_project = [coldwarm_env.active_project[1], prewarm_env.active_project[1]],
        profile_mode = [coldwarm_env.profile_mode[1], prewarm_env.profile_mode[1]],
        offline_prewarm_enabled = [coldwarm_env.offline_prewarm_enabled[1], prewarm_env.offline_prewarm_enabled[1]],
        run_count = [coldwarm_env.run_count[1], prewarm_env.run_count[1]],
    )
    rename!(env_view, Dict(
        :run_group => Symbol("run group"),
        :timestamp => Symbol("timestamp"),
        :julia_version => Symbol("Julia version"),
        :julia_threads => Symbol("Julia threads"),
        :cpu_threads => Symbol("CPU threads"),
        :kernel => Symbol("kernel"),
        :machine => Symbol("machine"),
        :architecture => Symbol("architecture"),
        :active_project => Symbol("active project"),
        :profile_mode => Symbol("profile mode"),
        :offline_prewarm_enabled => Symbol("offline prewarm enabled"),
        :run_count => Symbol("run count"),
    ))
    return env_view
end

function _ndmc_save_workflow_timing_views!(report; output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR)
    mkpath(output_dir)
    paths = (
        summary_view = joinpath(output_dir, "ndmc_conductivity_eopt_profile_summary_readable.csv"),
        runtime_split_view = joinpath(output_dir, "ndmc_conductivity_eopt_profile_runtime_split_readable.csv"),
        callback_breakdown_view = joinpath(output_dir, "ndmc_conductivity_eopt_profile_callback_breakdown_readable.csv"),
        comparison_view = joinpath(output_dir, "ndmc_conductivity_eopt_profile_comparison_readable.csv"),
        budget_view = joinpath(output_dir, "ndmc_conductivity_eopt_profile_budget_readable.csv"),
        online_step_stats = joinpath(output_dir, "ndmc_conductivity_eopt_profile_online_step_stats_readable.csv"),
        env_view = joinpath(output_dir, "ndmc_conductivity_eopt_profile_env_readable.csv"),
    )
    CSV.write(paths.summary_view, report.summary_view)
    CSV.write(paths.runtime_split_view, report.runtime_split_view)
    CSV.write(paths.callback_breakdown_view, report.callback_breakdown_view)
    CSV.write(paths.comparison_view, report.comparison_view)
    CSV.write(paths.budget_view, report.budget_view)
    CSV.write(paths.online_step_stats, report.online_step_stats)
    CSV.write(paths.env_view, report.env_view)
    return paths
end

function _ndmc_workflow_timing_plots(profile_summary::DataFrame,
                                     profile_steps::DataFrame,
                                     budget_view::DataFrame)
    cold_summary = _ndmc_subset_run(profile_summary, "cold")
    warm_summary = _ndmc_subset_run(profile_summary, "warm")
    prewarmed_summary = _ndmc_subset_run(profile_summary, "prewarmed_main")

    stages = [
        "offline\nwarm-up",
        "MTK model\nbuild",
        "controller\nsetup",
        "optimization model\nassembly",
        "initial NMPC solve\nat t = 0",
        "closed-loop solve wall\ntime incl callbacks",
        "total measured\nruntime",
    ]
    stage_matrix = hcat(
        _ndmc_stage_values(cold_summary[1, :]),
        _ndmc_stage_values(warm_summary[1, :]),
        _ndmc_stage_values(prewarmed_summary[1, :]),
    )

    stage_plot = bar(
        stages,
        stage_matrix;
        xlabel = "timing stage",
        ylabel = "seconds",
        label = [
            "cold / first call",
            "warm / second run",
            "main run after offline prewarm",
        ],
        title = "NDMC Timing Across Cold, Warm, and Offline Prewarm Modes",
    )
    budget_plot = bar(
        [
            "cold\nrun total",
            "warm\nrun total",
            "main run after\noffline prewarm",
            "offline warm-up\n+ main run",
        ],
        budget_view[!, Symbol("seconds")];
        xlabel = "timing budget",
        ylabel = "seconds",
        label = false,
        color = [:steelblue, :darkorange, :forestgreen, :firebrick],
        title = "End-to-End Timing Budgets",
    )

    run_colors = Dict(
        "cold" => :steelblue,
        "warm" => :darkorange,
        "prewarmed_main" => :forestgreen,
    )
    run_names = Dict(
        "cold" => "cold / first call",
        "warm" => "warm / second run",
        "prewarmed_main" => "main run after offline prewarm",
    )
    profile_plot = plot(
        xlabel = "MPC step index",
        ylabel = "MPC solve-step wall time (s)",
        title = "Per-Step MPC Solve-Step Wall Times",
    )
    for label in ["cold", "warm", "prewarmed_main"]
        run_steps = filter(:run_label => ==(label), profile_steps)
        sort!(run_steps, :step_index)
        plot!(
            profile_plot,
            run_steps.step_index,
            run_steps.solve_tracking_mpc_s;
            label = run_names[label],
            lw = 2,
            marker = :circle,
            ms = 4,
            color = run_colors[label],
        )
    end
    dashboard = plot(
        stage_plot,
        budget_plot,
        profile_plot;
        layout = (3, 1),
        size = (1100, 1200),
    )
    return (
        stage_plot = stage_plot,
        budget_plot = budget_plot,
        profile_plot = profile_plot,
        dashboard = dashboard,
    )
end

"""
    run_ndmc_workflow_timing_report(cfg=NDMCMPCConfig(); kwargs...)

Run the NDMC cold/warm/offline-prewarm timing profile and return the tables
used by the notebook. The raw timing CSVs and the readable reporting CSVs are
written when `save_outputs=true`.
"""
function run_ndmc_workflow_timing_report(cfg::NDMCMPCConfig=NDMCMPCConfig();
                                         output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR,
                                         save_outputs::Bool=true,
                                         announce_outputs::Bool=true)
    coldwarm_outputs = run_ndmc_profile_cases(
        cfg;
        profile_mode = :cold_and_warm,
        offline_prewarm = false,
        save_outputs = save_outputs,
        announce_outputs = announce_outputs,
        output_dir = output_dir,
    )
    prewarm_outputs = run_ndmc_profile_cases(
        cfg;
        profile_mode = :single,
        offline_prewarm = true,
        save_outputs = save_outputs,
        announce_outputs = announce_outputs,
        output_dir = output_dir,
    )

    prewarm_summary = copy(prewarm_outputs.summary_df)
    prewarm_steps = copy(prewarm_outputs.steps_df)
    prewarm_env = copy(prewarm_outputs.env_df)
    prewarm_summary.run_label .= "prewarmed_main"
    prewarm_steps.run_label .= "prewarmed_main"

    profile_summary = _ndmc_sort_run_rows(vcat(coldwarm_outputs.summary_df, prewarm_summary; cols=:union))
    profile_steps = _ndmc_sort_run_rows(vcat(coldwarm_outputs.steps_df, prewarm_steps; cols=:union))
    profile_env = vcat(coldwarm_outputs.env_df, prewarm_env; cols=:union)

    report_without_paths = (
        coldwarm_outputs = coldwarm_outputs,
        prewarm_outputs = prewarm_outputs,
        profile_summary = profile_summary,
        profile_steps = profile_steps,
        profile_env = profile_env,
        summary_view = _ndmc_profile_summary_view(profile_summary),
        runtime_split_view = _ndmc_runtime_split_view(profile_summary),
        callback_breakdown_view = _ndmc_callback_breakdown_view(profile_summary),
        comparison_view = _ndmc_profile_comparison_view(profile_summary),
        budget_view = _ndmc_profile_budget_view(profile_summary),
        online_step_stats = _ndmc_online_step_stats(profile_steps),
        env_view = _ndmc_profile_env_view(coldwarm_outputs.env_df, prewarm_env),
    )
    readable_paths = save_outputs ?
        _ndmc_save_workflow_timing_views!(report_without_paths; output_dir=output_dir) :
        nothing
    plots = _ndmc_workflow_timing_plots(
        report_without_paths.profile_summary,
        report_without_paths.profile_steps,
        report_without_paths.budget_view,
    )

    if announce_outputs && save_outputs
        println("Saved NDMC readable workflow timing tables to:")
        for path in readable_paths
            println("  " * path)
        end
    end

    return (
        report_without_paths...,
        readable_paths = readable_paths,
        plots = plots,
        profile_dashboard = plots.dashboard,
    )
end

function _ndmc_readable_btime_tables(btime_profile)
    btime_summary = select(
        btime_profile.summary_df,
        :stage,
        :includes,
        :excludes,
        :closest_profile_quantity,
        :samples,
        :evals,
        :minimum_s,
        :median_s,
        :mean_s,
        :p95_s,
        :memory_mib,
        :allocs,
    )
    rename!(btime_summary, Dict(
        :stage => Symbol("isolated stage"),
        :includes => Symbol("what this isolated timing includes in concrete NDMC terms"),
        :excludes => Symbol("what this isolated timing excludes from the NDMC run"),
        :closest_profile_quantity => Symbol("closest workflow-profile quantity"),
        :samples => Symbol("samples"),
        :evals => Symbol("evals per sample"),
        :minimum_s => Symbol("minimum time (s)"),
        :median_s => Symbol("median time (s)"),
        :mean_s => Symbol("mean time (s)"),
        :p95_s => Symbol("p95 time (s)"),
        :memory_mib => Symbol("memory (MiB)"),
        :allocs => Symbol("allocations"),
    ))

    btime_snapshot = copy(btime_profile.snapshot_df)
    rename!(btime_snapshot, Dict(
        :snapshot_time => Symbol("shared fixed snapshot time used for all isolated timings (s)"),
        :Q_air => Symbol("current control level at snapshot used as previous-control input"),
        :seconds_budget => Symbol("BenchmarkTools seconds budget"),
    ))
    return (btime_summary=btime_summary, btime_snapshot=btime_snapshot)
end

function _ndmc_workbook_pairs(workflow_report, btime_summary::DataFrame, btime_snapshot::DataFrame)
    pairs = Pair{String,DataFrame}[]
    if workflow_report !== nothing
        append!(
            pairs,
            [
                "workflow_summary_view" => workflow_report.summary_view,
                "workflow_runtime_view" => workflow_report.runtime_split_view,
                "workflow_callback_view" => workflow_report.callback_breakdown_view,
                "workflow_compare_view" => workflow_report.comparison_view,
                "workflow_budget_view" => workflow_report.budget_view,
                "workflow_online_steps" => workflow_report.online_step_stats,
                "workflow_env_view" => workflow_report.env_view,
                "workflow_raw_summary" => workflow_report.profile_summary,
                "workflow_raw_steps" => workflow_report.profile_steps,
                "workflow_raw_env" => workflow_report.profile_env,
            ],
        )
    end
    push!(pairs, "btime_summary" => btime_summary)
    push!(pairs, "btime_snapshot" => btime_snapshot)
    return pairs
end

"""
    run_ndmc_isolated_btime_report(cfg=NDMCMPCConfig(); kwargs...)

Run the isolated BenchmarkTools timings for separable NDMC stages. If a
workflow timing report is passed as `workflow_report`, the Excel workbook also
includes those workflow-profile tables.
"""
function run_ndmc_isolated_btime_report(cfg::NDMCMPCConfig=NDMCMPCConfig();
                                        snapshot_time::Real=ndmc_default_btime_snapshot_time(cfg),
                                        samples::Integer=10,
                                        evals::Integer=1,
                                        seconds::Real=2.0,
                                        output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR,
                                        workflow_report=nothing,
                                        save_outputs::Bool=true,
                                        save_workbook::Bool=true,
                                        announce_outputs::Bool=true)
    btime_profile = ndmc_btime_profile_tables(
        cfg;
        snapshot_time = snapshot_time,
        samples = samples,
        evals = evals,
        seconds = seconds,
    )
    readable = _ndmc_readable_btime_tables(btime_profile)

    btime_summary_path = joinpath(output_dir, "ndmc_conductivity_eopt_btime_summary.csv")
    btime_snapshot_path = joinpath(output_dir, "ndmc_conductivity_eopt_btime_snapshot.csv")
    if save_outputs
        mkpath(output_dir)
        CSV.write(btime_summary_path, readable.btime_summary)
        CSV.write(btime_snapshot_path, readable.btime_snapshot)
    end

    timing_workbook_path = nothing
    if save_outputs && save_workbook
        workbook_pairs = _ndmc_workbook_pairs(
            workflow_report,
            readable.btime_summary,
            readable.btime_snapshot,
        )
        timing_workbook_path = save_ndmc_timing_workbook!(
            workbook_pairs...;
            output_dir = output_dir,
            announce = false,
        )
    end

    if announce_outputs && save_outputs
        println("Saved NDMC isolated BenchmarkTools tables to:")
        println("  " * btime_summary_path)
        println("  " * btime_snapshot_path)
        if timing_workbook_path !== nothing
            println("Saved combined NDMC timing workbook to:")
            println("  " * timing_workbook_path)
        end
    end

    return (
        raw = btime_profile,
        summary = readable.btime_summary,
        snapshot = readable.btime_snapshot,
        summary_path = btime_summary_path,
        snapshot_path = btime_snapshot_path,
        workbook_path = timing_workbook_path,
    )
end

"""
    run_ndmc_whole_workflow_btime(cfg=NDMCMPCConfig(); samples=3, evals=1, seconds=30.0)

Benchmark the whole NDMC run as one warm BenchmarkTools target. This includes
model build, MPC problem build, initial MPC solve, closed-loop ODE solve, and
the MPC callback that writes the computed `Q_air` back into the plant state.
"""
function run_ndmc_whole_workflow_btime(cfg::NDMCMPCConfig=NDMCMPCConfig();
                                       samples::Integer=3,
                                       evals::Integer=1,
                                       seconds::Real=30.0)
    bench_cfg = ndmc_config_with(
        cfg;
        show_detailed_status = false,
        show_generic_status = false,
    )
    whole_workflow_bench = BenchmarkTools.@benchmarkable begin
        run_ndmc_case(
            $bench_cfg;
            save_outputs = false,
            announce_outputs = false,
            perform_offline_prewarm = false,
        )
    end
    whole_workflow_bench.params.samples = Int(samples)
    whole_workflow_bench.params.evals = Int(evals)
    whole_workflow_bench.params.seconds = Float64(seconds)
    whole_workflow_trial = BenchmarkTools.run(whole_workflow_bench)

    whole_workflow_times_s = whole_workflow_trial.times ./ 1.0e9
    return DataFrame(
        benchmark_scope = ["whole NDMC workflow as one warm run"],
        includes = ["Build MTK model + build MPC problem + initial MPC solve + closed-loop ODE solve + MPC callback, including writing computed Q_air back into the plant state."],
        excludes = ["First-run compilation overhead, offline prewarm, file writing, and console logging."],
        samples = [Int(samples)],
        evals_per_sample = [Int(evals)],
        minimum_time_s = [minimum(whole_workflow_times_s)],
        median_time_s = [median(whole_workflow_times_s)],
        mean_time_s = [mean(whole_workflow_times_s)],
        p95_time_s = [quantile(whole_workflow_times_s, 0.95)],
        memory_mib = [whole_workflow_trial.memory / 2.0^20],
        allocations = [whole_workflow_trial.allocs],
    )
end

"""
    main(; output_dir=NDMC_DEFAULT_OUTPUT_DIR)

Run the NDMC example from the script entry point.

Environment variables only change run settings such as horizon length,
profiling, and live status printing.
"""
function main(; output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR)
    cfg = load_config_from_env()
    profiling_enabled = ndmc_profiling_enabled()
    profile_mode_value = profiling_enabled ? ndmc_profile_mode() : :single
    offline_prewarm = ndmc_offline_prewarm_enabled()

    if profiling_enabled
        return run_ndmc_profile_cases(
            cfg;
            profile_mode = profile_mode_value,
            offline_prewarm = offline_prewarm,
            save_outputs = true,
            announce_outputs = true,
            output_dir = output_dir,
        )
    end

    return run_ndmc_case(
        cfg;
        save_outputs = true,
        announce_outputs = true,
        perform_offline_prewarm = offline_prewarm,
        output_dir = output_dir,
    )
end
