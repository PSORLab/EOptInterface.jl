# Optional timing, profiling, and reporting for the NDMC experiment.
# The plant model and closed-loop controller are defined in ndmc_case.jl.

using BenchmarkTools
using Statistics
using XLSX

const NDMC_PROFILE_SUMMARY_FILENAME = "ndmc_conductivity_eopt_profile_summary.csv"
const NDMC_PROFILE_STEPS_FILENAME = "ndmc_conductivity_eopt_profile_steps.csv"
const NDMC_PROFILE_ENV_FILENAME = "ndmc_conductivity_eopt_profile_env.csv"
const NDMC_TIMING_WORKBOOK_FILENAME = "ndmc_conductivity_eopt_timing_tables.xlsx"

function _ndmc_profile_mode(profile_mode)::Symbol
    mode = profile_mode isa Symbol ? profile_mode : Symbol(lowercase(strip(String(profile_mode))))
    mode in (:single, :cold_and_warm) ||
        error("NDMC profile mode must be `:single` or `:cold_and_warm`.")
    return mode
end
function ndmc_profile_mode()::Symbol
    return _ndmc_profile_mode(get(ENV, "NDMC_PROFILE_MODE", "single"))
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
