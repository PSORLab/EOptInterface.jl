# NDMC conductivity MPC case.
#
# Reading path for the scientific example:
# 1. constants and rate expressions define the plant;
# 2. `NDMCPlant` gives the ModelingToolkit equations;
# 3. `build_ndmc_controller(...)` defines the MPC problem;
# 4. `run_ndmc_case(...)` simulates the closed loop.
#
# Optional profiling and BenchmarkTools reports are kept separately in
# `ndmc_profile.jl`.


# Values carried over from the legacy NDMC script so the comparison uses the
# same sample time and horizons.
const NDMC_LEGACY_DT = 20.0
const NDMC_LEGACY_PSPAN = 400.0
const NDMC_LEGACY_MSPAN = 60.0
const NDMC_LEGACY_PH = Int(round(NDMC_LEGACY_PSPAN / NDMC_LEGACY_DT))
const NDMC_LEGACY_CH = Int(round(NDMC_LEGACY_MSPAN / NDMC_LEGACY_DT))
const NDMC_DEFAULT_SAVE_DT = 10.0

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


function ndmc_case_output_paths(output_dir::AbstractString=NDMC_DEFAULT_OUTPUT_DIR)
    return (
        closed_loop = joinpath(output_dir, NDMC_CLOSED_LOOP_FILENAME),
        applied = joinpath(output_dir, NDMC_APPLIED_CONTROL_FILENAME),
        name_map = joinpath(output_dir, NDMC_NAME_MAP_FILENAME),
    )
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
    save_dt::Float64 = NDMC_DEFAULT_SAVE_DT
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
        save_dt = parse(Float64, get(ENV, "NDMC_SAVE_DT", string(NDMC_DEFAULT_SAVE_DT))),
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
