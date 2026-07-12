const EOI_COLORS = Dict(
    :signal => :steelblue3,
    :setpoint => :darkorange3,
    :error => :firebrick3,
    :mv1 => :forestgreen,
    :mv2 => :mediumpurple3,
    :guide => :gray30,
    :band => :gray55,
    :fill => :gray85,
)

function apply_eoi_publication_style!()
    default(
        fmt = :png,
        dpi = 220,
        size = (1000, 720),
        thickness_scaling = 1.15,
        framestyle = :box,
        grid = true,
        gridalpha = 0.16,
        gridlinewidth = 0.5,
        minorgrid = false,
        legend = :topright,
        legendfont = font(10, "Helvetica"),
        guidefont = font(12, "Helvetica"),
        tickfont = font(10, "Helvetica"),
        titlefont = font(13, "Helvetica"),
        background_color = :white,
        foreground_color_legend = nothing,
        background_color_legend = :transparent,
        palette = [
            EOI_COLORS[:signal],
            EOI_COLORS[:setpoint],
            EOI_COLORS[:error],
            EOI_COLORS[:mv1],
            EOI_COLORS[:mv2],
        ],
    )
    return nothing
end

const NDMC_PLOT_CONFIG = (
    eoi_color = :steelblue3,
    legacy_color = :black,
    c1_color = :steelblue3,
    c2_color = :darkorange3,
    c3_color = :firebrick3,
    cmix_color = :indianred3,
    co_color = :deepskyblue4,
    qair_color = :forestgreen,
    setpoint_color = :gray30,
    guide_color = :gray35,
    line_width = 2.6,
    legacy_width = 2.2,
    setpoint_width = 1.4,
    guide_width = 1.2,
    legend = :topright,
    shock_window = (2000.0, 4000.0),
    disturbance_window = (2100.0, 2250.0),
    shock_pair_size = (860, 640),
    comparison_pair_size = (980, 760),
    outputs_pair_size = (980, 760),
    support_pair_size = (980, 760),
    diff_pair_size = (980, 800),
)

const NDMC_NOTEBOOK_PLOT_CONFIG = NDMC_PLOT_CONFIG

_ndmc_var_title(var::AbstractString, suffix::AbstractString) = "$(var): $(suffix)"

function _ndmc_windowed(df::DataFrame, time_col::Symbol, tspan)
    tspan === nothing && return df
    tlo, thi = float(first(tspan)), float(last(tspan))
    mask = (Float64.(df[!, time_col]) .>= tlo) .& (Float64.(df[!, time_col]) .<= thi)
    return df[mask, :]
end

function _ndmc_add_disturbance_guides!(p, disturbance_window=NDMC_PLOT_CONFIG.disturbance_window)
    disturbance_window === nothing && return p
    tlo, thi = disturbance_window
    vline!(p, [tlo]; color=NDMC_PLOT_CONFIG.guide_color, linestyle=:dot, linewidth=NDMC_PLOT_CONFIG.guide_width, label="")
    vline!(p, [thi]; color=NDMC_PLOT_CONFIG.guide_color, linestyle=:dash, linewidth=NDMC_PLOT_CONFIG.guide_width, label="")
    return p
end

function _ndmc_finalize_plot!(p; xlabel="Time (s)", ylabel="", title="", legend=NDMC_PLOT_CONFIG.legend, size=nothing, grid=true, xlims=nothing, ylims=nothing, xticks=nothing, yticks=nothing)
    plot!(
        p;
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        legend = legend,
        grid = grid,
        gridalpha = 0.10,
        gridlinewidth = 0.45,
        minorgrid = false,
    )
    xlims === nothing || plot!(p; xlims=xlims)
    ylims === nothing || plot!(p; ylims=ylims)
    xticks === nothing || plot!(p; xticks=xticks)
    yticks === nothing || plot!(p; yticks=yticks)
    size === nothing || plot!(p; size=size)
    return p
end

function ndmc_plot_shock_window_c3(
    current_df::DataFrame,
    legacy_df::DataFrame;
    current_label="EOI-based NMPC",
    legacy_label="Legacy MPC",
    setpoint=280.0,
    time_col::Symbol=:time,
    current_col::Symbol=:C3,
    legacy_col::Symbol=:C3,
    tspan=NDMC_PLOT_CONFIG.shock_window,
    disturbance_window=NDMC_PLOT_CONFIG.disturbance_window,
    title="c3: shock window",
    legend=NDMC_PLOT_CONFIG.legend,
)
    current_view = _ndmc_windowed(current_df, time_col, tspan)
    legacy_view = _ndmc_windowed(legacy_df, time_col, tspan)
    p = plot(
        current_view[!, time_col],
        current_view[!, current_col];
        label = current_label,
        color = NDMC_PLOT_CONFIG.eoi_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
    )
    plot!(
        p,
        legacy_view[!, time_col],
        legacy_view[!, legacy_col];
        label = legacy_label,
        color = NDMC_PLOT_CONFIG.legacy_color,
        linewidth = NDMC_PLOT_CONFIG.legacy_width,
        linestyle = :dash,
    )
    hline!(
        p,
        [setpoint];
        label = "setpoint",
        color = NDMC_PLOT_CONFIG.setpoint_color,
        linewidth = NDMC_PLOT_CONFIG.setpoint_width,
        linestyle = :dash,
    )
    _ndmc_add_disturbance_guides!(p, disturbance_window)
    c3_vals = vcat(Float64.(current_view[!, current_col]), Float64.(legacy_view[!, legacy_col]), [float(setpoint)])
    c3_lo = floor((minimum(c3_vals) - 0.5) * 2) / 2
    c3_hi = ceil((maximum(c3_vals) + 0.5) * 2) / 2
    return _ndmc_finalize_plot!(
        p;
        ylabel = "c3 (uS/cm)",
        title = title,
        legend = legend,
        grid = false,
        xlims = tspan,
        ylims = (c3_lo, c3_hi),
        xticks = collect(2000:250:4000),
    )
end

function ndmc_plot_shock_window_qair(
    current_df::DataFrame,
    legacy_df::DataFrame;
    current_label="EOI-based NMPC",
    legacy_label="Legacy MPC",
    time_col::Symbol=:time,
    current_col::Symbol=:Q_air,
    legacy_col::Symbol=:Q_air,
    tspan=NDMC_PLOT_CONFIG.shock_window,
    disturbance_window=NDMC_PLOT_CONFIG.disturbance_window,
    title="Aeration input: shock window",
    legend=false,
)
    current_view = _ndmc_windowed(current_df, time_col, tspan)
    legacy_view = _ndmc_windowed(legacy_df, time_col, tspan)
    p = plot(
        current_view[!, time_col],
        current_view[!, current_col];
        label = current_label,
        color = NDMC_PLOT_CONFIG.eoi_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
        seriestype = :steppost,
    )
    plot!(
        p,
        legacy_view[!, time_col],
        legacy_view[!, legacy_col];
        label = legacy_label,
        color = NDMC_PLOT_CONFIG.legacy_color,
        linewidth = NDMC_PLOT_CONFIG.legacy_width,
        linestyle = :dash,
        seriestype = :steppost,
    )
    _ndmc_add_disturbance_guides!(p, disturbance_window)
    return _ndmc_finalize_plot!(
        p;
        ylabel = "Aeration input (mg/s)",
        title = title,
        legend = legend,
        grid = false,
        xlims = tspan,
        ylims = (-20, 840),
        xticks = collect(2000:250:4000),
        yticks = collect(0:200:800),
    )
end

function ndmc_plot_shock_window_pair(
    current_df::DataFrame,
    legacy_df::DataFrame;
    tspan=NDMC_PLOT_CONFIG.shock_window,
    disturbance_window=NDMC_PLOT_CONFIG.disturbance_window,
    current_label="EOI-based NMPC",
    legacy_label="Legacy MPC",
)
    p_c3 = ndmc_plot_shock_window_c3(
        current_df,
        legacy_df;
        tspan = tspan,
        disturbance_window = disturbance_window,
        current_label = current_label,
        legacy_label = legacy_label,
        legend = NDMC_PLOT_CONFIG.legend,
    )
    p_q = ndmc_plot_shock_window_qair(
        current_df,
        legacy_df;
        tspan = tspan,
        disturbance_window = disturbance_window,
        current_label = current_label,
        legacy_label = legacy_label,
        legend = false,
    )
    return plot(p_c3, p_q; layout=(2, 1), size=NDMC_PLOT_CONFIG.shock_pair_size)
end

function ndmc_plot_full_horizon_pair(
    df::DataFrame;
    time_col::Symbol=:time,
    disturbance_window=nothing,
    setpoint=280.0,
)
    p_outputs = ndmc_plot_closed_loop_outputs(
        df;
        time_col = time_col,
        disturbance_window = disturbance_window,
        setpoint = setpoint,
        title = "Conductivity outputs: full horizon",
    )
    p_input = ndmc_plot_closed_loop_input(
        df;
        time_col = time_col,
        disturbance_window = disturbance_window,
        title = "Aeration input: full horizon",
    )
    return plot(p_outputs, p_input; layout=(2, 1), size=NDMC_PLOT_CONFIG.outputs_pair_size)
end

function ndmc_plot_closed_loop_outputs(
    df::DataFrame;
    time_col::Symbol=:time,
    tspan=nothing,
    disturbance_window=nothing,
    setpoint=280.0,
    title="Closed-loop outputs",
)
    view_df = _ndmc_windowed(df, time_col, tspan)
    p = plot(
        view_df[!, time_col],
        view_df.C1;
        label = "c1",
        color = NDMC_PLOT_CONFIG.c1_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
    )
    plot!(p, view_df[!, time_col], view_df.C2; label="c2", color=NDMC_PLOT_CONFIG.c2_color, linewidth=NDMC_PLOT_CONFIG.line_width)
    plot!(p, view_df[!, time_col], view_df.C3; label="c3", color=NDMC_PLOT_CONFIG.c3_color, linewidth=NDMC_PLOT_CONFIG.line_width)
    hline!(p, [setpoint]; label="setpoint", color=NDMC_PLOT_CONFIG.setpoint_color, linewidth=NDMC_PLOT_CONFIG.setpoint_width, linestyle=:dash)
    _ndmc_add_disturbance_guides!(p, disturbance_window)
    return _ndmc_finalize_plot!(p; ylabel="Conductivity (uS/cm)", title=title, legend=NDMC_PLOT_CONFIG.legend)
end

function ndmc_plot_closed_loop_input(
    df::DataFrame;
    time_col::Symbol=:time,
    tspan=nothing,
    disturbance_window=nothing,
    title="Applied aeration input",
)
    view_df = _ndmc_windowed(df, time_col, tspan)
    p = plot(
        view_df[!, time_col],
        view_df.Q_air;
        label = "Aeration input",
        color = NDMC_PLOT_CONFIG.qair_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
        seriestype = :steppost,
    )
    _ndmc_add_disturbance_guides!(p, disturbance_window)
    return _ndmc_finalize_plot!(p; ylabel="Aeration input (mg/s)", title=title, legend=NDMC_PLOT_CONFIG.legend)
end

function ndmc_plot_supporting_states(
    df::DataFrame;
    time_col::Symbol=:time,
    tspan=nothing,
    disturbance_window=nothing,
)
    view_df = _ndmc_windowed(df, time_col, tspan)
    p_mix = plot(
        view_df[!, time_col],
        view_df.Cmix;
        label = "cmix",
        color = NDMC_PLOT_CONFIG.cmix_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
    )
    _ndmc_add_disturbance_guides!(p_mix, disturbance_window)
    _ndmc_finalize_plot!(p_mix; ylabel="Mixing-zone conductivity (uS/cm)", title="Mixing-zone conductivity state", legend=NDMC_PLOT_CONFIG.legend)

    p_o = plot(
        view_df[!, time_col],
        view_df.cO;
        label = "cO",
        color = NDMC_PLOT_CONFIG.co_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
    )
    _ndmc_add_disturbance_guides!(p_o, disturbance_window)
    _ndmc_finalize_plot!(p_o; ylabel="Dissolved oxygen (mg/L)", title="Dissolved oxygen state", legend=NDMC_PLOT_CONFIG.legend)

    return plot(p_mix, p_o; layout=(2, 1), size=NDMC_PLOT_CONFIG.support_pair_size)
end

function ndmc_plot_legacy_comparison_pair(
    current_df::DataFrame,
    legacy_df::DataFrame;
    tspan=nothing,
    disturbance_window=nothing,
    title_suffix="full horizon",
    current_label="EOI-based MPC",
    legacy_label="Legacy DMC",
)
    p_c3 = ndmc_plot_shock_window_c3(
        current_df,
        legacy_df;
        tspan = tspan,
        disturbance_window = disturbance_window,
        current_label = current_label,
        legacy_label = legacy_label,
        title = _ndmc_var_title("c3", title_suffix),
        legend = NDMC_PLOT_CONFIG.legend,
    )
    p_q = ndmc_plot_shock_window_qair(
        current_df,
        legacy_df;
        tspan = tspan,
        disturbance_window = disturbance_window,
        current_label = current_label,
        legacy_label = legacy_label,
        title = _ndmc_var_title("Aeration input", title_suffix),
        legend = false,
    )
    return plot(p_c3, p_q; layout=(2, 1), size=NDMC_PLOT_CONFIG.comparison_pair_size)
end

function ndmc_plot_difference_pair(
    comparison_df::DataFrame;
    time_col::Symbol=:time,
    disturbance_window=nothing,
)
    p_states = plot(
        comparison_df[!, time_col],
        comparison_df.diff_C1;
        label = "c1 diff.",
        color = NDMC_PLOT_CONFIG.c1_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
    )
    plot!(p_states, comparison_df[!, time_col], comparison_df.diff_C2; label="c2 diff.", color=NDMC_PLOT_CONFIG.c2_color, linewidth=NDMC_PLOT_CONFIG.line_width)
    plot!(p_states, comparison_df[!, time_col], comparison_df.diff_C3; label="c3 diff.", color=NDMC_PLOT_CONFIG.c3_color, linewidth=NDMC_PLOT_CONFIG.line_width)
    plot!(p_states, comparison_df[!, time_col], comparison_df.diff_Cmix; label="cmix diff.", color=NDMC_PLOT_CONFIG.cmix_color, linewidth=NDMC_PLOT_CONFIG.line_width)
    _ndmc_add_disturbance_guides!(p_states, disturbance_window)
    _ndmc_finalize_plot!(p_states; ylabel="Conductivity diff. (uS/cm)", title="Conductivity-state differences", legend=NDMC_PLOT_CONFIG.legend)

    p_q = plot(
        comparison_df[!, time_col],
        comparison_df.diff_Q;
        label = "Aeration input diff.",
        color = NDMC_PLOT_CONFIG.qair_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
        seriestype = :steppost,
    )
    _ndmc_add_disturbance_guides!(p_q, disturbance_window)
    _ndmc_finalize_plot!(p_q; ylabel="Aeration-input diff. (mg/s)", title="Aeration-input difference", legend=:topright)

    p_o = plot(
        comparison_df[!, time_col],
        comparison_df.diff_cO;
        label = "cO diff.",
        color = NDMC_PLOT_CONFIG.co_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
    )
    _ndmc_add_disturbance_guides!(p_o, disturbance_window)
    _ndmc_finalize_plot!(p_o; ylabel="Dissolved-oxygen diff. (mg/L)", title="Dissolved-oxygen difference", legend=:topright)

    return plot(p_states, p_q, p_o; layout=(3, 1), size=(980, 1120))
end

function ndmc_plot_penalty_terms(
    df::DataFrame;
    time_col::Symbol=:time,
    disturbance_window=nothing,
)
    p_track = plot(
        df[!, time_col],
        df.track_C1;
        label = "c1 tracking",
        color = NDMC_PLOT_CONFIG.c1_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
        seriestype = :steppost,
    )
    plot!(p_track, df[!, time_col], df.track_C2; label="c2 tracking", color=NDMC_PLOT_CONFIG.c2_color, linewidth=NDMC_PLOT_CONFIG.line_width, seriestype=:steppost)
    plot!(p_track, df[!, time_col], df.track_C3; label="c3 tracking", color=NDMC_PLOT_CONFIG.c3_color, linewidth=NDMC_PLOT_CONFIG.line_width, seriestype=:steppost)
    _ndmc_add_disturbance_guides!(p_track, disturbance_window)
    _ndmc_finalize_plot!(p_track; ylabel="Tracking term", title="Tracking penalty terms", legend=NDMC_PLOT_CONFIG.legend)

    p_move = plot(
        df[!, time_col],
        df.move_Q_air;
        label = "move penalty",
        color = NDMC_PLOT_CONFIG.qair_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
        seriestype = :steppost,
    )
    plot!(p_move, df[!, time_col], df.move1_Q_air; label="first-move penalty", color=NDMC_PLOT_CONFIG.guide_color, linewidth=NDMC_PLOT_CONFIG.legacy_width, linestyle=:dash, seriestype=:steppost)
    _ndmc_add_disturbance_guides!(p_move, disturbance_window)
    _ndmc_finalize_plot!(p_move; ylabel="Move term", title="Move penalty terms", legend=NDMC_PLOT_CONFIG.legend)

    p_obj = plot(
        df[!, time_col],
        df.objective;
        label = "objective",
        color = NDMC_PLOT_CONFIG.legacy_color,
        linewidth = NDMC_PLOT_CONFIG.line_width,
        seriestype = :steppost,
    )
    _ndmc_add_disturbance_guides!(p_obj, disturbance_window)
    _ndmc_finalize_plot!(p_obj; ylabel="Objective value", title="Total objective", legend=:topright)

    return plot(p_track, p_move, p_obj; layout=(3, 1), size=(980, 1080))
end
