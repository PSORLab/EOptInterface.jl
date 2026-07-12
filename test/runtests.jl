using EOptInterface
using Test
using JuMP
using Ipopt
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

struct FakeNode
    name::Symbol
    children::Dict{Symbol, Any}
end

Base.getproperty(node::FakeNode, sym::Symbol) =
    sym === :name || sym === :children ? getfield(node, sym) : node.children[sym]

struct FakeIntegrator
    t::Float64
    u::Vector{Float64}
    ps::Dict{Any, Float64}
end

@testset "EOptInterface.jl" begin
    @testset "generic MPC helpers" begin
        outlet = FakeNode(Symbol("clarifier₊outlet_stream"), Dict(:S_N2O => :clf_n2o))
        reactor1 = FakeNode(:reactor1, Dict(:S_O => :r1_so, :S_NH => :r1_snh))
        clarifier = FakeNode(:clarifier, Dict(:outlet_stream => outlet))
        influent = FakeNode(:Influent, Dict(:COD => :influent_cod))
        sys = FakeNode(:sys, Dict(:reactor1 => reactor1, :clarifier => clarifier, :Influent => influent))

        vars = [
            "reactor1₊S_O(t)",
            "clarifier₊outlet_stream₊S_N2O(t)",
            "reactor1₊S_NH(t)",
            "Influent₊COD(t)",
        ]

        grouped = group_state_symbols(vars)
        @test grouped[:reactor1] == (:S_O, :S_NH)
        @test grouped[:clarifier] == (:S_N2O,)
        @test grouped[:Influent] == (:COD,)

        subsystems = collect_subsystems(sys, vars)
        @test subsystems == Any[reactor1, outlet, influent]
        @test pretty_subsystems(sys, subsystems) == [
            "sys.reactor1",
            "sys.clarifier.outlet_stream",
            "sys.Influent",
        ]

        model = Model()
        x_vars, c_ic = build_state_trajs_from_vars!(model, sys, vars, 4; lb=-1.0, ub=2.0, rhs0=0.5)
        @test length(x_vars) == 4
        @test haskey(c_ic, :r1_so)
        @test JuMP.lower_bound(x_vars[:r1_so][2]) == -1.0
        @test JuMP.upper_bound(x_vars[:r1_so][4]) == 2.0
        @test JuMP.normalized_rhs(c_ic[:r1_so]) == 0.5
        @test model.ext[:x_vars] === x_vars

        current_values = Dict(:r1_so => 0.5, :clf_n2o => 0.5)
        sync = check_ic_sync(current_values, c_ic; atol=1e-8)
        @test sync.ok
        @test sync.max_err == 0.0

        raw_ic = Dict(:r1_so => -1.0e-9, :r1_snh => 1.2)
        zero_small_values!(raw_ic; atol=1e-6)
        @test raw_ic[:r1_so] == 0.0

        @test zero_small_values!([1.0e-9, -2.0e-7, 0.3]; atol=1e-6) == [0.0, 0.0, 0.3]

        sync_state_trajectories!(c_ic, x_vars, Dict(:r1_so => -0.2, :r1_snh => 1.5, :clf_n2o => 0.9, :influent_cod => 0.1); lower_clip=0.0)
        @test JuMP.normalized_rhs(c_ic[:r1_so]) == 0.0
        @test JuMP.start_value(x_vars[:r1_snh][3]) == 1.5

        control_model = Model()
        @variable(control_model, u[1:4])
        warm_start_with_constant!(u, 3.0)
        @test all(JuMP.start_value(ui) == 3.0 for ui in u)
        shift_warm_start!(u, [10.0, 20.0, 30.0, 40.0])
        @test JuMP.start_value(u[1]) == 20.0
        @test JuMP.start_value(u[4]) == 40.0
        bounds = set_first_move_bounds!(u, 50.0; lower=0.0, upper=100.0, delta_max=15.0)
        @test bounds == (35.0, 65.0)
        @test JuMP.lower_bound(u[1]) == 35.0
        @test JuMP.upper_bound(u[1]) == 65.0

        idx_map = Dict(:x1 => 1, :x2 => 2)
        idx_integ = FakeIntegrator(0.0, [3.0, 4.0], Dict{Any, Float64}())
        @test current_state_map(idx_integ, idx_map) == Dict(:x1 => 3.0, :x2 => 4.0)

        @test default_mpc_accepted_statuses() == (
            JuMP.MOI.OPTIMAL,
            JuMP.MOI.LOCALLY_SOLVED,
            JuMP.MOI.ALMOST_OPTIMAL,
            JuMP.MOI.ALMOST_LOCALLY_SOLVED,
        )
        @test is_accepted_mpc_status(JuMP.MOI.LOCALLY_SOLVED)
        @test !is_accepted_mpc_status(JuMP.MOI.NUMERICAL_ERROR)

        summary_model = Model()
        @variable(summary_model, z >= 1.0)
        @objective(summary_model, Min, z^2)
        @test isnan(objective_value_or_nan(summary_model))
        unsolved = summarize_mpc_solve(summary_model)
        @test unsolved.status == JuMP.MOI.OPTIMIZE_NOT_CALLED
        @test !unsolved.accepted
        @test isnan(unsolved.objective)
    end

    @testset "constraint inspection helpers" begin
        model = Model()
        @variable(model, x)
        @variable(model, y)
        c_eq = @constraint(model, x + y == 3.0)
        c_ub = @constraint(model, x <= 1.0)
        @constraint(model, y >= 0.0)

        eq_hits = find_constraints_with_vars(model, [x]; verbose=false)
        @test eq_hits == JuMP.ConstraintRef[c_eq]

        mixed_hits = find_constraints_with_vars(
            model,
            [x];
            constraint_sets=(JuMP.MOI.EqualTo{Float64}, JuMP.MOI.LessThan{Float64}),
            verbose=false,
        )
        @test c_eq in mixed_hits
        @test c_ub in mixed_hits
    end

    @testset "conflict helpers" begin
        model = Model()
        @variable(model, 0.0 <= z <= 1.0)
        @constraint(model, z == 2.0)

        problems = find_conflicts(model)
        @test any(occursin("violates ub", problem) for problem in problems)

        stats = summarize_ic_uniqueness(model)
        @test stats.total_ic_vars == 0

        ic_model = Model()
        @variable(ic_model, q[1:2])
        get_ic_constraint!(ic_model, q; idx=1, rhs_if_missing=1.25)
        ic_stats = summarize_ic_uniqueness(ic_model)
        @test ic_stats.total_ic_vars == 1
        @test isempty(ic_stats.non_unique)
    end

    @testset "MPC logging helpers" begin
        state_vars = [:x1, :x2]
        logctx = make_mpc_log(state_vars; control_keys=[:u], predicted_keys=[:x1], metric_keys=[:obj])
        @test make_state_index_map(state_vars) == Dict(:x1 => 1, :x2 => 2)

        integ0 = FakeIntegrator(0.0, [1.0, 2.0], Dict(:u => 10.0))
        log_mpc_state!(logctx, integ0; control_values=Dict(:u => 10.0))
        record_mpc_prediction!(logctx, 0.0, Dict(:x1 => [1.0, 1.2]))
        record_mpc_metrics!(logctx, Dict(:obj => 5.0))

        integ1 = FakeIntegrator(1.0, [1.25, 2.2], Dict(:u => 11.0))
        log_mpc_state!(logctx, integ1; control_values=Dict(:u => 11.0))
        record_mpc_prediction!(logctx, 1.0, Dict(:x1 => [1.25, 1.3]))
        record_mpc_metrics!(logctx, Dict(:obj => 4.0))

        errs = compute_step_prediction_errors(logctx, :x1)
        @test errs.times == [1.0]
        @test errs.predicted == [1.2]
        @test errs.actual == [1.25]
        @test errs.errors[1] ≈ 0.05
        @test logctx.Controlhist[:u] == [10.0, 11.0]
        @test logctx.Metrichist[:obj] == [5.0, 4.0]
    end

    @testset "DMC helpers" begin
        model = Model()
        step_response = [0.35, 0.62, 0.8, 0.9, 0.95]
        ctx = register_dmcsystem(
            model,
            step_response,
            4;
            M = 2,
            y_meas0 = [0.4],
            u_hist0 = fill(0.3, length(step_response) - 1),
            u_bounds = (0.0, 1.0),
            base_name = "demo_dmc",
        )

        summary = summarize_dmc_registration(ctx)
        @test summary.N == 5
        @test summary.P == 4
        @test summary.M == 2
        @test summary.Ny == 1
        @test summary.u_count == 2
        @test summary.du_count == 1
    end

    @testset "tracking MPC builder helpers" begin
        ModelingToolkit.@parameters u = 0.0 d = 1.0
        ModelingToolkit.@variables x(t)
        @named sys = ODESystem([D(x) ~ -x + u + d], t, [x], [u, d])
        x_sym = only(collect(unknowns(sys)))

        model = Model()
        cfg = TrackingMPCConfig(
            PH = 3,
            CH = 2,
            dt = 1.0,
            integrator = "IE",
            system_kind = :ode,
            state_lower = -2.0,
            state_upper = 5.0,
            rhs0 = 0.0,
        )

        ctrl = build_tracking_mpc(
            model,
            sys;
            control_specs = [
                MPCControlSpec(
                    sym = sys.u,
                    lower = 0.0,
                    upper = 2.0,
                    delta_max = 0.4,
                    move_weight = 1.5,
                    first_move_weight = 2.0,
                ),
            ],
            output_specs = [
                MPCOutputSpec(
                    sym = x_sym,
                    setpoint = 1.2,
                    track_weight = 3.0,
                    terminal_weight = 4.0,
                    lower_soft = 0.8,
                    upper_soft = 1.4,
                    slack_weight = 100.0,
                ),
            ],
            stage_param_defaults = Dict(sys.d => [1.0, 1.1, 1.2, 1.3]),
            config = cfg,
        )

        u_key = only(collect(keys(ctrl.control_vars)))
        d_key = only(collect(keys(ctrl.stage_param_vars)))

        @test ctrl.N == 4
        @test ctrl.Nc == 2
        @test length(ctrl.control_vars[u_key]) == 4
        @test JuMP.fix_value(ctrl.stage_param_vars[d_key][2]) == 1.1
        @test JuMP.fix_value(ctrl.target_vars[x_sym]) == 1.2
        @test haskey(model.ext, :tracking_mpc)
        @test model.ext[:tracking_mpc] === ctrl

        update_stage_parameter!(ctrl, d_key, 2.0)
        @test all(JuMP.fix_value(v) == 2.0 for v in ctrl.stage_param_vars[d_key])

        update_tracking_targets!(
            ctrl;
            setpoints = Dict(x_sym => 1.5),
            lower_soft = Dict(x_sym => 1.0),
            upper_soft = Dict(x_sym => 1.8),
        )
        @test JuMP.fix_value(ctrl.target_vars[x_sym]) == 1.5
        @test JuMP.fix_value(ctrl.lower_ref_vars[x_sym]) == 1.0
        @test JuMP.fix_value(ctrl.upper_ref_vars[x_sym]) == 1.8

        prepare_tracking_mpc_step!(ctrl, Dict(x_sym => 0.25), Dict(u_key => 0.9))
        @test JuMP.normalized_rhs(ctrl.c_ic[x_sym]) == 0.25
        @test JuMP.normalized_rhs(ctrl.first_move_anchors[u_key]) == 0.9
        @test JuMP.lower_bound(ctrl.control_vars[u_key][1]) == 0.5
        @test JuMP.upper_bound(ctrl.control_vars[u_key][1]) == 1.3
        @test JuMP.start_value(ctrl.control_vars[u_key][1]) == 0.9

        integ = FakeIntegrator(0.0, [0.33], Dict(u_key => 0.9))
        cmap = current_state_map(integ, sys)
        @test cmap[x_sym] == 0.33

        logctx = make_mpc_log([x_sym]; control_keys=[u_key], predicted_keys=[x_sym], metric_keys=[:objective])
        seed_mpc_log!(logctx, Dict(x_sym => 0.33), 0.0; control_values=Dict(u_key => 0.9))
        @test logctx.ts == [0.0]
        @test logctx.Xhist[x_sym] == [0.33]
        @test logctx.Controlhist[u_key] == [0.9]
    end

    @testset "tracking MPC ODE integrators" begin
        ModelingToolkit.@parameters u_integrator = 0.0
        ModelingToolkit.@variables x_integrator(t) = 0.0
        ModelingToolkit.@named integrator_sys = ODESystem(
            [D(x_integrator) ~ -x_integrator + u_integrator],
            t,
            [x_integrator],
            [u_integrator],
        )
        x_key = only(collect(unknowns(integrator_sys)))

        for integrator in ("EE", "IE", "BDF1", "RK4", "IRK4")
            model = Model(Ipopt.Optimizer)
            set_silent(model)
            cfg = TrackingMPCConfig(
                PH = 2,
                CH = 2,
                dt = 0.2,
                integrator = integrator,
                system_kind = :ode,
                state_lower = -2.0,
                state_upper = 2.0,
                rhs0 = 0.0,
            )
            ctrl = build_tracking_mpc(
                model,
                integrator_sys;
                control_specs = [
                    MPCControlSpec(
                        sym = integrator_sys.u_integrator,
                        lower = -1.0,
                        upper = 1.0,
                        delta_max = 0.5,
                        move_weight = 0.1,
                    ),
                ],
                output_specs = [
                    MPCOutputSpec(
                        sym = x_key,
                        setpoint = 0.5,
                        track_weight = 1.0,
                    ),
                ],
                config = cfg,
            )
            result = solve_tracking_mpc!(
                ctrl,
                Dict(x_key => 0.0),
                Dict(integrator_sys.u_integrator => 0.0),
            )
            @test is_accepted_mpc_status(result.status)
            @test all(isfinite, result.controls[integrator_sys.u_integrator])
            @test all(isfinite, result.predictions[x_key])
        end
    end

    @testset "indexed MTK state helpers" begin
        ModelingToolkit.@parameters u_idx = 0.0
        ModelingToolkit.@variables x_idx(t)[1:2]
        @named idx_sys = ODESystem(
            [
                D(x_idx[1]) ~ -x_idx[1] + u_idx,
                D(x_idx[2]) ~ x_idx[1] - x_idx[2],
            ],
            t,
            collect(x_idx),
            [u_idx],
        )

        path_syms, state_sym = EOptInterface.split_mtk_state_path(idx_sys.x_idx[2])
        @test state_sym == :x_idx_2
        @test last(path_syms) == :idx_sys

        model = Model()
        x_vars_idx, c_ic_idx = EOptInterface.build_state_trajs_from_vars!(model, idx_sys, [idx_sys.x_idx[2]], 4)
        @test haskey(x_vars_idx, idx_sys.x_idx[2])
        @test haskey(c_ic_idx, idx_sys.x_idx[2])
        @test length(x_vars_idx[idx_sys.x_idx[2]]) == 4
    end

    @testset "legacy registration compatibility" begin
        @testset "legacy ode_model no-keyword registration" begin
            @mtkmodel KineticParameterEstimation begin
                @parameters begin
                    T = 273
                    K_2 = 46 * exp(6500 / T - 18)
                    K_3 = 2 * K_2
                    k_1 = 53
                    k_1s = k_1 * 1e-6
                    k_5 = 1.2e-3
                    c_O2 = 2e-3
                    k_2f
                    k_3f
                    k_4
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
                    D(x_A) ~ k_1 * x_Z * x_Y - c_O2 * (k_2f + k_3f) * x_A + k_2f / K_2 * x_D + k_3f / K_3 * x_B - k_5 * x_A^2
                    D(x_B) ~ c_O2 * k_3f * x_A - (k_3f / K_3 + k_4) * x_B
                    D(x_D) ~ c_O2 * k_2f * x_A - k_2f / K_2 * x_D
                    D(x_Y) ~ -k_1s * x_Z * x_Y
                    D(x_Z) ~ -k_1 * x_Z * x_Y
                    I ~ x_A + 2 / 21 * x_B + 2 / 21 * x_D
                end
            end

            @mtkcompile o = KineticParameterEstimation()
            tspan = (0.0, 0.1)
            tstep = 0.01
            N = Int(floor((tspan[2] - tspan[1]) / tstep)) + 1
            V = length(unknowns(o))

            model = Model()
            @variable(model, -75 <= z[1:V, 1:N] <= 150.0)
            pL = [10.0, 10.0, 0.001]
            pU = [1200.0, 1200.0, 40.0]
            @variable(model, pL[i] <= p[i = 1:3] <= pU[i])

            register_odesystem(model, o, tspan, tstep, "EE")

            @test JuMP.num_variables(model) == V * N + 3
            @test JuMP.num_constraints(model; count_variable_in_set_constraints = false) > 0
        end

        @testset "low-level dae registration" begin
            @parameters u = 0.0
            ModelingToolkit.@variables x_dae(t) = 0.0 z_dae(t) = 0.0
            @named dae_keep = System([
                D(x_dae) ~ -x_dae + z_dae + u,
                0 ~ z_dae + sin(z_dae) - x_dae,
            ], t)

            model = Model()
            dt_dae = 0.2
            tspan_dae = (0.0, 0.4)
            horizon = Int(floor((tspan_dae[2] - tspan_dae[1]) / dt_dae)) + 1
            ctx = decision_vars(dae_keep, Num[]; model = model, horizon = horizon, build_state_trajs = true, lb = -4.0, ub = 4.0)

            x_key = only(filter(sym -> endswith(string(sym), "x_dae(t)"), collect(keys(ctx.x_vars))))
            z_key = only(filter(sym -> endswith(string(sym), "z_dae(t)"), collect(keys(ctx.x_vars))))
            u_key = Num(only(filter(sym -> endswith(string(sym), "u"), collect(ModelingToolkit.parameters(dae_keep)))))

            x_traj = ctx.x_vars[x_key]
            z_traj = ctx.x_vars[z_key]
            @variable(model, -2.0 <= u_traj[1:horizon] <= 2.0)

            for var in x_traj
                JuMP.set_start_value(var, 0.0)
            end
            for var in z_traj
                JuMP.set_start_value(var, 0.0)
            end

            register_daesystem(
                model,
                dae_keep,
                tspan_dae,
                dt_dae,
                "IE";
                p_disc = [u_key],
                p_disc_vars = Dict(u_key => u_traj),
                x_vars = ctx.x_vars,
            )

            @test JuMP.num_constraints(model; count_variable_in_set_constraints = false) > 0
            @test length(x_traj) == horizon
            @test length(z_traj) == horizon
        end
    end
end
