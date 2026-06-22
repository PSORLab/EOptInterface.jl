using EOptInterface
using Ipopt
using JuMP
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SciCompDSL
using Test

@testset "Algebraic Model" begin

    @connector Stream begin
        @variables begin
            F(t),   [input=true]
            y_A(t), [input=true]
            y_B(t), [input=true]
            y_C(t), [input=true]
        end
        @parameters begin
            V_A = 8.937e-2
            V_B = 1.018e-1
            V_C = 1.13e-1
        end
    end
    @mtkmodel Influent begin
        @components begin
            out = Stream()
        end
        @parameters begin
            F
            y_A = 1.0
            y_B = 0.0
            y_C = 0.0
        end
        @equations begin
            out.F ~ F
            out.y_A ~ y_A
            out.y_B ~ y_B
            out.y_C ~ y_C
        end
    end
    @mtkmodel Mixer begin
        @components begin
            in1 = Stream()
            in2 = Stream()
            out = Stream()
        end
        @equations begin
            out.F ~ in1.F + in2.F
            out.y_A ~ (in1.y_A*in1.F + in2.y_A*in2.F)/(in1.F + in2.F)
            out.y_B ~ (in1.y_B*in1.F + in2.y_B*in2.F)/(in1.F + in2.F)
            out.y_C ~ (in1.y_C*in1.F + in2.y_C*in2.F)/(in1.F + in2.F)
        end
    end
    @mtkmodel CSTR begin
        @components begin
            in = Stream()
            out = Stream()
        end
        @parameters begin
            V
            k_1 = 0.4
            k_2 = 0.055
        end
        begin
            r_1 = k_1*out.y_A/(out.y_A*in.V_A + out.y_B*in.V_B + out.y_C*in.V_C)
            r_2 = k_2*out.y_B/(out.y_A*in.V_A + out.y_B*in.V_B + out.y_C*in.V_C)
        end
        @equations begin
            out.F ~ in.F
            out.y_A + out.y_B + out.y_C ~ 1.0
            out.y_B*out.F ~ in.y_B*in.F + (r_1 - r_2)*V
            out.y_C*out.F ~ in.y_C*in.F + r_2*V
        end
    end
    @mtkmodel Separator1 begin
        @components begin
            in = Stream()
            outV = Stream()
            outL = Stream()
        end
        @equations begin
            in.F ~ outV.F + outL.F
            in.y_B*in.F ~ outL.y_B*outL.F
            in.y_C*in.F ~ outL.y_C*outL.F
            outV.y_A + outV.y_B + outV.y_C ~ 1.0
            outV.y_C ~ 0.0
            outV.y_B ~ 0.0
            outL.y_A + outL.y_B + outL.y_C ~ 1.0
            outL.y_A ~ 0.0
        end
    end
    @mtkmodel Separator2 begin
        @components begin
            in = Stream()
            outV = Stream()
            outL = Stream()
        end
        @equations begin
            in.F ~ outV.F + outL.F
            in.y_B*in.F ~ outV.F
            outV.y_A + outV.y_B + outV.y_C ~ 1.0
            outV.y_A ~ 0.0
            outV.y_C ~ 0.0
            outL.y_A + outL.y_B + outL.y_C ~ 1.0
            outL.y_A ~ 0.0
            outL.y_B ~ 0.0
        end
    end
    @mtkmodel ReactorSeparatorRecycle begin
        @components begin
            influent = Influent()
            mixer = Mixer()
            cstr = CSTR()
            sep1 = Separator1()
            sep2 = Separator2()
        end
        @equations begin
            connect(influent.out, mixer.in1)
            connect(mixer.out, cstr.in)
            connect(cstr.out, sep1.in)
            connect(sep1.outV, mixer.in2)
            connect(sep1.outL, sep2.in)
        end
    end

    @mtkcompile system = ReactorSeparatorRecycle()

    exprF5 = system.sep2.outV.F
    exprTau = system.cstr.V/(system.cstr.out.F*(system.cstr.out.y_A*system.cstr.in.V_A + system.cstr.out.y_B*system.cstr.in.V_B + system.cstr.out.y_C*system.cstr.in.V_C))
    f_CSTR = (25764.0 + 8178.0*system.cstr.V)/2.5
    s1cap = 132718.0 + system.cstr.out.F*(369.0*system.cstr.out.y_A - 1113.9*system.cstr.out.y_B)
    s2cap = 25000.0 + system.sep1.outL.F*(6984.5*system.sep1.outL.y_B - 3869.53*system.sep1.outL.y_C^2)
    s1op = system.cstr.out.F*(3.0 + 36.11*system.cstr.out.y_A + 7.71*system.cstr.out.y_B)*26.32e-3
    s2op = system.sep1.outL.F*(26.21 + 29.45*system.sep1.outL.y_B)*26.32e-3
    f_Sep = (s1cap + s2cap)/2.5 + 0.52*(s1op + s2op)
    g1 = 25.0 - exprF5
    g2 = 475/3600 - exprTau
    obj = f_CSTR + f_Sep

    model = JuMP.Model(Ipopt.Optimizer)
    EOptInterface.decision_vars(system)
    xL = zeros(6)
    xU = [100.0, 1.0, 1.0, 1.0, 100.0, 10.0]
    JuMP.@variable(model, xL[i] <= x[i=1:6] <= xU[i])
    EOptInterface.register_nlsystem(model, system, obj, [g1, g2])
    JuMP.optimize!(model)
    soln_dict = EOptInterface.full_solution(model, system)

    @test JuMP.termination_status(model) == JuMP.LOCALLY_SOLVED
    @test JuMP.primal_status(model) == JuMP.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(model), 169869.99931631665, atol=1e-3)

    model = JuMP.Model(Ipopt.Optimizer)
    JuMP.@variable(model, xL[i] <= x[i=1:6] <= xU[i])
    EOptInterface.register_nlsystem(model, system, obj)

end

@testset "ODE Model" begin

    @mtkmodel KineticParameterEstimation begin
        @parameters begin
            T = 273.0
            K_2 = 46.0*exp(6500.0/T - 18.0)
            K_3 = 2.0*K_2
            k_1 = 53.0
            k_1s = k_1*1e-6
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
            D(x_A) ~ k_1*x_Z*x_Y - c_O2*(k_2f + k_3f)*x_A + k_2f/K_2*x_D + k_3f/K_3*x_B - k_5*x_A^2
            D(x_B) ~ c_O2*k_3f*x_A - (k_3f/K_3 + k_4)*x_B
            D(x_D) ~ c_O2*k_2f*x_A - k_2f/K_2*x_D
            D(x_Y) ~ -k_1s*x_Z*x_Y
            D(x_Z) ~ -k_1*x_Z*x_Y
            I ~ x_A + 2/21*x_B + 2/21*x_D
        end
    end

    @mtkcompile system = KineticParameterEstimation()

    tspan = (0.0, 2.0)
    tstep = 0.01
    N = Int(floor((tspan[2] - tspan[1])/tstep)) + 1

    include("kinetic_intensity_data.jl")
    intensity(x_A, x_B, x_D) = x_A + 2/21*x_B + 2/21*x_D

    model = JuMP.Model(Ipopt.Optimizer)
    V = length(unknowns(system))
    JuMP.@variable(model, -75.0 <= z[1:V,1:N] <= 150.0)
    pL = [10.0, 10.0, 0.001]
    pU = [1200.0, 1200.0, 40.0]
    JuMP.@variable(model, pL[i] <= p[i=1:3] <= pU[i])

    @test_throws ErrorException("Available integrators: EE, IE") EOptInterface.register_odesystem(model, system, tspan, tstep, "RK4")

    EOptInterface.register_odesystem(model, system, tspan, tstep, "EE")

    JuMP.@objective(model, Min, sum((intensity(z[5,i], z[4,i], z[3,i]) - data[i-1])^2 for i in 2:N))
    JuMP.optimize!(model)

    @test JuMP.termination_status(model) == JuMP.LOCALLY_SOLVED
    @test JuMP.primal_status(model) == JuMP.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(model), 9622.762852574022, atol=1e-3)

    model = JuMP.Model(Ipopt.Optimizer)
    V = length(unknowns(system))
    JuMP.@variable(model, -75.0 <= z[1:V,1:N] <= 150.0)
    pL = [10.0, 10.0, 0.001]
    pU = [1200.0, 1200.0, 40.0]
    JuMP.@variable(model, pL[i] <= p[i=1:3] <= pU[i])
    EOptInterface.register_odesystem(model, system, tspan, tstep, "IE")
    JuMP.@objective(model, Min, sum((intensity(z[5,i], z[4,i], z[3,i]) - data[i-1])^2 for i in 2:N))
    JuMP.optimize!(model)

    @test JuMP.termination_status(model) == JuMP.LOCALLY_SOLVED
    @test JuMP.primal_status(model) == JuMP.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(model), 16796.032234817612, atol=1e-3)

end
