using Test, EOptInterface, JuMP, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "EOptInterface.jl" begin
    @mtkmodel AlgebraicTest begin
        @parameters begin
            k₁ = 0.40
            k₂ = 0.055
            V_A = 8.937e-2
            V_B = 1.018e-1
            V_C = 1.130e-1
            V
            F₁
        end
        @variables begin
            F₂(t); F₃(t); F₄(t); F₅(t); F₆(t); F₇(t)
            y_3A(t); y_3B(t); y_3C(t)
            y_4B(t); y_4C(t)
        end
        begin
            r₁ = (k₁*y_3A)/(y_3A*V_A + y_3B*V_B + y_3C*V_C)
            r₂ = (k₂*y_3B)/(y_3A*V_A + y_3B*V_B + y_3C*V_C)
        end
        @equations begin
            F₅ ~ (y_4B * F₄)
            F₁ + F₇ ~ F₂
            (y_3A * F₃) ~ F₂ - r₁*V
            (y_3B * F₃) ~ (r₁ - r₂)*V
            (y_3C * F₃) ~ r₂*V
            F₃ ~ F₄ + F₇
            (y_3B * F₃) ~ (y_4B * F₄)
            (y_3C * F₃) ~ (y_4C * F₄)
            F₄ ~ F₅ + F₆
            y_3A + y_3B + y_3C ~ 1
            y_4B + y_4C ~ 1
        end
    end
    @mtkcompile n = AlgebraicTest()
    g1 = 25 - n.F₅
    g2 = 475/3600 - n.V/(n.F₃*(n.y_3A*n.V_A + n.y_3B*n.V_B + n.y_3C*n.V_C))
    f_CSTR = (25764 + 8178*n.V)/2.5
    s1cap = 132718 + n.F₃*(369*n.y_3A - 1113.9*n.y_3B)
    s2cap = 25000 + n.F₄*(6984.5*n.y_4B - 3869.53*n.y_4C^2)
    s1op = n.F₃*(3+36.11*n.y_3A + 7.71*n.y_3B)*26.32e-3
    s2op = n.F₄*(26.21 + 29.45*n.y_4B)*26.32e-3;
    f_Sep = (s1cap+s2cap)/2.5 + 0.52*(s1op+s2op)
    obj = f_CSTR + f_Sep
    using Ipopt
    model = JuMP.Model(Ipopt.Optimizer)
    EOptInterface.decision_vars(n)
    xL = zeros(6)
    xU = [1, 100, 1, 1, 10, 100]
    @variable(model, xL[i] <= xvar[i in 1:6] <= xU[i])
    EOptInterface.register_nlsystem(model, n, obj, [g1, g2])
    JuMP.optimize!(model)
    @test abs(JuMP.value.(xvar)[6] - 26.316700011101723)/26.316700011101723 < 1e-6
    @test JuMP.termination_status(model) == JuMP.LOCALLY_SOLVED
    @test abs(EOptInterface.full_solutions(model, n)[n.y_3C] - 0.01385685924868818)/0.01385685924868818 < 1e-6

    @mtkmodel KineticParameterEstimation begin
        @parameters begin
            T = 273
            K_2 = 46*exp(6500/T-18)
            K_3 = 2*K_2
            k_1 = 53
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
    @mtkcompile o = KineticParameterEstimation()
    tspan = (0.0,2.0)
    tstep = 0.01
    data = [
        66.0952
        104.762
        110.333
        114.905
        122.238
        125.429
        125.429
        123.476
        121.286
        118.857
        117.667
        116.143
        113.857
        111.571
        108.81
        105.952
        104.048
        102.048
        100.143
        98.5238
        96.2381
        94.381
        91.6667
        89.5714
        87.1429
        84.8571
        83.4286
        81.1905
        78.9048
        77.0476
        75.4762
        73.4762
        71.8095
        70.6667
        68.381
        67.3333
        65.0952
        63.7143
        62.0476
        60.8571
        59.619
        58.2857
        57.4762
        56.4762
        55.8095
        54.5238
        53
        51.8571
        50.4286
        49.381
        47.9524
        47.3714
        46.8952
        46.4857
        45.9048
        45.0762
        44.3238
        43.4143
        43.5429
        42.3619
        41.8381
        40.2381
        39.1286
        38.7857
        37.081
        36.9524
        36.581
        36.281
        35.3476
        34.8905
        34.1667
        33.6714
        32.9667
        31.8429
        31.5429
        31.1476
        30.9905
        29.9571
        29.1333
        28.7857
        28.4429
        28.3476
        27.5429
        27.4333
        27.6048
        27.1762
        27.2
        26.4333
        25.7619
        24.8095
        24.7429
        24.2857
        24.1714
        23.5667
        23.5476
        23.3952
        22.919
        22.3095
        21.8048
        21.2857
        21.2048
        20.8429
        20.4429
        20.0048
        19.9381
        19.5
        19.8667
        18.9333
        19.1381
        18.9619
        18.5476
        17.9048
        17.7571
        18.5333
        18.3762
        18.3571
        18.3286
        18.2762
        18.3952
        17.5952
        18.1524
        18.1952
        17.8476
        17.9095
        17.5048
        17.5
        15.9619
        16.2095
        16.181
        15.6952
        15.7095
        15.4619
        15.9476
        16
        16.1952
        16.1143
        15.7429
        15.5762
        15.7048
        15.8095
        15.6667
        14.9048
        14.5857
        14.7524
        14.7571
        14.9762
        14.5333
        14.5524
        14.0143
        13.6286
        13.4429
        13.4667
        13.319
        12.9333
        13.1238
        12.7476
        12.9333
        13.0714
        13.0714
        12.7619
        12.4238
        12.5143
        12.9143
        12.5714
        13.3667
        13.2286
        13.7905
        13.7571
        13.5905
        12.9667
        12.981
        12.8857
        12.919
        13.0143
        13.0095
        12.3857
        12.5571
        12.3429
        12.7571
        12.681
        12.5429
        12.1857
        12.7905
        12.5571
        12.8429
        12.5476
        12.5714
        12.3762
        11.9952
        11.4571
        11.3
        11.1524
        11.681
        11.619
        11.9048
        12
        12.0762
        11.9143
        11.7619
        11.5333
    ]
    intensity(x_A,x_B,x_D) = x_A + 2/21*x_B + 2/21*x_D
    model = Model(Ipopt.Optimizer)
    N = Int(floor((tspan[2] - tspan[1])/tstep))+1
    V = length(unknowns(o))
    zL = zeros(V)
    zU = [140.0, 0.4, 140.0, 140.0, 140.0]
    @variable(model, zL[i] <= z[i in 1:V,1:N] <= zU[i])
    pL = [10, 10, 0.001]
    pU = [1200, 1200, 40]
    @variable(model, pL[i] <= p[i=1:3] <= pU[i])
    register_odesystem(model, o, tspan, tstep, "EE")
    @objective(model, Min, sum((intensity(z[5,i],z[4,i],z[3,i]) - data[i-1])^2 for i in 2:N))
    JuMP.optimize!(model)
    @test JuMP.termination_status(model) == JuMP.LOCALLY_SOLVED
    @test abs(JuMP.objective_value(model) - 9622.762852574022)/9622.762852574022 < 1e-6
end