using EAGO
using EOptInterface
using JuMP
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# Create algebraic ModelingToolkit system
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
        # Unknown parameters (free design variables)
        F

        # Known parameters
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
        # Unknown parameters (free design variables)
        V

        # Known parameters
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

# Compile system
@mtkcompile system = ReactorSeparatorRecycle()

# Define symbolic expressions of constraints and objective
# Use syntax System.Component.Connector.Variable, System.Component.Component.Parameter, or System.Component.Parameter
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

# Solve reduced-space model
# Create JuMP model
model = Model(EAGO.Optimizer)

# Retrieve decision variables from ModelingToolkit system
# Returns [sep1.in.F(t), sep1.in.y_B(t), sep1.in.y_C(t), sep1.outL.y_C(t), influent₊F, cstr₊V]
decision_vars(system)

# Create decision variables
# z = [sep1.in.F(t), sep1.in.y_B(t), sep1.in.y_C(t), sep1.outL.y_C(t)]
# p = [influent₊F, cstr₊V]
xL = zeros(6)
xU = [100.0, 1.0, 1.0, 1.0, 100.0, 10.0]
@variable(model, xL[i] <= x[i=1:6] <= xU[i])

# Register ModelingToolkit nonlinear system as constraints and register objective
register_nlsystem(model, system, obj, [g1, g2])

# Optimize model
optimize!(model)

# Display results
println("Termination Status: $(termination_status(model))")
println("Primal Status: $(primal_status(model))")
println("Solve Time: $(round.(solve_time(model), digits=5))")
println("f^* = $(round(objective_value(model), digits=5))")
println("x* = $(round.(value.(x), digits=3))")

# Retrieve full-space solution
full_solution(model, system)

# Solve full-space model
# Compile system
@named full_system = ReactorSeparatorRecycle()

# Create JuMP model
full_model = Model(EAGO.Optimizer)

# Create decision variables
xL = zeros(50)
xU = vcat(repeat([100.0, 1.0, 1.0, 1.0], 12), 100.0, 10.0)
@variable(full_model, xL[i] <= x[i=1:50] <= xU[i])

# Register ModelingToolkit nonlinear system as constraints and objective
register_nlsystem(full_model, full_system, obj, [g1, g2])

# Optimize model and retrieve results
optimize!(full_model)
println("Termination Status: $(termination_status(full_model))")
println("Primal Status: $(primal_status(full_model))")
println("Solve Time: $(round.(solve_time(full_model), digits=5))")
println("f^* = $(round(objective_value(full_model), digits=5))")
println("x* = $(round.(value.(x), digits=3))")
