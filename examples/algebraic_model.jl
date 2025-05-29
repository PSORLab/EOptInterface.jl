using ModelingToolkit, JuMP, EOptInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

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
        y_A = 1
        y_B = 0
        y_C = 0
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
        out.y_A + out.y_B + out.y_C ~ 1
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
        
        outV.y_A + outV.y_B + outV.y_C ~ 1
        outV.y_C ~ 0
        outV.y_B ~ 0

        outL.y_A + outL.y_B + outL.y_C ~ 1
        outL.y_A ~ 0
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

        outV.y_A + outV.y_B + outV.y_C ~ 1
        outV.y_A ~ 0
        outV.y_C ~ 0

        outL.y_A + outL.y_B + outL.y_C ~ 1
        outL.y_A ~ 0
        outL.y_B ~ 0
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

@mtkcompile s = ReactorSeparatorRecycle()

exprF5 = s.sep2.outV.F
exprTau = s.cstr.V/(s.cstr.out.F*(s.cstr.out.y_A*s.cstr.in.V_A + s.cstr.out.y_B*s.cstr.in.V_B + s.cstr.out.y_C*s.cstr.in.V_C))
f_CSTR = (25764 + 8178*s.cstr.V)/2.5
s1cap = 132718 + s.cstr.out.F*(369*s.cstr.out.y_A - 1113.9*s.cstr.out.y_B)
s2cap = 25000 + s.sep1.outL.F*(6984.5*s.sep1.outL.y_B - 3869.53*s.sep1.outL.y_C^2)
s1op = s.cstr.out.F*(3+36.11*s.cstr.out.y_A + 7.71*s.cstr.out.y_B)*26.32e-3
s2op = s.sep1.outL.F*(26.21 + 29.45*s.sep1.outL.y_B)*26.32e-3;
f_Sep = (s1cap+s2cap)/2.5 + 0.52*(s1op+s2op)
g1 = 25 - exprF5
g2 = 475/3600 - exprTau
obj = f_CSTR + f_Sep

using EAGO
model = Model(EAGO.Optimizer)
decision_vars(s)
xL = zeros(6)
xU = [100, 1, 1, 1, 100, 10]
@variable(model, xL[i] <= x[i=1:6] <= xU[i])
register_nlsystem(model, s, obj, [g1, g2])
JuMP.optimize!(model)
JuMP.value.(x)
full_solutions(model, s)