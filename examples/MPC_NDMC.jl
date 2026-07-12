using JuMP, Ipopt
using DelimitedFiles
using Statistics
using DataFrames, CSV

"""
    MPC(t0, h, i_out, s, u0, Pspan, M, y0, Cin, Cs, k, ht, tspan, w)

Run the legacy NDMC DMC-style controller over the full simulation horizon.

This is a stand-alone historical example that predates the reusable
`build_tracking_mpc(...)` API. It keeps the original logic for reference:
solve a finite-horizon DMC problem, apply the first move, simulate the plant
forward, and repeat.
"""
function MPC(t0,h,i_out,s,u0,Pspan,M,y0,Cin,Cs,k,ht,tspan,w)
    Nt = convert(Int,round(tspan/ht))
    uopt = zeros(Float64,Nt)
    t0_Data = 0.0
    yk = y0[1:3]
    uold = u0
    for i = 1:(Nt)
        # Receding-horizon loop:
        # 1. solve for a future move sequence,
        # 2. apply only the first move,
        # 3. simulate the plant up to the next trigger,
        # 4. repeat with the updated measurement.
        t0 = ht*(i-1)
        usol = Optim(t0,h,i_out,yk,s,u0,Pspan,M,w)
        uopt[i] = usol[1]
        t_re = t0
        uold = [uold;usol[1]]
        Data_span = t_re - t0_Data
        n_out_Data = convert(Int,round(Data_span/h/i_out))
        yk = Get_Data(t0_Data,y0,h,n_out_Data,i_out,Cin,Cs,k,ht,t_re,uold)
        u0 = [u0;usol[1]]
        u0 = u0[2:end]
    end
    return uopt
end

"""
    Optim(t0, h, i_out, yk, s, u0, Pspan, M, w)

Solve one finite-horizon DMC optimization problem.

The decision variables are the next `M` control moves. The objective is
provided by `obj_DMC(...)`.
"""
function Optim(t0,h,i_out,yk,s,u0,Pspan,M,w)
    u_lo = zeros(Float64,M)
    u_hi = 800.0*ones(Float64,M)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol",1e-3)
    # Keep the historical solver setting deterministic so repeated benchmark
    # runs produce the same reference trajectory.
    set_optimizer_attribute(model, "bound_frac",0.25)
    set_optimizer_attribute(model, "print_level",0)
    function obj(u...)
        return obj_DMC(t0,h,i_out,yk,s,u0,Pspan,M,w,u...)
    end
    JuMP.register(model,:obj,M,obj,autodiff=true)
    #JuMP.register(model,:DMC,M,DMC,autodiff=true)
    @variable(model, u_lo[i] <= u[i = 1:M] <= u_hi[i])
    @NLobjective(model, Min, obj(u...))
    time = @elapsed JuMP.optimize!(model)
    fval = JuMP.objective_value(model)
    usol = JuMP.value.(u)
    return usol
end

"""
    obj_DMC(t0, h, i_out, yk, s, u0, Pspan, M, w, u...)

Return the scalar DMC objective value for one candidate control sequence.

The cost is:
- squared tracking error for the three conductivity outputs,
- plus a move penalty on changes between future control moves.
"""
function obj_DMC(t0,h,i_out,yk,s,u0,Pspan,M,w,u...)
    yp = DMC(t0,h,i_out,yk,s,u0,u,Pspan,ht)
    P = length(yp[:,1])

    SP = 280.0
    e = zeros(typeof(u[1]),P,3)
    for j = 1:3
        for i = 1:P
            e[i,j] = SP - yp[i,j]
        end
    end

    du = zeros(typeof(u[1]),M-1)
    for i = 1:(M-1)
        du[i] = u[i+1] - u[i]
    end

    SSE = 0
    for j = 1:3
        for i = 1:P
            SSE = SSE + e[i,j]^2
        end
    end

    for i = 1:(M-1)
        SSE = SSE + w*du[i]^2
    end

    f = SSE
    return f
end

"""
    Get_Data(t0, y0, h, n_out, i_out, Cin, Cs, k, ht, t_re, uold)

Simulate the plant up to the current review time and return the measured
conductivity states used by the controller.
"""
function Get_Data(t0,y0,h,n_out::Int,i_out::Int,Cin,Cs,k,ht,t_re,uold)
    tout,yout = EE_DataModel(t0,y0,h,n_out,i_out,Cin,Cs,k,ht,t_re,uold)
    yk = yout[end,1:3]
    return yk
end

"""
    EE_DataModel(t0, y0, h, n_out, i_out, Cin, Cs, k, ht, t_re, uold)

Explicit-Euler simulation of the plant over one data-acquisition interval.

This is the "plant side" integrator used to update the current measurement
between MPC triggers.
"""
function EE_DataModel(t0,y0,h,n_out::Int,i_out::Int,Cin,Cs,k,ht,t_re,uold)
    tout = zeros(Float64,n_out+1)
    yout = zeros(Float64,n_out+1,length(y0))
    tout[1] = copy(t0)
    yout[1,:] = copy(y0)
    t = t0
    y = y0
    for j = 2:n_out+1
        for l = 1:i_out
            y = y + h*DataModel(t,y,Cin,Cs,k,ht,t_re,uold)
            t = t + h
        end
        tout[j] = t
        yout[j,:] = y
    end
    return tout, yout
end

"""
    DataModel(t, x, Cin, Cs, k, ht, t_re, uold)

Right-hand side of the legacy NDMC plant model used during closed-loop
simulation.

The state vector is:
- zone-1 conductivity,
- zone-2 conductivity,
- zone-3 conductivity,
- mixed-zone conductivity,
- dissolved oxygen.
"""
function DataModel(t,x,Cin,Cs,k,ht,t_re,uold)
    if t >= 2100.0 && t <= 2250.0
        Cin1 = Cin[1]
        Cin2 = Cin[2]
        Cin3 = Cin[3]
    else
        Cin1 = Cs
        Cin2 = Cs
        Cin3 = Cs
    end

    # Reconstruct the piecewise-constant applied control from the historical
    # control sequence `uold`.
    Q = 0.0
    idx_t_re = convert(Int,round(t_re/ht))
    for i = 1:idx_t_re
        if t >= ht*(i-1) && t < ht*i
            Q = uold[end - idx_t_re + i]
        end
    end

    # Oxygen flow
    cO = x[5]
    SOTE = 0.1

    # Reaction term
    M_N = 14.0067 # g/mol
    r_AOmax = 0.67 # mg N-NH+ /g VSS_AO min
    K_OAO = 0.3 # mg O2/L
    r_AO = r_AOmax*(cO/(K_OAO + cO)) # mg N-NH+ /g VSS_AO min
    X_AO = 0.505 # g VSS/L
    R = -r_AO*X_AO # mg N-NH+/(L*min)
    R = R*0.001/60/M_N # mol/(L*s)
    R_C = -c_to_C(abs(R)) # uS/cm/s

    # Volumetric flow
    Vl = 1000.0
    V = Vl/4 # L
    min = Vl/240/4 # L/s
    mout = Vl/240/4 # L/s

    h = zeros(Float64,5)
    h[1] = 1/V*(k[1]*(x[4] - x[1]) + min*Cin1 - mout*x[1]) + R_C  # zone 1
    h[2] = 1/V*(k[2]*(x[4] - x[2]) + min*Cin2 - mout*x[2]) + R_C  # zone 2
    h[3] = 1/V*(k[3]*(x[4] - x[3]) + min*Cin3 - mout*x[3]) + R_C  # zone 3
    h[4] = 1/V*k[4]*(x[1] + x[2] + x[3] - 3*x[4]) + R_C           # mixing zone
    h[5] = Otrans_model(t,cO,Q,SOTE)
    return h
end

"""
    DMC(t0, h, i_out, yk, s, u0, u, Pspan, ht)

Predict future outputs from the step response `s` and the proposed move
sequence `u`.

This is the classical DMC prediction step: extend the manipulated-variable
trajectory, form move increments, then convolve them with the step response.
"""
function DMC(t0,h,i_out,yk,s,u0,u,Pspan,ht)
    P = convert(Int,Pspan/ht)
    M = length(u)
    N = length(s)
    ut = zeros(typeof(u[1]),N+P)
    for i = 1:N-1
        ut[i] = u0[i]
    end

    for i = 1:M
        ut[i+N-1] = u[i]
    end

    for i = (N+M):(N+P)
        ut[i] = u[M]
    end
    #ut = [u0,u,u[Nc]*ones(typeof(u[1]),Np-Nc+1)]
    dut = zeros(typeof(u[1]),length(ut)-1)
    for i = 1:length(ut)-1
        dut[i] = ut[i+1] - ut[i]
    end

    # First reconstruct the current model output implied by the previous move
    # history, then use the measurement mismatch as a bias term.
    yp0 = 0.0
    for i = 1:N-1
        yp0 = yp0 + s[i]*dut[N-i]
    end
    yp0 = yp0 + s[N]*ut[1]

    dk1 = yk[1] - yp0
    dk2 = yk[2] - yp0
    dk3 = yk[3] - yp0
    dk = [dk1,dk2,dk3]

    yp = zeros(typeof(u[1]),P,3)
    for l = 1:3
        for j = 1:P
            for i = 1:N-1
                yp[j,l] = yp[j,l] + s[i]*dut[N-i+j]
            end
            yp[j,l] = yp[j,l] + s[N]*ut[j+1] + dk[l]
        end
    end

    return yp
end

"""
    Step(t0, y0, h, n_out, i_out, Cs, k, ht, u, tspan, ts, Nspan)

Generate the step response used by the legacy DMC controller.

The function simulates one input step, samples the output at the controller
sample time, and returns the resulting response coefficients.
"""
function Step(t0,y0,h,n_out,i_out,Cs,k,ht,u,tspan,ts,Nspan)
    idx_ts = convert(Int,round(ts/ht))
    tout, yout = EE_StepModel(t0,y0,h,n_out,i_out,Cs,k,ht,u,tspan)
    y1 = yout[:,1]
    y2 = yout[:,2]
    y3 = yout[:,3]
    dt = convert(Int,round(ht/h/i_out))
    y1con = y1[1:dt:length(y1)]
    y2con = y2[1:dt:length(y2)]
    y3con = y3[1:dt:length(y3)]
    N = convert(Int,round(Nspan/ht))
    s = zeros(Float64,N)
    for i = 1:N
        s[i] = y1con[idx_ts+i] - y1con[idx_ts]
        #s[i,1] = y1con[idx_ts+i] - y1con[idx_ts]
        #s[i,2] = y2con[idx_ts+i] - y2con[idx_ts]
        #s[i,3] = y3con[idx_ts+i] - y3con[idx_ts]
    end
    return s
end

"""
    EE_StepModel(t0, y0, h, n_out, i_out, Cs, k, ht, u, tspan)

Explicit-Euler simulation used to build the DMC step response.
"""
function EE_StepModel(t0,y0,h,n_out::Int,i_out::Int,Cs,k,ht,u,tspan)

    tout = zeros(typeof(u[1]),n_out+1)
    yout = zeros(typeof(u[1]),n_out+1,length(y0))
    tout[1] = copy(t0)
    yout[1,:] = copy(y0)
    t = t0
    y = y0
    for j = 2:n_out+1
        for l = 1:i_out
            y = y + h*StepModel(t,y,Cs,k,ht,u,tspan)
            t = t + h
        end
        tout[j] = t
        yout[j,:] = y
    end
    return tout, yout
end

"""
    StepModel(t, x, Cs, k, ht, u, tspan)

Right-hand side used while generating the DMC step response.

This uses the same plant physics as `DataModel(...)`, but with a fixed clean
influent and a known input step profile.
"""
function StepModel(t,x,Cs,k,ht,u,tspan)
    Cin1 = Cs
    Cin2 = Cs
    Cin3 = Cs

    # Control
    Q = 0.0
    for i = 1:convert(Int,round(tspan/ht))
        if t>= ht*(i-1) && t <= ht*i
            Q = u[i]
        end
    end

    # Oxygen flow
    cO = x[5]
    SOTE = 0.1

    # Reaction term
    M_N = 14.0067; # g/mol
    r_AOmax = 0.67; # mg N-NH+ /g VSS_AO min
    K_OAO = 0.3; # mg O2/L
    r_AO = r_AOmax*(cO/(K_OAO + cO)); # mg N-NH+ /g VSS_AO min
    X_AO = 0.505; # g VSS/L
    R = -r_AO*X_AO; # mg N-NH+/(L*min)
    R = R*0.001/60/M_N; # mol/(L*s)
    R_C = -c_to_C(abs(R)); # uS/cm/s

    # Volumetric flow
    Vl = 1000.0
    V = Vl/4; # L
    min = Vl/240/4; # L/s
    mout = Vl/240/4; # L/s

    h = zeros(typeof(u[1]),5)
    h[1] = 1/V*(k[1]*(x[4] - x[1]) + min*Cin1 - mout*x[1]) + R_C  # zone 1
    h[2] = 1/V*(k[2]*(x[4] - x[2]) + min*Cin2 - mout*x[2]) + R_C  # zone 2
    h[3] = 1/V*(k[3]*(x[4] - x[3]) + min*Cin3 - mout*x[3]) + R_C  # zone 3
    h[4] = 1/V*k[4]*(x[1] + x[2] + x[3] - 3*x[4]) + R_C           # mixing zone
    h[5] = Otrans_model(t,cO,Q,SOTE)
    return h
end

"""
    Otrans_model(t, cO, Q, SOTE)

Return the dissolved-oxygen derivative caused by aeration and oxygen
consumption.

This is shared by both the step-response generation and the closed-loop plant
simulation.
"""
function Otrans_model(t,cO,Q,SOTE)
    W = 0.2967*Q
    SOTR = SOTE*W; # mg/s
    cO_s = 9.1; # mg/L
    Vl = 1000.0; # L
    kla = SOTR/cO_s/Vl; # s^(-1)
    r_AOmax = 0.67; # mg N-NH+ /g VSS_AO min
    K_OAO = 0.3; # mg O2/L
    r_AO = r_AOmax*(cO/(K_OAO + cO)); # mg N-NH+ /g VSS_AO min
    Phi_OAO = 2.5; # mg O/ mg N-NH+
    X_AO = 0.505 # g VSS/L
    dcO = kla*(cO_s - cO) - r_AO*Phi_OAO*X_AO/60; # mg O/(L*s)
    return dcO
end

"""
    c_to_C(c)

Convert concentration-like reaction rate information into conductivity units.

The piecewise formula matches the original NDMC legacy implementation and is
kept here for reproducibility of the old example.
"""
function c_to_C(c)
    A = 60.2;
    B = 0.229;
    La_0 = 149.6; # cm^2*S/mol
    l = 10.7828767123288;
    k = 132323.287671233;
    if c <= 1e-3
        La = La_0 - (A + B*La_0)*c^0.5; # cm^2*S/mol
        C = La*c*1e3; # uS/mol
    else
        C = l + k*c;
    end
    return C
end

"""
    EE_SimModel(t0, y0, h, n_out, i_out, Cin, Cs, k, ht, tspan, uopt)

Explicit-Euler simulation of the final closed-loop trajectory after the MPC
move sequence `uopt` has been computed.
"""
function EE_SimModel(t0,y0,h,n_out::Int,i_out::Int,Cin,Cs,k,ht,tspan,uopt)
    tout = zeros(Float64,n_out+1)
    yout = zeros(Float64,n_out+1,length(y0))
    tout[1] = copy(t0)
    yout[1,:] = copy(y0)
    t = t0
    y = y0
    for j = 2:n_out+1
        for l = 1:i_out
            y = y + h*DataModel(t,y,Cin,Cs,k,ht,tspan,uopt)
            t = t + h
        end
        tout[j] = t
        yout[j,:] = y
    end
    return tout, yout
end

"""
    SimModel(t, x, Cin, Cs, k, ht, tspan, uopt)

Right-hand side for the post-optimization closed-loop simulation.

This is similar to `DataModel(...)`, but reads the final optimized move
sequence `uopt` instead of the historical online sequence `uold`.
"""
function SimModel(t,x,Cin,Cs,k,ht,tspan,uopt)
    if t >= 2100.0 && t <= 2250.0
        Cin1 = Cin[1]
        Cin2 = Cin[2]
        Cin3 = Cin[3]
    else
        Cin1 = Cs
        Cin2 = Cs
        Cin3 = Cs
    end

    # Control
    Q = 0.0
    idx_t = convert(Int,round(tspan/ht))
    for i = 1:idx_t
        if t >= ht*(i-1) && t < ht*i
            Q = uopt[i]
        end
    end

    # Oxygen flow
    cO = x[5]
    SOTE = 0.1

    # Reaction term
    M_N = 14.0067 # g/mol
    r_AOmax = 0.67 # mg N-NH+ /g VSS_AO min
    K_OAO = 0.3 # mg O2/L
    r_AO = r_AOmax*(cO/(K_OAO + cO)) # mg N-NH+ /g VSS_AO min
    X_AO = 0.505 # g VSS/L
    R = -r_AO*X_AO # mg N-NH+/(L*min)
    R = R*0.001/60/M_N # mol/(L*s)
    R_C = -c_to_C(abs(R)) # uS/cm/s

    # Volumetric flow
    Vl = 1000.0
    V = Vl/4 # L
    min = Vl/240/4 # L/s
    mout = Vl/240/4 # L/s

    h = zeros(Float64,5)
    h[1] = 1/V*(k[1]*(x[4] - x[1]) + min*Cin1 - mout*x[1]) + R_C  # zone 1
    h[2] = 1/V*(k[2]*(x[4] - x[2]) + min*Cin2 - mout*x[2]) + R_C  # zone 2
    h[3] = 1/V*(k[3]*(x[4] - x[3]) + min*Cin3 - mout*x[3]) + R_C  # zone 3
    h[4] = 1/V*k[4]*(x[1] + x[2] + x[3] - 3*x[4]) + R_C           # mixing zone
    h[5] = Otrans_model(t,cO,Q,SOTE)
    return h
end

## Optimization
# The remainder of this file is the runnable legacy example:
# 1. define calibrated parameters,
# 2. generate the step response,
# 3. solve the DMC problem,
# 4. simulate the resulting closed loop,
# 5. export the applied moves.

# k
k1_hi = 2.37570423107917E-03
k2_hi = 1.41065048786101E-03
k3_hi = 1.50410563455869E-03
k4_hi = 9.47661458622374E-01

k1_me = 1.38465080722995E-03
k2_me = 2.91428968468164E-03
k3_me = 2.57662674484102E-03
k4_me = 1.89896205118280E+00

k1 = (k1_hi + k1_me)/2
k2 = (k2_hi + k2_me)/2
k3 = (k3_hi + k3_me)/2
k4 = (k4_hi + k4_me)/2
k = (1000/0.38).*[k1;k2;k3;k4]

# Step
# Build the step response coefficients used by the DMC predictor.
t0 = 0.0
y0 = [280.0; 280.0; 280.0; 280.0; 0.0]
step_span = 4000.0
h = 1.0e-2
i_out = 1000
n_out = convert(Int,round(step_span/h/i_out))
Cs = 285.0
ht = 20.0
Nt = convert(Int,round(step_span/ht))
u_set = 0.0*ones(Float64,Nt)
ts = 2000.0
idx_ts = convert(Int,round(ts/ht))
for i = idx_ts+1:Nt
    u_set[i] = 1.0
end
Nspan = 800.0
s = Step(t0,y0,h,n_out,i_out,Cs,k,ht,u_set,step_span,ts,Nspan)
#s = readdlm("s.csv", ',')

# MPC using DMC
# Solve for the full sequence of applied aeration moves.
t0 = 0.0
h = 1.0e-2
i_out = 1000
tspan = 4000.0
Pspan = 400.0 # Prediction Horizon
Mspan = 60.0 # Control Horizon
ht = 20.0
M = convert(Int,round(Mspan/ht))
u0 = 168.0*ones(Float64,length(s)-1)
y0 = [280.0; 280.0; 280.0; 280.0; 0.0]
Cs = 285.0
#Cin = [320.0; Cs; Cs]
#Cin = [Cs; 320.0; Cs]
Cin = [Cs; Cs; 320.0]
#Cin = [300.0;350.0;310.0]
w = 0.0
uopt = MPC(t0,h,i_out,s,u0,Pspan,M,y0,Cin,Cs,k,ht,tspan,w)

# Simulation Results
# Re-simulate the plant with the optimized move sequence and save the control
# history for later plotting and notebook use.
i_out_sim = 100
n_out_sim =  convert(Int,round(tspan/h/i_out_sim))
time, yout = EE_SimModel(t0,y0,h,n_out_sim,i_out_sim,Cin,Cs,k,ht,tspan,uopt)
C1 = Vector(yout[:,1])
C2 = Vector(yout[:,2])
C3 = Vector(yout[:,3])
writedlm("usol_MPC.csv",uopt,',')


## Test
#=
t0_Data = 0.0
t_re = 3000.0
Simu_span = t_re - t0_Data
y0 = [280.0; 280.0; 280.0; 280.0; 0.0]
h = 1.0e-2
i_out = 1000
n_out = convert(Int,round(Simu_span/h/i_out))
Cs = 285.0
Cin = [320.0; Cs; Cs]
u0 = 168.0*ones(Float64,length(s)-1)
uold = 168.0*ones(Float64,400)
yk = Get_Data(t0_Data,y0,h,n_out,i_out,Cin,Cs,k,ht,t_re,uold)

Pspan = 200.0
t0 = t_re
i_out = 1000
u0 = 168.0*ones(Float64,length(s)-1)
u = 169.0*ones(Float64,3)
yp = DMC(t0,h,i_out,yk,s,u0,u,Pspan,ht)
=#
