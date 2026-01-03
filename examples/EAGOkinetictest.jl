using DataFrames, CSV

# CODE FROM THE EAGO NOTEBOOK
data = CSV.read("kinetic_intensity_data.csv", DataFrame)
function explicit_euler_integration(p::T...) where {T}
    x = zeros(T, 1005)
    x[4] = 0.4
    x[5] = 140

    # sets known parameter values
    Temp = 273.0
    K2 = 46.0*exp(6500.0/Temp - 18.0)
    K3 = 2.0*K2
    k1 = 53.0
    k1s = k1*10^(-6)
    k5 = 0.0012
    cO2 = 0.002

    h = 0.01
    # offset by 1, since the initial condition is x[1:5]
    for i = 1:200
        term1 = k1*x[5i-1]*x[5i] - cO2*(p[1] + p[2])*x[5i-4]
        term2 = p[1]*x[5i-2]/K2 + p[2]*x[5i-3]/K3 - k5*x[5i-4]^2
        x[5i+1] = x[5i-4] + h*(term1 + term2)
        x[5i+2] = x[5i-3] + h*(p[2]*cO2*x[5i-4] - (p[2]/K3 + p[3])*x[5i-3])
        x[5i+3] = x[5i-2] + h*(p[1]*cO2*x[5i-4] - p[1]*x[5i-2]/K2)
        x[5i+4] = x[5i-1] + h*(-k1s*x[5i-1]*x[5i])
        x[5i+5] = x[5i] + h*(-k1*x[5i-1]*x[5i])
    end
    return x
end
intensity(xA, xB, xD) = xA + (2/21)*xB + (2/21)*xD
function objective(p::T...) where {T}
    x = explicit_euler_integration(p...)
    SSE = zero(T)
    for i = 1:200
        SSE += (intensity(x[5i+1], x[5i+2], x[5i+3]) - data[!,:intensity][i])^2
    end
    return SSE
end

# USING EAGO NOTEBOOK p* (EE)
pEAGO = [
    828.0651634543585
    385.732149484276
    14.567016208483889
    ]
objective(pEAGO...) # 9622.76285258131 (MATCHES EOI REUSLTS)
# EAGO NOTEBOOK REPORTED UPPER BOUND: 9627.574244007586 (DOES NOT MATCH ITS OWN RESULT)

# USING EOI IPOPT p* (EE) (USING 0.05 OR 1E-4 TOLERANCE GIVES THE EXACT SAME RESULTS)
pEOI = [
    828.0653808971397
    385.7322187557441
    14.56636196685152
    ]
objective(pEOI...) # 9622.762852633188
# EOI IPOPT RESULT: 9622.76259884641 (MATCHES)