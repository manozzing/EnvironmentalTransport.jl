export advect1d_vertical_L94, advect1d_vertical_ppm

"""
L94 advection in 1-D
"""
function advect1d_vertical_L94(s1d_in, vel_1d, nz, Δt, p)
    ϕ1 = zeros(Float64, nz + 5)
    ϕ2 = zeros(Float64, nz + 5)
    δϕ1 = zeros(Float64, nz + 5)
    Δϕ1_avg = zeros(Float64, nz + 4)
    Δϕ1_mono = zeros(Float64, nz + 3)
    δϕ2 = zeros(Float64, nz + 4)
    Δϕ2_avg = zeros(Float64, nz + 3)
    Δϕ2_mono = zeros(Float64, nz + 3)
    FLUX = zeros(Float64, nz + 1)
    U = zeros(Float64, nz + 1)
    courant = zeros(Float64, nz + 1)
    ΔP = zeros(Float64, nz + 2)
    s1d_out = zeros(Float64, nz + 5)

    #U[1] = vel_1d[1]; U[nz+1] = vel_1d[nz];  ## In case we want Neumann BC.
    U[1] = 0
    U[nz + 1] = 0  ## Dirichlet BC with zero on top and bottom of the column.

    for i in 2:nz
        U[i] = (vel_1d[i - 1] + vel_1d[i]) / 2.0
    end

    for i in 2:(nz + 1)
        ΔP[i] = (p[i - 1] - p[i]) * 100
    end
    ΔP[1] = ΔP[2]
    ΔP[nz + 2] = ΔP[nz + 1]

    ϕ1[3:(nz + 2)] = s1d_in[:] #.* p[1:nz]
    ϕ1[1] = ϕ1[3]
    ϕ1[2] = ϕ1[3]
    ϕ1[nz + 3] = ϕ1[nz + 2]
    ϕ1[nz + 4] = ϕ1[nz + 2]
    ϕ1[nz + 5] = ϕ1[nz + 2]

    for i in 2:(nz + 5)
        δϕ1[i] = ϕ1[i] - ϕ1[i - 1]
    end

    for i in 2:(nz + 4)
        Δϕ1_avg[i] = (δϕ1[i] + δϕ1[i + 1]) / 2.0
    end

    for i in 2:(nz + 3) ## Monotonicity slope limiter
        ϕ1_min = minimum((ϕ1[i - 1], ϕ1[i], ϕ1[i + 1]))
        ϕ1_max = maximum((ϕ1[i - 1], ϕ1[i], ϕ1[i + 1]))
        Δϕ1_mono[i] = sign(Δϕ1_avg[i]) *
                      minimum((abs(Δϕ1_avg[i]), 2 * (ϕ1[i] - ϕ1_min), 2 * (ϕ1_max - ϕ1[i])))
    end

    for i in 1:(nz + 1)
        if U[i] >= 0
            courant = U[i] * Δt / (ΔP[i + 1])
            FLUX[i] = U[i] * (ϕ1[i + 1] + Δϕ1_mono[i + 1] * (1 - courant) / 2.0)
        elseif U[i] < 0
            courant = U[i] * Δt / (ΔP[i])
            FLUX[i] = U[i] * (ϕ1[i + 2] - Δϕ1_mono[i + 2] * (1 + courant) / 2.0)
        end
    end

    for i in 3:(nz + 2)
        ϕ2[i] = ϕ1[i] - Δt / ΔP[i - 1] * (FLUX[i - 1] - FLUX[i - 2])
    end

    for i in 3:(nz + 2)
        s1d_out[i] = ϕ2[i] #./ p[i-2]
        if s1d_out[i] < 0
            s1d_out[i] = 0
        end
    end

    #s1d_out .= s1d_out .- s1d_in

    return s1d_out[3:(nz + 2)]#, FLUX*Δt/Δz
end

"""
PPM advection in 1-D vertical domain
"""
function advect1d_vertical_ppm(s1d_in::Vector{Float64}, vel_1d::Vector{Float64}, nz, Δt, p)
    ϵ = 0.01
    η⁽¹⁾ = 20
    η⁽²⁾ = 0.05

    ΔP = zeros(Float64, nz + 7)

    ϕ1 = zeros(Float64, nz + 6)
    ϕ2 = zeros(Float64, nz + 6)
    δϕ = zeros(Float64, nz + 6)
    δₘϕ = zeros(Float64, nz + 6)
    ϕ₊½ = zeros(Float64, nz + 6)
    ϕL = zeros(Float64, nz + 6)
    ϕR = zeros(Float64, nz + 6)
    ϕLᵈ = zeros(Float64, nz + 6)
    ϕRᵈ = zeros(Float64, nz + 6)
    Δϕ = zeros(Float64, nz + 6)
    δ²ϕ = zeros(Float64, nz + 6)
    η̃ = zeros(Float64, nz + 6)
    η = zeros(Float64, nz + 6)
    ϕ₆ = zeros(Float64, nz + 6)

    FLUX = zeros(Float64, nz + 1)
    courant = zeros(Float64, nz + 1)
    U = zeros(Float64, nz + 1)
    s1d_out = zeros(Float64, nz)

    for i in 1:nz
        ΔP[i + 3] = (p[i] - p[i + 1]) * 100
    end
    ΔP[1] = ΔP[4]
    ΔP[2] = ΔP[4]
    ΔP[3] = ΔP[4]
    ΔP[nz + 4] = ΔP[nz + 3]
    ΔP[nz + 5] = ΔP[nz + 3]
    ΔP[nz + 6] = ΔP[nz + 3]
    ΔP[nz + 7] = ΔP[nz + 3]

    #U[1] = vel_1d[1]; U[nz+1] = vel_1d[nz];
    U[1] = 0
    U[nz + 1] = 0
    for i in 2:nz
        U[i] = (vel_1d[i - 1] + vel_1d[i]) / 2.0
    end

    ## Initialization -- I think it is reasonable to assume boundary values = 0 in the vertical layer.
    ϕ1[4:(nz + 3)] = s1d_in[1:nz]
    #ϕ1[1] = 0; ϕ1[2] = 0; ϕ1[3] = 0
    #ϕ1[nz+4] = 0; ϕ1[nz+5] = 0; ϕ1[nz+6] = 0
    ϕ1[1] = ϕ1[4]
    ϕ1[2] = ϕ1[4]
    ϕ1[3] = ϕ1[4]
    ϕ1[nz + 4] = ϕ1[nz + 3]
    ϕ1[nz + 5] = ϕ1[nz + 3]
    ϕ1[nz + 6] = ϕ1[nz + 3]

    ## Edge value calculation
    for i in 2:(nz + 5)
        δϕ[i] = (ΔP[i] / (ΔP[i - 1] + ΔP[i] + ΔP[i + 1])) *
                (((2 * ΔP[i - 1] + ΔP[i]) / (ΔP[i + 1] + ΔP[i])) * (ϕ1[i + 1] - ϕ1[i]) +
                 ((ΔP[i] + 2 * ΔP[i + 1]) / (ΔP[i - 1] + ΔP[i])) * (ϕ1[i] - ϕ1[i - 1]))
    end

    for i in 2:(nz + 5)
        if (ϕ1[i + 1] - ϕ1[i]) * (ϕ1[i] - ϕ1[i - 1]) > 0
            δₘϕ[i] = min(
                abs(δϕ[i]), 2 * abs(ϕ1[i] - ϕ1[i - 1]), 2 * abs(ϕ1[i + 1] - ϕ1[i])) *
                     sign(δϕ[i])
            #δₘϕ[i] = min(abs(δϕ[i-1]), 2*abs(ϕ1[i]-ϕ1[i-1]), 2*abs(ϕ1[i+1]-ϕ1[i])) * sign(δϕ[i-1])
        else
            δₘϕ[i] = 0
        end
    end

    for i in 2:(nz + 5)
        #ϕ₊½[i] = 1/2*(ϕ1[i+1]+ϕ1[i]) - 1/6*(δₘϕ[i+1]-δₘϕ[i])
        ϕ₊½[i] = ϕ1[i] +
                 (ΔP[i]) / (ΔP[i] + ΔP[i + 1]) * (ϕ1[i + 1] - ϕ1[i]) +
                 1 / (ΔP[i - 1] + ΔP[i] + ΔP[i + 1] + ΔP[i + 2]) * (
                     (2 * ΔP[i + 1] * ΔP[i]) / (ΔP[i] + ΔP[i + 1]) *
                     ((ΔP[i - 1] + ΔP[i]) / (2 * ΔP[i] + ΔP[i + 1]) -
                      (ΔP[i + 2] + ΔP[i + 1]) / (2 * ΔP[i + 1] + ΔP[i])) *
                     (ϕ1[i + 1] - ϕ1[i]) -
                     ΔP[i] * (ΔP[i - 1] + ΔP[i]) / (2 * ΔP[i] + ΔP[i + 1]) * δₘϕ[i + 1] +
                     ΔP[i + 1] * (ΔP[i + 1] + ΔP[i + 2]) / (ΔP[i] + 2 * ΔP[i + 1]) * δₘϕ[i])
    end

    for i in 2:(nz + 4)
        ϕL[i] = ϕ₊½[i]
        ϕR[i] = ϕ₊½[i + 1]
    end

    ## Discontinuity detection
    for i in 2:(nz + 5)
        #δ²ϕ[i] = 1/(6*ΔP[i-1]^2)*(ϕ1[i+1]-2*ϕ1[i]+ϕ1[i-1])
        δ²ϕ[i] = 1 / (ΔP[i - 1] + ΔP[i] + ΔP[i + 1]) *
                 ((ϕ1[i + 1] - ϕ1[i]) / (ΔP[i + 1] + ΔP[i]) -
                  (ϕ1[i] - ϕ1[i - 1]) / (ΔP[i] + ΔP[i - 1]))
    end

    for i in 2:(nz + 5)
        if -δ²ϕ[i + 1] * δ²ϕ[i - 1] * abs(ϕ1[i + 1] - ϕ1[i - 1]) -
           ϵ * min(abs(ϕ1[i + 1]), abs(ϕ1[i - 1])) > 0
            #η_tilde[i] = -(δ²ϕ[i+1]-δ²ϕ[i-1]) * (ΔP[i-1]^2) / (ϕ1[i+1]-ϕ1[i-1])
            η̃[i] = -(δ²ϕ[i + 1] - δ²ϕ[i - 1]) / (ΔP[i + 1] + ΔP[i]) *
                    (ΔP[i]^3 + ΔP[i - 1]^3) / (ϕ1[i + 1] - ϕ1[i - 1])
        else
            η̃[i] = 0
        end
    end

    for i in 2:(nz + 5)
        η[i] = max(0, min(η⁽¹⁾ * (η̃[i] - η⁽²⁾), 1))
    end

    for i in 2:(nz + 4)
        ϕLᵈ[i] = ϕ1[i] + 1 / 2 * δₘϕ[i]
        ϕRᵈ[i] = ϕ1[i + 1] + 1 / 2 * δₘϕ[i + 1]
    end

    for i in 2:(nz + 4)
        ϕL[i] = ϕL[i] * (1 - η[i]) + ϕLᵈ[i] * η[i]
        ϕR[i] = ϕR[i] * (1 - η[i]) + ϕRᵈ[i] * η[i]
    end

    ## Monotonicity examination
    for i in 2:(nz + 4)
        if (ϕR[i] - ϕ1[i]) * (ϕ1[i] - ϕL[i]) <= 0
            ϕL[i] = ϕ1[i]
            ϕR[i] = ϕ1[i]
        elseif (ϕR[i] - ϕL[i]) * (ϕ1[i] - 1 / 2 * (ϕL[i] + ϕR[i])) >= (ϕR[i] - ϕL[i])^2 / 6
            ϕL[i] = 3 * ϕ1[i] - 2 * ϕR[i]
        elseif -(ϕR[i] - ϕL[i])^2 / 6 > (ϕR[i] - ϕL[i]) * (ϕ1[i] - 1 / 2 * (ϕR[i] + ϕL[i]))
            ϕR[i] = 3 * ϕ1[i] - 2 * ϕL[i]
        end
    end

    ## Compute flux
    for i in 2:(nz + 4)
        Δϕ[i] = ϕR[i] - ϕL[i]
    end

    for i in 2:(nz + 4)
        ϕ₆[i] = 6 * (ϕ1[i] - 1 / 2 * (ϕL[i] + ϕR[i]))
    end

    for i in 3:(nz + 3)
        if U[i - 2] >= 0
            courant = U[i - 2] * Δt / ΔP[i + 1]
            FLUX[i - 2] = courant * (ϕR[i] -
                           1 / 2 * courant * (Δϕ[i] - (1 - 2 / 3 * courant) * ϕ₆[i]))
        else
            courant = U[i - 2] * Δt / ΔP[i]
            FLUX[i - 2] = courant * (ϕL[i + 1] -
                           1 / 2 * courant *
                           (Δϕ[i + 1] - (1 + 2 / 3 * courant) * ϕ₆[i + 1]))
        end
    end

    for i in 4:(nz + 3)
        ϕ2[i] = ϕ1[i] + (FLUX[i - 3] - FLUX[i - 2])
    end

    for i in 1:nz
        s1d_out[i] = ϕ2[i + 3] * (ϕ2[i + 3] >= 0)
    end

    return s1d_out
end
