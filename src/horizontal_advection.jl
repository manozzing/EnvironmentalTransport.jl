export advect1d_l94, advect1d_ppm, advect2D

"""
L94 advection in 1-D (Lin et al., 1994)
"""
function advect1d_l94(s1d_in::Vector{Float64}, U::Vector{Float64}, nz, Δt, Δz)
    s1d_out = zeros(Float64, nz)

    ## Initialization
    ϕ1(i) = i < 3 ? s1d_in[1] : i > nz + 2 ? s1d_in[end] : s1d_in[i - 2]

    δϕ1(i) = ϕ1(i) - ϕ1(i - 1)

    Δϕ1_avg(i) = (δϕ1(i) + δϕ1(i + 1)) / 2.0

    ## Monotonicity slope limiter
    ϕ1_min(i) = minimum((ϕ1(i - 1), ϕ1(i), ϕ1(i + 1)))
    ϕ1_max(i) = maximum((ϕ1(i - 1), ϕ1(i), ϕ1(i + 1)))
    function Δϕ1_mono(i)
        sign(Δϕ1_avg(i)) *
        minimum((abs(Δϕ1_avg(i)), 2 * (ϕ1(i) - ϕ1_min(i)),
            2 * (ϕ1_max(i) - ϕ1(i))))
    end

    courant(i) = U[i] * Δt / (Δz)

    function FLUX(i)
        ifelse(U[i] >= 0,
            (U[i] * (ϕ1(i + 1) + Δϕ1_mono(i + 1) * (1 - courant(i)) / 2.0)),
            (U[i] * (ϕ1(i + 2) - Δϕ1_mono(i + 2) * (1 + courant(i)) / 2.0)))
    end

    ϕ2(i) = ϕ1(i) - Δt / Δz * (FLUX(i - 1) - FLUX(i - 2))

    for i in 1:nz
        s1d_out[i] = max(zero(eltype(s1d_in)), ϕ2(i + 2))
    end

    return s1d_out
end

"""
PPM advection in 1-D (Collela and Woodward, 1984)
"""
function advect1d_ppm(s1d_in::Vector{Float64}, U::Vector{Float64}, nz, Δt, Δz)
    ϵ = 0.01
    η⁽¹⁾ = 20
    η⁽²⁾ = 0.05

    s1d_out = zeros(Float64, nz)

    ## Initialization
    ϕ1(i) = i < 4 ? s1d_in[1] : i > nz + 3 ? s1d_in[end] : s1d_in[i - 3]

    ## Edge value calculation
    δϕ(i) = 1 / 2 * (ϕ1(i + 1) - ϕ1(i - 1))

    function δₘϕ(i)
        ifelse((ϕ1(i + 1) - ϕ1(i)) * (ϕ1(i) - ϕ1(i - 1)) > 0,
            min(
                abs(δϕ(i - 1)), 2 * abs(ϕ1(i) - ϕ1(i - 1)), 2 * abs(ϕ1(i + 1) - ϕ1(i))) *
            sign(δϕ(i - 1)),
            zero(eltype(s1d_in))
        )
    end

    ϕ₊½(i) = 1 / 2 * (ϕ1(i + 1) + ϕ1(i)) - 1 / 6 * (δₘϕ(i + 1) - δₘϕ(i))

    ## Discontinuity detection
    δ²ϕ(i) = 1 / (6 * Δz^2) * (ϕ1(i + 1) - 2 * ϕ1(i) + ϕ1(i - 1))

    function η_tilde(i)
        ifelse(
            -δ²ϕ(i + 1) * δ²ϕ(i - 1) * abs(ϕ1(i + 1) - ϕ1(i - 1)) -
            ϵ * min(abs(ϕ1(i + 1)), abs(ϕ1(i - 1))) > 0,
            -(δ²ϕ(i + 1) - δ²ϕ(i - 1)) * (Δz^2) / (ϕ1(i + 1) - ϕ1(i - 1)),
            zero(eltype(s1d_in))
        )
    end

    η(i) = clamp(η⁽¹⁾ * (η_tilde(i) - η⁽²⁾), 0, 1)

    ϕLᵈ(i) = ϕ1(i) + 1 / 2 * δₘϕ(i)
    ϕRᵈ(i) = ϕ1(i + 1) + 1 / 2 * δₘϕ(i + 1)

    ϕL₀(i) = ϕ₊½(i) * (1 - η(i)) + ϕLᵈ(i) * η(i)
    ϕR₀(i) = ϕ₊½(i + 1) * (1 - η(i)) + ϕRᵈ(i) * η(i)

    ## Monotonicity examination
    function ϕL(i)
        ifelse((ϕR₀(i) - ϕ1(i)) * (ϕ1(i) - ϕL₀(i)) <= 0,
            ϕ1(i),
            ifelse(
                (ϕR₀(i) - ϕL₀(i)) * (ϕ1(i) - 1 / 2 * (ϕL₀(i) + ϕR₀(i))) >=
                (ϕR₀(i) - ϕL₀(i))^2 / 6,
                3 * ϕ1(i) - 2 * ϕR₀(i),
                ϕL₀(i)
            ))
    end

    function ϕR(i)
        ifelse((ϕR₀(i) - ϕ1(i)) * (ϕ1(i) - ϕL₀(i)) <= 0,
            ϕ1(i),
            ifelse(
                -(ϕR₀(i) - ϕL₀(i))^2 / 6 >
                (ϕR₀(i) - ϕL₀(i)) * (ϕ1(i) - 1 / 2 * (ϕR₀(i) + ϕL₀(i))),
                3 * ϕ1(i) - 2 * ϕL₀(i),
                ϕR₀(i)
            ))
    end

    ## Compute flux
    courant(i) = U[i] * Δt / Δz

    Δϕ(i) = ϕR(i) - ϕL(i)

    ϕ₆(i) = 6 * (ϕ1(i) - 1 / 2 * (ϕL(i) + ϕR(i)))

    function FLUX(i)
        ifelse(U[i] >= 0,
            courant(i) * (ϕR(i + 2) -
             1 / 2 * courant(i) *
             (Δϕ(i + 2) - (1 - 2 / 3 * courant(i)) * ϕ₆(i + 2))),
            courant(i) * (ϕL(i + 3) -
             1 / 2 * courant(i) *
             (Δϕ(i + 3) - (1 + 2 / 3 * courant(i)) * ϕ₆(i + 3)))
        )
    end

    ϕ2(i) = ϕ1(i) + (FLUX(i - 3) - FLUX(i - 2))

    for i in 1:nz
        s1d_out[i] = ϕ2(i + 3)
    end

    return s1d_out
end

"""
Advection in 2-D
Using the output of 1-D advection, we implemented the non-directional splitting 2-D advection. -- LR96 approach (Lin and Rood 1996)
"""
function advect2D(dt, nx, ny, s1, u, v, lons, lats, advect1d)

    ### BC implementation -- Need to revisit here later
    s2d = zeros(nx + 4, ny + 4)
    s2d[3:(nx + 2), 3:(ny + 2)] = s1'

    #bc!(nx, ny, s2d)

    u_edge = zeros(nx + 1, ny)
    v_edge = zeros(nx, ny + 1)

    u_edge[1, :] = u'[1, :]
    u_edge[nx + 1, :] = u'[nx, :]
    v_edge[:, 1] = v'[:, 1]
    v_edge[:, ny + 1] = v'[:, ny]

    for j in 1:ny
        for i in 2:nx
            u_edge[i, j] = (u'[i - 1, j] + u'[i, j]) / 2.0
        end
    end

    for i in 1:nx
        for j in 2:ny
            v_edge[i, j] = (v'[i, j - 1] + v'[i, j]) / 2.0
        end
    end

    ### Populate the arrays
    s1d_xin = zeros(nx + 4)
    s1d_xout = zeros(nx + 4)
    s1d_yin = zeros(ny + 4)
    s1d_yout = zeros(ny + 4)
    vel_x1d = zeros(nx + 1)
    vel_y1d = zeros(ny + 1)

    F_Q = zeros(nx + 4, ny + 4)
    G_Q = zeros(nx + 4, ny + 4)

    F_GQ = zeros(nx + 4, ny + 4)
    G_FQ = zeros(nx + 4, ny + 4)

    s2d_xout = zeros(nx + 4, ny + 4)
    s2d_yout = zeros(nx + 4, ny + 4)

    Re = 6378000
    dy = Re * (lats[2] - lats[1]) * π / 180

    for j in 1:ny
        for i in 1:(nx + 4)
            s1d_xin[i] = s2d[i, j + 2]
        end
        for i in 1:(nx + 1)
            vel_x1d[i] = u_edge[i, j]
        end
        dx = Re * cos(lats[j] * π / 180) * (lons[2] - lons[1]) * π / 180
        s1d_xout = advect1d(s1d_xin, vel_x1d, nx, dt, dx)

        s1d_xout = [
            s1d_xout[begin], s1d_xout[begin], s1d_xout..., s1d_xout[end], s1d_xout[end]]

        F_Q[:, j + 2] = s1d_xout[:]
    end

    for i in 1:nx
        for j in 1:(ny + 4)
            s1d_yin[j] = s2d[i + 2, j]
        end
        for j in 1:(ny + 1)
            vel_y1d[j] = v_edge[i, j]
        end
        s1d_yout = advect1d(s1d_yin, vel_y1d, ny, dt, dy)
        s1d_yout = [
            s1d_yout[begin], s1d_yout[begin], s1d_yout..., s1d_yout[end], s1d_yout[end]]

        G_Q[i + 2, :] = s1d_yout[:]
    end

    for j in 1:ny
        for i in 1:(nx + 4)
            s1d_xin[i] = s2d[i, j + 2] + 0.5 * G_Q[i, j + 2]
        end
        for i in 1:(nx + 1)
            vel_x1d[i] = u_edge[i, j]
        end
        dx = Re * cos(lats[j] * π / 180) * (lons[2] - lons[1]) * π / 180
        s1d_xout = advect1d(s1d_xin, vel_x1d, nx, dt, dx)
        s1d_xout = [
            s1d_xout[begin], s1d_xout[begin], s1d_xout..., s1d_xout[end], s1d_xout[end]]

        F_GQ[:, j + 2] = s1d_xout[:]
    end

    for i in 1:nx
        for j in 1:(ny + 4)
            s1d_yin[j] = s2d[i + 2, j] + 0.5 * F_Q[i + 2, j]
        end
        for j in 1:(ny + 1)
            vel_y1d[j] = v_edge[i, j]
        end
        s1d_yout = advect1d(s1d_yin, vel_y1d, ny, dt, dy)
        s1d_yout = [
            s1d_yout[begin], s1d_yout[begin], s1d_yout..., s1d_yout[end], s1d_yout[end]]

        G_FQ[i + 2, :] = s1d_yout[:]
    end

    for i in 3:(nx + 2)
        for j in 3:(ny + 2)
            s2d[i, j] = s2d[i, j] + F_GQ[i, j] + G_FQ[i, j]
            s1[j - 2, i - 2] = s2d[i, j]
        end
    end

    return s1
end
