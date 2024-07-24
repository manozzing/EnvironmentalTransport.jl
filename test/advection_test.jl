using EnvironmentalTransport: advect_1d_op, tensor_advection_op
using EnvironmentalTransport
using Test
using LinearAlgebra

c = [0.0, 1, 2, 3, 4, 5]
v = [10.0, 8, 6, 4, 2, 0, 1]
Δt = 0.05
Δz = 0.5

@testset "1D advection" begin
    c2 = [c[1], c[1], c[1], c..., c[end], c[end], c[end], c[end]]
    adv_op = advect_1d_op(
        Float64, (length(c),), ppm_stencil, (i, t) -> v[i], (i, t) -> Δz, Δt)
    result = adv_op(c2, nothing, nothing)
    @test result ≈ [0.0, -12, -4, 4, 12, -10.0]

    result .= 0
    adv_op(result, c2, nothing, 1.0)
    @test result ≈ [0.0, -12, -4, 4, 12, -10.0]

    c3 = repeat(c2, 1, 3)
    adv_op = advect_1d_op(Float64, (length(c), 3), ppm_stencil,
        (i, j, t) -> v[i], (i, j, t) -> Δz, Δt)
    result = adv_op(c3, nothing, nothing)
    @test result ≈ repeat([0.0, -12, -4, 4, 12, -10.0], 1, 3)

    result .= 0
    adv_op(result, c3, nothing, 1.0)
    @test result ≈ repeat([0.0, -12, -4, 4, 12, -10.0], 1, 3)
end

@testset "Tensor advection" begin
    # We have a 4D grid, 1st dimension is state variables, 2=x direction, 3=y direction, 4=vertical levels
    nv = 2
    nx = 5
    ny = 4
    nz = 3
    c = reshape(1.0:*(nv, nx, ny, nz), nv, nx, ny, nz)
    u = (1.0:(nx + 1)) .- 3

    tensor_adv_x, _ = tensor_advection_op(Float64, (nv, nx, ny, nz), 2, ppm_stencil,
        (i, j, t) -> u[i], (i, j, t) -> Δz, Δt
    )

    runf(c) = c[:] + (tensor_adv_x * c[:]) .* Δt

    @test norm(runf(c[:])) ≈ 665.961920833316
    @test norm(runf(runf(c[:]))) ≈ 582.6715980721902
    @test norm(runf(runf(runf(c[:])))) ≈ 511.3511651771215
end