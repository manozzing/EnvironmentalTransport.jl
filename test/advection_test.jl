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

    op_f, _ = tensor_advection_op(Float64, (nv, nx, ny, nz), 2, ppm_stencil)
    tensor_adv_x = op_f((i, j, t) -> u[i], (i, j, t) -> Δz, Δt)

    runf(c) = c[:] + (tensor_adv_x * c[:]) .* Δt

    # Norm decreases as pollution spreads out.
    @test norm(runf(c[:])) ≈ 686.4747336938192
    @test norm(runf(runf(c[:]))) ≈ 617.1569463272692
    @test norm(runf(runf(runf(c[:])))) ≈ 554.9089673126576
end

nv = 2
nx = 4
ny = 5
nz = 6
c = zeros(nv, nx, ny, nz)
c[2, 2, 3, 4] = 1.0

@testset "Advection X" begin
    op_f, _ = tensor_advection_op(Float64, (nv, nx, ny, nz), 2, ppm_stencil)
    tensor_adv_x = op_f((i, j, t) -> 1, (i, j, t) -> Δz, Δt)

    c1 = tensor_adv_x * c[:]
    c1_want = zeros(nv, nx, ny, nz)
    c1_want[2, 2, 3, 4] = -2
    c1_want[2, 3, 3, 4] = 2
    @test c1 ≈ c1_want[:]
end

@testset "Advection Y" begin
    op_f, _ = tensor_advection_op(Float64, (nv, nx, ny, nz), 3, ppm_stencil)
    tensor_adv_y = op_f((i, j, t) -> 1, (i, j, t) -> Δz, Δt)

    c1 = tensor_adv_y * c[:]

    c1_want = zeros(nv, nx, ny, nz)
    c1_want[2, 2, 3, 4] = -2
    c1_want[2, 2, 4, 4] = 2

    @test c1 ≈ c1_want[:]
end

@testset "Advection Z" begin
    op_f, _ = tensor_advection_op(Float64, (nv, nx, ny, nz), 4, ppm_stencil)
    tensor_adv_z = op_f((i, j, t) -> 1, (i, j, t) -> Δz, Δt)

    c1 = tensor_adv_z * c[:]
    c1_want = zeros(nv, nx, ny, nz)
    c1_want[2, 2, 3, 4] = -2
    c1_want[2, 2, 3, 5] = 2
    @test c1 ≈ c1_want[:]
end