using Test
using AllocCheck

"""
$(SIGNATURES)

Add zero gradient Newmann boundary conditions to the 1D array `u` with the given stencil size.
Results are optionally stored in the preallocated array `du`.
"""
function zerograd_bcf(u, p, t; stencil_size = (3, 4))
    lpad, rpad = stencil_size
    [repeat([u[begin]], lpad); u; repeat([u[end]], rpad)]
end
function zerograd_bcf(du::T, u::T, p, t; stencil_size = (3, 4))::T where {T}
    lpad, rpad = stencil_size
    @views begin
        du[begin:(begin + lpad - 1)] .= u[begin]
        du[(end - rpad + 1):end] .= u[end]
        #du[(begin + lpad):(end - rpad)] .= u # Allocates memory for some reason.
        @inbounds for i in eachindex(u)
            du[i + lpad] = u[i]
        end
    end
    du
end

@code_warntype zerograd_bcf(zeros(10), [1.0, 2, 3], [1.0], 1.0)

@test zerograd_bcf([1, 2, 3], [1.0], 1.0) == [1, 1, 1, 1, 2, 3, 3, 3, 3, 3]
du = zeros(Int, 10)
zerograd_bcf(du, [1, 2, 3], [1.0], 1.0)
@test du == [1, 1, 1, 1, 2, 3, 3, 3, 3, 3]

@testset "Allocations" begin
    @check_allocs ff(dx, x, p) = zerograd_bcf(dx, x, p, 1.0)
    zerograd_bcf(zeros(11), zeros(3), zeros(1), 1.0)
    @test_nowarn ff(zeros(11), zeros(3), zeros(1))
end

u = [1.0, 2, 3]
u2 = zeros(10)

function zerograd_bc(dtype, vec_length, stencil)
    lpad, rpad = stencil_size(stencil)
    FunctionOperator(
        zerograd_bcf, zeros(dtype, vec_length), zeros(dtype, vec_length + lpad + rpad);
        p = [1.0], t = 1.0, accepted_kwargs = (:stencil_size,))
end

zbc_op = zerograd_bc(Float64, 3, ppm_stencil)
@test zbc_op(u, [1.0], 1.0; stencil_size = (3, 4)) ≈ [1.0, 1, 1, 1, 2, 3, 3, 3, 3, 3]
@test begin
    zbc_op(u2, u, [1.0], 1.0; stencil_size = (3, 4))
    u2 ≈ [1.0, 1, 1, 1, 2, 3, 3, 3, 3, 3]
end

"""
$(SIGNATURES)

Advect the 1D array `u` with the given velocity field `v` and timestep `Δt`, and 
grid spacing `Δz`, which should be stored in the parameter array `p` as `p = (v, Δt, Δz)`.
The stencil algorithm used for the advection can be specified with the `stencil` keyword argument.

Results are optionally stored in the preallocated array `du`.
"""
function advect_1d(u, p, t; stencil = ppm_stencil)
    v, Δt, Δz = p
    lpad, rpad = stencil_size(stencil)
    [stencil(u[(i - lpad):(i + rpad)], v[(i - lpad):(i - lpad + 1)], Δt, Δz)
     for i in (firstindex(u) + lpad):(lastindex(u) - rpad)]
end
function advect_1d(du, u, p, t; stencil = ppm_stencil)
    v, Δt, Δz = p
    lpad, rpad = stencil_size(stencil)
    for i in (firstindex(u) + lpad):(lastindex(u) - rpad)
        du[i - lpad] = stencil(
            u[(i - lpad):(i + rpad)], v[(i - lpad):(i - lpad + 1)], Δt, Δz)
    end
    du
end

c = [0.0, 1, 2, 3, 4, 5]
v = [10.0, 8, 6, 4, 2, 0, 1]
Δt = 0.05
Δz = 0.5

@testset "ppm advection" begin
    c2 = [c[1], c[1], c[1], c..., c[end], c[end], c[end], c[end]]
    result = advect_1d(c2, (v, Δt, Δz), 1.0)
    @test result ≈ [0.0, 0.3999999999999999, 1.8, 3.2, 4.6, 4.5]

    result .= 0
    advect_1d(result, c2, (v, Δt, Δz), 1.0)
    @test result ≈ [0.0, 0.3999999999999999, 1.8, 3.2, 4.6, 4.5]
end

using Plots

size(zbc_op)
u = collect(1.0:9)

b = [1 0 0
     1.0 0 0
     0 1 0
     0 0 1
     0 0 1]
d = DiagonalOperator(ones(3))


function permute_op(dtype, shape, perm)
    f(u, p, t) = permutedims(u, perm)
    f(du, u, p, t) = permutedims!(du, u, perm)
    x = zeros(dtype, shape...)
    y = f(x, nothing, nothing)
    FunctionOperator(f, x, y; islinear=true)
end

@testset "permute" begin
    p_op = permute_op(Float64, (3,3), (2,1))
    x = 1:9
    @test reshape(p_op(x, nothing, nothing),3,3) == transpose(reshape(x,3,3))
end

permutedims(rand(2,2), (1,2))

@testset "permutation_matrix" begin
    for i in 1:10
        x = rand(i, i)
        xx = reshape(permutation_matrix(i) * reshape(x, :), i, i)
        @test xx ≈ transpose(x)
        @test xx ≈ permutedims(x, (2,1))
        xxx = reshape(permutation_matrix(i) * reshape(xx, :), i, i)
        @test xxx ≈ x
    end
end

plot(
    heatmap(reshape(u, 3,3)),
    heatmap(reshape(permutation_matrix(3) * u, 3,3))
)


xx = TensorProductOperator(b, d) * u
yy = permute_op(Float64, (3,5), (2,1)) * TensorProductOperator(b, d) * permute_op(Float64, (3,3), (2,1)) *  u

plot(
    heatmap(reshape(xx, 3, 5)),
    heatmap(reshape(yy, 5, 3))
)

a = [1.0, 2, 3]

using LinearAlgebra

b * a