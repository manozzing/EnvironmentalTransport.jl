using Test
using AllocCheck

"""
$(SIGNATURES)

Add zero gradient Newmann boundary conditions to the 1D array `u` with the given stencil size.
Results are optionally stored in the preallocated array `du`.
"""
function zerograd_bc(vec_length, stencil)
    lpad, rpad = stencil_size(stencil)
    padmat = vcat(
        repeat(hcat(1, repeat([0], vec_length - 1)...), lpad),
        I(vec_length),
        repeat(hcat(repeat([0], vec_length - 1)..., 1), rpad)
    )
    MatrixOperator(padmat)
end

u = [1.0, 2, 3]
u2 = zeros(10)

zbc_op = zerograd_bc(3, ppm_stencil)
@test zbc_op(u, [1.0], 1.0) ≈ [1.0, 1, 1, 1, 2, 3, 3, 3, 3, 3]
@test begin
    zbc_op(u2, u, [1.0], 1.0)
    u2 ≈ [1.0, 1, 1, 1, 2, 3, 3, 3, 3, 3]
end

"""
$(SIGNATURES)

Create an operator to perform 1D advection.

Arguments:
    * `dtype`: The data type of the input and output arrays, e.g. `Float64` or `Float32`.

Advect the 1D array `u` with the given velocity field `v` and timestep `Δt`, and 
grid spacing `Δz`, which should be stored in the parameter array `p` as `p = (v, Δt, Δz)`.
The stencil algorithm used for the advection can be specified with the `stencil` keyword argument.

Results are optionally stored in the preallocated array `du`.
"""
function advect_1d_op(dtype, shape, stencil, v_f, Δt_f, Δz_f)
    lpad, rpad = stencil_size(stencil)
    function f(u, p, t) # Out-of-place, 1D
        [stencil(u[(i - lpad):(i + rpad)],
             (v_f(i - lpad, t), v_f(i - lpad + 1, t)), Δt_f(i, t), Δz_f(i, t))
         for i in (firstindex(u) + lpad):(lastindex(u) - rpad)]
    end
    function f(du, u, p, t) # In-place, 1D
        for i in (firstindex(u) + lpad):(lastindex(u) - rpad)
            du[i - lpad] = stencil(
                u[(i - lpad):(i + rpad)], (v_f(i - lpad, t), v_f(i - lpad + 1, t)), Δt_f(
                    i, t), Δz_f(i, t))
        end
        du
    end
    function f(u::Matrix, p, t) # Out-of-place, 2D
        hcat([[stencil(col[(i - lpad):(i + rpad)],
                   (v_f(i - lpad, j, t), v_f(i - lpad + 1, j, t)), Δt_f(i, j, t), Δz_f(
                       i, j, t))
               for i in (firstindex(col) + lpad):(lastindex(col) - rpad)]
              for (j, col) in enumerate(eachcol(u))]...)
    end
    function f(du::Matrix, u::Matrix, p, t) # In-place, 2D
        @views begin
            for j in 1:size(u, 2)
                ddu = du[:, j]
                uu = u[:, j]
                for i in (firstindex(uu) + lpad):(lastindex(uu) - rpad)
                    ddu[i - lpad] = stencil(
                        uu[(i - lpad):(i + rpad)],
                        (v_f(i - lpad, j, t), v_f(i - lpad + 1, j, t)),
                        Δt_f(i, j, t), Δz_f(i, j, t)
                    )
                end
            end
        end
        du
    end
    x = zeros(dtype, shape[1] + lpad + rpad, shape[2:end]...)
    y = f(x, nothing, nothing)
    FunctionOperator(f, x, y, batch = true)
end

c = [0.0, 1, 2, 3, 4, 5]
v = [10.0, 8, 6, 4, 2, 0, 1]
Δt = 0.05
Δz = 0.5

@testset "ppm advection" begin
    c2 = [c[1], c[1], c[1], c..., c[end], c[end], c[end], c[end]]
    adv_op = advect_1d_op(
        Float64, (length(c),), ppm_stencil, (i, t) -> v[i], (i, t) -> Δt, (i, t) -> Δz)
    result = adv_op(c2, nothing, nothing)
    @test result ≈ [0.0, 0.3999999999999999, 1.8, 3.2, 4.6, 4.5]

    result .= 0
    adv_op(result, c2, nothing, 1.0)
    @test result ≈ [0.0, 0.3999999999999999, 1.8, 3.2, 4.6, 4.5]

    c3 = repeat(c2, 1, 3)
    adv_op = advect_1d_op(Float64, (length(c), 3), ppm_stencil,
        (i, j, t) -> v[i], (i, j, t) -> Δt, (i, j, t) -> Δz)
    result = adv_op(c3, nothing, nothing)
    @test result ≈ repeat([0.0, 0.3999999999999999, 1.8, 3.2, 4.6, 4.5], 1, 3)

    result .= 0
    adv_op(result, c3, nothing, 1.0)
    @test result ≈ repeat([0.0, 0.3999999999999999, 1.8, 3.2, 4.6, 4.5], 1, 3)
end

function permute_op(dtype, shape, perm)
    f(u, p, t) = permutedims(u, perm)
    f(du, u, p, t) = permutedims!(du, u, perm)
    x = zeros(dtype, shape...)
    y = f(x, nothing, nothing)
    FunctionOperator(f, x, y; islinear = true)
end

@testset "permute" begin
    p_op = permute_op(Float64, (3, 3), (2, 1))
    x = 1:9
    @test reshape(p_op(x, nothing, nothing), 3, 3) == transpose(reshape(x, 3, 3))
end

using Plots

size(zbc_op)
u = collect(1.0:9)

d = DiagonalOperator(ones(3))

xx = TensorProductOperator(zbc_op, d) * u
yy = permute_op(Float64, (3, 10), (2, 1)) * TensorProductOperator(zbc_op, d) *
     permute_op(Float64, (3, 3), (2, 1)) * u

@test reshape(xx, 3, 10) ≈ [1.0 1.0 1.0 1.0 4.0 7.0 7.0 7.0 7.0 7.0
       2.0 2.0 2.0 2.0 5.0 8.0 8.0 8.0 8.0 8.0
       3.0 3.0 3.0 3.0 6.0 9.0 9.0 9.0 9.0 9.0]

@test reshape(yy, 10, 3) ≈ [1.0 4.0 7.0
       1.0 4.0 7.0
       1.0 4.0 7.0
       1.0 4.0 7.0
       2.0 5.0 8.0
       3.0 6.0 9.0
       3.0 6.0 9.0
       3.0 6.0 9.0
       3.0 6.0 9.0
       3.0 6.0 9.0]

plot(
    heatmap(reshape(xx, 3, 10)),
    heatmap(reshape(yy, 10, 3))
)

uu = 1.0:3

zbc_op * uu

function orderby_op(dtype, shape::AbstractVector, index::Int)
    vec_length = *(shape...)
    ii = 1:length(shape)
    iiremainder = setdiff(ii, index)
    fwd_idx = tuple(index, iiremainder...)
    function fwd(u, p, t)
        reshape(permutedims(reshape(u, shape...), fwd_idx), shape[index], :)
    end
    newshape = tuple(shape[[index, iiremainder...]]...)
    rev_idx = tuple(insert!(collect(2:length(shape)), index, 1)...)
    function rev(u, p, t)
        permutedims(reshape(u, newshape...), rev_idx)[:]
    end
    x = zeros(dtype, shape[1], Int(vec_length / shape[1]))
    y = zeros(dtype, shape[index], Int(vec_length / shape[index]))
    FunctionOperator(fwd, x, x; op_adjoint = rev, op_inverse = rev,
        op_adjoint_inverse = fwd, islinear = true)
end

q = collect(1:12)

@testset "dim 1" begin
    oop = orderby_op(Float64, [2, 3, 2], 1)
    oop_inv = inv(oop)
    r = oop * q
    @test r[:] ≈ q
    @test oop_inv * r ≈ q
end

@testset "dim 2" begin
    oop = orderby_op(Float64, [2, 3, 2], 2)
    oop_inv = inv(oop)
    r = oop * q
    @test r[:] ≈ reshape(permutedims(reshape(q, 2, 3, 2), (2, 1, 3)), 3, :)[:]
    @test oop_inv * r ≈ q
end

@testset "dim 3" begin
    oop = orderby_op(Float64, [2, 3, 2], 3)
    oop_inv = inv(oop)
    r = oop * q
    @test r[:] ≈ reshape(permutedims(reshape(q, 2, 3, 2), (3, 1, 2)), 2, :)[:]
    @test oop_inv * r ≈ q
end

function tensor_advection_op(dtype, shape, index, stencil; bc_opf=zerograd_bc)
    shape = [shape...]
    order = orderby_op(dtype, shape, index)
    xbc = bc_opf(shape[index], stencil)
    ncols = *(shape...) ÷ shape[index]
    adv_op = advect_1d_op(dtype, (shape[index], ncols), stencil,
        (i, j, t) -> u[i], (i, j, t) -> Δt, (i, j, t) -> Δz)
    
    inv(order) * TensorProductOperator((adv_op * xbc), I(ncols)) * order
end

# We have a 4D grid, 1st dimension is state variables, 2=x direction, 3=y direction, 4=vertical levels
nv = 2
nx = 5
ny = 4
nz = 3
c = reshape(1.0:*(nv, nx, ny, nz), nv, nx, ny, nz)
u = (1.0:(nx + 1)) .- 3

tensor_adv_x = tensor_advection_op(Float64, (nv, nx, ny, nz), 2, ppm_stencil)

@test norm(tensor_adv_x * c[:]) ≈ 665.961920833316
@test norm(tensor_adv_x * tensor_adv_x * c[:]) ≈ 582.6715980721902
@test norm(tensor_adv_x * tensor_adv_x *tensor_adv_x * c[:]) ≈ 511.3511651771215
