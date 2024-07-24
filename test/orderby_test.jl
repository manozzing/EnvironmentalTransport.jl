using EnvironmentalTransport: orderby_op
using Test

q = collect(1:12)

@testset "dim 1" begin
    oop, idx_f = orderby_op(Float64, [2, 3, 2], 1)
    oop_inv = inv(oop)
    r = oop * q
    @test r[:] ≈ q
    @test oop_inv * r ≈ q
    @test idx_f(3) ==  [CartesianIndex(1, 3, 1), CartesianIndex(2, 3, 1)]
end

@testset "dim 2" begin
    oop, idx_f = orderby_op(Float64, [2, 3, 2], 2)
    oop_inv = inv(oop)
    r = oop * q
    @test r[:] ≈ reshape(permutedims(reshape(q, 2, 3, 2), (2, 1, 3)), 3, :)[:]
    @test oop_inv * r ≈ q
    @test idx_f(3) == [CartesianIndex(1, 1, 2), CartesianIndex(1, 2, 2), CartesianIndex(1, 3, 2)]
end

@testset "dim 3" begin
    oop, idx_f = orderby_op(Float64, [2, 3, 2], 3)
    oop_inv = inv(oop)
    r = oop * q
    @test r[:] ≈ reshape(permutedims(reshape(q, 2, 3, 2), (3, 1, 2)), 2, :)[:]
    @test oop_inv * r ≈ q
    @test idx_f(4) == [CartesianIndex(2, 2, 1), CartesianIndex(2, 2, 2)]
end