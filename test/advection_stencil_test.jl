using Main.EnvironmentalTransport

using Test

c = [0.0, 1, 2, 3, 4, 5]
v = [10.0, 8, 6, 4, 2, 0, 1]
Δt = 0.05
Δz = 0.5

@testset "l94 1" begin
    c2 = [c[1], c[1], c..., c[end], c[end]]
    result = [l94_stencil(c2[(i - 2):(i + 2)], v[(i - 2):(i - 1)], Δt, Δz) for i in 3:8]
    @test max.((0,), c .+ result .* Δt) ≈  [0.0, 0.28, 1.8, 3.24, 4.68, 4.5]
end

@testset "ppm 1" begin
    c2 = [c[1], c[1], c[1], c..., c[end], c[end], c[end], c[end]]
    result = [ppm_stencil(c2[(i - 3):(i + 4)], v[(i - 3):(i - 2)], Δt, Δz) for i in 4:9]
    @test c .+ result .* Δt ≈ [0.0, 0.3999999999999999, 1.8, 3.2, 4.6, 4.5]
end

@testset "upwind1 1" begin
    c2 = [c[1], c..., c[end]]
    result = [upwind1_stencil(c2[(i - 1):(i + 1)], v[(i - 1):i], Δt, Δz) for i in 2:7]
    @test c .+ result .* Δt ≈ [0.0, 1.8, 2.6, 3.4, 4.2, 5.0]
end

@testset "upwind2 1" begin
    c2 = [c[1], c[1], c..., c[end], c[end]]
    result = [upwind2_stencil(c2[(i - 2):(i + 2)], v[(i - 2):(i - 1)], Δt, Δz) for i in 3:8]
    @test max.((0,), c .+ result .* Δt) ≈  [0.0, 2.2, 2.6, 3.4, 4.2, 5.0]
end


c = [6.0, 6, 5, 5, 6, 6]
v = [2.0, 2, 2, 2, 2, 2, 2]

@testset "l94 2" begin
    c2 = [c[1], c[1], c..., c[end], c[end]]
    result = [l94_stencil(c2[(i - 2):(i + 2)], v[(i - 2):(i - 1)], Δt, Δz) for i in 3:8]
    @test max.((0,), c .+ result .* Δt) ≈ [6.0, 6.0, 5.2, 5.0, 5.8, 6.0]
end

@testset "ppm 2" begin
    c2 = [c[1], c[1], c[1], c..., c[end], c[end], c[end], c[end]]
    result = [ppm_stencil(c2[(i - 3):(i + 4)], v[(i - 3):(i - 2)], Δt, Δz) for i in 4:9]
    @test c .+ result .* Δt ≈ [6.0, 6.0, 5.2, 5.0, 5.8, 6.0]
end

@testset "upwind1 2" begin
    c2 = [c[1], c..., c[end]]
    result = [upwind1_stencil(c2[(i - 1):(i + 1)], v[(i - 1):i], Δt, Δz) for i in 2:7]
    @test c .+ result .* Δt ≈ [6.0, 6.0, 4.8, 5.0, 6.2, 6.0]
end

@testset "upwind2 2" begin
    c2 = [c[1], c[1], c..., c[end], c[end]]
    result = [upwind2_stencil(c2[(i - 2):(i + 2)], v[(i - 2):(i - 1)], Δt, Δz) for i in 3:8]
    @test max.((0,), c .+ result .* Δt) ≈  [6.0, 6.0, 4.7, 5.1, 6.3, 5.9]
end
