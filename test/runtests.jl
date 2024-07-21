using EnvironmentalTransport
using Test, SafeTestsets

@testset "EnvironmentalTransport.jl" begin
    @safetestset "Horizontal Advection" begin include("horizontal_advection_test.jl") end
    @safetestset "Vertical Advection" begin include("advect1d_vertical_test.jl") end
end