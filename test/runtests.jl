using EnvironmentalTransport
using Test, SafeTestsets

@testset "EnvironmentalTransport.jl" begin
    @safetestset "Advection Stencils" begin include("advection_stencil_test.jl") end
    @safetestset "Boundary Conditions" begin include("boundary_conditions_test.jl") end
    @safetestset "Order By" begin include("orderby_test.jl") end
    @safetestset "Advection" begin include("advection_test.jl") end
    @safetestset "Vertical Advection" begin include("advect1d_vertical_test.jl") end
end