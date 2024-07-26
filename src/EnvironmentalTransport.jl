module EnvironmentalTransport

using DocStringExtensions
using Tullio
using SciMLOperators
using LinearAlgebra
using SciMLBase: NullParameters
using EarthSciMLBase, EarthSciData

include("advection_stencils.jl")
include("boundary_conditions.jl")
include("orderby.jl")
include("advection.jl")
include("advect1d_vertical.jl")

end
