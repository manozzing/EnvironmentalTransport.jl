module EnvironmentalTransport

using DocStringExtensions
using SciMLOperators, OrdinaryDiffEq
using LinearAlgebra
using SciMLBase: NullParameters
using EarthSciMLBase, EarthSciData

include("advection_stencils.jl")
include("boundary_conditions.jl")
include("orderby.jl")
include("advection.jl")

end
