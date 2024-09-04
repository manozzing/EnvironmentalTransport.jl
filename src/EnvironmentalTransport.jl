module EnvironmentalTransport

using DocStringExtensions
using SciMLOperators
using LinearAlgebra
using SciMLBase: NullParameters
using EarthSciMLBase, EarthSciData

include("advection_stencils.jl")
include("boundary_conditions.jl")
include("advection.jl")

end
