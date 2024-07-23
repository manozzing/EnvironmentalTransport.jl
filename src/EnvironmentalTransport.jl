module EnvironmentalTransport

using DocStringExtensions
using Tullio
using SciMLOperators

include("horizontal_advection.jl")
include("advect1d_vertical.jl")

end
