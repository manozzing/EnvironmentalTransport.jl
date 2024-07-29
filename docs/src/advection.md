# Numerical Advection Operator

We have two ways to represent phenomena that occur across space such as advection: through symbolically-defined partial differential equation systems, which will be covered elsewhere in
documentation, and through numerically-implemented algorithms.
This is an example of the latter. (Currently, symbolically defined PDEs are too slow to be
used in large-scale simulations.)

To demonstrate how it works, let's first set up our environment:

```@example adv
using EnvironmentalTransport
using EarthSciMLBase, EarthSciData
using ModelingToolkit, DomainSets, DifferentialEquations
using Distributions, LinearAlgebra
using Dates
using NCDatasets, Plots
nothing #hide
```

## Emissions

Next, let's set up an emissions scenario to advect.
We have some emissions centered around Portland, starting at the beginning of the simulation and then tapering off:

```@example adv
@parameters lon=0.0 lat=0.0 lev=1.0 t

function emissions(t, μ_lon, μ_lat, σ)
    @variables c(t) = 0.0
    dist =MvNormal([starttime, μ_lon, μ_lat, 1], Diagonal(map(abs2, [3600.0, σ, σ, 1])))
    D = Differential(t)
    ODESystem([D(c) ~ pdf(dist, [t, lon, lat, lev]) * 50],
        t, name = :Test₊emissions)
end

emis = emissions(t, deg2rad(-122.6), deg2rad(45.5), 0.1)
```

## Domain

Next, let's set up a spatial and temporal domain for our simulation.

```@example adv
starttime = datetime2unix(DateTime(2022, 5, 1, 0, 0))
endtime = datetime2unix(DateTime(2022, 5, 1, 0, 30))

domain = DomainInfo(
    [partialderivatives_δxyδlonlat,
        partialderivatives_δPδlev_geosfp(geosfp)],
    constIC(16.0, t ∈ Interval(starttime, endtime)),
    constBC(16.0,
        lon ∈ Interval(deg2rad(-130.0), deg2rad(-60.0)),
        lat ∈ Interval(deg2rad(9.75), deg2rad(60.0)),
        lev ∈ Interval(1, 15)))
nothing # hide
```

## Coupled System

Now, let's set up some input data from GEOS-FP to get wind fields for our advection.
We need to use `coord_defaults` in this case to get it to work correctly, but 
it doesn't matter what the defaults are.
We also set up an [outputter](https://data.earthsci.dev/stable/api/#EarthSciData.NetCDFOutputter) to save the results of our simulation, and couple the components we've created so far into a 
single system.

```@example adv
geosfp = GEOSFP("4x5", t; dtype = Float64,
    coord_defaults = Dict(:lon => 0.0, :lat => 0.0, :lev => 1.0))


outfile = "out.nc"
output = NetCDFOutputter(outfile, 3600.0)

csys = couple(emis, domain, geosfp, output) 
```
## Advection Operator

Next, we create an [`AdvectionOperator`](@ref) to perform advection. 
We need to specify a time step (600 s in this case), as stencil algorithm to do the advection (current options are [`l94_stencil`](@ref) and [`ppm_stencil`](@ref)), and a time integration scheme (`SSPRK22` in this case).
Refer [here](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) for the available time integrator choices.

Then, we couple the advection operator to the rest of the system.

!!! warning
    The advection operator will automatically couple itself to available wind fields such as those from GEOS-FP, but the wind-field component (e.g.. `geosfp`) must already be present
    in the coupled system for this to work correctly.

```@example adv
adv = AdvectionOperator(600.0, l94_stencil, SSPRK22())

csys = couple(csys, adv)
```
Now, we initialize a [`Simulator`](https://base.earthsci.dev/dev/simulator/) to run our demonstration. 
We specify a horizontal resolution of 4 degrees and a vertical resolution of 1 level, and use the `Tsit5` time integrator for our emissions system of equations.
Then, we run the simulation.

```@example adv
sim = Simulator(csys, [deg2rad(4), deg2rad(4), 1], Tsit5())

@time run!(sim)
```

## Visualization

Finally, we can visualize the results of our simulation:

```@example adv
ds = NCDataset(outfile, "r")

anim = @animate for i ∈ 1:size(ds["Test₊emissions₊c"])[4]
    plot(
        heatmap(ds["Test₊emissions₊c"][:, :, 1, i]', title="Ground-Level"),
        heatmap(ds["Test₊emissions₊c"][:, 10, :, i]', title="Vertical Cross-Section"),
    )
end
gif(anim, fps = 15)
```

```@setup adv
rm(outfile, force=true)
```