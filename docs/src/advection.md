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
using ModelingToolkit: t, D
using DynamicQuantities
using Distributions, LinearAlgebra
using Dates
using NCDatasets, Plots
nothing #hide
```

## Emissions

Next, let's set up an emissions scenario to advect.
We have some emissions centered around Portland, starting at the beginning of the simulation and then tapering off:

```@example adv
starttime = datetime2unix(DateTime(2022, 5, 1, 0, 0))
endtime = datetime2unix(DateTime(2022, 6, 1, 0, 0))

@parameters(
    lon=-97.0, [unit=u"rad"],
    lat=30.0, [unit=u"rad"],
    lev=1.0,
)

function emissions(μ_lon, μ_lat, σ)
    @variables c(t) = 0.0 [unit=u"kg"]
    @constants v_emis = 50.0 [unit=u"kg/s"]
    @constants t_unit = 1.0 [unit=u"s"] # Needed so that arguments to `pdf` are unitless.
    dist = MvNormal([starttime, μ_lon, μ_lat, 1], Diagonal(map(abs2, [3600.0*24*3, σ, σ, 1])))
    ODESystem([D(c) ~ pdf(dist, [t/t_unit, lon, lat, lev]) * v_emis],
        t, name = :emissions)
end

emis = emissions(deg2rad(-122.6), deg2rad(45.5), deg2rad(1))
```

## Coupled System

Next, let's set up a spatial and temporal domain for our simulation, and
some input data from GEOS-FP to get wind fields for our advection.
We need to use `coord_defaults` in this case to get the GEOS-FP data to work correctly, but 
it doesn't matter what the defaults are.
We also set up an [outputter](https://data.earthsci.dev/stable/api/#EarthSciData.NetCDFOutputter) to save the results of our simulation, and couple the components we've created so far into a 
single system.

```@example adv
geosfp, geosfp_updater = GEOSFP("0.5x0.625_NA"; dtype = Float64,
    coord_defaults = Dict(:lon => -97.0, :lat => 30.0, :lev => 1.0))

domain = DomainInfo(
    [partialderivatives_δxyδlonlat,
        partialderivatives_δPδlev_geosfp(geosfp)],
    constIC(16.0, t ∈ Interval(starttime, endtime)),
    constBC(16.0,
        lon ∈ Interval(deg2rad(-129), deg2rad(-61)),
        lat ∈ Interval(deg2rad(11), deg2rad(59)),
        lev ∈ Interval(1, 30)),
    dtype = Float64)

outfile = ("RUNNER_TEMP" ∈ keys(ENV) ? ENV["RUNNER_TEMP"] : tempname()) * "out.nc" # This is just a location to save the output.
output = NetCDFOutputter(outfile, 3600.0)

csys = couple(emis, domain, geosfp, geosfp_updater, output) 
```
## Advection Operator

Next, we create an [`AdvectionOperator`](@ref) to perform advection. 
We need to specify a time step (600 s in this case), as stencil algorithm to do the advection (current options are [`l94_stencil`](@ref) and [`ppm_stencil`](@ref)).

Then, we couple the advection operator to the rest of the system.

!!! warning
    The advection operator will automatically couple itself to available wind fields such as those from GEOS-FP, but the wind-field component (e.g.. `geosfp`) must already be present
    in the coupled system for this to work correctly.

```@example adv
adv = AdvectionOperator(300.0, upwind1_stencil, ZeroGradBCArray)

csys = couple(csys, adv)
```
Now, we initialize a [`Simulator`](https://base.earthsci.dev/dev/simulator/) to run our demonstration. 
We specify a horizontal resolution of 4 degrees and a vertical resolution of 1 level, and use the `Tsit5` time integrator for our emissions system of equations, and a time integration scheme for our advection operator (`SSPRK22` in this case).
Refer [here](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) for the available time integrator choices.
We also choose a operator splitting interval of 600 seconds.
Then, we run the simulation.

```@example adv
sim = Simulator(csys, [deg2rad(1), deg2rad(1), 1])
st = SimulatorStrangThreads(Tsit5(), SSPRK22(), 300.0)

@time run!(sim, st, save_on=false, save_start=false, save_end=false, 
    initialize_save=false)
```

## Visualization

Finally, we can visualize the results of our simulation:

```@example adv
ds = NCDataset(outfile, "r")

imax = argmax(reshape(maximum(ds["emissions₊c"][:, :, :, :], dims=(1, 3, 4)), :))
anim = @animate for i ∈ 1:size(ds["emissions₊c"])[4]
    plot(
        heatmap(rad2deg.(sim.grid[1]), rad2deg.(sim.grid[2]), 
            ds["emissions₊c"][:, :, 1, i]', title="Ground-Level", xlabel="Longitude", ylabel="Latitude"),
        heatmap(rad2deg.(sim.grid[1]), sim.grid[3], ds["emissions₊c"][:, imax, :, i]', 
            title="Vertical Cross-Section (lat=$(round(rad2deg(sim.grid[2][imax]), digits=1)))", 
            xlabel="Longitude", ylabel="Vertical Level"),
    )
end
gif(anim, fps = 15)
```

```@setup adv
rm(outfile, force=true)
```