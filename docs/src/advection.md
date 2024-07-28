# Numerical Advection Operator



```@example adv
using EnvironmentalTransport
using EarthSciMLBase, EarthSciData
using ModelingToolkit, DomainSets, DifferentialEquations
using Distributions, LinearAlgebra
using Dates
using NCDatasets, Plots

@parameters lon=0.0 lat=0.0 lev=1.0 t
starttime = datetime2unix(DateTime(2022, 5, 1, 0, 0))
endtime = datetime2unix(DateTime(2022, 5, 1, 3, 0))

# Need to add coord_defaults or it won't work.
geosfp = GEOSFP("4x5", t; dtype = Float64,
    coord_defaults = Dict(:lon => 0.0, :lat => 0.0, :lev => 1.0))

domain = DomainInfo(
    [partialderivatives_δxyδlonlat,
        partialderivatives_δPδlev_geosfp(geosfp)],
    constIC(16.0, t ∈ Interval(starttime, endtime)),
    constBC(16.0,
        lon ∈ Interval(deg2rad(-130.0), deg2rad(-60.0)),
        lat ∈ Interval(deg2rad(9.75), deg2rad(60.0)),
        lev ∈ Interval(1, 15)))

function emissions(t, μ_lon, μ_lat, σ)
    @variables c(t) = 0.0
    dist =MvNormal([starttime, μ_lon, μ_lat, 1], Diagonal(map(abs2, [3600.0, σ, σ, 1])))
    D = Differential(t)
    ODESystem([D(c) ~ pdf(dist, [t, lon, lat, lev]) * 50],
        t, name = :Test₊emissions)
end

emis = emissions(t, deg2rad(-122.6), deg2rad(45.5), 0.1)
outfile = tempname() * ".nc"
output = NetCDFOutputter(outfile, 3600.0)

adv = AdvectionOperator(600.0, l94_stencil, SSPRK22())

# Make sure that geos-fp is added first
csys = couple(emis, domain, geosfp, output)

csys = couple(csys, adv)

sim = Simulator(csys, [deg2rad(4), deg2rad(4), 1], Tsit5())

run!(sim)

ds = NCDataset(outfile, "r")


anim = @animate for i ∈ 1:size(ds["Test₊emissions₊c"])[4]
    plot(
        heatmap(ds["Test₊emissions₊c"][:, :, 1, i]', title="Ground-Level"),
        heatmap(ds["Test₊emissions₊c"][:, 10, :, i]', title="Vertical Cross-Section"),
    )
end
gif(anim, fps = 15)

rm(outfile, force=true)
```