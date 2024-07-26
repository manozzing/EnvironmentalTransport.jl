
using EnvironmentalTransport
using EarthSciMLBase, EarthSciData
using ModelingToolkit, DomainSets, OrdinaryDiffEq
using Distributions, LinearAlgebra
using Dates

@parameters lon=0.0 lat=0.0 lev=1.0 t
lat = GlobalScope(lat)
lon = GlobalScope(lon)
lev = GlobalScope(lev)
starttime = datetime2unix(DateTime(2022, 5, 1))
endtime = datetime2unix(DateTime(2022, 5, 1, 1, 0, 5))

geosfp = GEOSFP("4x5", t; dtype = Float64)

domain = DomainInfo(
    [partialderivatives_δxyδlonlat,
        partialderivatives_δPδlev_geosfp(geosfp)],
    constIC(16.0, t ∈ Interval(starttime, endtime)),
    constBC(16.0,
        lon ∈ Interval(deg2rad(-130.0), deg2rad(-60.0)),
        lat ∈ Interval(deg2rad(9.75), deg2rad(60.0)),
        lev ∈ Interval(1, 3)))

function emissions(t, μ_lon, μ_lat, σ)
    @variables c(t) = 0.0
    dist =MvNormal([starttime, μ_lon, μ_lat, 1], Diagonal(map(abs2, [3600.0, σ, σ, 1])))
    D = Differential(t)
    ODESystem([D(c) ~ pdf(dist, [t, lon, lat, lev]) * 50],
        t, name = :Test₊emissions)
end

emis = emissions(t, deg2rad(-122.6), deg2rad(45.5), 0.1)
output = NetCDFOutputter("out.nc", 3600.0)

csys = couple(emis, domain, geosfp, output)

sim = Simulator(csys, [deg2rad(4), deg2rad(4), 1], Tsit5())

run!(sim)

@test norm(sim.u) ≈ 451.09230204187736

op = AdvectionOperator(100.0, l94_stencil, SSPRK22())

@test isnothing(op.vardict) # Before coupling, there shouldn't be anything here.

csys = couple(csys, op)

@test !isnothing(op.vardict) # after coupling, there should be something here.

run!(sim)
