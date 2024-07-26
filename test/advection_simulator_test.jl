using EnvironmentalTransport
using EnvironmentalTransport: get_vf, get_Δ, orderby_op

using Test
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

# With advection, the norm should be lower because the pollution is more spread out.
@test norm(sim.u) ≈ 330.19092435836507 

@testset "get_vf lon" begin
    pvaridx = findfirst(
        isequal("lon"), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))
    _, idx_f = orderby_op(
        EarthSciMLBase.utype(sim.domaininfo), [size(sim.u)...], 1 + pvaridx)

    @test get_vf(sim, "lon", sim.obs_fs[sim.obs_fs_idx[op.vardict["lon"]]], idx_f)(
        2, 3, starttime) ≈-0.8448177656085027
end

@testset "get_vf lat" begin
    pvaridx = findfirst(
        isequal("lat"), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))
    _, idx_f = orderby_op(
        EarthSciMLBase.utype(sim.domaininfo), [size(sim.u)...], 1 + pvaridx)

    @test get_vf(sim, "lat", sim.obs_fs[sim.obs_fs_idx[op.vardict["lat"]]], idx_f)(
        2, 3, starttime) ≈ 3.8061048027295152
end

@testset "get_vf lev" begin
    pvaridx = findfirst(
        isequal("lev"), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))
    _, idx_f = orderby_op(
        EarthSciMLBase.utype(sim.domaininfo), [size(sim.u)...], 1 + pvaridx)

    @test get_vf(sim, "lev", sim.obs_fs[sim.obs_fs_idx[op.vardict["lev"]]], idx_f)(
        2, 3, starttime) ≈ -0.006329571396093452
end

@testset "get_Δ" begin
    @test get_Δ(sim, "lat")(2, 3, starttime) ≈ 445280.0
    @test get_Δ(sim, "lon")(2, 3, starttime) ≈ 424080.6852300487
    @test get_Δ(sim, "lev")(2, 3, starttime) ≈ -1526.0725231324905
end