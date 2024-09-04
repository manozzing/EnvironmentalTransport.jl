using EnvironmentalTransport
using EnvironmentalTransport: get_vf, get_Δ

using Test
using EarthSciMLBase, EarthSciData
using ModelingToolkit, DomainSets, OrdinaryDiffEq
using ModelingToolkit: t, D
using Distributions, LinearAlgebra
using DynamicQuantities
using Dates

@parameters(
    lon=0.0, [unit=u"rad"],
    lat=0.0, [unit=u"rad"],
    lev=1.0,
)
starttime = datetime2unix(DateTime(2022, 5, 1))
endtime = datetime2unix(DateTime(2022, 5, 1, 1, 0, 5))

geosfp, geosfp_updater = GEOSFP("4x5"; dtype = Float64,
    coord_defaults = Dict(:lon => 0.0, :lat => 0.0, :lev => 1.0))

domain = DomainInfo(
    [partialderivatives_δxyδlonlat,
        partialderivatives_δPδlev_geosfp(geosfp)],
    constIC(16.0, t ∈ Interval(starttime, endtime)),
    constBC(16.0,
        lon ∈ Interval(deg2rad(-130.0), deg2rad(-60.0)),
        lat ∈ Interval(deg2rad(9.75), deg2rad(60.0)),
        lev ∈ Interval(1, 3)))

function emissions(μ_lon, μ_lat, σ)
    @variables c(t) = 0.0 [unit=u"kg"]
    @constants v_emis = 50.0 [unit=u"kg/s"]
    @constants t_unit = 1.0 [unit=u"s"] # Needed so that arguments to `pdf` are unitless.
    dist = MvNormal([starttime, μ_lon, μ_lat, 1], Diagonal(map(abs2, [3600.0, σ, σ, 1])))
    ODESystem([D(c) ~ pdf(dist, [t/t_unit, lon, lat, lev]) * v_emis],
        t, name = :Test₊emissions)
end

emis = emissions(deg2rad(-122.6), deg2rad(45.5), 0.1)

csys = couple(emis, domain, geosfp, geosfp_updater)

sim = Simulator(csys, [deg2rad(4), deg2rad(4), 1])
st = SimulatorStrangThreads(Tsit5(), SSPRK22(), 1.0)

sol = run!(sim, st)

@test 310 < norm(sol.u[end]) < 330

op = AdvectionOperator(100.0, l94_stencil, ZeroGradBC())

@test isnothing(op.vardict) # Before coupling, there shouldn't be anything here.

csys = couple(csys, op)

@test !isnothing(op.vardict) # after coupling, there should be something here.

sol = run!(sim, st)

# With advection, the norm should be lower because the pollution is more spread out.
@test 310 < norm(sol.u[end]) < 350

@testset "get_vf lon" begin
    f = sim.obs_fs[sim.obs_fs_idx[op.vardict["lon"]]]
    @test get_vf(sim, "lon", f)(2, 3, 1, starttime) ≈ -6.816295428727573
end

@testset "get_vf lat" begin
f = sim.obs_fs[sim.obs_fs_idx[op.vardict["lat"]]]
    @test get_vf(sim, "lat", f)(3, 2, 1, starttime) ≈ -5.443038969820774
end

@testset "get_vf lev" begin
    f = sim.obs_fs[sim.obs_fs_idx[op.vardict["lev"]]]
    @test get_vf(sim, "lev", f)(3, 1, 2, starttime) ≈ -0.019995461793337128
end

@testset "get_Δ" begin
    @test get_Δ(sim, "lat")(2, 3, 1, starttime) ≈ 445280.0
    @test get_Δ(sim, "lon")(3, 2, 1, starttime) ≈ 432517.0383085161
    @test get_Δ(sim, "lev")(3, 1, 2, starttime) ≈ -1516.7789198950632
end
