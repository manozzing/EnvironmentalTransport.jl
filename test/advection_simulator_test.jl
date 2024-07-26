using EnvironmentalTransport
using EnvironmentalTransport: get_vf, get_Δ

using Test
using ModelingToolkit, DomainSets, OrdinaryDiffEq
using Distributions
using Dates
using Plots


@parameters lon=0.0 lat=0.0 lev=1.0 t
lat = GlobalScope(lat)
lon = GlobalScope(lon)
lev = GlobalScope(lev)
starttime = datetime2unix(DateTime(2022, 5, 1))
endtime = datetime2unix(DateTime(2022, 5, 1, 3))

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
    dist = MvNormal([starttime, μ_lon, μ_lat, 1], [3600.0, σ, σ, 1])
    D = Differential(t)
    ODESystem([D(c) ~ pdf(dist, [t, lon, lat, lev]) * 50],
        t, name = :Test₊emissions)
end

emis = emissions(t, deg2rad(-122.6), deg2rad(45.5), 0.1)
output = NetCDFOutputter("out.nc", 3600.0)

csys = couple(emis, domain, geosfp, output)

sim = Simulator(csys, [deg2rad(2), deg2rad(2), 1], Tsit5())

run!(sim)


pvaridx = findfirst(
    isequal("lon"), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))
_, idx_f = orderby_op(
    EarthSciMLBase.utype(sim.domaininfo), [size(sim.u)...], 1 + pvaridx)

@test get_vf(sim, "lon", gfp_vars, idx_f)(2, 3, starttime) ≈ -0.8947638543146277

pvaridx = findfirst(
    isequal("lat"), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))
_, idx_f = orderby_op(
    EarthSciMLBase.utype(sim.domaininfo), [size(sim.u)...], 1 + pvaridx)

@test get_vf(sim, "lat", gfp_vars, idx_f)(2, 3, starttime) ≈ 3.820172123212213

pvaridx = findfirst(
    isequal("lev"), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))
_, idx_f = orderby_op(
    EarthSciMLBase.utype(sim.domaininfo), [size(sim.u)...], 1 + pvaridx)

@test get_vf(sim, "lev", gfp_vars, idx_f)(2, 3, starttime) ≈ -0.006279307244207394



@test get_Δ(sim, "lat")(2, 3, starttime) ≈ 222640.0
@test get_Δ(sim, "lon")(2, 3, starttime) ≈ 216258.51915425804
@test get_Δ(sim, "lev")(2, 3, starttime) ≈ -1526.087459321192





op = AdvectionOperator(100.0, SSPRK22())

EarthSciMLBase.initialize!(op, sim)

csys = couple(csys, op)

sim = Simulator(csys, [deg2rad(2), deg2rad(2), 1], Tsit5())

run!(sim)

@profview EarthSciMLBase.run!(op, sim, starttime, op.Δt)

du = op.op(sim.u[:], nothing, starttime)
heatmap(reshape(du, size(sim.u)...)[1, :, :, 1]')
op.op = cache_operator(op.op, sim.du[:])
op.op(sim.du[:], sim.u[:], NullParameters(), starttime)

prob = ODEProblem(op.op, sim.u[:], (starttime, endtime))
solve(prob, SSPRK22(), dt = op.dt)

lat_adv = simulator_advection_1d(sim, "lat", gfp_vars)
lon_adv = simulator_advection_1d(sim, "lon", gfp_vars)
lev_adv = simulator_advection_1d(sim, "lev", gfp_vars)

Δt = 300.0
du_lat = lat_adv(sim.u[:], nothing, starttime) .* Δt
du_lon = lon_adv(sim.u[:], nothing, starttime) .* Δt
du_lev = lev_adv(sim.u[:], nothing, starttime) .* Δt
plot(
    heatmap(sim.u[1, :, :, 1]', title = "U"),
    heatmap(reshape(du_lon, size(sim.u)...)[1, :, :, 1]', title = "lon"),
    heatmap(reshape(du_lat, size(sim.u)...)[1, :, :, 1]', title = "lat"),
    heatmap(reshape(du_lev, size(sim.u)...)[1, :, :, 2]', title = "lev")
)

heatmap(sim.u[1, :, :, 1])