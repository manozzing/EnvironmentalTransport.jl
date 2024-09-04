using BenchmarkTools
using EnvironmentalTransport

using EarthSciMLBase, EarthSciData
using ModelingToolkit, DomainSets
using Dates

starttime = datetime2unix(DateTime(2022, 5, 1))

function setup_advection_simulator(lonres, latres, stencil)
    @parameters lon=0.0 lat=0.0 lev=1.0 t
    endtime = datetime2unix(DateTime(2022, 5, 1, 1, 0, 5))

    geosfp, updater = GEOSFP("4x5"; dtype = Float64,
        coord_defaults = Dict(:lon => 0.0, :lat => 0.0, :lev => 1.0))

    domain = DomainInfo(
        [partialderivatives_δxyδlonlat,
            partialderivatives_δPδlev_geosfp(geosfp)],
        constIC(16.0, t ∈ Interval(starttime, endtime)),
        constBC(16.0,
            lon ∈ Interval(deg2rad(-129), deg2rad(-61)),
            lat ∈ Interval(deg2rad(11), deg2rad(59)),
            lev ∈ Interval(1, 3)))

    function emissions(t)
        @variables c(t) = 1.0
        D = Differential(t)
        ODESystem([D(c) ~ lat + lon + lev], t, name = :emissions)
    end

    emis = emissions(t)

    csys = couple(emis, domain, geosfp, updater)
    op = AdvectionOperator(100.0, stencil, ZeroGradBCArray)
    csys = couple(csys, op)
    sim = Simulator(csys, [deg2rad(lonres), deg2rad(latres), 1])
    scimlop = EarthSciMLBase.get_scimlop(only(csys.ops), sim)
    scimlop, init_u(sim)
end

suite = BenchmarkGroup()
suite["Advection Simulator"] = BenchmarkGroup(["advection", "simulator"])
suite["Advection Simulator"]["in-place"] = BenchmarkGroup()
suite["Advection Simulator"]["out-of-place"] = BenchmarkGroup()

for stencil in [l94_stencil, ppm_stencil]
    suite["Advection Simulator"]["in-place"][stencil] = BenchmarkGroup()
    suite["Advection Simulator"]["out-of-place"][stencil] = BenchmarkGroup()
    for lonres in [4, 2]
        for latres in [5, 3]
            @info "setting up $lonres x $latres with $stencil"
            op, u = setup_advection_simulator(lonres, latres, stencil)
            suite["Advection Simulator"]["in-of-place"][stencil][length(u)] = @benchmarkable $(op)(
                $(u[:]), $(u[:]), [0.0], $starttime)
            suite["Advection Simulator"]["out-of-place"][stencil][length(u)] = @benchmarkable $(op)(
                $(u[:]), [0.0], $starttime)
        end
    end
end

#op, u = setup_advection_simulator(0.1, 0.1, l94_stencil)
#@profview op(u[:], u[:], [0.0], starttime)

tune!(suite)
results = run(suite, verbose = true)

BenchmarkTools.save("output.json", median(results))
