"""
$(SIGNATURES)

Create a SciMLOperator to perform 1D advection.

Arguments:
    * `dtype`: The data type of the input and output arrays, e.g. `Float64` or `Float32`.
    * `shape`: The shape of the input vector or matrix, e.g. (8,) or (8,10). 
            If the input is a matrix, 1D advection will be applied to each of the
            columns in the matrix.
    * `stencil`: The stencil operator, e.g. `l94_stencil` or `ppm_stencil`.
    *  `v_f`: A function to get the wind velocity at a given place and time. For vector inputs
            the function signature should be `v_f(i, t)`, where `i` is the staggered-grid 
            index and `t` is time. For matrix, inputs, the signature should be `v_f(i,j,t)`,
            where `i` is a staggered grid index, `j` is the non-staggered grid column index
            (where `max(i) == shape(2)`)`, and `t` is time.
    * `Δx_f`: A function to get the grid spacing at a given place and time. For vector inputs
            the function signature should be `Δx_f(i, t)`, where `i` is the grid 
            index and `t` is time. For matrix, inputs, the signature should be `v_f(i,j,t)`,
            where `i` is a grid index, `j` is the column index
            (where `max(i) == shape(2)`)`, and `t` is time.
    * `Δt`: The time step size, which is assumed to be fixed.
"""
function advect_1d_op(dtype, shape, stencil, v_f, Δx_f, Δt)
    lpad, rpad = stencil_size(stencil)
    function f(u::AbstractVector, p, t) # Out-of-place, vector
        [stencil(u[(i - lpad):(i + rpad)],
             (v_f(i - lpad, t), v_f(i - lpad + 1, t)), Δt, Δx_f(i - lpad, t))
         for i in (firstindex(u) + lpad):(lastindex(u) - rpad)]
    end
    function f(du::AbstractVector, u::AbstractVector, p, t) # In-place, vector
        for i in (firstindex(u) + lpad):(lastindex(u) - rpad)
            du[i - lpad] = stencil(
                u[(i - lpad):(i + rpad)], (v_f(i - lpad, t), v_f(i - lpad + 1, t)), Δt, Δx_f(
                    i - lpad, t))
        end
        du
    end
    function f(u::AbstractMatrix, p, t) # Out-of-place, matrix
        hcat([[stencil(col[(i - lpad):(i + rpad)],
                   (v_f(i - lpad, j, t), v_f(i - lpad + 1, j, t)), Δt, Δx_f(i - lpad, j, t))
               for i in (firstindex(col) + lpad):(lastindex(col) - rpad)]
              for (j, col) in enumerate(eachcol(u))]...)
    end
    function f(du::AbstractMatrix, u::AbstractMatrix, p, t) # In-place, matrix
        @views begin
            for j in 1:size(u, 2)
                ddu = du[:, j]
                uu = u[:, j]
                for i in (firstindex(uu) + lpad):(lastindex(uu) - rpad)
                    ddu[i - lpad] = stencil(
                        uu[(i - lpad):(i + rpad)],
                        (v_f(i - lpad, j, t), v_f(i - lpad + 1, j, t)),
                        Δt, Δx_f(i - lpad, j, t)
                    )
                end
            end
        end
        du
    end
    indata = zeros(dtype, shape[1] + lpad + rpad, shape[2:end]...)
    outdata = zeros(dtype, shape[1], shape[2:end]...)
    FunctionOperator(f, indata, outdata, batch = true)
end

"""
$(SIGNATURES)

Return a SciMLOperator that performs 1D advection on the dimension of a 
tensor given by `index`, where the original tensor has the given data type `dtype` 
(e.g. `Float32` or `Float64`), and given shape (e.g. `(10, 20, 30)`).
Advection is performed using the given `stencil` operator 
(e.g. `l94_stencil` or `ppm_stencil`). 

Additional arguments:
    *  `v_f`: A function to get the wind velocity at a given place and time. For vector inputs
            the function signature should be `v_f(i, t)`, where `i` is the staggered-grid 
            index and `t` is time. For matrix, inputs, the signature should be `v_f(i,j,t)`,
            where `i` is a staggered grid index, `j` is the non-staggered grid column index
            (where `max(i) == shape(2)`)`, and `t` is time.
    * `Δx_f`: A function to get the grid spacing at a given place and time. For vector inputs
            the function signature should be `Δx_f(i, t)`, where `i` is the grid 
            index and `t` is time. For matrix, inputs, the signature should be `v_f(i,j,t)`,
            where `i` is a grid index, `j` is the column index
            (where `max(i) == shape(2)`)`, and `t` is time.
    * `Δt`: The time step size, which is assumed to be fixed.

Optionally, a function can be 
specified to create boundary conditions, where the function should have the signature
`bc_opf(vector_length, stencil)`. See the default boundary condition operator 
[`EnvironmentalTransport.zerograd_bc_op`](@ref) for more information.

This function returns the operator and also a function `idx_f` that takes a 
column number of the transformed matrix as an input and returns a vector of 
`CartesianIndex`es in the original tensor that make up that transformed matrix
column.
"""
function tensor_advection_op(
        dtype, shape, index, stencil, v_f, Δx_f, Δt; bc_opf = zerograd_bc_op)
    shape = [shape...]
    order, idx_f = orderby_op(dtype, shape, index)
    bc = bc_opf(shape[index], stencil)
    ncols = *(shape...) ÷ shape[index]
    adv_op = advect_1d_op(dtype, (shape[index], ncols), stencil, v_f, Δx_f, Δt)

    # How this operator works:
    # Starting from the right side and moving toward the left, first we 
    # reorder the tensor into a matrix where each column is an element
    # in the requested index of the tensor ('order').
    # Next we apply the boundary conditions ('bc') and advection ('adv_op')
    # to each column in the matrix, reorder the tensor back into the 
    # original configuration ('inv(order)').
    inv(order) * TensorProductOperator((adv_op * bc), I(ncols)) * order, idx_f
end

using ModelingToolkit, DomainSets, OrdinaryDiffEq
using Distributions
using Dates
using Plots

@parameters lon=0.0 lat=0.0 lev=1.0 t
lat = GlobalScope(lat)
lon = GlobalScope(lon)
lev = GlobalScope(lev)
starttime = datetime2unix(DateTime(2022, 5, 1))
endtime = datetime2unix(DateTime(2022, 5, 2))

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

gfp_vars = Dict(
    "lon" => sim.obs_fs[sim.obs_fs_idx[geosfp.A3dyn₊U]],
    "lat" => sim.obs_fs[sim.obs_fs_idx[geosfp.A3dyn₊V]],
    "lev" => sim.obs_fs[sim.obs_fs_idx[geosfp.A3dyn₊OMEGA]]
)

function get_vf(sim, varname::AbstractString, vardict)
    pvaridx = findfirst(
        isequal(varname), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))

    _, idx_f = orderby_op(
        EarthSciMLBase.utype(sim.domaininfo), [size(sim.u)...], 1 + pvaridx)
    vf = vardict[varname]
    if varname ∈ ("lon", "x")
        return (i, j, t) -> begin
            idx = idx_f(min(i, size(sim.u, pvaridx + 1)), j) # Avoid out-of-bounds because CartesianIndex isn't on staggered grid.
            vf(t,
                sim.grid[1][idx[2]] - sim.Δs[1] / 2, # Staggered grid 
                sim.grid[2][idx[3]],
                sim.grid[3][idx[4]])
        end
    elseif varname ∈ ("lat", "y")
        return (i, j, t) -> begin
            idx = idx_f(min(i, size(sim.u, pvaridx + 1)), j) # Avoid out-of-bounds because CartesianIndex isn't on staggered grid.
            vf(t,
                sim.grid[1][idx[2]],
                sim.grid[2][idx[3]] - sim.Δs[2] / 2, # Staggered grid 
                sim.grid[3][idx[4]])
        end
    elseif varname == "lev"
        return (i, j, t) -> begin
            idx = idx_f(min(i, size(sim.u, pvaridx + 1)), j) # Avoid out-of-bounds because CartesianIndex isn't on staggered grid.
            vf(t,
                sim.grid[1][idx[2]],
                sim.grid[2][idx[3]],
                idx[4] > 1 ? sim.grid[3][idx[4]] - sim.Δs[3] / 2 : sim.grid[3][idx[4]] # Staggered grid 
            )
        end
    else
        error("Invalid variable name $(varname).")
    end
end

using Test

@test get_vf(sim, "lon", gfp_vars)(2, 3, starttime) ≈ -0.8947638543146277
@test get_vf(sim, "lat", gfp_vars)(2, 3, starttime) ≈ 3.820172123212213
@test get_vf(sim, "lev", gfp_vars)(2, 3, starttime) ≈ -0.006279307244207394

function get_Δ(sim, varname::AbstractString)
    pvaridx = findfirst(
        isequal(varname), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))
    tff = sim.tf_fs[pvaridx]

    _, idx_f = orderby_op(
        EarthSciMLBase.utype(sim.domaininfo), [size(sim.u)...], 1 + pvaridx)

    return (i, j, t) -> begin
        idx = idx_f(i, j)
        sim.Δs[pvaridx] / tff(t,
            sim.grid[1][idx[2]],
            sim.grid[2][idx[3]],
            sim.grid[3][idx[4]])
    end
end
@test get_Δ(sim, "lat")(2, 3, starttime) ≈ 222640.0
@test get_Δ(sim, "lon")(2, 3, starttime) ≈ 216258.51915425804
@test get_Δ(sim, "lev")(2, 3, starttime) ≈ -1526.087459321192

function simulator_advection_1d(sim, varname, vardict)
    pvaridx = findfirst(
        isequal(varname), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))

    tensor_advection_op(
        EarthSciMLBase.utype(sim.domaininfo),
        size(sim.u),
        1 + pvaridx,
        ppm_stencil,
        (i, j, t) -> 1.0, #get_vf(sim, varname, vardict),
        get_Δ(sim, varname),
        100.0
    )
end

lat_adv, _ = simulator_advection_1d(sim, "lat", gfp_vars)
lon_adv, _ = simulator_advection_1d(sim, "lon", gfp_vars)
lev_adv, _ = simulator_advection_1d(sim, "lev", gfp_vars)

Δt = 300.0
du_lat = lat_adv(sim.u[:], nothing, starttime) .* Δt
du_lon = lon_adv(sim.u[:], nothing, starttime) .* Δt
du_lev = lev_adv(sim.u[:], nothing, starttime) .* Δt
plot(
    heatmap(sim.u[1, :, :, 1]', title="U"),
    heatmap(reshape(du_lon, size(sim.u)...)[1, :, :, 1]', title="lon"),
    heatmap(reshape(du_lat, size(sim.u)...)[1, :, :, 1]', title="lat"),
    heatmap(reshape(du_lev, size(sim.u)...)[1, :, :, 1]', title="lev"),
)