export AdvectionOperator

#=
An advection kernel for a 4D array, where the first dimension is the state variables
and the next three dimensions are the spatial dimensions.
=#
function advection_kernel_4d(u, stencil, vs, Δs, Δt, idx, p = NullParameters())
    lpad, rpad = stencil_size(stencil)
    offsets = ((CartesianIndex(0, lpad, 0, 0), CartesianIndex(0, rpad, 0, 0)),
        (CartesianIndex(0, 0, lpad, 0), CartesianIndex(0, 0, rpad, 0)),
        (CartesianIndex(0, 0, 0, lpad), CartesianIndex(0, 0, 0, rpad))
    )
    du = zero(eltype(u))
    @inbounds for i in eachindex(vs, Δs, offsets)
        v, Δ, (l, r) = vs[i], Δs[i], offsets[i]
        uu = @view u[(idx - l):(idx + r)]
        du += stencil(uu, v, Δt, Δ; p)
    end
    du
end

function get_vs(v_fs, i, j, k, t)
    (
        (v_fs[1](i, j, k, t), v_fs[1](i + 1, j, k, t)),
        (v_fs[2](i, j, k, t), v_fs[2](i, j + 1, k, t)),
        (v_fs[3](i, j, k, t), v_fs[3](i, j, k + 1, t))
    )
end
get_vs(v_fs, idx::CartesianIndex{4}, t) = get_vs(v_fs, idx[2], idx[3], idx[4], t)

get_Δs(Δ_fs, i, j, k, t) = (Δ_fs[1](i, j, k, t), Δ_fs[2](i, j, k, t), Δ_fs[3](i, j, k, t))
get_Δs(Δ_fs, idx::CartesianIndex{4}, t) = get_Δs(Δ_fs, idx[2], idx[3], idx[4], t)

#=
A function to create an advection operator for a 4D array,

Arguments:
    * `u_prototype`: A prototype array of the same size and type as the input array.
    * `stencil`: The stencil operator, e.g. `l94_stencil` or `ppm_stencil`.
    * `v_fs`: A vector of functions to get the wind velocity at a given place and time.
            The function signature should be `v_fs(i, j, k, t)`.
    * `Δ_fs`: A vector of functions to get the grid spacing at a given place and time.
            The function signature should be `Δ_fs(i, j, k, t)`.
    * `Δt`: The time step size, which is assumed to be fixed.
    * `bc_arraytype`: The boundary condition array type, e.g. `ZeroGradBCArray`.
=#
function advection_op(u_prototype, stencil, v_fs, Δ_fs, Δt, bc_arraytype;
        p = NullParameters())
    sz = size(u_prototype)
    v_fs = tuple(v_fs...)
    Δ_fs = tuple(Δ_fs...)
    function advection(u, p, t) # Out-of-place
        u = bc_arraytype(reshape(u, sz...))
        du = [advection_kernel_4d(u, stencil, get_vs(v_fs, idx, t),
                  get_Δs(Δ_fs, idx, t), Δt, idx, p)
              for idx in CartesianIndices(u)]
        reshape(du, :)
    end
    function advection(du, u, p, t) # In-place
        u = bc_arraytype(reshape(u, sz...))
        du = reshape(du, sz...)
        for i in 1:sz[2], j in 1:sz[3], k in 1:sz[4]
            vs = get_vs(v_fs, i, j, k, t)
            Δs = get_Δs(Δ_fs, i, j, k, t)
            for s in 1:sz[1]
                idx = CartesianIndex(s, i, j, k)
                du[s, i, j, k] = advection_kernel_4d(u, stencil, vs, Δs, Δt, idx, p)
            end
        end
    end
    FunctionOperator(advection, reshape(u_prototype, :), p = p)
end

"Get a value from the x-direction velocity field."
function vf_x(args1, args2)
    i, j, k, t = args1
    data_f, grid1, grid2, grid3, Δ = args2
    x1 = grid1[min(i, length(grid1))] - Δ / 2 # Staggered grid
    x2 = grid2[j]
    x3 = grid3[k]
    data_f(t, x1, x2, x3)
end

"Get a value from the y-direction velocity field."
function vf_y(args1, args2)
    i, j, k, t = args1
    data_f, grid1, grid2, grid3, Δ = args2
    x1 = grid1[i]
    x2 = grid2[min(j, length(grid2))] - Δ / 2 # Staggered grid
    x3 = grid3[k]
    data_f(t, x1, x2, x3)
end

"Get a value from the z-direction velocity field."
function vf_z(args1, args2)
    i, j, k, t = args1
    data_f, grid1, grid2, grid3, Δ = args2
    x1 = grid1[i]
    x2 = grid2[j]
    x3 = k > 1 ? grid3[min(k, length(grid3))] - Δ / 2 : grid3[k]
    data_f(t, x1, x2, x3) # Staggered grid
end
tuplefunc(vf) = (i, j, k, t) -> vf((i, j, k, t))

"""
$(SIGNATURES)

Return a function that gets the wind velocity at a given place and time for the given `varname`.
`vardict` should be a dictionary with keys that are strings with the
possible `varname`s and values that are the corresponding variables in the
ODESystem that should be used to get the wind velocity values.
`idx_f` should be an index function of the type returned by [`EnvironmentalTransport.orderby_op`](@ref).
`data_f` should be a function that takes a time and three spatial coordinates and returns the value of
the wind speed in the direction indicated by `varname`.
"""
function get_vf(sim, varname::AbstractString, data_f)
    if varname ∈ ("lon", "x")
        vf = Base.Fix2(
            vf_x, (data_f, sim.grid[1], sim.grid[2], sim.grid[3], sim.Δs[1]))
        return tuplefunc(vf)
    elseif varname ∈ ("lat", "y")
        vf = Base.Fix2(
            vf_y, (data_f, sim.grid[1], sim.grid[2], sim.grid[3], sim.Δs[2]))
        return tuplefunc(vf)
    elseif varname == "lev"
        vf = Base.Fix2(
            vf_z, (data_f, sim.grid[1], sim.grid[2], sim.grid[3], sim.Δs[3]))
        return tuplefunc(vf)
    else
        error("Invalid variable name $(varname).")
    end
end

"function to get grid deltas."
function Δf(args1, args2)
    i, j, k, t = args1
    tff, Δ, grid1, grid2, grid3 = args2
    c1, c2, c3 = grid1[i], grid2[j], grid3[k]
    Δ / tff(t, c1, c2, c3)
end

"""
$(SIGNATURES)

Return a function that gets the grid spacing at a given place and time for the given `varname`.
"""
function get_Δ(sim::EarthSciMLBase.Simulator, varname::AbstractString)
    pvaridx = findfirst(
        isequal(varname), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))
    tff = sim.tf_fs[pvaridx]

    tuplefunc(Base.Fix2(
        Δf, (tff, sim.Δs[pvaridx], sim.grid[1], sim.grid[2], sim.grid[3])))
end

"""
$(SIGNATURES)

Create an `EarthSciMLBase.Operator` that performs advection.
Advection is performed using the given `stencil` operator
(e.g. `l94_stencil` or `ppm_stencil`).
`p` is an optional parameter set to be used by the stencil operator.
`bc_arraytype` is the boundary condition array type, e.g. `ZeroGradBCArray`.
"""
mutable struct AdvectionOperator <: EarthSciMLBase.Operator
    Δt::Any
    stencil::Any
    bc_arraytype::Any
    vardict::Any

    function AdvectionOperator(Δt, stencil, bc_arraytype)
        new(Δt, stencil, bc_arraytype, nothing)
    end
end

function EarthSciMLBase.get_scimlop(op::AdvectionOperator, sim::Simulator, u = nothing)
    u = isnothing(u) ? init_u(sim) : u
    pvars = EarthSciMLBase.pvars(sim.domaininfo)
    pvarstrs = [String(Symbol(pv)) for pv in pvars]

    v_fs = []
    Δ_fs = []
    for varname in pvarstrs
        data_f = sim.obs_fs[sim.obs_fs_idx[op.vardict[varname]]]
        push!(v_fs, get_vf(sim, varname, data_f))
        push!(Δ_fs, get_Δ(sim, varname))
    end
    scimlop = advection_op(u, op.stencil, v_fs, Δ_fs, op.Δt, op.bc_arraytype, p = sim.p)
    cache_operator(scimlop, u[:])
end

"""
$(SIGNATURES)

Couple the advection operator into the CoupledSystem.
This function mutates the operator to add the windfield variables.
There must already be a source of wind data in the coupled system for this to work.
Currently the only valid source of wind data is `EarthSciData.GEOSFP`.
"""
function EarthSciMLBase.couple(c::CoupledSystem, op::AdvectionOperator)::CoupledSystem
    found = 0
    for sys in c.systems
        if EarthSciMLBase.get_coupletype(sys) == EarthSciData.GEOSFPCoupler
            found += 1
            op.vardict = Dict(
                "lon" => sys.A3dyn₊U,
                "lat" => sys.A3dyn₊V,
                "lev" => sys.A3dyn₊OMEGA
            )
        end
    end
    if found == 0
        error("Could not find a source of wind data in the coupled system. Valid sources are currently {EarthSciData.GEOSFP}.")
    elseif found > 1
        error("Found multiple sources of wind data in the coupled system. Valid sources are currently {EarthSciData.GEOSFP}")
    end
    push!(c.ops, op)
    c
end
