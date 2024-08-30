export AdvectionOperator

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
    * `p`: Optional parameters to pass to the stencil function.
"""
function advect_1d_op(dtype, shape, stencil, v_f, Δx_f, Δt; p = NullParameters())
    lpad, rpad = stencil_size(stencil)
    function f(u::AbstractVector, p, t) # Out-of-place, vector
        [stencil(u[(i - lpad):(i + rpad)],
             (v_f(i - lpad, t), v_f(i - lpad + 1, t)), Δt, Δx_f(i - lpad, t), p = p)
         for i in (firstindex(u) + lpad):(lastindex(u) - rpad)]
    end
    function f(du::AbstractVector, u::AbstractVector, p, t) # In-place, vector
        for i in (firstindex(u) + lpad):(lastindex(u) - rpad)
            du[i - lpad] = stencil(
                u[(i - lpad):(i + rpad)], (v_f(i - lpad, t), v_f(i - lpad + 1, t)), Δt, Δx_f(
                    i - lpad, t), p = p)
        end
        du
    end
    function f(u::AbstractMatrix, p, t) # Out-of-place, matrix
        hcat([[stencil(col[(i - lpad):(i + rpad)],
                   (v_f(i - lpad, j, t), v_f(i - lpad + 1, j, t)), Δt, Δx_f(i - lpad, j, t), p = p)
               for i in (firstindex(col) + lpad):(lastindex(col) - rpad)]
              for (j, col) in enumerate(eachcol(u))]...)
    end
    function f(du::AbstractMatrix, u::AbstractMatrix, p, t) # In-place, matrix
        begin
            Threads.@threads for j in 1:size(u, 2)
                ddu = view(du, :, j)
                uu = view(u, :, j)
                for i in (firstindex(uu) + lpad):(lastindex(uu) - rpad)
                    uuslice = view(uu, (i - lpad):(i + rpad))
                    ddu[i - lpad] = stencil(uuslice,
                        (v_f(i - lpad, j, t), v_f(i - lpad + 1, j, t)),
                        Δt, Δx_f(i - lpad, j, t), p = p
                    )
                end
            end
        end
        du
    end
    indata = zeros(dtype, shape[1] + lpad + rpad, shape[2:end]...)
    outdata = zeros(dtype, shape[1], shape[2:end]...)
    FunctionOperator(f, indata, outdata, batch = true, p = p)
end

"""
$(SIGNATURES)

Return a SciMLOperator that performs 1D advection on the dimension of a
tensor given by `index`, where the original tensor has the given data type `dtype`
(e.g. `Float32` or `Float64`), and given shape (e.g. `(10, 20, 30)`).
Advection is performed using the given `stencil` operator
(e.g. `l94_stencil` or `ppm_stencil`).

Optionally, a function can be
specified to create boundary conditions, where the function should have the signature
`bc_opf(vector_length, stencil)`. See the default boundary condition operator
[`EnvironmentalTransport.zerograd_bc_op`](@ref) for more information.

This function returns a function that creates the the operator and also a
function `idx_f` that takes a
column number of the transformed matrix as an input and returns a vector of
`CartesianIndex`es in the original tensor that make up that transformed matrix
column.

The first returned function (the function that creates the operator) can be called with the
following arguments to create the operator:

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
    * `p`: Optional parameters to pass to the stencil function.

The reason for the two-step process is that `Δv_f` and `Δx_f`` may be dependent on
`idx_f`, so we need to return `idx_f` before we create the operator.
"""
function tensor_advection_op(dtype, shape, index, stencil; bc_opf = zerograd_bc_op,
        p = NullParameters(), kwargs...)
    shape = [shape...]
    order, idx_f = orderby_op(dtype, shape, index, p = p)
    bc = bc_opf(zeros(dtype, shape[index]), stencil, p = p)
    ncols = *(shape...) ÷ shape[index]
    return (v_f, Δx_f, Δt) -> begin
        adv_op = advect_1d_op(dtype, (shape[index], ncols), stencil, v_f, Δx_f, Δt; p = p)

        # How this operator works:
        # Starting from the right side and moving toward the left, first we
        # reorder the tensor into a matrix where each column is an element
        # in the requested index of the tensor ('order').
        # Next we apply the boundary conditions ('bc') and advection ('adv_op')
        # to each column in the matrix, reorder the tensor back into the
        # original configuration ('inv(order)').
        A = I(ncols)
        B = adv_op * bc
        inv(order) * TensorProductOperator(A, B) * order
    end,
    idx_f
end

"Get a value from the x-direction velocity field."
function vf_x(args1, args2)
    i, j, t = args1
    data_f, d, grid1, grid2, grid3, idx_f, Δ = args2
    idx = idx_f(min(i, d), j) # Avoid out-of-bounds because CartesianIndex isn't on staggered grid.
    x1 = grid1[idx[2]] - Δ / 2 # Staggered grid
    x2 = grid2[idx[3]]
    x3 = grid3[idx[4]]
    data_f(t, x1, x2, x3)
end

"Get a value from the y-direction velocity field."
function vf_y(args1, args2)
    i, j, t = args1
    data_f, d, grid1, grid2, grid3, idx_f, Δ = args2
    idx = idx_f(min(i, d), j) # Avoid out-of-bounds because CartesianIndex isn't on staggered grid.
    x1 = grid1[idx[2]]
    x2 = grid2[idx[3]] - Δ / 2 # Staggered grid
    x3 = grid3[idx[4]]
    data_f(t, x1, x2, x3)
end

"Get a value from the z-direction velocity field."
function vf_z(args1, args2)
    i, j, t = args1
    data_f, d, grid1, grid2, grid3, idx_f, Δ = args2
    idx = idx_f(min(i, d), j) # Avoid out-of-bounds because CartesianIndex isn't on staggered grid.
    x1 = grid1[idx[2]]
    x2 = grid2[idx[3]]
    x3 = idx[4] > 1 ? grid3[idx[4]] - Δ / 2 : grid3[idx[4]]
    data_f(t, x1, x2, x3) # Staggered grid
end
tuplefunc(vf) = (i, j, t) -> vf((i, j, t))

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
function get_vf(sim, varname::AbstractString, data_f, idx_f)
    pvaridx = findfirst(
        isequal(varname), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))
    d = size(sim)[pvaridx + 1]

    if varname ∈ ("lon", "x")
        vf = Base.Fix2(
            vf_x, (data_f, d, sim.grid[1], sim.grid[2], sim.grid[3], idx_f, sim.Δs[1]))
        return tuplefunc(vf)
    elseif varname ∈ ("lat", "y")
        vf = Base.Fix2(
            vf_y, (data_f, d, sim.grid[1], sim.grid[2], sim.grid[3], idx_f, sim.Δs[2]))
        return tuplefunc(vf)
    elseif varname == "lev"
        vf = Base.Fix2(
            vf_z, (data_f, d, sim.grid[1], sim.grid[2], sim.grid[3], idx_f, sim.Δs[3]))
        return tuplefunc(vf)
    else
        error("Invalid variable name $(varname).")
    end
end

"function to get grid deltas."
function Δf(args1, args2)
    i, j, t = args1
    idx_f, tff, Δ, grid1, grid2, grid3 = args2
    idx = idx_f(i, j)
    c1 = grid1[idx[2]]
    c2 = grid2[idx[3]]
    c3 = grid3[idx[4]]
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

    _, idx_f = orderby_op(
        EarthSciMLBase.utype(sim.domaininfo), [size(sim)...], 1 + pvaridx)

    tuplefunc(Base.Fix2(
        Δf, (idx_f, tff, sim.Δs[pvaridx], sim.grid[1], sim.grid[2], sim.grid[3])))
end

"""
$(SIGNATURES)

Create an `EarthSciMLBase.Operator` that performs advection.
Advection is performed using the given `stencil` operator
(e.g. `l94_stencil` or `ppm_stencil`).
`p` is an optional parameter set to be used by the stencil operator.
"""
mutable struct AdvectionOperator <: EarthSciMLBase.Operator
    Δt::Any
    stencil::Any
    vardict::Any

    function AdvectionOperator(Δt, stencil)
        new(Δt, stencil, nothing)
    end
end

"""
$(SIGNATURES)

Create a 1D advection SciMLOperator for the given variable name `varname`.
`vardict` should be a dictionary with keys that are strings with the
possible `varname`s and values that are the corresponding variables in the
ODESystem that should be used to get the wind velocity values.
`p` is an optional parameter set that can be passed to the stencil function.
"""
function simulator_advection_1d(sim::EarthSciMLBase.Simulator, op::AdvectionOperator,
                                varname; p = NullParameters())
    pvaridx = findfirst( # Get the index of the variable in the domaininfo
        isequal(varname), String.(Symbol.(EarthSciMLBase.pvars(sim.domaininfo))))

    op_f, idx_f = tensor_advection_op(
        EarthSciMLBase.utype(sim.domaininfo),
        size(sim),
        1 + pvaridx,
        op.stencil;
        p = p
    )

    data_f = sim.obs_fs[sim.obs_fs_idx[op.vardict[varname]]]
    op_f(
        get_vf(sim, varname, data_f, idx_f),
        get_Δ(sim, varname),
        100.0
    )
end

function EarthSciMLBase.get_scimlop(op::AdvectionOperator, s::Simulator)
    pvars = EarthSciMLBase.pvars(s.domaininfo)
    pvarstrs = [String(Symbol(pv)) for pv in pvars]
    # Create advection operators in each of the three dimensions and add them together.
    # TODO: Turn vertical advection back on.
    #op = +([simulator_advection_1d(s, op, pv, p = s.p) for pv in pvarstrs]...)
    op = +([simulator_advection_1d(s, op, pv, p = s.p) for pv in pvarstrs[1:2]]...)
    u = zeros(EarthSciMLBase.utype(s.domaininfo), size(s)...)
    cache_operator(op, u[:])
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
