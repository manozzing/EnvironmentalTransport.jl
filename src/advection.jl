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
             (v_f(i - lpad, t), v_f(i - lpad + 1, t)), Δt, Δx_f(i, t))
         for i in (firstindex(u) + lpad):(lastindex(u) - rpad)]
    end
    function f(du::AbstractVector, u::AbstractVector, p, t) # In-place, vector
        for i in (firstindex(u) + lpad):(lastindex(u) - rpad)
            du[i - lpad] = stencil(
                u[(i - lpad):(i + rpad)], (v_f(i - lpad, t), v_f(i - lpad + 1, t)), Δt, Δx_f(
                    i, t))
        end
        du
    end
    function f(u::AbstractMatrix, p, t) # Out-of-place, matrix
        hcat([[stencil(col[(i - lpad):(i + rpad)],
                   (v_f(i - lpad, j, t), v_f(i - lpad + 1, j, t)), Δt, Δx_f(i, j, t))
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
                        Δt, Δx_f(i, j, t)
                    )
                end
            end
        end
        du
    end
    indata = zeros(dtype, shape[1] + lpad + rpad, shape[2:end]...)
    outdata = f(indata, nothing, nothing)
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
"""
function tensor_advection_op(
        dtype, shape, index, stencil, v_f, Δx_f, Δt; bc_opf = zerograd_bc_op)
    shape = [shape...]
    order = orderby_op(dtype, shape, index)
    xbc = bc_opf(shape[index], stencil)
    ncols = *(shape...) ÷ shape[index]
    adv_op = advect_1d_op(dtype, (shape[index], ncols), stencil, v_f, Δx_f, Δt)

    # How this operator works:
    # Starting from the right side and moving toward the left, first we 
    # reorder the tensor into a matrix where each column is an element
    # in the requested index of the tensor.
    # Next we apply the boundary conditions (xbc) and advection (adv_op)
    # to each column in the matrix, reorder the tensor back into the 
    # original configuration.
    inv(order) * TensorProductOperator((adv_op * xbc), I(ncols)) * order
end