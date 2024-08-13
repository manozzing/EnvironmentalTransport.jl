"""
$(SIGNATURES)

Create a SciMLOperator to reorder a tensor into a matrix, 
where each column in the matrix is an element in the given 
`index` of the original tensor. `dtype` is the type of 
the input and output data, e.g. `Float64` or `Float32`,
and `shape` is the shape of the tensor, e.g. [2,3,5].

The inverse of the operator can be used to put the tensor
back into the orinal order. So if `op` is the operator and
 `u` is the original tensor, then: `x = op * u[:]` returns
the reordered tensor `x`, and `u_prime = inv(op) * x` returns
a vector version of the orginal tensor, so `uprime == u[:]`.

This function returns the operator and also a function `idx_f` that takes a 
column number of the transformed matrix as an input and returns a vector of 
`CartesianIndex`es in the original tensfor that make up that transformed matrix
column.
"""
function orderby_op(dtype, shape::AbstractVector, index::Int; p=NullParameters())
    vec_length = *(shape...)
    ii = 1:length(shape)
    iiremainder = setdiff(ii, index)
    fwd_idx = tuple(index, iiremainder...)
    function fwd(u, p, t)
        reshape(permutedims(reshape(u, shape...), fwd_idx), shape[index], :)
    end
    fwd(du, u, p, t) = du[:] .= fwd(u, p, t)[:]
    newshape = tuple(shape[[index, iiremainder...]]...)
    rev_idx = tuple(insert!(collect(2:length(shape)), index, 1)...)
    function rev(u, p, t)
        permutedims(reshape(u, newshape...), rev_idx)[:]
    end
    rev(du, u, p, t) = du[:] .= rev(u, p, t)[:]
    idx_all = reshape(1:vec_length, shape...)
    c_ii = CartesianIndices(idx_all)
    idx_reshaped = fwd(idx_all, nothing, nothing)
    idx_f(col) = view(c_ii, view(idx_reshaped, :, col)) # function to transform the indexes
    idx_f(row, col) = view(c_ii, view(idx_reshaped, :, col))[row] # function to transform the indexes
    x = zeros(dtype, shape[1], Int(vec_length / shape[1]))
    FunctionOperator(fwd, x, x; op_adjoint = rev, op_inverse = rev,
        op_adjoint_inverse = fwd, islinear = true, p = p),
    idx_f
end
