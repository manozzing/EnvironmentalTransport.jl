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
"""
function orderby_op(dtype, shape::AbstractVector, index::Int)
    vec_length = *(shape...)
    ii = 1:length(shape)
    iiremainder = setdiff(ii, index)
    fwd_idx = tuple(index, iiremainder...)
    function fwd(u, p, t)
        reshape(permutedims(reshape(u, shape...), fwd_idx), shape[index], :)
    end
    newshape = tuple(shape[[index, iiremainder...]]...)
    rev_idx = tuple(insert!(collect(2:length(shape)), index, 1)...)
    function rev(u, p, t)
        permutedims(reshape(u, newshape...), rev_idx)[:]
    end
    x = zeros(dtype, shape[1], Int(vec_length / shape[1]))
    y = zeros(dtype, shape[index], Int(vec_length / shape[index]))
    FunctionOperator(fwd, x, x; op_adjoint = rev, op_inverse = rev,
        op_adjoint_inverse = fwd, islinear = true)
end
