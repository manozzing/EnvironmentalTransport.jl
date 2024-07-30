"""
$(SIGNATURES)

Create an SciMLOperator which adds zero gradient Newmann boundary conditions to a 1D array of
length `vec_length`, which work with the given `stencil`.
"""
function zerograd_bc_op(vec_length, stencil)
    lpad, rpad = stencil_size(stencil)
    padmat = vcat(
        repeat(hcat(1, repeat([0], vec_length - 1)...), lpad),
        I(vec_length),
        repeat(hcat(repeat([0], vec_length - 1)..., 1), rpad)
    )
    MatrixOperator(padmat)
end
