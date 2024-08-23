"""
$(SIGNATURES)

Create an SciMLOperator which adds zero gradient Newmann boundary conditions to a 1D array of
length `vec_length`, which work with the given `stencil`.
"""
function zerograd_bc_op(u_prototype, stencil; p = NullParameters())
    lpad, rpad = stencil_size(stencil)
    function f(u::AbstractVector, p, t)
        [repeat([u[begin]], lpad); u; repeat([u[end]], rpad)]
    end
    function f(v, u::AbstractVector, p, t)
        @views begin
            v[begin:(begin + lpad - 1)] .= u[begin]
            v[(begin + lpad):(end - rpad)] .= u
            v[(end - rpad + 1):end] .= u[end]
        end
    end
    function f(u::AbstractMatrix, p, t)
        hcat([f(view(u, :, i), p, t) for i in 1:size(u, 2)]...)
    end
    function f(v, u::AbstractMatrix, p, t)
        for j in 1:size(u, 2)
            @views begin
                v[begin:(begin + lpad - 1), j] .= u[begin, j]
                v[(begin + lpad):(end - rpad), j] .= u[:, j]
                v[(end - rpad + 1):end, j] .= u[end, j]
            end
        end
    end
    v_prototype = f(u_prototype, nothing, nothing)
    FunctionOperator(f, u_prototype, v_prototype; p = p, islinear = true, batch = true)
end
