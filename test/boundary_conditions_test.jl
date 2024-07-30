using EnvironmentalTransport: zerograd_bc_op
using EnvironmentalTransport
using Test

u = [1.0, 2, 3]
u2 = zeros(10)

zbc_op = zerograd_bc_op(3, ppm_stencil)
@test zbc_op(u, nothing, 1.0) ≈ [1.0, 1, 1, 1, 2, 3, 3, 3, 3, 3]
@test begin
    zbc_op(u2, u, nothing, 1.0)
    u2 ≈ [1.0, 1, 1, 1, 2, 3, 3, 3, 3, 3]
end
