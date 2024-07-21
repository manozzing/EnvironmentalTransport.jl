using EnvironmentalTransport

using Test

c = [0.0, 1, 2, 3, 4, 5]
v = [10.0, 8, 6, 4, 2, 0, 1]
nz = 6
Δt = 0.05
Δz = 0.5

result = advect1d_l94(c, v, nz, Δt, Δz)

@test result ≈ [0.0, 0.28, 1.8, 3.24, 4.68, 4.5]

result = advect1d_ppm(c, v, nz, Δt, Δz)

@test result ≈ [0.0, 0.3999999999999999, 1.8, 3.2, 4.6, 4.5]


c = [6.0, 6, 5, 5, 6, 6]
v = [2.0, 2, 2, 2, 2, 2, 2]

result = advect1d_l94(c, v, nz, Δt, Δz)

@test result ≈ [6.0, 6.0, 5.2, 5.0, 5.8, 6.0]

result = advect1d_ppm(c, v, nz, Δt, Δz)

@test result ≈ [6.0, 6.0, 5.2, 5.0, 5.8, 6.0]

u = repeat([10.0, 8, 6, 4, 5, 7]', 6, 1)
v = repeat([0.0, 2, 4, 6, 4, 2], 1, 6)
c = reshape(1.0:36, 6, 6)
lons = -90:30:60
lats = -70:30:80

@test_broken advect2D(Δt, 6, 6, c, u, v, lons, lats, advect1d_l94)