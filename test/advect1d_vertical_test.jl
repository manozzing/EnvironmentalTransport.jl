using EnvironmentalTransport

using Test

c = [0.0, 1, 2, 3, 4, 5]
v = [10.0, 8, 6, 4, 2, 0, 1]
nz = 6
Δt = 50
p = [1000.0, 900, 850, 700, 500, 300, 250]

result = advect1d_vertical_L94(c, v, nz, Δt, p)

@test result ≈ [0.0, 0.8958166666666667, 1.9931652777777777, 3.00495, 4.014984375, 5.04495]

result = advect1d_vertical_ppm(c, v, nz, Δt, p)

@test result ≈ [0.0, 0.9766666666666667, 1.9983333333333333, 3.0025, 3.9825, 5.04]
