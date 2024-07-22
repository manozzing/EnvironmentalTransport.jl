using EnvironmentalTransport

using Test

c = [0.0, 1, 2, 3, 4, 5]
v = [10.0, 8, 6, 4, 2, 0, 1]
nz = 6
Δt = 0.05
Δz = 0.5

@testset "l94 1" begin
    c2 = [c[1], c[1], c..., c[end], c[end]]
    result = [l94_stencil(c2[(i - 2):(i + 2)], v[(i - 2):(i - 1)], Δt, Δz) for i in 3:8]
    @test result ≈ [0.0, 0.28, 1.8, 3.24, 4.68, 4.5]
end

@testset "ppm 1" begin
    c2 = [c[1], c[1], c[1], c..., c[end], c[end], c[end], c[end]]
    result = [ppm_stencil(c2[(i - 3):(i + 4)], v[(i - 3):(i - 2)], Δt, Δz) for i in 4:9]
    @test result ≈ [0.0, 0.3999999999999999, 1.8, 3.2, 4.6, 4.5]
end

c = [6.0, 6, 5, 5, 6, 6]
v = [2.0, 2, 2, 2, 2, 2, 2]

@testset "l94 2" begin
    c2 = [c[1], c[1], c..., c[end], c[end]]
    result = [l94_stencil(c2[(i - 2):(i + 2)], v[(i - 2):(i - 1)], Δt, Δz) for i in 3:8]
    @test result ≈ [6.0, 6.0, 5.2, 5.0, 5.8, 6.0]
end

@testset "ppm 2" begin
    c2 = [c[1], c[1], c[1], c..., c[end], c[end], c[end], c[end]]
    result = [ppm_stencil(c2[(i - 3):(i + 4)], v[(i - 3):(i - 2)], Δt, Δz) for i in 4:9]
    @test result ≈ [6.0, 6.0, 5.2, 5.0, 5.8, 6.0]
end

u = repeat([10.0, 8, 6, 4, 5, 7]', 6, 1)
v = repeat([0.0, 2, 4, 6, 4, 2], 1, 6)
c = collect(reshape(1.0:36, 6, 6))
Δt = 10
lons = -90:30:60
lats = -70:30:80

c_out = advect2D(Δt, 6, 6, c, u, v, lons, lats, l94_stencil)

@test c_out ≈
      [3.9996758431717745 27.999041560938153 51.999741040919446 75.99954752332437 99.9986360821562 123.99895033186321;
       7.999682626397016 31.99957465837063 55.99976078189158 79.99960682862074 103.99912843125838 127.99921682736436;
       11.999652774976862 35.99964475716122 59.99975988865124 83.99961504780215 107.99921479385594 131.9992706307947;
       15.999597225268698 39.99976292370086 63.99994456097572 87.99985852339972 111.99950162297623 135.99963864458493;
       19.999386077695178 43.99984458910325 68.00016102937246 92.00009314262749 115.9996246218866 139.9999014575152;
       23.997233510747368 47.99923509733475 72.00009877379722 95.9995656781324 119.99753232647797 143.99834416314692]

c = collect(reshape(1.0:36, 6, 6))

c_out = advect2D(Δt, 6, 6, c, u, v, lons, lats, l94_stencil)

@test c_out ≈
      [3.9996758431717745 27.999041560938153 51.999741040919446 75.99954752332437 99.9986360821562 123.99895033186321;
       7.999682626397016 31.99957465837063 55.99976078189158 79.99960682862074 103.99912843125838 127.99921682736436;
       11.999652774976862 35.99964475716122 59.99975988865124 83.99961504780215 107.99921479385594 131.9992706307947;
       15.999597225268698 39.99976292370086 63.99994456097572 87.99985852339972 111.99950162297623 135.99963864458493;
       19.999386077695178 43.99984458910325 68.00016102937246 92.00009314262749 115.9996246218866 139.9999014575152;
       23.997233510747368 47.99923509733475 72.00009877379722 95.9995656781324 119.99753232647797 143.99834416314692]
