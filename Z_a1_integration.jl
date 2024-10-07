using Combinatorics
using LinearAlgebra
using Cuba
using DelimitedFiles
include("/ssd/ri47hud/Projects/3d Cosmology/codes/definitions.jl")
include("/ssd/ri47hud/Projects/3d Cosmology/codes/1slice.jl")
include("/ssd/ri47hud/Projects/3d Cosmology/codes/wynn.jl")

ϕ0 = 2.0
ϕ2 = 4.0
m = 0.05
a0 = 10.0
a2 = 30.0

ref_step = 0.1
b_range = 0.5:ref_step:50.0
N = length(b_range)
Z_Re = zeros(Float64, N-1, N-1)
Z_Im = zeros(Float64, N-1, N-1)

Threads.@threads for i in 2:N

    for j in 2:N

        Z_Re[i-1, j-1] = vegas((x,f) -> f[1] = real(1/((1 - x[1])^2) * Ampl_1slice(a0,  1/(1 - x[1]) - 1/2, a2, b_range[i], b_range[j], ϕ0, ϕ2, m)), 1, 1, minevals=1e5, maxevals=1e6)[1][1]
        Z_Im[i-1, j-1] = vegas((x,f) -> f[1] = imag(1/((1 - x[1])^2) * Ampl_1slice(a0,  1/(1 - x[1]) - 1/2, a2, b_range[i], b_range[j], ϕ0, ϕ2, m)), 1, 1, minevals=1e5, maxevals=1e6)[1][1]

    end

end

writedlm("/ssd/ri47hud/Projects/3d Cosmology/codes/data/Z_Re_a0=10_a2=30_phi0=2_phi2=4_m=005_step=01.txt", Z_Re)
writedlm("/ssd/ri47hud/Projects/3d Cosmology/codes/data/Z_Im_a0=10_a2=30_phi0=2_phi2=4_m=005_step=01.txt", Z_Im)