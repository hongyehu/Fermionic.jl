using LinearAlgebra
using SparseArrays
using Fermionic
using Arpack
using Optimization
include("Hamiltonian.jl")

# J = 2.0
# U = 2.0
# t = 1.0
# δ1up = 0.0
# δ1down = 0.0
# δ2up = 0.0
# δ2down = 0.0
# θ = 0.4
particle_sector = 1
num_ham = 4
Js = rand(num_ham)
Us = rand(num_ham)
θs = rand(num_ham)
t = 1.0
δ1ups = rand(num_ham)
δ1downs = rand(num_ham)
δ2ups = rand(num_ham)
δ2downs = rand(num_ham)
map = zeros(ComplexF64, 16,16)
for i in 1:num_ham
    J = Js[i]
    U = Us[i]
    θ = θs[i]
    δ1up = δ1ups[i]
    δ1down = δ1downs[i]
    δ2up = δ2ups[i]
    δ2down = δ2downs[i]
    map += ShadowMap(J,U,t,δ1up,δ1down,δ2up,δ2down,θ,particle_sector)/4.0
end
eigen_res = eigen(map)
eigen_res.values