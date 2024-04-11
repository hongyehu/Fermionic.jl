using Fermionic
using LinearAlgebra
using SparseArrays
using Arpack
using PythonCall
using KrylovKit
using LinearSolve
using BenchmarkTools
using ProgressMeter
plt = pyimport("matplotlib.pyplot")

binomial(2*rows*cols,particle_number)
rows = 2
cols = 1
particle_number = 2
print("local Hilbert space dimension: ", binomial(2*rows*cols,particle_number),"\n")
num_of_Hamiltonians = 1000
shadow_map = zeros(ComplexF64, binomial(2*rows*cols,particle_number)^2,binomial(2*rows*cols,particle_number)^2)
@showprogress for i in 1:num_of_Hamiltonians
    δ = random_δ([rows,cols],1.0,0.0)
    t = 10.0
    Jv = 5.0*rand()+1
    Jh = 0.0
    U = 5.0+1
    zfield = rand()*10
    θ = rand()*π/2
    tmp_basis = ShadowBasis_Global_Xrot([rows,cols],Jh,Jv,U,t,zfield,δ,particle_number,θ)
    @show size(tmp_basis[1])
    for j in 1:binomial(2*rows*cols,particle_number)
        flat_basis = reshape(tmp_basis[j],binomial(2*rows*cols,particle_number)^2)
        shadow_map += flat_basis*flat_basis'
    end
end
shadow_map = shadow_map/num_of_Hamiltonians
eigen_res = eigen(shadow_map)
eigen_res.values
plt.plot(abs.(eigen_res.values), "s-")
plt.show()
print(eigen_res.values)

rows = 2
cols = 1
all_eigen_values = []
for pn in 1:3
    particle_number = pn
    num_of_Hamiltonians = 50
    shadow_map = zeros(ComplexF64, binomial(2*rows*cols,particle_number)^2,binomial(2*rows*cols,particle_number)^2)
    @showprogress for i in 1:num_of_Hamiltonians
        δ = random_δ([rows,cols],0.0,0.0)
        t = 1.0
        Jv = rand()*10.0
        Jh = 0.0
        U = 10.0 #rand()*1.0+1
        zfield = rand()*10
        θ = rand()*π/2
        tmp_basis = ShadowBasis_Global_Xrot([rows,cols],Jh,Jv,U,t,zfield,δ,particle_number,θ)
        @show size(tmp_basis[1])
        for j in 1:binomial(2*rows*cols,particle_number)
            flat_basis = reshape(tmp_basis[j],binomial(2*rows*cols,particle_number)^2)
            shadow_map += flat_basis*flat_basis'
        end
    end
    shadow_map = shadow_map/num_of_Hamiltonians
    eigen_res = eigen(shadow_map)
    push!(all_eigen_values,eigen_res.values)
end

plt.plot(all_eigen_values[1], "s-", label="U(1) number = 1")
plt.plot(all_eigen_values[2], "s-", label="U(1) number = 2")
plt.plot(all_eigen_values[3], "s-", label="U(1) number = 3")
plt.yscale("log")
plt.xlabel("Eigenvalue index")
plt.ylabel("Eigenvalue")
plt.title("Eigenvalues of the shadow map")
plt.legend()
plt.show()

v = [1 + 2im, 3 - 4im, 5 + 6im]
v'
v*v'
tmp = [0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]
exp(Matrix(-1im*0.3*tmp))