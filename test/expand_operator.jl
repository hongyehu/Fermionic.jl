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
rows = 2
cols = 2
particle_number = 2#rows*cols
num_of_Hamiltonians = 1000
shadow_map = zeros(ComplexF64, binomial(2*rows*cols,particle_number)^2,binomial(2*rows*cols,particle_number)^2)
@showprogress for i in 1:num_of_Hamiltonians
    δ = random_δ([rows,cols],0.0,2.0)
    t = 200.0
    Jv = rand()*5.0+3
    Jh = rand()*5.0+2
    U = rand()*5.0+1
    tmp_basis = ShadowBasis([rows,cols],Jh,Jv,U,t,δ,particle_number)
    @show size(tmp_basis[1])
    for j in 1:binomial(2*rows*cols,particle_number)
        flat_basis = reshape(tmp_basis[j],binomial(2*rows*cols,particle_number)^2)
        shadow_map += flat_basis*flat_basis'
    end
end
shadow_map = shadow_map/num_of_Hamiltonians
eigen_res = eigen(shadow_map)
plt.plot(eigen_res.values, "o")
plt.yscale("log")
plt.ylim(1e-4,1e0)
plt.show()
eigen_res.values
first_ind = findfirst(x -> x > 1e-4, eigen_res.values)
plt.plot(eigen_res.values[first_ind:end], "o")
plt.yscale("log")
plt.show()
target_operator  = ada(Op_fixed(2*rows*cols, particle_number),Cartesian2Index([1,1],[rows,cols],1),Cartesian2Index([1,2],[rows,cols],1))
rank(eigen_res.vectors[:,first_ind:end])
rank(hcat([eigen_res.vectors[:,first_ind:end],reshape(target_operator,binomial(2*rows*cols,particle_number)^2)]...))