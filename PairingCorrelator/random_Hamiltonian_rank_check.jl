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
# The following test shows 1000 Hamiltonian is tomography complete
rows = 2
cols = 2
particle_number = 2#rows*cols
num_of_Hamiltonians = 1000
shadow_map = zeros(ComplexF64, binomial(2*rows*cols,particle_number)^2,binomial(2*rows*cols,particle_number)^2)
@showprogress for i in 1:num_of_Hamiltonians
    δ = random_δ([rows,cols],0.0,2.0)
    t = 1.0
    Jv = rand()*5.0+1
    Jh = 0.0#rand()*5.0+1
    U = rand()*1.0+1
    tmp_basis = ShadowBasis([rows,cols],Jh,Jv,U,t,δ,particle_number)
    @show size(tmp_basis[1])
    for j in 1:binomial(2*rows*cols,particle_number)
        flat_basis = reshape(tmp_basis[j],binomial(2*rows*cols,particle_number)^2)
        shadow_map += flat_basis*flat_basis'
    end
end
shadow_map = shadow_map/num_of_Hamiltonians
eigen_res = eigen(shadow_map)
minimum(eigen_res.values)
# target_operator  = ada(Op_fixed(2*rows*cols, particle_number),Cartesian2Index([1,1],[rows,cols],1),Cartesian2Index([1,2],[rows,cols],1))
target_operator = Dwave([rows,cols],1,2,particle_number)
first_ind = findfirst(x -> x > 1e-4, eigen_res.values)
rank(eigen_res.vectors[:,first_ind:end])
rank(hcat([eigen_res.vectors[:,first_ind:end],reshape(target_operator,binomial(2*rows*cols,particle_number)^2)]...))

## Global X rotation
U = exp(Matrix(-1im*π/2*global_spin_X_rot([rows,cols],particle_number)))
plt.imshow(real(U))
plt.colorbar()
plt.show()

# test global X rotation
rows = 2
cols = 2
particle_number = 4#rows*cols
num_of_Hamiltonians = 500
shadow_map = zeros(ComplexF64, binomial(2*rows*cols,particle_number)^2,binomial(2*rows*cols,particle_number)^2)
@showprogress for i in 1:num_of_Hamiltonians
    δ = random_δ([rows,cols],0.0,0.0)
    t = 1.0
    Jv = rand()*5.0+1
    Jh = 0.0
    U = 5.0 
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
minimum(eigen_res.values)

# target_operator  = ada(Op_fixed(2*rows*cols, particle_number),Cartesian2Index([1,1],[rows,cols],1),Cartesian2Index([1,2],[rows,cols],1))
target_operator = transformed_Dwave([rows,cols],1,2,particle_number)
print((target_operator))
first_ind = findfirst(x -> x > 1e-4, eigen_res.values)
old_rank = rank(eigen_res.vectors[:,first_ind:end])
new_rank = rank(hcat([eigen_res.vectors[:,first_ind:end],reshape(target_operator,binomial(2*rows*cols,particle_number)^2)]...))
if old_rank == new_rank
    println("It is in the visible space!")
else
    println("Not in the visible space :(")
end

particle_number = 4
target_operator = transformed_Dwave([rows,cols],1,2,particle_number)
Matrix(target_operator)
sum(abs.(Matrix(target_operator)))
plt.imshow(real(Matrix(target_operator)))
plt.show()
plt.plot(diag(adjoint(target_operator)*target_operator))
plt.show()



# test local X rotation
rows = 2
cols = 2
particle_number = 2#rows*cols
num_of_Hamiltonians = 50
shadow_map = zeros(ComplexF64, binomial(2*rows*cols,particle_number)^2,binomial(2*rows*cols,particle_number)^2)
@showprogress for i in 1:num_of_Hamiltonians
    δ = random_δ([rows,cols],0.0,2.0)
    t = 1.0
    Jv = rand()*5.0+1
    Jh = 0.0#rand()*5.0+1
    U = rand()*1.0+1
    θ = random_θ([rows,cols],0.0,π/2)
    tmp_basis = ShadowBasis_Local_Xrot([rows,cols],Jh,Jv,U,t,δ,particle_number,θ)
    @show size(tmp_basis[1])
    for j in 1:binomial(2*rows*cols,particle_number)
        flat_basis = reshape(tmp_basis[j],binomial(2*rows*cols,particle_number)^2)
        shadow_map += flat_basis*flat_basis'
    end
end
shadow_map = shadow_map/num_of_Hamiltonians
eigen_res = eigen(shadow_map)
minimum(eigen_res.values)
# target_operator  = ada(Op_fixed(2*rows*cols, particle_number),Cartesian2Index([1,1],[rows,cols],1),Cartesian2Index([1,2],[rows,cols],1))
target_operator = transformed_Dwave([rows,cols],1,2,particle_number)
first_ind = findfirst(x -> x > 1e-4, eigen_res.values)
old_rank = rank(eigen_res.vectors[:,first_ind:end])
new_rank = rank(hcat([eigen_res.vectors[:,first_ind:end],reshape(target_operator,binomial(2*rows*cols,particle_number)^2)]...))
if old_rank == new_rank
    println("It is in the visible space!")
else
    println("Not in the visible space :(")
end


# test local X rotation (different columns share θ)
rows = 2
cols = 2
particle_number = 2#rows*cols
num_of_Hamiltonians = 50
shadow_map = zeros(ComplexF64, binomial(2*rows*cols,particle_number)^2,binomial(2*rows*cols,particle_number)^2)
@showprogress for i in 1:num_of_Hamiltonians
    δ = random_δ([rows,cols],0.0,2.0)
    t = 1.0
    Jv = rand()*5.0+1
    Jh = 0.0#rand()*5.0+1
    U = rand()*1.0+1
    θ = random_horizontal_δ([rows,cols],0.0,π/2)
    tmp_basis = ShadowBasis_Local_Xrot([rows,cols],Jh,Jv,U,t,δ,particle_number,θ)
    @show size(tmp_basis[1])
    for j in 1:binomial(2*rows*cols,particle_number)
        flat_basis = reshape(tmp_basis[j],binomial(2*rows*cols,particle_number)^2)
        shadow_map += flat_basis*flat_basis'
    end
end
shadow_map = shadow_map/num_of_Hamiltonians
eigen_res = eigen(shadow_map)
minimum(eigen_res.values)
# target_operator  = ada(Op_fixed(2*rows*cols, particle_number),Cartesian2Index([1,1],[rows,cols],1),Cartesian2Index([1,2],[rows,cols],1))
target_operator = transformed_Dwave([rows,cols],1,2,particle_number)
first_ind = findfirst(x -> x > 1e-4, eigen_res.values)
old_rank = rank(eigen_res.vectors[:,first_ind:end])
new_rank = rank(hcat([eigen_res.vectors[:,first_ind:end],reshape(target_operator,binomial(2*rows*cols,particle_number)^2)]...))
if old_rank == new_rank
    println("It is in the visible space!")
else
    println("Not in the visible space :(")
end

# test local X rotation (different rows share θ)
rows = 2
cols = 2
particle_number = 2#rows*cols
num_of_Hamiltonians = 50
shadow_map = zeros(ComplexF64, binomial(2*rows*cols,particle_number)^2,binomial(2*rows*cols,particle_number)^2)
@showprogress for i in 1:num_of_Hamiltonians
    δ = random_δ([rows,cols],0.0,2.0)
    t = 1.0
    Jv = rand()*5.0+1
    Jh = 0.0#rand()*5.0+1
    U = rand()*1.0+1
    θ = random_vertical_δ([rows,cols],0.0,π/2)
    tmp_basis = ShadowBasis_Local_Xrot([rows,cols],Jh,Jv,U,t,δ,particle_number,θ)
    @show size(tmp_basis[1])
    for j in 1:binomial(2*rows*cols,particle_number)
        flat_basis = reshape(tmp_basis[j],binomial(2*rows*cols,particle_number)^2)
        shadow_map += flat_basis*flat_basis'
    end
end
shadow_map = shadow_map/num_of_Hamiltonians
eigen_res = eigen(shadow_map)
minimum(eigen_res.values)
# target_operator  = ada(Op_fixed(2*rows*cols, particle_number),Cartesian2Index([1,1],[rows,cols],1),Cartesian2Index([1,2],[rows,cols],1))
target_operator = transformed_Dwave([rows,cols],1,2,particle_number)
first_ind = findfirst(x -> x > 1e-4, eigen_res.values)
old_rank = rank(eigen_res.vectors[:,first_ind:end])
new_rank = rank(hcat([eigen_res.vectors[:,first_ind:end],reshape(target_operator,binomial(2*rows*cols,particle_number)^2)]...))
if old_rank == new_rank
    println("It is in the visible space!")
else
    println("Not in the visible space :(")
end