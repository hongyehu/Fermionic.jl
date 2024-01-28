using Fermionic
using LinearAlgebra
using SparseArrays
using Arpack
using PythonCall
using KrylovKit
using LinearSolve
# using Distributions
plt = pyimport("matplotlib.pyplot")
rows = 2
cols = 3
J = Dict{Any,Any}()
for j in 1:cols-1
    J[((1,j),(1,j+1))] = 3.0
    J[((2,j),(2,j+1))] = 3.0
end
for j in 1:cols
    J[((1,j),(2,j))] = 2.0
end
δ = Dict{Any,Any}()
for i in 1:rows
    for j in 1:cols
        δ[((i,j),(i,j))] = 0.0
    end
end

function Hopping(Lattice::AbstractArray{<:Int}, Jh::Float64, Jv::Float64)
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        J = Dict{Any,Any}()
        # Horizontal hopping
        for j in 1:cols-1
            for i in 1:rows
                J[((i,j),(i,j+1))] = Jh
            end
        end
        # Vertical hopping
        for j in 1:cols
            for i in 1:rows-1
                J[((i,j),(i+1,j))] = Jv
            end
        end
        return J
    else
        error("Dimension not implemented")
    end
end
test = randn(1000)
plt.hist(test)
plt.show()
function random_δ(Lattice::AbstractArray{<:Int}, μ::Float64, var::Float64)
    δ = Dict{Any,Any}()
    for i in 1:rows
        for j in 1:cols
            δ[((i,j),(i,j))] = randn()*sqrt(var) + μ
        end
    end
    return δ
end

J = Hopping([rows,cols],3.0,2.0)
δ = random_δ([rows,cols],0.0,0.1)
Ham_sub = SpinFullFermiHubbardSubspace([rows,cols],J,7.2,δ,6)
test_dict = Dict{Any,Any}()
test_dict[1] = Dict{Any,Any}("Jh"=>2.0,"Jv"=>1.5,"U"=>2.0,"δ"=>random_δ([rows,cols],0.0,0.1))
test_dict[2] = Dict{Any,Any}("Jh"=>2.0,"Jv"=>1.5,"U"=>2.0,"δ"=>random_δ([rows,cols],0.0,0.1))
test_dict[1]["δ"]


### Test for states ###
st = GHZSubspace([rows,cols],6)
length(st)
up_st, info = exponentiate(Ham_sub,3.0*im,st; ishermitian = true)
up_st
print("norm of updated state: ", up_st'*up_st,"\n")
test_up_st = exp(Matrix(1im*Ham_sub*3.0))*st
print("norm of updated state: ", test_up_st'*test_up_st,"\n")
norm(Vector(up_st-test_up_st),1)

_, st = measure!(st)
st
dist = Categorical(abs.(st).^2)
plt.hist(rand(dist,10000))
plt.show()

b, ind = basis_m(5,2)
b[1,:]
st = SpinFullRandomStateSubspace([rows,cols],true,6)

"""
=======
Input: 
    - Lattice: AbstractArray{<:AbstractFloat}
    - Params: Dict{Any,Any}. Example: Dict(1=>Dict("Jh"=>2.0,"Jv"=>1.5,"U"=>2.0,"δ"=>Dict{Any,Any}(),"t"=>0.4))
    - N: Int
=======
Output:
    - basis: AbstractArray{<:AbstractFloat}
"""
function ShadowBasis(
    Lattice::AbstractArray{<:Int},
    Params::Dict{Any,Any},
    N::Int,
    )
    """
    Row first ordering
    """
    dimension = length(Lattice)
    if dimension == 1
        error("Dimension not implemented")
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        subspace_dim = binomial(2*rows*cols,N)
        # sub_basis, _ = basis_m(2*rows*cols,N)
        Results = Dict{Any,Any}()
        for (key,_) in Params
            print("key: ", key,"\n")
            tmp = Dict{Any,Any}()
            Ham = SpinFullFermiHubbardSubspace([rows,cols],Hopping([rows,cols],Params[key]["Jh"],Params[key]["Jv"]),Params[key]["U"],Params[key]["δ"],N)
            ones_vector = ones(ComplexF64, binomial(2*rows*cols,particle_number))
            sparse_identity_matrix = spdiagm(0 => ones_vector)
            batch_evolve_state = sparse(exp(Matrix(1im*Ham*Params[key]["t"]))*sparse_identity_matrix)
            for k in 1:subspace_dim
                print("k: ", k,"\n")
                tmp[k] = batch_evolve_state[:,k]*batch_evolve_state[:,k]'
            end
            # for k in 1:subspace_dim
            #     print("k: ", k,"\n")
            #     # evolve_state = exp(Matrix(1im*Ham*Params[key]["t"]))*bitstring_state
            #     sparse_identity_matrix[:,k],info = exponentiate(Ham,Params[key]["t"]*1im,sparse_identity_matrix[:,k]; ishermitian = true,eager=true)
            #     tmp[k] = sparse_identity_matrix[:,k]*sparse_identity_matrix[:,k]'
            # end
            Results[key] = tmp
        end
        return Results
    else
        error("Dimension not implemented")
    end
    return basis
end
"""
Input: 
    - Lattice: AbstractArray{<:AbstractFloat}
    - N: Int
    - l: Int Distance of the pairing operator
"""
function PairingOperator(Lattice::AbstractArray{<:Int}, N::Int, I::Int, J::Int, l::Int)
    dimension = length(Lattice)
    if dimension == 1
        error("Dimension not implemented")
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        o = Op_fixed(2*rows*cols, N)
        pairing = ada(o,Cartesian2Index([I,J],[rows,cols],1),Cartesian2Index([I,J+l],[rows,cols],1))*ada(o,Cartesian2Index([I+1,J],[rows,cols],2),Cartesian2Index([I+1,J+l],[rows,cols],2))
        return pairing
    else
        error("Dimension not implemented")
    end
end


## Test for shadow basis
rows = 2
cols = 2
Params_test = Dict{Any,Any}()
num_of_Ham = 1000
for i in 1:num_of_Ham
    δ = random_δ([rows,cols],0.0,2.0)
    Params_test[i] = Dict{Any,Any}("Jh"=>2.0,"Jv"=>1.5,"U"=>2.0,"δ"=>δ,"t"=>1.0)
end

# δ = random_δ([rows,cols],0.0,0.1)
# ind = 1
# for t in range(1.0,stop=100.0,length=10)
#     Params_test[ind] = Dict{Any,Any}("Jh"=>3.0,"Jv"=>2.0,"U"=>5.4,"δ"=>δ,"t"=>t)
#     ind += 1
# end
# δ = random_δ([rows,cols],0.0,0.1)
# for t in range(1.0,stop=100.0,length=10)
#     Params_test[ind] = Dict{Any,Any}("Jh"=>2.7,"Jv"=>4.3,"U"=>1.4,"δ"=>δ,"t"=>t)
#     ind += 1
# end

particle_number = 4
# target_operator = PairingOperator([rows,cols],particle_number,1,1,1)
o = Op_fixed(2*rows*cols, particle_number)
target_operator  = ada(o,Cartesian2Index([1,1],[rows,cols],1),Cartesian2Index([1,1],[rows,cols],1))
Lattice = [rows,cols]
shadowbasis_test = ShadowBasis(Lattice,Params_test,particle_number)
reshaped_basis = spzeros(ComplexF64,binomial(2*rows*cols,particle_number)^2,length(shadowbasis_test)*length(shadowbasis_test[1]))
shadowbasis_test[1]

for i in 1:length(shadowbasis_test)
    for j in 1:length(shadowbasis_test[i])
        print("i: ", i, " j: ", j,"\n")
        reshaped_basis[:,(i-1)*length(shadowbasis_test[i])+j] = reshape(shadowbasis_test[i][j],binomial(2*rows*cols,particle_number)^2)
    end
end
reshape(shadowbasis_test[21][1],binomial(2*rows*cols,particle_number)^2)
rank(reshaped_basis)
rank(hcat([reshaped_basis,reshape(target_operator,binomial(2*rows*cols,particle_number)^2)]...))



## Test for grouping bitstring state
@time Ham_sub = SpinFullFermiHubbardSubspace([rows,cols],Hopping([rows,cols],3.0,2.0),7.2,δ,particle_number);
@time exp(Matrix(1im*Ham_sub*3.0));
ones_vector = ones(ComplexF64, binomial(2*rows*cols,particle_number))
sparse_identity_matrix = spdiagm(0 => ones_vector)

up12 = exp(Matrix(1im*Ham_sub*3.0))*sparse_identity_matrix
up12

@time test = exponentiate(Ham_sub,3.0*im,sparse_identity_matrix[:,1]; ishermitian = true);
@time exp(Matrix(1im*Ham_sub*3.0))*sparse_identity_matrix[:,1];

binomial(2*2*4,3)
## Test for linear solution
sigma_x = [0 1.0; 1.0 0]
sigma_y = [0 -1.0im; 1im 0]
sigma_z = [1 0; 0 -1]
identity = [1 0; 0 1]
target_rho = 0.5*(identity + 0.3*sigma_x+0.3*sigma_y+0.6*sigma_z)
vec_basis = hcat([reshape(identity,4),reshape(sigma_y,4),reshape(sigma_z,4)]...)
prob = LinearProblem(vec_basis, reshape(target_rho,4))
sol = solve(prob)
sol.u
rank(vec_basis)
rank([vec_basis reshape(target_rho,4)])