using Fermionic
using LinearAlgebra
using SparseArrays
using Arpack
using PythonCall
using KrylovKit
using LinearSolve
using BenchmarkTools


rows = 2
cols = 2
δ = random_δ([rows,cols],0.0,2.0)
δ
@benchmark ShadowBasis_krylov([rows,cols],2.0,3.0,2.3,1.0,δ,rows*cols)
@profview ShadowBasis_krylov([2,3],2.0,3.0,2.3,1.0,δ,6)
@profview_allocs ShadowBasis_krylov([2,3],2.0,3.0,2.3,1.0,δ,6)

@benchmark ShadowBasis([rows,cols],2.0,3.0,2.3,1.0,δ,rows*cols)
@profview ShadowBasis([2,3],2.0,3.0,2.3,1.0,δ,6)

Params_test = Dict{Any,Any}()
num_of_Ham = 3
for i in 1:num_of_Ham
    δ = random_δ([rows,cols],0.0,2.0)
    Params_test[i] = Dict{Any,Any}("Jh"=>2.0,"Jv"=>1.5,"U"=>2.0,"δ"=>δ,"t"=>1.0)
end
particle_number = rows*cols
# target_operator = PairingOperator([rows,cols],particle_number,1,1,1)
# o = Op_fixed(2*rows*cols, particle_number)
# target_operator  = ada(o,Cartesian2Index([1,1],[rows,cols],1),Cartesian2Index([1,1],[rows,cols],1))
Lattice = [rows,cols]
binomial(2*rows*cols,particle_number)
Params_test[1]
@benchmark shadowbasis_test = ShadowBasis_krylov(Lattice,Params_test,particle_number)
@benchmark shadowbasis_test = ShadowBasis(Lattice,Params_test,particle_number)

ones_vector = ones(ComplexF64, binomial(2*rows*cols,particle_number))
sparse_identity_matrix = spdiagm(0 => ones_vector)
key = 1
Ham = SpinFullFermiHubbardSubspace([rows,cols],Hopping([rows,cols],Params_test[key]["Jh"],Params_test[key]["Jv"]),Params_test[key]["U"],Params_test[key]["δ"],particle_number)
@benchmark exponentiate(Ham,ComplexF64(1.0*1im),sparse_identity_matrix[:,1]; ishermitian = true,eager=true)
@benchmark exp(Matrix(1im*Ham*1.0))*sparse_identity_matrix[:,1]

st = @view sparse_identity_matrix[:, 1]