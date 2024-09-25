using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using Statistics
using ArgParse
using JSON3
using JLD2
using PythonCall
include("./analytical_2d_bond.jl")
include("./PHTransformed_2d_bond.jl")
op = Op(4)
basis(op)
Ob = ada(op,1,4)-ada(op,2,3)+(ada(op,1,4)-ada(op,2,3))'
us, vs = eigen(Matrix(Ob))
begin
    plt.imshow(Matrix(Ob))
    plt.show()
end
begin
    plt.plot(us,"o")
    plt.show()
end
op = Op_fixed(4,1)
Ob = ada(op,1,4)-ada(op,2,3)+(ada(op,1,4)-ada(op,2,3))'
eigen(Matrix(Ob))
