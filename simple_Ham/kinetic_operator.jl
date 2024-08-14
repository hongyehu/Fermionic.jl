using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using Statistics
using ArgParse
using JSON3
using JLD2
using PythonCall
using CSV
using DataFrames
using Distributions
using JSON
plt = pyimport("matplotlib.pyplot")
o = Op_fixed(4,2); 
basis(o)
obs = ada(o,1,3)+ada(o,3,1)+ada(o,2,4)+ada(o,4,2)
u,v=eigen(Matrix(obs))
obs2 = [0 1 0 0 1 0;1 0 0 0 0 1;0 0 0 0 0 0;0 0 0 0 0 0;1 0 0 0 0 1;0 1 0 0 1 0]
u2,v2 = eigen(Matrix(obs2))
v
test = Matrix([1 1 0 0 1 1;1 -1 0 0 -1 1;1 0 1 -1 0 1;1 0 -1 1 0 1;1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0;0 0 0 0 0 1;1 0 1 1 0 -1;])'
cv = [v v2]
M_basis = zeros(36,size(cv)[2])
for i in 1:size(cv)[2]
    M_basis[:,i]=reshape(cv[:,i]*cv[:,i]',36)
end
rank(M_basis)