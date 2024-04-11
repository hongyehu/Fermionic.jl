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

row = 2; col = 1
particle_number = 2
op_size = 2*rows*cols
o = Op_fixed(op_size,particle_number); # fixed particle number subspace
opt1 = ada(o,Cartesian2Index([1,2],[rows,cols],1),Cartesian2Index([1,1],[rows,cols],2))