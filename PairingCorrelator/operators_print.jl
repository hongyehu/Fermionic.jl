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
cols = 1
U1number = 1
op_size = 2*rows*cols
o = Op_fixed(op_size,U1number);
i=2
j=1
global_spin_X_rot([rows,cols],U1number)
n = ada(o,Cartesian2Index([i,j],[rows,cols],2),Cartesian2Index([i,j],[rows,cols],2))
Int = sum([ada(o,i,i)*ada(o,i+1,i+1) for i in 1:2:2*rows*cols-1])
Hopping = (ada(o,Cartesian2Index([i,j],[rows,cols],1),Cartesian2Index([i+1,j],[rows,cols],1))+
                    ada(o,Cartesian2Index([i+1,j],[rows,cols],1),Cartesian2Index([i,j],[rows,cols],1))+
                    ada(o,Cartesian2Index([i,j],[rows,cols],2),Cartesian2Index([i+1,j],[rows,cols],2))+
                    ada(o,Cartesian2Index([i+1,j],[rows,cols],2),Cartesian2Index([i,j],[rows,cols],2))
                )
