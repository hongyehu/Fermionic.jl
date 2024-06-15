using Fermionic
using PythonCall
using LinearAlgebra
plt = pyimport("matplotlib.pyplot")

particle_sector = 2
o = Op_fixed(4,particle_sector)
basis(o)
obs = ada(o,1,4)+ada(o,4,1)
eigen(Matrix(obs))

particle_sector = 2
o = Op_fixed(2,particle_sector)
basis(o)
obs = ada(o,1,2)+ada(o,2,1)
eigen(Matrix(obs))
15*14/2