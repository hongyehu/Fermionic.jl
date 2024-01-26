using Fermionic
using LinearAlgebra
using SparseArrays
using Arpack
using PythonCall
using KrylovKit
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

Ham_sub = SpinFullFermiHubbardSubspace([rows,cols],J,7.2,δ,6)


### Test for states ###
st = GHZSubspace([rows,cols],6)
length(st)
up_st, info = exponentiate(-1im*Ham_sub,3.0,st)
print("norm of updated state: ", st'*st,"\n")
test_up_st = exp(Matrix(-1im*Ham_sub*3.0))*st
Vector(up_st-test_up_st)

_, st = measure!(st)
st
dist = Categorical(abs.(st).^2)
plt.hist(rand(dist,10000))
plt.show()
