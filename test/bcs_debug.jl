using Fermionic
using PythonCall
using LinearAlgebra
plt = pyimport("matplotlib.pyplot")

# Minh's code
using LinearAlgebra, SparseArrays
o=Op(3);
basis(o)
ad(o,1)
# G = BCS_G(3, :d)
G = getBCSaij(3,:s)
G
L = Int64(sqrt(size(G)[1]/2))
Gtmp,cs = shuffle_G(G,[1,2,9,10],2*L^2);
Gtmp[1:4,1:4]

new_G,cs = shuffle_G(G,[1,2,9,10],18)
new_G[1:4,1:4]
round.(RDM(G, [1,2,9,10]),digits=4);