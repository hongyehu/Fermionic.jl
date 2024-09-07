using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using Statistics
using ArgParse
using JSON3
using JLD2
using PythonCall
plt = pyimport("matplotlib.pyplot")
cm = pyimport("matplotlib.cm")
patches = pyimport("matplotlib.patches")

function DoubleWell(x,a,b,c,shift,Δ)
    return a*((x-shift)^2-Δ^2)^2+b*((x-shift))+c
end
x = -5:0.01:5
y = [DoubleWell(x,1,-6,0,0,3) for x in x]
begin
    f = plt.figure(figsize = (2,2))
    plt.plot(x,y,linewidth=5)
    # plt.xlabel("x",fontsize = 14)
    # plt.ylabel("V(x)",fontsize = 14)
    plt.axis("off")
    plt.show()
    f.savefig("./double_well_potential_-6.pdf",bbox_inches = "tight")
end