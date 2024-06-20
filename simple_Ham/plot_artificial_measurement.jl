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

L = 40
Δ=0.2
μ=0.5
# d-wave
type = 1
d_true = []
d_mean = []
d_std = []
for dis in 1:7
    filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/L$(L)_type$(type)_dis$(dis)_Delta$(Δ)_mu$(μ)_samples500.json"
    data = JSON3.read(filename)
    append!(d_true,data["true_mean"])
    append!(d_mean,data["mean"])
    append!(d_std,data["std"])
end
#s-wave
type = 0
s_true = []
s_mean = []
s_std = []
for dis in 1:7
    filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/L$(L)_type$(type)_dis$(dis)_Delta$(Δ)_mu$(μ)_samples500.json"
    data = JSON3.read(filename)
    append!(s_true,data["true_mean"])
    append!(s_mean,data["mean"])
    append!(s_std,data["std"])
end
begin
    f = plt.figure(figsize = (6,4))
    plt.plot(collect(1:7),d_true,color="C0",linestyle="--")
    plt.errorbar(collect(1:7),d_mean,yerr=2*d_std,marker="o",color="C0")
    plt.plot(collect(1:7),s_true,color="C1",linestyle="--")
    plt.errorbar(collect(1:7),s_mean,yerr=2*d_std,marker="o",color="C1")
    plt.ylim([-3,3])
    plt.legend(["dwave theory","swave theory","dwave measure","swave measure"],fontsize = 9.5)
    plt.xlabel("distance |i-j|",fontsize = 14)
    plt.ylabel(raw"$\langle (\Delta^{\dagger}_i+\Delta_i)(\Delta^\dagger_j+\Delta_j)\rangle$",fontsize = 14)
    plt.show()
    # f.savefig("./direct_measurement_500shots.pdf",bbox_inches = "tight")
end