using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using PythonCall
using Statistics
plt = pyimport("matplotlib.pyplot")
obs = two_paddle_obs(left=1,right=1)
u, v = eigen(Matrix(obs))
tmp = v'*(Matrix(obs))*v
L = 10
type = :d
G = BCS_G(L,type)
ρ = RDM_PH(G,[[2,1],[2,2],[2,5],[2,6]])
tr(ρ*obs)

measure_times = 10000
rotate_ρ = v'*ρ*v
measure_res = []
for i in 1:measure_times
    ind,_ = ρmeasure!(rotate_ρ)
    push!(measure_res,u[ind])
end
mean(measure_res)