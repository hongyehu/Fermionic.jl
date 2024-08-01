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
plt = pyimport("matplotlib.pyplot")
function create_antidiagonal_matrix(n::Int)
    mat = zeros(Float64, n, n)
    for i in 1:n
        mat[i, n - i + 1] = 1.0
    end
    return mat
end
df = CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData/rhoA_1_3.csv", DataFrame; header=false)
Matrix(df)

obs = Matrix{Float64}(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData/PHTransformedDWaveOrderParameter.csv", DataFrame; header=false))
expectations = []
for j in 2:2:50
    println(j)
    ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData/rhoA_1_$j.csv", DataFrame; header=false))
    push!(expectations,(-1)^mod(j-1,2)*tr(ρ*obs))
end

begin
    plt.plot(collect(2:2:50),(expectations),"o-")
    plt.yscale("log")
    plt.xscale("log") 
    plt.show()
end
## Simulating measurement
obs = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData/PHTransformedDWaveOrderParameter.csv", DataFrame; header=false));
ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData/rhoA_1_2.csv", DataFrame; header=false));
u, v = eigen(obs)
rotate_ρ = transpose(conj(v))*ρ*v;
measure_res = []
for _ in 1:100
    ind,_ = ρmeasure!(rotate_ρ)
    push!(measure_res,u[ind])
end
u
begin
    plt.plot(u)
    plt.show()
end
tr(obs)
transpose(obs)≈obs

o = Op(8);
Matrix(basis(o))
Utransform = create_antidiagonal_matrix(256)
begin
    plt.imshow(Utransform'*ρ*Utransform)
    plt.show()
end