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
o = Op(8);
my_obs = Matrix(((-ada(o,1,4)+ada(o,3,2))+(-ada(o,1,4)+ada(o,3,2))')*((ada(o,5,8)-ada(o,7,6))+(ada(o,5,8)-ada(o,7,6))'))
expectations = []
my_expectations = []
for j in 2:2:50
    println(j)
    ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData/rhoA_1_$j.csv", DataFrame; header=false))
    push!(expectations,(-1)^mod(j-1,2)*tr(ρ*obs))
    push!(my_expectations,tr(ρ*my_obs))
end
begin
    plt.plot(collect(2:2:50),(expectations.*2.0),"-")
    plt.plot(collect(2:2:50),(my_expectations),"o")
    plt.yscale("log")
    plt.xscale("log") 
    plt.show()
end
## Simulating measurement
function new_ρmeasure!(state::Matrix{Float64})
    prob = real(diag(state))
    prob = prob./sum(prob)
    dist = Categorical(prob)
    idn = rand(dist,1)[1]
    state = spzeros(Complex{Float64},size(state))
    state[idn,idn] = 1.0
    return idn,state
end
o = Op(8);
my_obs = Matrix(((-ada(o,1,4)+ada(o,3,2))+(-ada(o,1,4)+ada(o,3,2))')*((ada(o,5,8)-ada(o,7,6))+(ada(o,5,8)-ada(o,7,6))'))
ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData/rhoA_1_2.csv", DataFrame; header=false));
u, v = eigen(my_obs)
rotate_ρ = transpose(conj(v))*ρ*v;
measure_res = []
for _ in 1:10000
    ind,_ = new_ρmeasure!(rotate_ρ)
    push!(measure_res,u[ind])
end
mean(measure_res)
std(measure_res)/sqrt(10000)
tr(ρ*my_obs)
#### Systematic study
means = []
stds = []
o = Op(8);
my_obs = Matrix(((-ada(o,1,4)+ada(o,3,2))+(-ada(o,1,4)+ada(o,3,2))')*((ada(o,5,8)-ada(o,7,6))+(ada(o,5,8)-ada(o,7,6))'))
u, v = eigen(my_obs);
for j in 2:2:50
    println(j)
    ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData/rhoA_1_$j.csv", DataFrame; header=false))
    rotate_ρ = transpose(conj(v))*ρ*v;
    measure_res = []
    for _ in 1:1000000
        ind,_ = new_ρmeasure!(rotate_ρ)
        push!(measure_res,u[ind])
    end
    push!(means,mean(measure_res))
    push!(stds,std(measure_res)/sqrt(1000000))
end


## save data
# data_saved = Dict("means"=>means,"stds"=>stds,"my_expectations"=>my_expectations)
# open("snr_plot_data.json", "w") do file
#     JSON.print(file, data_saved)
# end
# Read data
saved_data = JSON.parsefile("snr_plot_data.json")
means = saved_data["means"]
stds = saved_data["stds"]
my_expectations = saved_data["my_expectations"]
means
begin
    plt.errorbar(collect(2:2:50),abs.(means),yerr=stds)
    plt.plot(collect(2:2:50),(my_expectations),"o")
    plt.ylim([1e-4,0.6])
    plt.yscale("log")
    plt.xscale("log") 
    plt.show()
end
individual_stds = []
for k in 1:7
    j = collect(2:2:50)[k]
    println(j)
    push!(individual_stds,stds[j]*sqrt(1000000/((50-j-1)*1000)))
end
begin
    fig = plt.figure(figsize = (3.5,3.3))
    plt.errorbar(collect(2:2:50)[1:7],abs.(means)[1:7],yerr=individual_stds,ls="none",marker="o",markerfacecolor="#90BCD5", markeredgecolor="#5BAEFF", ecolor="#5BAEFF",markersize = 6)#"#90BCD5"
    plt.plot(collect(2:2:50)[1:13],(my_expectations)[1:13],linestyle="--",color="#FFD06F")
    plt.ylim([1e-3,1.0])
    plt.xlabel(raw"$r=|i-j|$",fontsize = 13)
    # plt.ylabel(raw"$d-wave pairing \n correlation$", fontsize=13)
    plt.ylabel(raw"Pairing correlation",fontsize = 13)
    plt.text(2, 0.0016, "1000 Meas.", fontsize=13, color="black")
    plt.yscale("log")
    plt.xscale("log") 
    plt.show()
    fig.savefig("./dmrg_1k_measure_2.svg",bbox_inches = "tight")
end
# data = Dict("x"=>collect(2:2:50),"exp_mean"=>means,"exp_std"=>stds,"theory"=>my_expectations)
# open("./dmrg_1kk_measurement.json", "w") do io
#     JSON.print(io, data)
# end