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

obs = Matrix{Float64}(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_exponential/PHTransformedDWaveOrderParameter.csv", DataFrame; header=false))
o = Op(8);
my_obs = Matrix(((-ada(o,1,4)+ada(o,3,2))+(-ada(o,1,4)+ada(o,3,2))')*((ada(o,5,8)-ada(o,7,6))+(ada(o,5,8)-ada(o,7,6))'))
expectations = []
my_expectations = []
for j in 2:2:50
    println(j)
    ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_exponential/rhoA_1_$j.csv", DataFrame; header=false))
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
ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_exponential/rhoA_1_2.csv", DataFrame; header=false));
u, v = eigen(my_obs)
rotate_ρ = transpose(conj(v))*ρ*v;
measure_res = []
for _ in 1:100000
    ind,_ = new_ρmeasure!(rotate_ρ)
    push!(measure_res,u[ind])
end
mean(measure_res)
std(measure_res)/sqrt(100000)
tr(ρ*my_obs)
#### Systematic study
means = []
stds = []
o = Op(8);
my_obs = Matrix(((-ada(o,1,4)+ada(o,3,2))+(-ada(o,1,4)+ada(o,3,2))')*((ada(o,5,8)-ada(o,7,6))+(ada(o,5,8)-ada(o,7,6))'))
u, v = eigen(my_obs);
for j in 2:2:8
    println(j)
    ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_exponential/rhoA_1_$j.csv", DataFrame; header=false))
    rotate_ρ = transpose(conj(v))*ρ*v;
    measure_res = []
    for _ in 1:100000000
        ind,_ = new_ρmeasure!(rotate_ρ)
        push!(measure_res,u[ind])
    end
    push!(means,mean(measure_res))
    push!(stds,std(measure_res)/sqrt(100000000))
end


## save data
# data_saved = Dict("means"=>means,"stds"=>stds,"my_expectations"=>my_expectations)
# open("snr_plot_exponential_data.json", "w") do file
#     JSON.print(file, data_saved)
# end
# Read data
saved_data_power_law = JSON.parsefile("snr_plot_data.json")
means_power_law = saved_data_power_law["means"]
stds_power_law = saved_data_power_law["stds"]
my_expectations_power_law = saved_data_power_law["my_expectations"]

saved_data = JSON.parsefile("snr_plot_exponential_data.json")
means = saved_data["means"]
stds = saved_data["stds"]
my_expectations = saved_data["my_expectations"]

means
stds
means[1:3]
collect(2:2:6)
begin
    plt.errorbar(collect(2:2:6),means[1:3],yerr=stds[1:3])
    plt.plot(collect(2:2:10),(my_expectations)[1:5],"o-")
    plt.ylim([1e-6,0.6])
    plt.yscale("log")
    plt.xscale("log") 
    plt.show()
end
individual_stds = []
for k in 1:4
    j = collect(2:2:8)[k]
    println(j)
    push!(individual_stds,stds[k]*sqrt(1000000/((50-j-1)*1000)))
end
individual_stds_power_law = []
for k in 1:7
    j = collect(2:2:50)[k]
    println(j)
    push!(individual_stds_power_law,stds_power_law[k]*sqrt(1000000/((50-j-1)*1000)))
end
individual_stds
np=pyimport("numpy")
individual_stds_power_law[1:7]
custom_ind_stds_power_law = [0.00728 0.0074 0.00758 0.00777 0 0 0;0.00728 0.0074 0.00758 0.00777 0.0079737 0.0081796 0.0084157]
individual_stds
custom_ind_stds = [0.0013123240346245792 0.0012226838673797674 0 0;0.0013123240346245792 0.0012226838673797674 0.0012460984939478839 0.0012756231672939892]
begin
    fig = plt.figure(figsize = (3.5,3.3))
    plt.errorbar(collect(2:2:50)[1:7],abs.(means_power_law)[1:7],yerr=custom_ind_stds_power_law,ls="none",marker="o",markerfacecolor="#90BCD5", markeredgecolor="#3a76ae", ecolor="#3a76ae",markersize = 6)
    plt.plot([collect(2:2:50)[5], collect(2:2:50)[5]], [1e-10, abs(means_power_law[5])], linestyle="--", color="#5BAEFF")
    plt.plot([collect(2:2:50)[6], collect(2:2:50)[6]], [1e-10, abs(means_power_law[6])], linestyle="--", color="#5BAEFF")
    plt.plot([collect(2:2:50)[7], collect(2:2:50)[7]], [1e-10, abs(means_power_law[7])], linestyle="--", color="#5BAEFF")
    plt.plot(collect(2:2:50)[1:13],(my_expectations_power_law)[1:13],linestyle="--",color="#9ed17b")

    plt.errorbar(collect(2:2:8),abs.(means),yerr=custom_ind_stds,ls="none",marker="o",markerfacecolor="#e3b0a3", markeredgecolor="#d65d48", ecolor="#d65d48",markersize = 6)
    plt.plot([collect(2:2:8)[3], collect(2:2:8)[3]], [1e-10, abs(means[3])], linestyle="--", color="#d65d48")
    plt.plot([collect(2:2:8)[4], collect(2:2:8)[4]], [1e-10, abs(means[4])], linestyle="--", color="#d65d48")
    plt.plot(collect(2:2:50)[1:13],(my_expectations)[1:13],linestyle="--",color="#FFD06F")
    plt.ylim([1e-5,1.0])
    plt.xlabel(raw"$r=|i-j|$",fontsize = 13)
    # plt.ylabel(raw"$d-wave pairing \n correlation$", fontsize=13)
    plt.ylabel(raw"Pairing correlation",fontsize = 13)
    plt.text(1.9, 0.000016, "1000 Meas.", fontsize=13, color="black")
    plt.yscale("log")
    plt.xscale("log") 
    plt.show()
    # fig.savefig("./dmrg_1k_measure_2curves.svg",bbox_inches = "tight")
end
# data = Dict("x"=>collect(2:2:50),"exp_mean"=>means,"exp_std"=>stds,"theory"=>my_expectations)
# open("./dmrg_1kk_measurement.json", "w") do io
#     JSON.print(io, data)
# end