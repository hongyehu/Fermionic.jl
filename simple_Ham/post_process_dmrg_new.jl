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

df = CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_0.875filling/rhoA_1_3.csv", DataFrame; header=false)
Matrix(df)
o = Op(8);
my_obs = Matrix(((-ada(o,1,4)+ada(o,3,2))+(-ada(o,1,4)+ada(o,3,2))')*((ada(o,5,8)-ada(o,7,6))+(ada(o,5,8)-ada(o,7,6))'))
my_expectations = []
for j in 2:2:40
    println(j)
    ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_0.875filling/rhoA_1_$j.csv", DataFrame; header=false))
    push!(my_expectations,tr(ρ*my_obs))
end
begin
    plt.plot(collect(2:2:40),(my_expectations),"-")
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
ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_0.875filling/rhoA_1_3.csv", DataFrame; header=false));
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


#### Systematic study ==== doped case
means_dope = []
stds_dope = []
o = Op(8);
my_obs = Matrix(((-ada(o,1,4)+ada(o,3,2))+(-ada(o,1,4)+ada(o,3,2))')*((ada(o,5,8)-ada(o,7,6))+(ada(o,5,8)-ada(o,7,6))'))
u, v = eigen(my_obs);
for j in 2:2:14
    println(j)
    ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_0.875filling/rhoA_1_$j.csv", DataFrame; header=false))
    rotate_ρ = transpose(conj(v))*ρ*v;
    measure_res = []
    for _ in 1:10000000
        ind,_ = new_ρmeasure!(rotate_ρ)
        push!(measure_res,u[ind])
    end
    push!(means_dope,mean(measure_res))
    push!(stds_dope,std(measure_res)/sqrt(10000000))
end

## theory
my_obs = Matrix(((-ada(o,1,4)+ada(o,3,2))+(-ada(o,1,4)+ada(o,3,2))')*((ada(o,5,8)-ada(o,7,6))+(ada(o,5,8)-ada(o,7,6))'))
theory_dope = []
for j in 2:2:14
    println(j)
    ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_0.875filling/rhoA_1_$j.csv", DataFrame; header=false))
    push!(theory_dope,tr(ρ*my_obs))
end

begin
    plt.errorbar(collect(2:2:14),(means_dope), yerr=stds_dope)
    plt.plot(collect(2:2:14),(theory_dope),"-")
    plt.yscale("log")
    plt.xscale("log")
    plt.show()
end

## Save data
# data_saved = Dict("means"=>means_dope,"stds"=>stds_dope,"my_expectations"=>theory_dope)
# open("dmrg_new_data_doping.json", "w") do file
#     JSON.print(file, data_saved)
# end

#### Systematic study ==== half case
means_half = []
stds_half = []
o = Op(8);
my_obs = Matrix(((-ada(o,1,4)+ada(o,3,2))+(-ada(o,1,4)+ada(o,3,2))')*((ada(o,5,8)-ada(o,7,6))+(ada(o,5,8)-ada(o,7,6))'))
u, v = eigen(my_obs);
for j in 2:2:8
    println(j)
    ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_0.5filling/rhoA_1_$j.csv", DataFrame; header=false))
    rotate_ρ = transpose(conj(v))*ρ*v;
    measure_res = []
    for _ in 1:100000000
        ind,_ = new_ρmeasure!(rotate_ρ)
        push!(measure_res,u[ind])
    end
    push!(means_half,mean(measure_res))
    push!(stds_half,std(measure_res)/sqrt(100000000))
end
## theory
my_obs = Matrix(((-ada(o,1,4)+ada(o,3,2))+(-ada(o,1,4)+ada(o,3,2))')*((ada(o,5,8)-ada(o,7,6))+(ada(o,5,8)-ada(o,7,6))'))
theory_half = []
for j in 2:2:8
    println(j)
    ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_0.5filling/rhoA_1_$j.csv", DataFrame; header=false))
    push!(theory_half,tr(ρ*my_obs))
end
## Save data
# data_saved = Dict("means"=>means_half,"stds"=>stds_half,"my_expectations"=>theory_half)
# open("dmrg_new_data_half_filling.json", "w") do file
#     JSON.print(file, data_saved)
# end

# begin
#     plt.errorbar(collect(2:2:8),(means_half), yerr=stds_half)
#     plt.plot(collect(2:2:8),(theory_half),"-")
#     plt.yscale("log")
#     plt.xscale("log")
#     plt.show()
# end


individual_stds = []
for k in 1:4
    j = collect(2:2:8)[k]
    println(j)
    push!(individual_stds,stds_half[k]*sqrt(100000000/((40-j-1)*5000)))
end
individual_stds_power_law = []
for k in 1:7
    j = collect(2:2:50)[k]
    println(j)
    push!(individual_stds_power_law,stds_dope[k]*sqrt(10000000/((40-j-1)*5000)))
end

individual_stds_power_law
# custom_ind_stds_power_law = [0.0036740181 0.00300515 0.002965821 0.003245168 0 0 0;0.0036740181 0.00300515 0.002965821 0.003245168 0.003351522 0.0032639 0.00335182]
custom_ind_stds_power_law = [0.00585604443166159 0.0048189677604732195 0.004787814685307327 0.005277926176996058 0 0 0;0.00585604443166159 0.0048189677604732195 0.004787814685307327 0.005277926176996058 0.0054965525684719345 0.005403534650254163 0.00560868339147644]
individual_stds
# custom_ind_stds = [0.0041496806 0.0033666 0 0;0.0041496806 0.0033666 0.003940  0.0040]
# custom_ind_stds = custom_ind_stds./sqrt(5)
custom_ind_stds = [0.006614206296788662 0.006200451931057928 0 0; 0.006614206296788662 0.006200451931057928 0.006360768172772285 0.006560780471703199]./1.5
## Plot Figure
means_half[2]
custom_ind_stds[2]
collect(2:2:50)
begin
    fig = plt.figure(figsize = (3.5,3.3))
    plt.errorbar(collect(2:2:50)[1:7],abs.(means_dope)[1:7],yerr=custom_ind_stds_power_law,ls="none",marker="o",markerfacecolor="#9ed17b", markeredgecolor="#50AA4B", ecolor="#50AA4B",markersize = 6)
    plt.plot([collect(2:2:50)[5], collect(2:2:50)[5]], [1e-10, abs(means_dope[5])], linestyle="--", color="#50AA4B")
    plt.plot([collect(2:2:50)[6], collect(2:2:50)[6]], [1e-10, abs(means_dope[6])], linestyle="--", color="#50AA4B")
    plt.plot([collect(2:2:50)[7], collect(2:2:50)[7]], [1e-10, abs(means_dope[7])], linestyle="--", color="#50AA4B")
    plt.plot(collect(2:2:50)[1:7],(theory_dope)[1:7],linestyle="--",color="#9ed17b")

    plt.errorbar(collect(2:2:8),abs.(means_half),yerr=custom_ind_stds,ls="none",marker="o",markerfacecolor="#FFD06F", markeredgecolor="#F19D54", ecolor="#F19D54",markersize = 6)
    plt.plot([collect(2:2:8)[3], collect(2:2:8)[3]], [1e-10, abs(means_half[3])], linestyle="--", color="#F19D54")
    plt.plot([collect(2:2:8)[4], collect(2:2:8)[4]], [1e-10, abs(means_half[4])], linestyle="--", color="#F19D54")
    plt.plot(collect(2:2:50)[1:4],(theory_half)[1:4],linestyle="--",color="#FFD06F")
    plt.yscale("log")
    plt.xscale("log") 
    plt.ylim([1e-5,1.0])
    plt.xlabel(raw"$r=|i-j|$",fontsize = 13)
    # plt.ylabel(raw"$d-wave pairing \n correlation$", fontsize=13)
    plt.ylabel(raw"Pairing correlation",fontsize = 13)
    plt.text(1.9, 0.000016, "1000 Meas.", fontsize=13, color="black")
    
    plt.show()
    fig.savefig("./dmrg_1k_measure_2curves_U8_data3.svg",bbox_inches = "tight")
end


# Test
means_dope_test = []
stds_dope_test = []
o = Op(8);
my_obs = Matrix(((-ada(o,1,4)+ada(o,3,2))+(-ada(o,1,4)+ada(o,3,2))')*((ada(o,5,8)-ada(o,7,6))+(ada(o,5,8)-ada(o,7,6))'))
u, v = eigen(my_obs);
for j in 2:2:14
    println(j)
    ρ = Matrix(CSV.read("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/dmrg_data/rhoAData_0.875filling/rhoA_1_$j.csv", DataFrame; header=false))
    rotate_ρ = transpose(conj(v))*ρ*v;
    measure_res = []
    for _ in 1:100000
        ind,_ = new_ρmeasure!(rotate_ρ)
        push!(measure_res,u[ind])
    end
    push!(means_dope_test,mean(measure_res))
    push!(stds_dope_test,std(measure_res)/sqrt(100000))
end
begin
    plt.errorbar(collect(2:2:14),(means_dope_test), yerr=stds_dope_test)
    plt.plot(collect(2:2:14),(theory_dope),"-")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylim([1e-5,1.0])
    plt.show()
end