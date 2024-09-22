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
std_means = []
std_stds = []
means = []
filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/hypothesis_data/L40_type1_dis5_Delta0.1_mu0.5_samples10000000_left1_right1.json"
data = JSON3.read(filename)
measure_data = data["measure_res"]
mean(measure_data)
begin
    plt.hist(measure_data,edgecolor="black")
    plt.show()
end
begin
    plt.plot(measure_data,"o")
    plt.show()
end

# tmp_mean = []
# tmp_std = []
# for bt in 1:100
#     append!(tmp_mean,data["mean_$bt"])
#     append!(tmp_std,data["std_$bt"])
# end
# tmp_mean
# begin
#     plt.hist(tmp_mean,bins=20,edgecolor="black")
#     plt.show()
# end
# begin
#     plt.plot(tmp_mean,"o")
#     plt.show()
# end