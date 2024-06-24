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
Δs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0]
for Δ in Δs
    filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/SNR_data/L56_type1_dis7_Delta$(Δ)_mu0.5_samples500_bootstrap10.json"
    data = JSON3.read(filename)
    tmp_mean = []
    tmp_std = []
    for bt in 1:10
        append!(tmp_mean,data["mean_$bt"])
        append!(tmp_std,data["std_$bt"])
    end
    append!(std_means,mean(tmp_std))
    append!(std_stds,std(tmp_std))
    append!(means,mean(tmp_mean))
end
SNR50 = √(50)*means./std_means
SNR100 = √(100)*means./std_means
SNR300 = √(300)*means./std_means

begin
    # Define the colors
    color_SNR50 = "#408fbf"
    color_SNR100 = "#89BED9"
    color_SNR300 = "#BFDAE9"
    f = plt.figure(figsize = (5.5,2.5))
    # plt.xlim([0.05,2.0])
    # plt.ylim([30,100])
    # plt.plot(Δs,20*log.(SNR),marker="o")
    plt.plot(Δs, SNR50, color=color_SNR50, label="50 Measurements")
    plt.plot(Δs, SNR100, color=color_SNR100, label="100 Measurements")
    plt.plot(Δs, SNR300, color=color_SNR300, label="300 Measurements")

    # Add labels and legend (optional)
    plt.xlabel("Δ")
    plt.ylabel("Linear SNR")
    plt.legend(frameon=false, ncol=1, loc="lower right", bbox_to_anchor=(1.0, -0.04))
    plt.show()
    # f.savefig("/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/SNR.pdf",bbox_inches = "tight")
end