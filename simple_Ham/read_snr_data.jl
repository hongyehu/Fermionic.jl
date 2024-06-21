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
filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/SNR_data/L56_type1_dis7_Delta0.1_mu0.5_samples500_bootstrap1.json"
data = JSON3.read(filename)
data["std_1"]