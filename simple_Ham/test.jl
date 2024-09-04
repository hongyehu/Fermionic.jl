using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using Statistics
using ArgParse
using JSON3
using JLD2
using PythonCall
include("./analytical_2d_bond.jl")
include("./PHTransformed_2d_bond.jl")
L=24
G_file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/G/L$(L)_type1_Delta0.3_mu0.5.jld2"
@load G_file_name G
bond_bond_correlation(G, 1,1,8)
PH_bond_bond_correlation(G, 1,1,8)