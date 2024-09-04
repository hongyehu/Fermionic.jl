using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using Statistics
using ArgParse
using JSON3
using JLD2
include("./PHTransformed_2d_bond.jl")
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--L"
            help = "lattice size"
            arg_type = Int
            default = 24
        "--type"
            help = "type of BCS state, 0 for s-wave, 1 for d-wave"
            arg_type = Int
            default = 1
        "--samples"
            help = "number of samples"
            arg_type = Int
            default = 1000
        "--super_x"
            help = "super lattice x coordinate of the second bond"
            arg_type = Int
            default = 1
        "--super_y"
            help = "super lattice y coordinate of the second bond"
            arg_type = Int
            default = 2
        "--index"
            help = "index of the second bond"
            arg_type = Int
            default = 1
        "--Delta"
            help = "Delta for BCS theory"
            arg_type = Float64
            default = 0.3
        "--mu"
            help = "mu for BCS theory"
            arg_type = Float64
            default = 0.5
    end
    return parse_args(s)
end
function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    sample_res = Dict{Any, Any}();
    G_file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/G/L$(parsed_args["L"])_type1_Delta$(parsed_args["Delta"])_mu$(parsed_args["mu"]).jld2"
    @load G_file_name G
    ρ, obs = PH_bond_bond_obs_and_rho(G, parsed_args["super_x"],parsed_args["super_y"],parsed_args["index"])
    u, v = eigen(Matrix(obs))
    rotate_ρ = v'*ρ*v
    measure_res = []
    for _ in 1:parsed_args["samples"]
        ind,_ = ρmeasure!(rotate_ρ)
        push!(measure_res,u[ind])
    end
    sample_res["mean"] = mean(measure_res)
    sample_res["std"] = std(measure_res)/sqrt(parsed_args["samples"])
    println("sample_res: ",sample_res)
    file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/bond_bond_data/L$(parsed_args["L"])_type$(parsed_args["type"])_superx$(parsed_args["super_x"])_supery$(parsed_args["super_y"])_index$(parsed_args["index"])_Delta$(parsed_args["Delta"])_mu$(parsed_args["mu"])_samples$(parsed_args["samples"]).json"
    open(file_name, "w") do io
        JSON3.write(io, sample_res)
    end
end
main()
