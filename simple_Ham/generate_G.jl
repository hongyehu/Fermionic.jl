using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using Statistics
using ArgParse
using JSON3
using JLD2

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
        "--Delta"
            help = "Delta for BCS theory"
            arg_type = Float64
            default = 0.1
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
    L = parsed_args["L"]
    if parsed_args["type"]==1
        type = :d
    else
        type = :s
    end
    G = BCS_G(L,type,Δ=parsed_args["Delta"],μ=parsed_args["mu"]);
    file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/G/L$(parsed_args["L"])_type$(parsed_args["type"])_Delta$(parsed_args["Delta"])_mu$(parsed_args["mu"]).jld2"
    @save file_name G
end
main()
# file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/G/L16_type1_Delta5.0_mu0.5.jld2"
# @load file_name G
# G