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
            default = 40
        "--type"
            help = "type of BCS state, 0 for s-wave, 1 for d-wave"
            arg_type = Int
            default = 1
        "--measurements"
            help = "number of measurement"
            arg_type = Int
            default = 300
        "--distance"
            help = "distance between two flower operators"
            arg_type = Int
            default = 5
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
    paddle_means = Dict{String,Vector{Float64}}()
    for left in 1:4
        for right in 1:4
            println("left: $left, right: $right")
            filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/hypothesis_data/L$(parsed_args["L"])_type$(parsed_args["type"])_dis$(parsed_args["distance"])_Delta$(parsed_args["Delta"])_mu$(parsed_args["mu"])_samples10000000_left$(left)_right$(right).json"
            data = JSON3.read(filename)
            tmp = []
            for bt in 1:10000
                if bt % 100 == 0
                    println("left: $left, right: $right, bt: $bt")
                end
                append!(tmp,mean(data["measure_res"][(bt-1)*parsed_args["measurements"]+1:bt*parsed_args["measurements"]]))
            end
            paddle_means["$(left)_$(right)"] = tmp
        end
    end
    save_file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/hypothesis_data/processed_measurement$(parsed_args["measurements"])_L$(parsed_args["L"])_type$(parsed_args["type"])_dis$(parsed_args["distance"])_Delta$(parsed_args["Delta"])_mu$(parsed_args["mu"])_samples10000000.json"
    open(save_file_name, "w") do io
        JSON3.write(io, paddle_means)
    end
end
main()



# # 300 measurements, 10000 bootstraps
# paddle_means = Dict{String,Vector{Float64}}()
# signs = [1.0,-1.0,-1.0,1.0];
# for left in 1:4
#     for right in 1:4
#         println("left: $left, right: $right")
#         filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/hypothesis_data/L40_type1_dis5_Delta0.1_mu0.5_samples10000000_left$(left)_right$(right).json"
#         data = JSON3.read(filename)
#         tmp = []
#         for bt in 1:10000
#             if bt % 100 == 0
#                 println("left: $left, right: $right, bt: $bt")
#             end
#             append!(tmp,mean(data["measure_res"][(bt-1)*300+1:bt*300]))
#         end
#         paddle_means["$(left)_$(right)"] = tmp
#     end
# end
# paddle_means
# file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/hypothesis_data/processed_L$(parsed_args["L"])_type$(parsed_args["type"])_dis$(parsed_args["distance"])_Delta$(parsed_args["Delta"])_mu$(parsed_args["mu"])_samples$(parsed_args["samples"])_left$(parsed_args["left"])_right$(parsed_args["right"]).json"
# open(file_name, "w") do io
#     JSON3.write(io, res)
# end
# results = []
# for bt in 1:10000
#     correlator = 0.0
#     for left in 1:4
#         for right in 1:4
#             correlator += signs[left]*signs[right]*paddle_means["$(left)_$(right)"][bt]
#         end
#     end
#     append!(results,correlator)
# end
# begin
#     plt.hist(results,bins = 40,edgecolor="black")
#     plt.show()
# end