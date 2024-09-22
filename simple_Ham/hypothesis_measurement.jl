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
        "--left"
            help = "left paddle"
            arg_type = Int
            default = 1
        "--right"
            help = "right paddle"
            arg_type = Int
            default = 1
        "--samples"
            help = "number of samples"
            arg_type = Int
            default = 1000000
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
    G_file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/G/L$(parsed_args["L"])_type$(parsed_args["type"])_Delta$(parsed_args["Delta"])_mu$(parsed_args["mu"]).jld2"
    @load G_file_name G
    dis = parsed_args["distance"];
    signs = [1.0,-1.0,-1.0,1.0];
    left = parsed_args["left"]
    right = parsed_args["right"]
    println("left: $left, right: $right")
    obs = two_paddle_obs(left=left,right=right)
    u, v = eigen(Matrix(obs))
    if left==1 && right==1
        ρ = RDM_PH(2.0.*G,[[2,1],[2,2],[2,1+4*dis],[2,2+4*dis]])
    elseif left==1 && right==2
        ρ = RDM_PH(2.0.*G,[[2,1],[2,2],[3,2+4*dis],[4,2+4*dis]])
    elseif left==1 && right==3
        ρ = RDM_PH(2.0.*G,[[2,1],[2,2],[1,3+4*dis],[2,3+4*dis]])
    elseif left==1 && right==4
        ρ = RDM_PH(2.0.*G,[[2,1],[2,2],[3,3+4*dis],[3,4+4*dis]])
    elseif left==2 && right==1
        ρ = RDM_PH(2.0.*G,[[3,2],[4,2],[2,1+4*dis],[2,2+4*dis]])
    elseif left==2 && right==2
        ρ = RDM_PH(2.0.*G,[[3,2],[4,2],[3,2+4*dis],[4,2+4*dis]])
    elseif left==2 && right==3
        ρ = RDM_PH(2.0.*G,[[3,2],[4,2],[1,3+4*dis],[2,3+4*dis]])
    elseif left==2 && right==4
        ρ = RDM_PH(2.0.*G,[[3,2],[4,2],[3,3+4*dis],[3,4+4*dis]])
    elseif left==3 && right==1
        ρ = RDM_PH(2.0.*G,[[1,3],[2,3],[2,1+4*dis],[2,2+4*dis]])
    elseif left==3 && right==2
        ρ = RDM_PH(2.0.*G,[[1,3],[2,3],[3,2+4*dis],[4,2+4*dis]])
    elseif left==3 && right==3
        ρ = RDM_PH(2.0.*G,[[1,3],[2,3],[1,3+4*dis],[2,3+4*dis]])
    elseif left==3 && right==4
        ρ = RDM_PH(2.0.*G,[[1,3],[2,3],[3,3+4*dis],[3,4+4*dis]])
    elseif left==4 && right==1
        ρ = RDM_PH(2.0.*G,[[3,3],[3,4],[2,1+4*dis],[2,2+4*dis]])
    elseif left==4 && right==2
        ρ = RDM_PH(2.0.*G,[[3,3],[3,4],[3,2+4*dis],[4,2+4*dis]])
    elseif left==4 && right==3
        ρ = RDM_PH(2.0.*G,[[3,3],[3,4],[1,3+4*dis],[2,3+4*dis]])
    elseif left==4 && right==4
        ρ = RDM_PH(2.0.*G,[[3,3],[3,4],[3,3+4*dis],[3,4+4*dis]])
    end
    rotate_ρ = v'*ρ*v
    measure_res = []
    for kt in 1:parsed_args["samples"]
        if kt%1000==0
            println("left($left), right($right)sample up to ($kt)")
        end
        ind,_ = ρmeasure!(rotate_ρ)
        push!(measure_res,u[ind])
    end
    println("true mean: $(tr(ρ*obs))")
    res = Dict("measure_res"=>measure_res)
    file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/hypothesis_data/L$(parsed_args["L"])_type$(parsed_args["type"])_dis$(parsed_args["distance"])_Delta$(parsed_args["Delta"])_mu$(parsed_args["mu"])_samples$(parsed_args["samples"])_left$(parsed_args["left"])_right$(parsed_args["right"]).json"
    open(file_name, "w") do io
        JSON3.write(io, res)
    end
end
main()
