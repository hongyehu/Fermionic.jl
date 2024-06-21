using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using Statistics
using ArgParse
using JSON3
using JLD2
## In the SNR experiment, the std we stored is the std of bare observables
## So the std of estimator should be std/sqrt(samples)
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--L"
            help = "lattice size"
            arg_type = Int
            default = 16
        "--type"
            help = "type of BCS state, 0 for s-wave, 1 for d-wave"
            arg_type = Int
            default = 1
        "--samples"
            help = "number of samples"
            arg_type = Int
            default = 500
        "--distance"
            help = "distance between two flower operators"
            arg_type = Int
            default = 1
        "--Delta"
            help = "Delta for BCS theory"
            arg_type = Float64
            default = 5.0
        "--mu"
            help = "mu for BCS theory"
            arg_type = Float64
            default = 0.5
        "--bootstrap"
            help = "number of bootstrap samples for estimating std of std estimation"
            arg_type = Int
            default = 10
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
    sample_res = Dict{Any, Any}();
    signs = [1.0,-1.0,-1.0,1.0];
    for bt in 1:parsed_args["bootstrap"]
        mean_value = 0.0;
        true_mean_value = 0.0;
        stds=[];
        for left in 1:4
            for right in 1:4
                println("bootstrap: $bt, left: $left, right: $right")
                obs = two_paddle_obs(left=left,right=right)
                u, v = eigen(Matrix(obs))
                # tmp = v'*(Matrix(obs))*v
                if left==1 && right==1
                    ρ = RDM_PH(G,[[2,1],[2,2],[2,1+4*dis],[2,2+4*dis]])
                elseif left==1 && right==2
                    ρ = RDM_PH(G,[[2,1],[2,2],[3,2+4*dis],[4,2+4*dis]])
                elseif left==1 && right==3
                    ρ = RDM_PH(G,[[2,1],[2,2],[1,3+4*dis],[2,3+4*dis]])
                elseif left==1 && right==4
                    ρ = RDM_PH(G,[[2,1],[2,2],[3,3+4*dis],[3,4+4*dis]])
                elseif left==2 && right==1
                    ρ = RDM_PH(G,[[3,2],[4,2],[2,1+4*dis],[2,2+4*dis]])
                elseif left==2 && right==2
                    ρ = RDM_PH(G,[[3,2],[4,2],[3,2+4*dis],[4,2+4*dis]])
                elseif left==2 && right==3
                    ρ = RDM_PH(G,[[3,2],[4,2],[1,3+4*dis],[2,3+4*dis]])
                elseif left==2 && right==4
                    ρ = RDM_PH(G,[[3,2],[4,2],[3,3+4*dis],[3,4+4*dis]])
                elseif left==3 && right==1
                    ρ = RDM_PH(G,[[1,3],[2,3],[2,1+4*dis],[2,2+4*dis]])
                elseif left==3 && right==2
                    ρ = RDM_PH(G,[[1,3],[2,3],[3,2+4*dis],[4,2+4*dis]])
                elseif left==3 && right==3
                    ρ = RDM_PH(G,[[1,3],[2,3],[1,3+4*dis],[2,3+4*dis]])
                elseif left==3 && right==4
                    ρ = RDM_PH(G,[[1,3],[2,3],[3,3+4*dis],[3,4+4*dis]])
                elseif left==4 && right==1
                    ρ = RDM_PH(G,[[3,3],[3,4],[2,1+4*dis],[2,2+4*dis]])
                elseif left==4 && right==2
                    ρ = RDM_PH(G,[[3,3],[3,4],[3,2+4*dis],[4,2+4*dis]])
                elseif left==4 && right==3
                    ρ = RDM_PH(G,[[3,3],[3,4],[1,3+4*dis],[2,3+4*dis]])
                elseif left==4 && right==4
                    ρ = RDM_PH(G,[[3,3],[3,4],[3,3+4*dis],[3,4+4*dis]])
                end
                rotate_ρ = v'*ρ*v
                measure_res = []
                for _ in 1:parsed_args["samples"]
                    ind,_ = ρmeasure!(rotate_ρ)
                    push!(measure_res,u[ind])
                end
                sample_res["$left$right"] = mean(measure_res)
                mean_value += signs[left]*signs[right]*mean(measure_res)
                true_mean_value += signs[left]*signs[right]*tr(ρ*obs)
                push!(stds,std(measure_res))
            end
        end
        total_std = sqrt(sum([x^2 for x in stds]))
        sample_res["mean_$bt"] = mean_value
        sample_res["true_mean_$bt"] = real(true_mean_value)
        sample_res["std_$bt"] = total_std
        if bt==1
            println("mean value: $mean_value")
            println("total std: $total_std")
            println("true mean value: $true_mean_value")
        end
    end
    file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/SNR_data/L$(parsed_args["L"])_type$(parsed_args["type"])_dis$(parsed_args["distance"])_Delta$(parsed_args["Delta"])_mu$(parsed_args["mu"])_samples$(parsed_args["samples"])_bootstrap$(parsed_args["bootstrap"]).json"
    open(file_name, "w") do io
        JSON3.write(io, sample_res)
    end
end
main()
