using LinearAlgebra
using SparseArrays
using Fermionic
using Arpack
using Optimization
using Optimization, OptimizationNOMAD
using ProgressMeter
using OptimizationOptimJL
using Optim
include("Hamiltonian.jl")
"""
p = [Int] is the number of rotation sequences
"""
function joint_constraint_objective(x)
    p=[1]
    target1 = target_op(1)+target_op(1)';
    target2 = target_op(2)+target_op(2)';
    @assert length(x)==3*p[1]
    for l in 1:p[1]
        Ul_particle_sector1 = Uevolve_uniform(x[(l-1)*3+1],x[(l-1)*3+2],1.0,x[(l-1)*3+3],1);
        Ul_particle_sector2 = Uevolve_uniform(x[(l-1)*3+1],x[(l-1)*3+2],1.0,x[(l-1)*3+3],2);
        target1 = Ul_particle_sector1*target1*Ul_particle_sector1';
        target2 = Ul_particle_sector2*target2*Ul_particle_sector2';
    end
    residual1 = target1 - Diagonal(diag(target1));
    residual2 = target2 - Diagonal(diag(target2));
    return sum(abs.(residual1))+sum(abs.(residual2))
end
function repeat_array(arr, L)
    return vcat([arr for _ in 1:L]...)
end
# sequence-1
seq_length = 1
constraint_x0 = rand(3*seq_length)
for l in 1:seq_length
    constraint_x0[(l-1)*2+1] += constraint_x0[(l-1)*2+1]+1
    constraint_x0[(l-1)*3+2] += constraint_x0[(l-1)*3+2]*0.01+0.1
end
p = [seq_length]
sol = optimize(joint_constraint_objective, [1.0, 0.1,0.0], [15.0,1.0,2*π], constraint_x0, )
sol.minimum
sol.lower_bound


sol.minimizer
for field in fieldnames(typeof(sol))
    println("$(field): $(getfield(sol, field))")
end
for field in fieldnames(typeof(sol))
    println("$(field)")
end