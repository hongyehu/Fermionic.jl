using LinearAlgebra
using SparseArrays
using Fermionic
using Arpack
using Optimization
using Optimization, OptimizationNOMAD
using ProgressMeter
using OptimizationOptimJL
using Optim
using QuadDIRECT
include("Hamiltonian.jl")
"""
p = [Int] is the number of rotation sequences
"""
function joint_constraint_objective(x,p)
    target1 = target_op(1)+target_op(1)';
    target2 = target_op(2)+target_op(2)';
    @assert length(x)==2*p[1]
    for l in 1:p[1]
        Ul_particle_sector1 = Uevolve_uniform(x[(l-1)*2+1],0.1,1.0,x[(l-1)*2+2],1);
        Ul_particle_sector2 = Uevolve_uniform(x[(l-1)*2+1],0.1,1.0,x[(l-1)*2+2],2);
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
constraint_x0 = rand(2*seq_length)
for l in 1:seq_length
    constraint_x0[(l-1)*2+1] += constraint_x0[(l-1)*2+1]+1
    # constraint_x0[(l-1)*3+2] += constraint_x0[(l-1)*3+2]*0.01+0.1
end
p = [seq_length]
constraint_bounded_prob = OptimizationProblem(joint_constraint_objective, constraint_x0,p,lb = repeat_array([1.0, 0.0],seq_length), ub = repeat_array([10.0,2*π],seq_length))
constraint_bounded_sol = Optimization.solve(constraint_bounded_prob, QuadDirect(),splits = ([1., 2.0, 4.0], [1.2,π, 0.8]))
constraint_bounded_sol.minimum
constraint_bounded_sol.u

# sequence-5
seq_length = 3
# constraint_x0 = rand(3*seq_length)
# for l in 1:seq_length
#     constraint_x0[(l-1)*3+1] += constraint_x0[(l-1)*3+1]*0.5+1
#     constraint_x0[(l-1)*3+2] += constraint_x0[(l-1)*3+2]*0.01+0.1
# end

constraint_x0 = [1.96349536863696,0.78539816449691,rand(),rand(),rand(),rand()]
p = [seq_length]
constraint_bounded_prob = OptimizationProblem(joint_constraint_objective, constraint_x0,p,lb = repeat_array([0.0, 0.0],seq_length), ub = repeat_array([10.0, 2*π],seq_length))
constraint_bounded_sol = Optimization.solve(constraint_bounded_prob, NOMADOpt())
constraint_bounded_sol.minimum
constraint_bounded_sol.u


rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
p = [1.0, 100.0]
f = OptimizationFunction(rosenbrock)
prob = Optimization.OptimizationProblem(f, x0, p, lb = [-1.0, -1.0], ub = [1.0, 1.0])
solve(prob, QuadDIRECT(), splits = ([-0.9, 0, 0.9], [-0.8, 0, 0.8]))


function camel(x)
    # 6-hump camel function. Typically evaluated over [-3,3] × [-2,2].
    x1, x2 = x[1], x[2]
    x1s = x1*x1
    x2s = x2*x2
    return (4 - 2.1*x1s + x1s*x1s/3)*x1s + x1*x2 + (-4 + 4*x2s)*x2s
end
lower, upper = [-3,-2], [3,2]   
splits = ([-2, 0, 2], [-1, 0, 1])   # 
root, x0 = analyze(camel, splits, lower, upper)
x0
minimum(root)