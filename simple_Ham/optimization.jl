using LinearAlgebra
using SparseArrays
using Fermionic
using Arpack
using Optimization
using Optimization, OptimizationNOMAD
using ProgressMeter
include("Hamiltonian.jl")

# function rosenbrock(u, p)
#     return (p[1] - u[1])^2 + p[2] * (u[2] - u[1]^2)^2
# end
# u0 = zeros(2)
# p = [1.0, 100.0]
# prob = OptimizationProblem(rosenbrock, u0, p, lb = [-1.0, -1.0], ub = [1.0, 1.0])
# sol = Optimization.solve(prob, NOMADOpt())
# print(sol)

"""
x = [J,U,t,δ1up,δ1down,δ2up,δ2down,θ]
"""
function objective(x,p)
    target = target_op(1)+target_op(1)';
    U = Uevolve(
    x[1],x[2],1.0,x[3],x[4],x[5],x[6],x[7],1
    );
    transformed = U*target*U';
    residual = transformed - Diagonal(diag(transformed));
    return sum(abs.(residual))
end
function transformed_obs(x)
    target = target_op(1)+target_op(1)';
    U = Uevolve(
    x[1],x[2],1.0,x[3],x[4],x[5],x[6],x[7],1
    );
    transformed = U*target*U';
    return transformed
end
x = [1.0,2.0,1.0,3.9,4.0,5.0,6.0]
objective(x,0)
x0 = rand(7)
prob = OptimizationProblem(objective, x0)
sol = Optimization.solve(prob, NOMADOpt())
sol.minimum
sol.u
transformed_obs(sol.u)

x0[2] = 1.2
x0
bounded_prob = OptimizationProblem(objective, x0,lb = [0.0, 1.0,-10.0,-10.0,-10.0,-10.0,0.0], ub = [10.0, 2.0,10.0,10.0,10.0,10.0,π*2])
bounded_sol = Optimization.solve(bounded_prob, NOMADOpt())