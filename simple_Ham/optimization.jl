using LinearAlgebra
using SparseArrays
using Fermionic
using Arpack
using Optimization
using Optimization, OptimizationNOMAD
using ProgressMeter
include("Hamiltonian.jl")


"""
x = [J,U,t,δ1up,δ1down,δ2up,δ2down,θ]
"""
function objective(x,p)
    target = target_op(p[1])+target_op(p[1])';
    U = Uevolve(
    x[1],x[2],1.0,x[3],x[4],x[5],x[6],x[7],p[1]
    );
    transformed = U*target*U';
    residual = transformed - Diagonal(diag(transformed));
    return sum(abs.(residual))
end
function transformed_obs(x,p)
    target = target_op(p[1])+target_op(p[1])';
    U = Uevolve(
    x[1],x[2],1.0,x[3],x[4],x[5],x[6],x[7],p[1]
    );
    transformed = U*target*U';
    return transformed
end
"""
x = [J,U,t,δ1up,δ1down,δ2up,δ2down,θ]
"""
function constraint_objective(x,p)
    target = target_op(p[1])+target_op(p[1])';
    U = Uevolve(
    x[1],x[2],1.0,0.0,0.0,0.0,0.0,x[3],p[1]
    );
    transformed = U*target*U';
    residual = transformed - Diagonal(diag(transformed));
    return sum(abs.(residual))
end
function joint_constraint_objective(x,p)
    target1 = target_op(1)+target_op(1)';
    U1 = Uevolve(
    x[1],x[2],1.0,0.0,0.0,0.0,0.0,x[3],1
    );
    transformed1 = U1*target1*U1';
    residual1 = transformed1 - Diagonal(diag(transformed1));
    target2 = target_op(2)+target_op(2)';
    U2 = Uevolve(
    x[1],x[2],1.0,0.0,0.0,0.0,0.0,x[3],2
    );
    transformed2 = U2*target2*U2';
    residual2 = transformed2 - Diagonal(diag(transformed2));
    return sum(abs.(residual1))+sum(abs.(residual2))
end
function constrait_transformed_obs(x,p)
    target = target_op(p[1])+target_op(p[1])';
    U = Uevolve(
    x[1],x[2],1.0,0.0,0.0,0.0,0.0,x[3],p[1]
    );
    transformed = U*target*U';
    return transformed
end
function residual_constrait_transformed_obs(x,p)
    target = target_op(p[1])+target_op(p[1])';
    U = Uevolve(
    x[1],x[2],1.0,0.0,0.0,0.0,0.0,x[3],p[1]
    );
    transformed = U*target*U';
    residual = transformed - Diagonal(diag(transformed));

    return U'*residual*U
end


x0 = rand(7)
x0[1]=x0[1]+1
x0[2]=x0[2] 
p=[2]
bounded_prob = OptimizationProblem(objective, x0,p,lb = [1.0, 0.1,-10.0,-10.0,-10.0,-10.0,0.0], ub = [10.0, 5.0,10.0,10.0,10.0,10.0,π*2])
bounded_sol = Optimization.solve(bounded_prob, NOMADOpt())
bounded_sol.minimum
bounded_sol.u



constraint_x0 = rand(3)
constraint_x0[1] += constraint_x0[1]+1
constraint_x0[2] += constraint_x0[2]+0.01
# constraint_x0[3]=0.005
p = [2]
constraint_bounded_prob = OptimizationProblem(constraint_objective, constraint_x0,p,lb = [1.0, 0.1,0.0], ub = [10.0, 2.0,2*π])
constraint_bounded_sol = Optimization.solve(constraint_bounded_prob, NOMADOpt())
constraint_bounded_sol.minimum
constraint_bounded_sol.u

constrait_transformed_obs(constraint_bounded_sol.u,p)


constraint_x0 = rand(3)
constraint_x0[1] += constraint_x0[1]+2
constraint_x0[2] += constraint_x0[2]+0.1
p = []
joint_constraint_bounded_prob = OptimizationProblem(joint_constraint_objective, constraint_x0,p,lb = [1.0, 0.0,0.0], ub = [6.0, 2.0,2*π])
joint_constraint_bounded_sol = Optimization.solve(joint_constraint_bounded_prob, NOMADOpt())
joint_constraint_bounded_sol.minimum
joint_constraint_bounded_sol.u
x_test = [4.32,0.1,0.7854]
x2_test = [3.53,0.1,0.7854]
round.(constrait_transformed_obs(x2_test,[1]),digits=3)
round.(constrait_transformed_obs(x2_test,[2]),digits=3)
round.(residual_constrait_transformed_obs(x2_test,[1]),digits=3)
round.(residual_constrait_transformed_obs(x2_test,[2]),digits=3)
target_op(2)+target_op(2)'



o = Op_fixed(4,2);
basis(o)
(ada(o,1,2))
(ada(o,4,3)+ada(o,3,4))*exp(im*π*Matrix(ada(o,3,3)))+(ada(o,2,1)+ada(o,1,2))*exp(im*π*Matrix(ada(o,2,2)))


o = Op_fixed(2*2,2)
basis(o)
target_op(2)+target_op(2)'

us, vs = eigen(Matrix(target_op(2)+target_op(2)'))
us
vs
Matrix(target_op(1)+target_op(1)')*vs[:,3]≈us[3]*vs[:,3]



v1 = [0;0;0;0;1/sqrt(2);0].+[0;1/sqrt(2);0;0;0;0]
o = Matrix(target_op(2)+target_op(2)')
o*v1
vs[:,1]
2*vs[:,1]*vs[:,1]'-2*vs[:,end]*vs[:,end]'
target_op(3)+target_op(3)'
ϕ1 = 0.5*([])


