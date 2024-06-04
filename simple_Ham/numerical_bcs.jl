using Fermionic
using PythonCall
using LinearAlgebra
include("./analytical_bcs.jl")

########### Compare with numerical result on <c*_{0,0,up}c*_{xj,yj,down}c_{xk,yk,up}c_{xl,yl,down}>
function numerical_four_fermion(xj::Int64,yj::Int64,xk::Int64,yk::Int64,xl::Int64,yl::Int64, L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    if type == :s
        G = BCS_G(L, :s;μ=μ,Δ=Δ);
    elseif type == :d
        G = BCS_G(L, :d;μ=μ,Δ=Δ);
    else
        error("type not implemented")
    end
    # i_site = Cartesian2Index([(L+1)÷2,(L+1)÷2],[L,L],1);
    # j_site = Cartesian2Index([xj+(L+1)÷2,yj+(L+1)÷2],[L,L],2);
    # k_site = Cartesian2Index([xk+(L+1)÷2,yk+(L+1)÷2],[L,L],1);
    # l_site = Cartesian2Index([xl+(L+1)÷2,yl+(L+1)÷2],[L,L],2);
    i_site = Cartesian2Index([1,1],[L,L],1);
    j_site = Cartesian2Index([mod(xj,L)+1,mod(yj,L)+1],[L,L],2);
    k_site = Cartesian2Index([mod(xk,L)+1,mod(yk,L)+1],[L,L],1);
    l_site = Cartesian2Index([mod(xl,L)+1,mod(yl,L)+1],[L,L],2);
    unsort = [i_site,j_site,k_site,l_site];
    Is = sort(unsort);
	ρs = RDM(G, Is);
    sort_indices = sortperm(unsort);
    o = Op(4);
    obs = ad(o,sort_indices[1])*ad(o,sort_indices[2])*a(o,sort_indices[3])*a(o,sort_indices[4]);
    return tr(ρs*obs)
end
