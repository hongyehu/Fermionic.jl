using Fermionic
using PythonCall
using LinearAlgebra
plt = pyimport("matplotlib.pyplot")
include("../simple_Ham/analytical_bcs.jl")

function generate_full_BCS_wavefunction(L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    o = Op(2*L^2)
    if type == :s
        G = BCS_G(L, :s;μ=μ,Δ=Δ);
    elseif type == :d
        G = BCS_G(L, :d;μ=μ,Δ=Δ);
    else
        error("type not implemented")
    end
    generator = zeros(ComplexF64,2^(2*L^2),2^(2*L^2))
    for i in 1:2*L^2
        for j in 1:2*L^2
            generator += G[i,j].*Matrix(ad(o,i)*ad(o,j))
        end
    end
    Projector = exp(1.0.*Matrix(generator))
    vaccum = zeros(ComplexF64,2^(2*L^2))
    vaccum[1] = 1.0
    groundstate = Projector*vaccum
    return groundstate./sqrt(groundstate'*groundstate)
end
gs = generate_full_BCS_wavefunction(2,:d;μ=0.5,Δ=5.0)
o = Op(8)
test_ob = ad(o,1)*ad(o,6)
num_exp = gs'*test_ob*gs
G = BCS_G(2, :d ;μ=0.5,Δ=5.0)
ρ=RDM(2*G,[1,6])
ro = Op(2)
tr(ρ*ad(ro,1)*ad(ro,2))

analytical_two_fermion(0,1,2,:d)