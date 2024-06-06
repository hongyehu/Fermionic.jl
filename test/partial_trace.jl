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
    Projector = exp(Matrix(generator))
    vaccum = zeros(ComplexF64,2^(2*L^2))
    vaccum[1] = 1.0
    groundstate = Projector*vaccum
    return G, groundstate./sqrt(groundstate'*groundstate)
end
function generate_full_wavefunction(L::Int64, G)
    o = Op(2*L^2)
    generator = zeros(ComplexF64,2^(2*L^2),2^(2*L^2))
    for i in 1:2*L^2
        for j in 1:2*L^2
            generator += G[i,j].*Matrix(ad(o,i)*ad(o,j))
        end
    end
    Projector = exp(Matrix(generator))
    vaccum = zeros(ComplexF64,2^(2*L^2))
    vaccum[1] = 1.0
    groundstate = Projector*vaccum
    return groundstate./sqrt(groundstate'*groundstate)
end
# Define the necessary functions
function two_site_test_adad(L, site1, site2, type)
    μ = rand()
    Δ = rand() * 5
    G, gs = generate_full_BCS_wavefunction(L, type; μ=μ, Δ=Δ)
    o = Op(2 * L^2)
    test_ob = ad(o, site1) * ad(o, site2)
    num_exp = gs' * test_ob * gs
    ρ = RDM(G, [site1, site2])
    ro = Op(2)
    return round(tr(ρ * ad(ro, 1) * ad(ro, 2)),digits=7) ≈ round(num_exp,digits=7)
end
function randG_two_site_test_adad(G, L, site1, site2)
    gs = generate_full_wavefunction(L, G)
    o = Op(2 * L^2)
    test_ob = ad(o, site1) * ad(o, site2)
    num_exp = gs' * test_ob * gs
    ρ = RDM(G, [site1, site2])
    ro = Op(2)
    println("expected: ",num_exp)
    return round(tr(ρ * ad(ro, 1) * ad(ro, 2)),digits=7) ≈ round(num_exp,digits=7)
end
function two_site_test_ada(L, site1, site2, type)
    μ = rand()
    Δ = rand() * 5
    G, gs = generate_full_BCS_wavefunction(L, type; μ=μ, Δ=Δ)
    o = Op(2 * L^2)
    test_ob = ad(o, site1) * a(o, site2)
    num_exp = gs' * test_ob * gs
    ρ = RDM(G, [site1, site2])
    ro = Op(2)
    return round(tr(ρ * ad(ro, 1) * a(ro, 2)),digits=7) ≈ round(num_exp,digits=7)
end
function randG_two_site_test_ada(G, L,site1,site2)
    gs = generate_full_wavefunction(L, G)
    o = Op(2 * L^2)
    test_ob = ad(o, site1) * a(o, site2)
    num_exp = gs' * test_ob * gs
    ρ = RDM(G, [site1, site2])
    ro = Op(2)
    println("expected: ",num_exp)
    return round(tr(ρ * ad(ro, 1) * a(ro, 2)),digits=7) ≈ round(num_exp,digits=7)
end
function four_site_test_adadaa(L,site1,site2,site3,site4,type)
    μ = rand()
    Δ = rand() * 5
    G, gs = generate_full_BCS_wavefunction(L, type; μ=μ, Δ=Δ)
    o = Op(2 * L^2)
    test_ob = ad(o, site1) * ad(o, site2) * a(o, site3) * a(o, site4)
    num_exp = gs' * test_ob * gs
    ρ = RDM(G, [site1, site2, site3, site4])
    ro = Op(4)
    println("expected: ",num_exp)
    return round(tr(ρ * ad(ro, 1) * ad(ro, 2)*a(ro,3)*a(ro,4)),digits=7) ≈ round(num_exp,digits=7)
end
function four_site_test_adadadad(L,site1,site2,site3,site4,type)
    μ = rand()
    Δ = rand() * 5
    G, gs = generate_full_BCS_wavefunction(L, type; μ=μ, Δ=Δ)
    o = Op(2 * L^2)
    test_ob = ad(o, site1) * ad(o, site2) * ad(o, site3) * ad(o, site4)
    num_exp = gs' * test_ob * gs
    ρ = RDM(G, [site1, site2, site3, site4])
    ro = Op(4)
    println("expected: ",num_exp)
    return round(tr(ρ * ad(ro, 1) * ad(ro, 2)*ad(ro,3)*ad(ro,4)),digits=7) ≈ round(num_exp,digits=7)
end

## Test two site partial trace
@assert two_site_test_adad(2, 3, 4, :s)
@assert two_site_test_ada(2, 3, 4, :d)

@assert four_site_test_adadaa(2, 1,4,5, 8, :s)
@assert four_site_test_adadaa(2, 1,4,5, 8, :d)

@assert four_site_test_adadadad(2, 1,4,7, 8, :s)
@assert four_site_test_adadadad(2, 1,4,5, 8, :d)

function generate_antisymmetric_matrix(L::Int)
    # Create a random LxL matrix
    A = rand(L, L)
    # Make the matrix anti-symmetric
    A_antisymmetric = A - transpose(A)
    return A_antisymmetric
end
G = Complex{Float64}.(generate_antisymmetric_matrix(8))
randG_two_site_test_adad(G, 2, 1, 2)
randG_two_site_test_ada(G, 2, 1, 2)

L = 2
μ = rand()
Δ = rand()*5
site1 = 1
site2 = 3
type = :s
G, gs = generate_full_BCS_wavefunction(L,type;μ=μ,Δ=Δ);
o = Op(2*L^2)
test_ob = ad(o,site1)*a(o,site2)
num_exp = round(gs'*test_ob*gs,digits=5)
ρ=RDM(G,[site1,site2])
ro = Op(2)
tr(ρ*ad(ro,1)*a(ro,2))
println(tr(ρ*ad(ro,1)*ad(ro,2))≈num_exp)

site1 = 3
site2 = 1
type = :s
G, gs = generate_full_BCS_wavefunction(L,type;μ=μ,Δ=Δ);
gs
gs'
o = Op(2*L^2)
test_ob = Matrix(ad(o,3)*a(o,1))
num_exp = round(gs'*test_ob*gs,digits=5)

ρ=RDM(G,[site1,site2])
ro = Op(2)
tr(ρ*ad(ro,1)*a(ro,2))
println(tr(ρ*ad(ro,1)*ad(ro,2))≈num_exp)

test_o  = Op(2)
basis(test_o)
ad(test_o,2)

begin
    plt.imshow(Matrix(a(o,3)*ad(o,1)))
    plt.show()
end


### zero momentum
L = 2
μ = rand()
Δ = rand()
type = :s
G, gs = generate_full_BCS_wavefunction(L,type;μ=μ,Δ=Δ);
o = Op(2*L^2)
test_ob = ad(o,1)*ad(o,6)
num_exp = round(gs'*test_ob*gs,digits=5)
analytical_two_fermion(0,1,L,type,μ=μ,Δ=Δ)
ρ=RDM(2*G,[1,6])
tr(ρ*ad(Op(2),1)*ad(Op(2),2))

