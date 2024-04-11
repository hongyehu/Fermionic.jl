using LinearAlgebra
using SparseArrays
using Fermionic
using Arpack

function Hopping(particle_sector::Int64)
    if particle_sector == 1
        Hop = ComplexF64[0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0; 1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0];
        return Hop
    elseif particle_sector == 2
        Hop = ComplexF64[0 0 -1 1 0 0;0 0 0 0 0 0;-1 0 0 0 0 -1;1 0 0 0 0 1;0 0 0 0 0 0; 0 0 -1 1 0 0];
        return Hop
    elseif particle_sector == 3
        Hop = ComplexF64[0 0 -1 0;0 0 0 -1;-1 0 0 0;0 -1 0 0];
        return Hop
    else
        error("particle_sector not implemented")
    end
end
function Int(particle_sector::Int64)
    if particle_sector == 1
        return 0
    elseif particle_sector == 2
        H = ComplexF64[1 0 0 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 1];
        return H
    elseif particle_sector == 3
        H = ComplexF64[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
        return H
    else
        error("particle_sector not implemented")
    end
end
function n(particle_sector::Int64;Site::Int64=1,Spin::Int64=1)
    #Spin = 1 for up, Spin = 2 for down
    if particle_sector == 1
        if Site==1&&Spin==1
            return Diagonal(ComplexF64[0,0,0,1])
        elseif Site==1&&Spin==2
            return Diagonal(ComplexF64[0,0,1,0])
        elseif Site==2&&Spin==1
            return Diagonal(ComplexF64[0,1,0,0])
        elseif Site==2&&Spin==2
            return Diagonal(ComplexF64[1,0,0,0])
        else
            error("Site or Spin not implemented (only 2 sites and 2 spins)")
        end
    elseif particle_sector == 2
        if Site==1&&Spin==1
            return Diagonal(ComplexF64[0,0,0,1,1,1])
        elseif Site==1&&Spin==2
            return Diagonal(ComplexF64[0,1,1,0,0,1])
        elseif Site==2&&Spin==1
            return Diagonal(ComplexF64[1,0,1,0,1,0])
        elseif Site==2&&Spin==2
            return Diagonal(ComplexF64[1,1,0,1,0,0])
        else
            error("Site or Spin not implemented (only 2 sites and 2 spins)")
        end
    elseif particle_sector == 3
        if Site==1&&Spin==1
            return Diagonal(ComplexF64[0,1,1,1])
        elseif Site==1&&Spin==2
            return Diagonal(ComplexF64[1,0,1,1])
        elseif Site==2&&Spin==1
            return Diagonal(ComplexF64[1,1,0,1])
        elseif Site==2&&Spin==2
            return Diagonal(ComplexF64[1,1,1,0])
        else
            error("Site or Spin not implemented (only 2 sites and 2 spins)")
        end
    else
        error("particle_sector not implemented")
    end
end
function Xrot(particle_sector::Int64)
    if particle_sector == 1
        return ComplexF64[0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0]
    elseif particle_sector == 2
        return ComplexF64[0 0 0 0 0 0;0 0 1 1 0 0;0 1 0 0 1 0;0 1 0 0 1 0;0 0 1 1 0 0;0 0 0 0 0 0]
    elseif particle_sector == 3
        return ComplexF64[0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0]
    else
        error("particle_sector not implemented")
    end
end
function ShadowMap(
    J::Float64,
    U::Float64,
    t::Float64,
    δ1up::Float64,
    δ1down::Float64,
    δ2up::Float64,
    δ2down::Float64,
    θ::Float64,
    particle_sector::Int64,
    )
# Attractive Hubbard model
    subspace_dim = binomial(4,particle_sector)
    ones_vector = ones(ComplexF64, subspace_dim)
    sparse_identity_matrix = spdiagm(0 => ones_vector)
    results = [spzeros(ComplexF64,subspace_dim,subspace_dim) for _ = 1:subspace_dim]
    Ham = -J*Hopping(particle_sector).-U*Int(particle_sector).+δ1up.*n(particle_sector,Site=1,Spin=1).+δ1down.*n(particle_sector,Site=1,Spin=2).+δ2up.*n(particle_sector,Site=2,Spin=1).+δ2down.*n(particle_sector,Site=2,Spin=2)
    sparse_identity_matrix = exp(Matrix(1im*t*Ham))*sparse_identity_matrix
    Urot = exp(Matrix(1im*θ*Xrot(particle_sector)))
    sparse_identity_matrix = Urot*sparse_identity_matrix
    for k in 1:subspace_dim
        results[k] = sparse_identity_matrix[:,k]*sparse_identity_matrix[:,k]'
    end
    shadow_map = zeros(ComplexF64,subspace_dim^2,subspace_dim^2)
    for j in 1:subspace_dim
        flat_basis = reshape(results[j],subspace_dim^2)
        shadow_map += flat_basis*flat_basis'
    end
    return shadow_map
end
function Uevolve(
    J::Float64,
    U::Float64,
    t::Float64,
    δ1up::Float64,
    δ1down::Float64,
    δ2up::Float64,
    δ2down::Float64,
    θ::Float64,
    particle_sector::Int64,
    )
    Ham = -J*Hopping(particle_sector).-U*Int(particle_sector).+δ1up.*n(particle_sector,Site=1,Spin=1).+δ1down.*n(particle_sector,Site=1,Spin=2).+δ2up.*n(particle_sector,Site=2,Spin=1).+δ2down.*n(particle_sector,Site=2,Spin=2)
    U = exp(Matrix(-1im*t*Ham))
    Urot = exp(Matrix(-1im*θ*Xrot(particle_sector)))
    return Urot*U
end
function target_op(particle_sector::Int64)
    o = Op_fixed(2*2,particle_sector); 
    return ada(o,3,2)-ada(o,1,4)
end

U = Uevolve(1.0,1.0,1.0,0.0,0.0,0.0,0.4,0.0,1)
U*target_op(1)*U'

