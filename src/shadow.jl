using KrylovKit
using SparseArrays
"""
=======
Input: 
    - Lattice: AbstractArray{<:AbstractFloat}
    - Params: Dict{Any,Any}. Example: Dict(1=>Dict("Jh"=>2.0,"Jv"=>1.5,"U"=>2.0,"δ"=>Dict{Any,Any}(),"t"=>0.4))
    - N: Int
=======
Output:
    - basis: AbstractArray{<:AbstractFloat}
"""
function ShadowBasis_krylov(
    Lattice::AbstractArray{<:Int},
    Params::Dict{Any,Any},
    N::Int,
    )
    """
    Row first ordering
    """
    dimension = length(Lattice)
    if dimension == 1
        error("Dimension not implemented")
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        subspace_dim = binomial(2*rows*cols,N)
        # sub_basis, _ = basis_m(2*rows*cols,N)
        Results = Dict{Any,Any}()
        for (key,_) in Params
            # print("key: ", key,"\n")
            tmp = Dict{Any,Any}()
            Ham = SpinFullFermiHubbardSubspace([rows,cols],Hopping([rows,cols],Params[key]["Jh"],Params[key]["Jv"]),Params[key]["U"],Params[key]["δ"],N)
            ones_vector = ones(ComplexF64, binomial(2*rows*cols,N))
            sparse_identity_matrix = spdiagm(0 => ones_vector)
            # batch_evolve_state = sparse(exp(Matrix(1im*Ham*Params[key]["t"]))*sparse_identity_matrix)
            # for k in 1:subspace_dim
            #     print("k: ", k,"\n")
            #     tmp[k] = batch_evolve_state[:,k]*batch_evolve_state[:,k]'
            # end
            for k in 1:subspace_dim
                # print("k: ", k,"\n")
                # sparse_identity_matrix[:,k] = exp(Matrix(1im*Ham*Params[key]["t"]))*sparse_identity_matrix[:,k]
                dt::ComplexF64 = ComplexF64(Params[key]["t"]*1im)
                sparse_identity_matrix[:,k], _ = exponentiate(Ham, dt, (@view sparse_identity_matrix[:,k]); ishermitian = true,eager=true)
                tmp[k] = sparse_identity_matrix[:,k]*sparse_identity_matrix[:,k]'
            end
            Results[key] = tmp
        end
        return Results
    else
        error("Dimension not implemented")
    end
    return basis
end

"""
=======
Input: 
    - Lattice: AbstractArray{<:AbstractFloat}
    - Params: Dict{Any,Any}. Example: Dict(1=>Dict("Jh"=>2.0,"Jv"=>1.5,"U"=>2.0,"δ"=>Dict{Any,Any}(),"t"=>0.4))
    - N: Int
=======
Output:
    - basis: AbstractArray{<:AbstractFloat}
"""
function ShadowBasis_krylov(
    Lattice::AbstractArray{<:Int},
    Jh::Float64,
    Jv::Float64,
    U::Float64,
    t::Float64,
    δ::Dict{Tuple{Int, Int}, Float64},
    N::Int,
    )
    """
    Row first ordering
    """
    dimension = length(Lattice)
    if dimension == 1
        error("Dimension not implemented")
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        subspace_dim = binomial(2*rows*cols,N)
        Ham = SpinFullFermiHubbardSubspace([rows,cols],Hopping([rows,cols],Jh,Jv),U,δ,N)
        ones_vector = ones(ComplexF64, binomial(2*rows*cols,N))
        sparse_identity_matrix = spdiagm(0 => ones_vector)
        results = [spzeros(ComplexF64,subspace_dim,subspace_dim) for _ = 1:subspace_dim]
        for k in 1:subspace_dim
            # print("k: ", k,"\n")
            # sparse_identity_matrix[:,k] = exp(Matrix(1im*Ham*Params[key]["t"]))*sparse_identity_matrix[:,k]
            dt = t*im
            sparse_identity_matrix[:,k], _ = exponentiate(Ham, dt, (@view sparse_identity_matrix[:,k]); ishermitian = true,eager=true)
            results[k] = sparse_identity_matrix[:,k]*sparse_identity_matrix[:,k]'
        end
        return results
    else
        error("Dimension not implemented")
    end
    return basis
end

function ShadowBasis(
    Lattice::AbstractArray{<:Int},
    Jh::Float64,
    Jv::Float64,
    U::Float64,
    t::Float64,
    δ::Dict{Tuple{Int, Int}, Float64},
    N::Int,
    )
    """
    Row first ordering
    """
    dimension = length(Lattice)
    if dimension == 1
        error("Dimension not implemented")
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        subspace_dim = binomial(2*rows*cols,N)
        # print("key: ", key,"\n")
        Ham = SpinFullFermiHubbardSubspace([rows,cols],Hopping([rows,cols],Jh,Jv),U,δ,N)
        ones_vector = ones(ComplexF64, binomial(2*rows*cols,N))
        sparse_identity_matrix = spdiagm(0 => ones_vector)
        results = [spzeros(ComplexF64,subspace_dim,subspace_dim) for _ = 1:subspace_dim]
        sparse_identity_matrix = exp(Matrix(im*t*Ham))*sparse_identity_matrix
        for k in 1:subspace_dim
            # print("k: ", k,"\n")
            results[k] = sparse_identity_matrix[:,k]*sparse_identity_matrix[:,k]'
        end
        return results
    else
        error("Dimension not implemented")
    end
    return basis
end