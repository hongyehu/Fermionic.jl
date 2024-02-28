function Cartesian2Index(C_coor::AbstractArray{<:Int64},Lattice::AbstractArray{<:Int64},Spin::Int)
    """
    Row first ordering. 
    Spin = 1 for up, Spin = 2 for down
    """
    dimension = length(Lattice)
    if dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        tmp = rows*(C_coor[2]-1)+C_coor[1]
        if Spin == 1
            return tmp*2-1
        elseif Spin == 2
            return tmp*2
        else
            error("Spin number should be 1(up) or 2(down)")
        end
    else
        error("Dimension not implemented")
    end
end

# Implement Dwave pairing 
function Dwave(
    Lattice::AbstractArray{<:Int64},
    i::Int,
    j::Int,
    N::Int,
)
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        op_size = 2*rows*cols
        o = Op_fixed(op_size,N); # fixed particle number subspace
        dwave = ada(o,Cartesian2Index([1,i],[rows,cols],1),Cartesian2Index([1,j],[rows,cols],1))*ada(o,Cartesian2Index([2,i],[rows,cols],2),Cartesian2Index([2,j],[rows,cols],2))
        dwave -= ada(o,Cartesian2Index([1,i],[rows,cols],1),Cartesian2Index([1,j],[rows,cols],2))*ada(o,Cartesian2Index([2,i],[rows,cols],2),Cartesian2Index([2,j],[rows,cols],1))
        dwave -= ada(o,Cartesian2Index([1,i],[rows,cols],2),Cartesian2Index([1,j],[rows,cols],1))*ada(o,Cartesian2Index([2,i],[rows,cols],1),Cartesian2Index([2,j],[rows,cols],2))
        dwave += ada(o,Cartesian2Index([1,i],[rows,cols],2),Cartesian2Index([1,j],[rows,cols],2))*ada(o,Cartesian2Index([2,i],[rows,cols],1),Cartesian2Index([2,j],[rows,cols],1))
        return dwave
    else
        error("Dimension not implemented")
    end
end
function transformed_Dwave(
    Lattice::AbstractArray{<:Int64},
    i::Int,
    j::Int,
    N::Int,
)
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        op_size = 2*rows*cols
        o = Op_fixed(op_size,N); # fixed particle number subspace
        if i % 2 == 1 && j % 2 == 1 # both i and j are odd numbers
            dwave = (ada(o,Cartesian2Index([2,i],[rows,cols],1),Cartesian2Index([1,i],[rows,cols],2))-ada(o,Cartesian2Index([1,i],[rows,cols],1),Cartesian2Index([2,i],[rows,cols],2)))*
            (ada(o,Cartesian2Index([1,j],[rows,cols],2),Cartesian2Index([2,j],[rows,cols],1))-ada(o,Cartesian2Index([2,j],[rows,cols],2),Cartesian2Index([1,j],[rows,cols],1)))
            return dwave
        elseif i % 2 == 1 && j % 2 == 0
            dwave = (ada(o,Cartesian2Index([2,i],[rows,cols],1),Cartesian2Index([1,i],[rows,cols],2))-ada(o,Cartesian2Index([1,i],[rows,cols],1),Cartesian2Index([2,i],[rows,cols],2)))*
            (ada(o,Cartesian2Index([2,j],[rows,cols],2),Cartesian2Index([1,j],[rows,cols],1))-ada(o,Cartesian2Index([1,j],[rows,cols],2),Cartesian2Index([2,j],[rows,cols],1)))
        else
            error("Not implemented")
        end
    else
        error("Dimension not implemented")
    end 
end

# Implement Spinfull Fermi-Hubbard model 
"""
`FermiHubbard(Dimension::Int, t::Float64, U::Float64, δ::AbstractArray{<:AbstractFloat})`

# Arguments
- `Dimension::Int`: Description of what this parameter represents.
- `t::Float64`: Description of what this parameter represents.
- `U::Float64`: Description of what this parameter represents.
- `δ::AbstractArray{<:AbstractFloat}`: Array of floating-point numbers.

# Returns
(Describe what the function returns, if anything)
"""
function SpinFullFermiHubbard(
    Lattice::AbstractArray{<:Int64}, 
    J::Dict{Any,Any}, 
    U::Float64, 
    δ::Dict{Any,Any}
    )
    # function body
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        op_size = 2*rows*cols
        # Repulsion
        Ham = U*sum([ada(Op(op_size),i,i)*ada(Op(op_size),i+1,i+1) for i in 1:2:2*rows*cols-1])
        # Kinetic term
        for i in 1:rows-1
            for j in 1:cols-1
                # Horizontal
                Ham += -J[((i,j),(i,j+1))]*(
                    ad(Op(op_size),Cartesian2Index([i,j],[rows,cols],1))*a(Op(op_size),Cartesian2Index([i,j+1],[rows,cols],1))+
                    ad(Op(op_size),Cartesian2Index([i,j+1],[rows,cols],1))*a(Op(op_size),Cartesian2Index([i,j],[rows,cols],1))+
                    ad(Op(op_size),Cartesian2Index([i,j],[rows,cols],2))*a(Op(op_size),Cartesian2Index([i,j+1],[rows,cols],2))+
                    ad(Op(op_size),Cartesian2Index([i,j+1],[rows,cols],2))*a(Op(op_size),Cartesian2Index([i,j],[rows,cols],2))
                )
                # Vertical
                Ham += -J[((i,j),(i+1,j))]*(
                    ad(Op(op_size),Cartesian2Index([i,j],[rows,cols],1))*a(Op(op_size),Cartesian2Index([i+1,j],[rows,cols],1))+
                    ad(Op(op_size),Cartesian2Index([i+1,j],[rows,cols],1))*a(Op(op_size),Cartesian2Index([i,j],[rows,cols],1))+
                    ad(Op(op_size),Cartesian2Index([i,j],[rows,cols],2))*a(Op(op_size),Cartesian2Index([i+1,j],[rows,cols],2))+
                    ad(Op(op_size),Cartesian2Index([i+1,j],[rows,cols],2))*a(Op(op_size),Cartesian2Index([i,j],[rows,cols],2))
                )
            end
        end
        # Horizontal last line  
        for i in 1:cols-1
            Ham += -J[((rows,i),(rows,i+1))]*(
                ad(Op(op_size),Cartesian2Index([rows,i],[rows,cols],1))*a(Op(op_size),Cartesian2Index([rows,i+1],[rows,cols],1))+
                ad(Op(op_size),Cartesian2Index([rows,i+1],[rows,cols],1))*a(Op(op_size),Cartesian2Index([rows,i],[rows,cols],1))+
                ad(Op(op_size),Cartesian2Index([rows,i],[rows,cols],2))*a(Op(op_size),Cartesian2Index([rows,i+1],[rows,cols],2))+
                ad(Op(op_size),Cartesian2Index([rows,i+1],[rows,cols],2))*a(Op(op_size),Cartesian2Index([rows,i],[rows,cols],2))
            )
        end
        # Vertical last line
        for i in 1:rows-1
            Ham += -J[((i,cols),(i+1,cols))]*(
                ad(Op(op_size),Cartesian2Index([i,cols],[rows,cols],1))*a(Op(op_size),Cartesian2Index([i+1,cols],[rows,cols],1))+
                ad(Op(op_size),Cartesian2Index([i+1,cols],[rows,cols],1))*a(Op(op_size),Cartesian2Index([i,cols],[rows,cols],1))+
                ad(Op(op_size),Cartesian2Index([i,cols],[rows,cols],2))*a(Op(op_size),Cartesian2Index([i+1,cols],[rows,cols],2))+
                ad(Op(op_size),Cartesian2Index([i+1,cols],[rows,cols],2))*a(Op(op_size),Cartesian2Index([i,cols],[rows,cols],2))
            )
        end
        # Random Potential
        for i in 1:rows
            for j in 1:cols
                Ham += δ[((i,j),(i,j))]*(
                    ad(Op(op_size),Cartesian2Index([i,j],[rows,cols],1))*a(Op(op_size),Cartesian2Index([i,j],[rows,cols],1))+
                    ad(Op(op_size),Cartesian2Index([i,j],[rows,cols],2))*a(Op(op_size),Cartesian2Index([i,j],[rows,cols],2))
                )
            end
        end
        return Ham
    else
        error("Dimension not implemented")
    end
end
function SpinFullFermiHubbardSubspace(
    Lattice::AbstractArray{<:Int64}, 
    J::Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int}},Float64}, 
    U::Float64, 
    δ::Dict{Tuple{Int,Int},Float64},
    N::Int
    )
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        op_size = 2*rows*cols
        o = Op_fixed(op_size,N); # fixed particle number subspace
        # Repulsion
        Ham = U*sum([ada(o,i,i)*ada(o,i+1,i+1) for i in 1:2:2*rows*cols-1])
        # Kinetic term
        for i in 1:rows-1
            for j in 1:cols-1
                # Horizontal
                Ham += -J[((i,j),(i,j+1))]*(
                    ada(o,Cartesian2Index([i,j],[rows,cols],1),Cartesian2Index([i,j+1],[rows,cols],1))+
                    ada(o,Cartesian2Index([i,j+1],[rows,cols],1),Cartesian2Index([i,j],[rows,cols],1))+
                    ada(o,Cartesian2Index([i,j],[rows,cols],2),Cartesian2Index([i,j+1],[rows,cols],2))+
                    ada(o,Cartesian2Index([i,j+1],[rows,cols],2),Cartesian2Index([i,j],[rows,cols],2))
                )
                # Vertical
                Ham += -J[((i,j),(i+1,j))]*(
                    ada(o,Cartesian2Index([i,j],[rows,cols],1),Cartesian2Index([i+1,j],[rows,cols],1))+
                    ada(o,Cartesian2Index([i+1,j],[rows,cols],1),Cartesian2Index([i,j],[rows,cols],1))+
                    ada(o,Cartesian2Index([i,j],[rows,cols],2),Cartesian2Index([i+1,j],[rows,cols],2))+
                    ada(o,Cartesian2Index([i+1,j],[rows,cols],2),Cartesian2Index([i,j],[rows,cols],2))
                )
            end
        end
        # Horizontal last line  
        for i in 1:cols-1
            Ham += -J[((rows,i),(rows,i+1))]*(
                ada(o,Cartesian2Index([rows,i],[rows,cols],1),Cartesian2Index([rows,i+1],[rows,cols],1))+
                ada(o,Cartesian2Index([rows,i+1],[rows,cols],1),Cartesian2Index([rows,i],[rows,cols],1))+
                ada(o,Cartesian2Index([rows,i],[rows,cols],2),Cartesian2Index([rows,i+1],[rows,cols],2))+
                ada(o,Cartesian2Index([rows,i+1],[rows,cols],2),Cartesian2Index([rows,i],[rows,cols],2))
            )
        end
        # Vertical last line
        for i in 1:rows-1
            Ham += -J[((i,cols),(i+1,cols))]*(
                ada(o,Cartesian2Index([i,cols],[rows,cols],1),Cartesian2Index([i+1,cols],[rows,cols],1))+
                ada(o,Cartesian2Index([i+1,cols],[rows,cols],1),Cartesian2Index([i,cols],[rows,cols],1))+
                ada(o,Cartesian2Index([i,cols],[rows,cols],2),Cartesian2Index([i+1,cols],[rows,cols],2))+
                ada(o,Cartesian2Index([i+1,cols],[rows,cols],2),Cartesian2Index([i,cols],[rows,cols],2))
            )
        end
        # Random Potential
        for i in 1:rows
            for j in 1:cols
                Ham += δ[(i,j)]*(
                    ada(o,Cartesian2Index([i,j],[rows,cols],1),Cartesian2Index([i,j],[rows,cols],1))+
                    ada(o,Cartesian2Index([i,j],[rows,cols],2),Cartesian2Index([i,j],[rows,cols],2))
                )
            end
        end
        return Ham
    else
        error("Dimension not implemented")
    end
end

function transformed_SpinFullFermiHubbardSubspace(
    Lattice::AbstractArray{<:Int64}, 
    J::Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int}},Float64}, 
    U::Float64, 
    δ::Dict{Tuple{Int,Int},Float64},
    N::Int
    )
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        op_size = 2*rows*cols
        o = Op_fixed(op_size,N); # fixed particle number subspace
        # Repulsion
        Ham = -U*sum([ada(o,i,i)*ada(o,i+1,i+1) for i in 1:2:2*rows*cols-1])
        for i in 1:rows
            for j in 1:cols
                Ham += U*ada(o,Cartesian2Index([i,j],[rows,cols],1),Cartesian2Index([i,j],[rows,cols],1))
            end
        end
        for i in 1:rows
            for j in 1:cols
                Ham += δ[(i,j)]*(
                    ada(o,Cartesian2Index([i,j],[rows,cols],1),Cartesian2Index([i,j],[rows,cols],1))-
                    ada(o,Cartesian2Index([i,j],[rows,cols],2),Cartesian2Index([i,j],[rows,cols],2))
                )
            end
        end
        # Kinetic term
        for i in 1:rows-1
            for j in 1:cols-1
                # Horizontal
                Ham += -J[((i,j),(i,j+1))]*(
                    ada(o,Cartesian2Index([i,j],[rows,cols],1),Cartesian2Index([i,j+1],[rows,cols],1))+
                    ada(o,Cartesian2Index([i,j+1],[rows,cols],1),Cartesian2Index([i,j],[rows,cols],1))+
                    ada(o,Cartesian2Index([i,j],[rows,cols],2),Cartesian2Index([i,j+1],[rows,cols],2))+
                    ada(o,Cartesian2Index([i,j+1],[rows,cols],2),Cartesian2Index([i,j],[rows,cols],2))
                )
                # Vertical
                Ham += -J[((i,j),(i+1,j))]*(
                    ada(o,Cartesian2Index([i,j],[rows,cols],1),Cartesian2Index([i+1,j],[rows,cols],1))+
                    ada(o,Cartesian2Index([i+1,j],[rows,cols],1),Cartesian2Index([i,j],[rows,cols],1))+
                    ada(o,Cartesian2Index([i,j],[rows,cols],2),Cartesian2Index([i+1,j],[rows,cols],2))+
                    ada(o,Cartesian2Index([i+1,j],[rows,cols],2),Cartesian2Index([i,j],[rows,cols],2))
                )
            end
        end
        # Horizontal last line  
        for i in 1:cols-1
            Ham += -J[((rows,i),(rows,i+1))]*(
                ada(o,Cartesian2Index([rows,i],[rows,cols],1),Cartesian2Index([rows,i+1],[rows,cols],1))+
                ada(o,Cartesian2Index([rows,i+1],[rows,cols],1),Cartesian2Index([rows,i],[rows,cols],1))+
                ada(o,Cartesian2Index([rows,i],[rows,cols],2),Cartesian2Index([rows,i+1],[rows,cols],2))+
                ada(o,Cartesian2Index([rows,i+1],[rows,cols],2),Cartesian2Index([rows,i],[rows,cols],2))
            )
        end
        # Vertical last line
        for i in 1:rows-1
            Ham += -J[((i,cols),(i+1,cols))]*(
                ada(o,Cartesian2Index([i,cols],[rows,cols],1),Cartesian2Index([i+1,cols],[rows,cols],1))+
                ada(o,Cartesian2Index([i+1,cols],[rows,cols],1),Cartesian2Index([i,cols],[rows,cols],1))+
                ada(o,Cartesian2Index([i,cols],[rows,cols],2),Cartesian2Index([i+1,cols],[rows,cols],2))+
                ada(o,Cartesian2Index([i+1,cols],[rows,cols],2),Cartesian2Index([i,cols],[rows,cols],2))
            )
        end
        return Ham
    end

end


function SpinFullRandomState(
    Lattice::AbstractArray{<:Int64}, 
    Sparse::Bool
    )
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        op_size = 2*rows*cols
        if Sparse
            state_ran_com = 2*rand(Complex{Float64},2^op_size).-1
            state_ran_com = state_ran_com/sqrt(state_ran_com'*state_ran_com)
            return sparse(state_ran_com)
        else
            state_ran_com = 2*rand(Complex{Float64},2^op_size).-1
            state_ran_com = state_ran_com/sqrt(state_ran_com'*state_ran_com)
            return state_ran_com
        end
    else
        error("Dimension not implemented")
    end

end

function SpinFullRandomState(
    Lattice::AbstractArray{<:Int64}, 
    Sparse::Bool,
    N::Int
    )
    # N is the number of particles working
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        op_size = 2*rows*cols
        if Sparse
            state_ran_com = spzeros(Complex{Float64},2^op_size)
            for i in 1:2^op_size
                if sum(basis(Op(2*rows*cols))[i,:]) == N
                    state_ran_com[i] = 2*rand(Complex{Float64},1)[1]-1
                end
            end
            state_ran_com = state_ran_com/sqrt(state_ran_com'*state_ran_com)
            return state_ran_com
        else
            state_ran_com = zeros(Complex{Float64},2^op_size)
            for i in 1:2^op_size
                if sum(basis(Op(2*rows*cols))[i,:]) == N
                    state_ran_com[i] = 2*rand(Complex{Float64},1)[1]-1
                end
            end
            state_ran_com = state_ran_com/sqrt(state_ran_com'*state_ran_com)
            return state_ran_com
        end
    else
        error("Dimension not implemented")
    end
end

function SpinFullRandomStateSubspace(
    Lattice::AbstractArray{<:Int64}, 
    Sparse::Bool,
    N::Int
    )
    # N is the number of particles working
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        op_size = 2*rows*cols
        if Sparse
            state_ran_com = 2*rand(Complex{Float64},binomial(2*rows*cols,N)).-1
            state_ran_com = state_ran_com/sqrt(state_ran_com'*state_ran_com)
            return sparse(state_ran_com)
        else
            state_ran_com = 2*rand(Complex{Float64},binomial(2*rows*cols,N)).-1
            state_ran_com = state_ran_com/sqrt(state_ran_com'*state_ran_com)
            return state_ran_com
        end
    else
        error("Dimension not implemented")
    end
end
function GHZSubspace(Lattice::AbstractArray{<:Int64}, 
    N::Int)
    # N is the number of particles working
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        state = spzeros(Complex{Float64},binomial(2*rows*cols,N))
        state[1]=1.0
        state[end]=1.0
        state = state/sqrt(state'*state)
        return state
    else
        error("Dimension not implemented")
    end
end
function measure!(state::AbstractArray{<:ComplexF64})
    prob = abs.(state).^2
    dist = Categorical(prob)
    idn = rand(dist,1)[1]
    state = spzeros(Complex{Float64},length(state))
    state[idn] = 1.0
    return idn,state
end

function evolve!(state::AbstractArray{<:ComplexF64},Ham::AbstractArray{<:ComplexF64},dt::Float64)
    state = exp(Matrix(-1im*Ham*dt))*state
    return state
end


function Hopping(Lattice::AbstractArray{<:Int}, Jh::Float64, Jv::Float64)
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        J = Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int}},Float64}()
        # Horizontal hopping
        for j in 1:cols-1
            for i in 1:rows
                J[((i,j),(i,j+1))] = Jh
            end
        end
        # Vertical hopping
        for j in 1:cols
            for i in 1:rows-1
                J[((i,j),(i+1,j))] = Jv
            end
        end
        return J
    else
        error("Dimension not implemented")
    end
end

function random_δ(Lattice::AbstractArray{<:Int}, μ::Float64, var::Float64)
    δ = Dict{Tuple{Int,Int},Float64}()
    dimension = length(Lattice)
    if dimension==2
        for i in 1:Lattice[1]
            for j in 1:Lattice[2]
                # δ[((i,j),(i,j))] = randn()*sqrt(var) + μ
                δ[(i,j)] = randn()*sqrt(var) + μ
            end
        end
        return δ
    else
        error("Dimension not implemented")
    end
end

function random_horizontal_δ(Lattice::AbstractArray{<:Int}, μ::Float64, var::Float64)
    δ = Dict{Tuple{Int,Int},Float64}()
    dimension = length(Lattice)
    if dimension==2
        for i in 1:Lattice[1]
            tmp = randn()*sqrt(var) + μ
            for j in 1:Lattice[2]
                δ[(i,j)] = tmp
            end
        end
        return δ
    else
        error("Dimension not implemented")
    end
end

function random_vertical_δ(Lattice::AbstractArray{<:Int}, μ::Float64, var::Float64)
    δ = Dict{Tuple{Int,Int},Float64}()
    dimension = length(Lattice)
    if dimension==2
        for j in 1:Lattice[2]
            tmp = randn()*sqrt(var) + μ
            for i in 1:Lattice[1]
                δ[(i,j)] = tmp
            end
        end
        return δ
    else
        error("Dimension not implemented")
    end
end

function random_θ(Lattice::AbstractArray{<:Int}, min::Float64, max::Float64)
    δ = Dict{Tuple{Int,Int},Float64}()
    dimension = length(Lattice)
    if dimension==2
        for i in 1:Lattice[1]
            for j in 1:Lattice[2]
                # δ[((i,j),(i,j))] = randn()*sqrt(var) + μ
                δ[(i,j)] = randn()*(max-min) + min
            end
        end
        return δ
    else
        error("Dimension not implemented")
    end
end

# """
# Input: 
#     - Lattice: AbstractArray{<:AbstractFloat}
#     - N: Int
#     - l: Int Distance of the pairing operator
# """
# function PairingOperator(Lattice::AbstractArray{<:Int}, N::Int, I::Int, J::Int, l::Int)
#     dimension = length(Lattice)
#     if dimension == 1
#         error("Dimension not implemented")
#     elseif dimension == 2
#         rows = Lattice[1]
#         cols = Lattice[2]
#         o = Op_fixed(2*rows*cols, N)
#         pairing = ada(o,Cartesian2Index([I,J],[rows,cols],1),Cartesian2Index([I,J+l],[rows,cols],1))*ada(o,Cartesian2Index([I+1,J],[rows,cols],2),Cartesian2Index([I+1,J+l],[rows,cols],2))
#         return pairing
#     else
#         error("Dimension not implemented")
#     end
# end
function global_spin_X_rot(
    Lattice::AbstractArray{<:Int},
    N::Int,
)
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        op_size = 2*rows*cols
        o = Op_fixed(op_size,N); # fixed particle number subspace
        Ham = spzeros(Complex{Float64},binomial(2*rows*cols,N),binomial(2*rows*cols,N))
        print("o dim: ", o.dim,"\n")
        for i in 1:rows
            for j in 1:cols
                Ham += (ada(o,Cartesian2Index([i,j],[rows,cols],1),Cartesian2Index([i,j],[rows,cols],2))+ada(o,Cartesian2Index([i,j],[rows,cols],2),Cartesian2Index([i,j],[rows,cols],1)))
            end
        end
        return Ham
    else
        error("Dimension not implemented")
    end
end

function local_spin_X_rot(
    Lattice::AbstractArray{<:Int},
    N::Int,
    θ::Dict{Tuple{Int,Int},Float64},
)
    dimension = length(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = Lattice[1]
        cols = Lattice[2]
        op_size = 2*rows*cols
        o = Op_fixed(op_size,N); # fixed particle number subspace
        Ham = spzeros(Complex{Float64},binomial(2*rows*cols,N),binomial(2*rows*cols,N))
        print("o dim: ", o.dim,"\n")
        for i in 1:rows
            for j in 1:cols
                Ham += θ[(i,j)]*(ada(o,Cartesian2Index([i,j],[rows,cols],1),Cartesian2Index([i,j],[rows,cols],2))+ada(o,Cartesian2Index([i,j],[rows,cols],2),Cartesian2Index([i,j],[rows,cols],1)))
            end
        end
        return Ham
    else
        error("Dimension not implemented")
    end
end

