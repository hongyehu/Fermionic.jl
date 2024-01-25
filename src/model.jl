using Fermionic
using LinearAlgebra
using SparseArrays
using Arpack
using PythonCall
plt = pyimport("matplotlib.pyplot")
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
Cartesian2Index([1,1],[2,2],1)
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
    J::Dict{Any,Any}, 
    U::Float64, 
    δ::Dict{Any,Any},
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
                Ham += δ[((i,j),(i,j))]*(
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


# function SpinFullRandomState(
#     Lattice::AbstractArray{<:Int64}, 
#     J::Dict{Any,Any}, 
#     U::Float64, 
#     δ::Dict{Any,Any}
#     )


rows = 2
cols = 3
J = Dict{Any,Any}()
for j in 1:cols-1
    J[((1,j),(1,j+1))] = 3.0
    J[((2,j),(2,j+1))] = 3.0
end
for j in 1:cols
    J[((1,j),(2,j))] = 2.0
end
J

δ = Dict{Any,Any}()
for i in 1:rows
    for j in 1:cols
        δ[((i,j),(i,j))] = 0.0
    end
end
Ham = SpinFullFermiHubbard([rows,cols],J,7.2,δ)
# Matrix(Ham)
vals, vecs = eigs(Ham, nev=30, which=:SR)
# vals2,vecs2 = eigen(Matrix(Ham))
vals[1:4]
# vals2[1:4]
plt.plot(vals, "o")
plt.show()
typeof(Ham)
Ham_sub = SpinFullFermiHubbardSubspace([rows,cols],J,7.2,δ,6)
vals, vecs = eigs(Ham_sub, nev=30, which=:SR)
plt.plot(vals, "o")
plt.show()
b = Op_fixed(12,6)
Matrix(b.basis)
vals[1:5]