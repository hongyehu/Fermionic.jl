using Fermionic
using LinearAlgebra
using SparseArrays
# Matrix(basis(Op(3)))
# Matrix(ad(Op(3),2))
# Matrix(ada(Op(3),1,2))
# ada(Op(3),1,2)
# d = 3; n = 2;
# b, index = basis_m(d,n);
# Matrix(b)
# index
# basis(Op(3))
# print(fixed_state(State([0,0,1,0],Op(2)),1))
# Array(basis_m(3,1)[1])
# Array(basis_m(4,2)[2])
# ad(3,1,1)

s = State([0.5,0,1,0.5],Op(2))
typeof(s)


sum([ada(Op(3),i,i) for i in 1:3])

# Implement coordinate transformation
function Cartesian2Index(C_coor::AbstractArray{<:AbstractFloat})
    """
    Row first ordering
    """
    
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
function SpinFullFermiHubbard(Lattice::AbstractArray{<:AbstractFloat}, J::Float64, U::Float64, δ::AbstractArray{<:AbstractFloat})
    # function body
    dimension = ndims(Lattice)
    if dimension == 1
        pass
    elseif dimension == 2
        rows = size(Lattice)[1]
        cols = size(Lattice)[2]
    else
        error("Dimension not implemented")
    end
end