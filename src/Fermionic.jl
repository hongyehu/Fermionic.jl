module Fermionic
using SparseArrays
using LinearAlgebra

include("base_functions.jl")
include("operators.jl")
include("states.jl")
include("logic_gates.jl")

export  Op, dim, basis, cm, cdm, cdcm, cmcd, cmcm, cdcd, vacuum
export State, State_sparse, State_complex, State_sparse_complex, st, ope, rhosp, eigensp, ssp, rhoqsp
export sigma_x, sigma_y, sigma_z, phase, hadamard, ucnot, swap
end # module
