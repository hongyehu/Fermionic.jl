module Fermionic
using SparseArrays
using LinearAlgebra

include("base_functions.jl")
include("operators.jl")
include("states.jl")
include("correlations.jl")
include("logic_gates.jl")
include("operators_fixed.jl")
include("mixed.jl")

export Op, dim, basis, cm, cdm, cdcm, cmcd, cmcm, cdcd, vacuum
export State, State2, st, ope, typ, rhosp, rhoqsp
export eigensp, ssp, eigenqsp, sqsp, majorization_sp, majorization_qsp, n_avg, rhom, rhomd, trp
export sigma_x, sigma_y, sigma_z, phase, hadamard, ucnot, swap
export fixed, basis_m, fixed_state, unfixed_state, cdc, ccd
export rhosp_mixed, eigensp_mixed
end # module
