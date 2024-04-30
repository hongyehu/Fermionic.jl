module Fermionic
using SparseArrays
using LinearAlgebra

include("base_functions.jl")
include("operators.jl")
include("operators_fixed.jl")
include("states.jl")
include("states_fixed.jl")
include("correlations.jl")
include("logic_gates.jl")
include("mixed.jl")


export Op, dim, basis, a, ad, ada, aad, aa, adad, adf, af, vacuum
export fixed, basis_m,  Op_fixed
export State, st, ope, typ, rhosp, rhoqsp, non_zero
export State_fixed, nume, fixed_state, unfixed_state
export eigensp, ssp, eigenqsp, sqsp, majorization_sp, majorization_qsp, n_avg, coef, rhom, rhomd, rhom2, trp
export sigma_x, sigma_y, sigma_z, phase, hadamard, ucnot, swap
export rhosp_mixed, eigensp_mixed

include("model.jl")
export Cartesian2Index, SpinFullFermiHubbard,SpinFullFermiHubbardSubspace,SpinFullRandomState
export SpinFullRandomStateSubspace,GHZSubspace,measure!
export Hopping, random_δ, random_θ, random_horizontal_δ, random_vertical_δ
export Dwave, global_spin_X_rot,local_spin_X_rot
export transformed_Dwave, transformed_SpinFullFermiHubbardSubspace

include("shadow.jl")
export  ShadowBasis, ShadowBasis_Global_Xrot,ShadowBasis_Local_Xrot

include("BCS.jl")
export createRowSwapMatrix

end # module
