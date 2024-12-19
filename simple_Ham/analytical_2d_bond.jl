using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using Statistics
using ArgParse
using JSON3
using JLD2
using PythonCall

# L=24
# type = :d
# G_file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/G/L$(L)_type1_Delta0.1_mu0.5.jld2"
# @load G_file_name G

# reduced_positions = [[1,1],[2,1],[1,2],[2,2]]
# ρ1 = RDM(2.0.*G,reduced_positions)
# o = Op(8)
# O1 = ad(o,1)*ad(o,6)-ad(o,2)*ad(o,5)+(ad(o,1)*ad(o,6)-ad(o,2)*ad(o,5))'
# O2 = ad(o,3)*ad(o,8)-ad(o,4)*ad(o,7)+(ad(o,3)*ad(o,8)-ad(o,4)*ad(o,7))'
# tr(ρ1*O1*O2)

# reduced_positions = [[1,1],[1,2],[3,2],[4,2]]
# ρ1 = RDM(2.0.*G,reduced_positions)
# o = Op(8)
# O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)+(ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3))'
# O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)+(ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7))'
# tr(ρ1*O1*O2)


# reduced_positions = [[1,1],[1,2],[1,3],[2,3]]
# ρ2 = RDM(2.0.*G,reduced_positions)
# O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)+(ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3))'
# O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)+(ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7))'
# tr(ρ2*O1*O2)

# reduced_positions = [[1,1],[1,2],[1,4],[2,4]]
# ρ2 = RDM(2.0.*G,reduced_positions)
# O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)+(ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3))'
# O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)+(ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7))'
# tr(ρ2*O1*O2)

# reduced_positions = [[1,1],[1,2],[1,3+4*3],[2,3+4*3]]
# ρ2 = RDM(2.0.*G,reduced_positions)
# O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)+(ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3))'
# O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)+(ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7))'
# tr(ρ2*O1*O2)

function bond_bond_correlation(G, super_x,super_y,index)
    super_x = super_x-1
    super_y = super_y-1
    if super_y>0
        reduced_positions = [[1,1],[1,2]]
        if index == 1
            append!(reduced_positions,[[4*super_x+1,4*super_y+1],[4*super_x+1,4*super_y+2]])
        elseif index == 2
            append!(reduced_positions,[[4*super_x+2,4*super_y+1],[4*super_x+2,4*super_y+2]])
        elseif index == 3
            append!(reduced_positions,[[4*super_x+3,4*super_y+1],[4*super_x+4,4*super_y+1]])
        elseif index == 4
            append!(reduced_positions,[[4*super_x+3,4*super_y+2],[4*super_x+4,4*super_y+2]])
        elseif index == 5
            append!(reduced_positions,[[4*super_x+1,4*super_y+3],[4*super_x+2,4*super_y+3]])
        elseif index == 6
            append!(reduced_positions,[[4*super_x+1,4*super_y+4],[4*super_x+2,4*super_y+4]])
        elseif index == 7
            append!(reduced_positions,[[4*super_x+3,4*super_y+3],[4*super_x+3,4*super_y+4]])
        elseif index == 8
            append!(reduced_positions,[[4*super_x+4,4*super_y+3],[4*super_x+4,4*super_y+4]])
        end
        ρ = RDM(2.0.*G,reduced_positions)
        o = Op(8)
        O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)+a(o,4)*a(o,1)-a(o,3)*a(o,2)
        O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)+a(o,8)*a(o,5)-a(o,7)*a(o,6)
        obs = O1*O2
        return tr(ρ*obs)
    else # in the same column
        if index == 5
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)+a(o,4)*a(o,1)-a(o,3)*a(o,2)
            O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)+a(o,8)*a(o,5)-a(o,7)*a(o,6)
            append!(reduced_positions,[[4*super_x+1,4*super_y+3],[4*super_x+2,4*super_y+3]])
        elseif index == 6
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)+a(o,4)*a(o,1)-a(o,3)*a(o,2)
            O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)+a(o,8)*a(o,5)-a(o,7)*a(o,6)
            append!(reduced_positions,[[4*super_x+1,4*super_y+4],[4*super_x+2,4*super_y+4]])
        elseif index == 7
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)+a(o,4)*a(o,1)-a(o,3)*a(o,2)
            O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)+a(o,8)*a(o,5)-a(o,7)*a(o,6)
            append!(reduced_positions,[[4*super_x+3,4*super_y+3],[4*super_x+3,4*super_y+4]])
        elseif index == 8
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)+a(o,4)*a(o,1)-a(o,3)*a(o,2)
            O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)+a(o,8)*a(o,5)-a(o,7)*a(o,6)
            append!(reduced_positions,[[4*super_x+4,4*super_y+3],[4*super_x+4,4*super_y+4]])
        elseif index == 1
            if super_x==0
                return 0.0
            else
                reduced_positions = [[1,1],[4*super_x+1,4*super_y+1],[1,2],[4*super_x+1,4*super_y+2]]
                o = Op(8)
                O1 = ad(o,1)*ad(o,6)-ad(o,2)*ad(o,5)+(ad(o,1)*ad(o,6)-ad(o,2)*ad(o,5))'
                O2 = ad(o,3)*ad(o,8)-ad(o,4)*ad(o,7)+(ad(o,3)*ad(o,8)-ad(o,4)*ad(o,7))'
            end
        elseif index == 2
            reduced_positions = [[1,1],[4*super_x+2,4*super_y+1],[1,2],[4*super_x+2,4*super_y+2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,6)-ad(o,2)*ad(o,5)+(ad(o,1)*ad(o,6)-ad(o,2)*ad(o,5))'
            O2 = ad(o,3)*ad(o,8)-ad(o,4)*ad(o,7)+(ad(o,3)*ad(o,8)-ad(o,4)*ad(o,7))'
        elseif index == 3
            reduced_positions = [[1,1],[4*super_x+3,4*super_y+1],[4*super_x+4,4*super_y+1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,8)-ad(o,2)*ad(o,7)+(ad(o,1)*ad(o,8)-ad(o,2)*ad(o,7))'
            O2 = ad(o,3)*ad(o,6)-ad(o,4)*ad(o,5)+(ad(o,3)*ad(o,6)-ad(o,4)*ad(o,5))'
        elseif index == 4
            reduced_positions = [[1,1],[1,2],[4*super_x+3,4*super_y+2],[4*super_x+4,4*super_y+2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)+(ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3))'
            O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)+(ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7))'
        end
        ρ = RDM(2.0.*G,reduced_positions)
        obs = O1*O2
        return tr(ρ*obs)
    end
end

function no_sym_bond_bond_correlation(G, super_x,super_y,index)
    super_x = super_x-1
    super_y = super_y-1
    if super_y>0
        reduced_positions = [[1,1],[1,2]]
        if index == 1
            append!(reduced_positions,[[4*super_x+1,4*super_y+1],[4*super_x+1,4*super_y+2]])
        elseif index == 2
            append!(reduced_positions,[[4*super_x+2,4*super_y+1],[4*super_x+2,4*super_y+2]])
        elseif index == 3
            append!(reduced_positions,[[4*super_x+3,4*super_y+1],[4*super_x+4,4*super_y+1]])
        elseif index == 4
            append!(reduced_positions,[[4*super_x+3,4*super_y+2],[4*super_x+4,4*super_y+2]])
        elseif index == 5
            append!(reduced_positions,[[4*super_x+1,4*super_y+3],[4*super_x+2,4*super_y+3]])
        elseif index == 6
            append!(reduced_positions,[[4*super_x+1,4*super_y+4],[4*super_x+2,4*super_y+4]])
        elseif index == 7
            append!(reduced_positions,[[4*super_x+3,4*super_y+3],[4*super_x+3,4*super_y+4]])
        elseif index == 8
            append!(reduced_positions,[[4*super_x+4,4*super_y+3],[4*super_x+4,4*super_y+4]])
        end
        ρ = RDM(2.0.*G,reduced_positions)
        o = Op(8)
        O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)
        O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)
        obs = O1*O2
        return tr(ρ*obs)
    else # in the same column
        if index == 5
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)
            O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)
            append!(reduced_positions,[[4*super_x+1,4*super_y+3],[4*super_x+2,4*super_y+3]])
        elseif index == 6
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)
            O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)
            append!(reduced_positions,[[4*super_x+1,4*super_y+4],[4*super_x+2,4*super_y+4]])
        elseif index == 7
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)
            O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)
            append!(reduced_positions,[[4*super_x+3,4*super_y+3],[4*super_x+3,4*super_y+4]])
        elseif index == 8
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)
            O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)
            append!(reduced_positions,[[4*super_x+4,4*super_y+3],[4*super_x+4,4*super_y+4]])
        elseif index == 1
            if super_x==0
                return 0.0
            else
                reduced_positions = [[1,1],[4*super_x+1,4*super_y+1],[1,2],[4*super_x+1,4*super_y+2]]
                o = Op(8)
                O1 = ad(o,1)*ad(o,6)-ad(o,2)*ad(o,5)
                O2 = ad(o,3)*ad(o,8)-ad(o,4)*ad(o,7)
            end
        elseif index == 2
            reduced_positions = [[1,1],[4*super_x+2,4*super_y+1],[1,2],[4*super_x+2,4*super_y+2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,6)-ad(o,2)*ad(o,5)
            O2 = ad(o,3)*ad(o,8)-ad(o,4)*ad(o,7)
        elseif index == 3
            reduced_positions = [[1,1],[4*super_x+3,4*super_y+1],[4*super_x+4,4*super_y+1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,8)-ad(o,2)*ad(o,7)
            O2 = ad(o,3)*ad(o,6)-ad(o,4)*ad(o,5)
        elseif index == 4
            reduced_positions = [[1,1],[1,2],[4*super_x+3,4*super_y+2],[4*super_x+4,4*super_y+2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)
            O2 = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7)
        end
        ρ = RDM(2.0.*G,reduced_positions)
        obs = O1*O2
        return tr(ρ*obs)
    end
end

function bare_bond_bond_correlation(G, super_x,super_y,index)
    super_x = super_x-1
    super_y = super_y-1
    if super_y>0
        reduced_positions = [[1,1],[1,2]]
        if index == 1
            append!(reduced_positions,[[4*super_x+1,4*super_y+1],[4*super_x+1,4*super_y+2]])
        elseif index == 2
            append!(reduced_positions,[[4*super_x+2,4*super_y+1],[4*super_x+2,4*super_y+2]])
        elseif index == 3
            append!(reduced_positions,[[4*super_x+3,4*super_y+1],[4*super_x+4,4*super_y+1]])
        elseif index == 4
            append!(reduced_positions,[[4*super_x+3,4*super_y+2],[4*super_x+4,4*super_y+2]])
        elseif index == 5
            append!(reduced_positions,[[4*super_x+1,4*super_y+3],[4*super_x+2,4*super_y+3]])
        elseif index == 6
            append!(reduced_positions,[[4*super_x+1,4*super_y+4],[4*super_x+2,4*super_y+4]])
        elseif index == 7
            append!(reduced_positions,[[4*super_x+3,4*super_y+3],[4*super_x+3,4*super_y+4]])
        elseif index == 8
            append!(reduced_positions,[[4*super_x+4,4*super_y+3],[4*super_x+4,4*super_y+4]])
        end
        ρ = RDM(2.0.*G,reduced_positions)
        o = Op(8)
        O1 = ad(o,1)*ad(o,4)
        O2 = ad(o,5)*ad(o,8)
        obs = O1*O2
        return tr(ρ*obs)
    else # in the same column
        if index == 5
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)
            O2 = ad(o,5)*ad(o,8)
            append!(reduced_positions,[[4*super_x+1,4*super_y+3],[4*super_x+2,4*super_y+3]])
        elseif index == 6
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)
            O2 = ad(o,5)*ad(o,8)
            append!(reduced_positions,[[4*super_x+1,4*super_y+4],[4*super_x+2,4*super_y+4]])
        elseif index == 7
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)
            O2 = ad(o,5)*ad(o,8)
            append!(reduced_positions,[[4*super_x+3,4*super_y+3],[4*super_x+3,4*super_y+4]])
        elseif index == 8
            reduced_positions = [[1,1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)
            O2 = ad(o,5)*ad(o,8)
            append!(reduced_positions,[[4*super_x+4,4*super_y+3],[4*super_x+4,4*super_y+4]])
        elseif index == 1
            if super_x==0
                return 0.0
            else
                reduced_positions = [[1,1],[4*super_x+1,4*super_y+1],[1,2],[4*super_x+1,4*super_y+2]]
                o = Op(8)
                O1 = ad(o,1)*ad(o,6)
                O2 = ad(o,3)*ad(o,8)
            end
        elseif index == 2
            reduced_positions = [[1,1],[4*super_x+2,4*super_y+1],[1,2],[4*super_x+2,4*super_y+2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,6)
            O2 = ad(o,3)*ad(o,8)
        elseif index == 3
            reduced_positions = [[1,1],[4*super_x+3,4*super_y+1],[4*super_x+4,4*super_y+1],[1,2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,8)
            O2 = ad(o,3)*ad(o,6)
        elseif index == 4
            reduced_positions = [[1,1],[1,2],[4*super_x+3,4*super_y+2],[4*super_x+4,4*super_y+2]]
            o = Op(8)
            O1 = ad(o,1)*ad(o,4)
            O2 = ad(o,5)*ad(o,8)
        end
        ρ = RDM(2.0.*G,reduced_positions)
        obs = O1*O2
        return tr(ρ*obs)
    end
end
G = BCS_G(24,:d,Δ=0.1,μ=0.5)
bond_bond_correlation(G, 1,1,1)
bond_bond_correlation(G, 1,1,2)
bond_bond_correlation(G, 1,1,3)
real(bond_bond_correlation(G, 3,6,3))/real(no_sym_bond_bond_correlation(G, 3,6,3))

begin
    indexI = 2
    indexJ = 5
    bond_ind = 8
    println("factor: ",real(bond_bond_correlation(G, indexI,indexJ,bond_ind))/real(bare_bond_bond_correlation(G, indexI,indexJ,bond_ind)))
    println("value: ", real(bare_bond_bond_correlation(G, indexI,indexJ,bond_ind)))
end



# bond_bond_correlation(G, 1,1,4)
# bond_bond_correlation(G, 1,1,5)
# bond_bond_correlation(G, 1,1,6)
# bond_bond_correlation(G, 1,1,7)
# bond_bond_correlation(G, 1,1,8)

