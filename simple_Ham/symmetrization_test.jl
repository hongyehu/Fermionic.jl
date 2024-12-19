using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using Statistics
using ArgParse
using JSON3
using JLD2
using PythonCall

L = 2
type = :s
G = BCS_G(L,type,Δ=0.1,μ=0.5)
ρ = RDM(1.0.*G,[[1,1],[2,1],[1,2],[2,2]])
o = Op(8)
op1 = ad(o,1)*ad(o,2)
op2 = ad(o,3)*ad(o,4)
op1_sym = (ad(o,1)*ad(o,2)-ad(o,2)*ad(o,1))+(ad(o,1)*ad(o,2)-ad(o,2)*ad(o,1))'
op2_sym = (ad(o,3)*ad(o,4)-ad(o,4)*ad(o,3))+(ad(o,3)*ad(o,4)-ad(o,4)*ad(o,3))'
tr(ρ*op2_sym)/tr(ρ*op2)

tr(ρ*op1_sym*op2_sym)

tr(ρ*op1*op2')

Cartesian2Index([2,1],[2,2],2)