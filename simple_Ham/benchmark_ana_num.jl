using Fermionic
using PythonCall
using LinearAlgebra
using Expokit
plt = pyimport("matplotlib.pyplot")
include("./analytical_bcs.jl")
include("./numerical_bcs.jl")

L = 13
analytical_four_fermion(0,-1,0,1,0,0,L,:d)
numerical_four_fermion(0,-1,0,1,0,0,L,:d)

"""
The following function calculate true pairing correlation using analytical result.
"""
function d_pairing_correlation(l::Int64, L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = 0.0+0.0im
    res += analytical_four_fermion(1,0,l%L,0,(l+1)%L,0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(1,0,l%L,0,(l-1)%L,0,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion(1,0,l%L,0,l%L,1,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion(1,0,l%L,0,l%L,(-1)%L,L,type;μ=μ,Δ=Δ)
    ##
    res += analytical_four_fermion((-1)%L,0,l%L,0,(l+1)%L,0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion((-1)%L,0,l%L,0,(l-1)%L,0,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion((-1)%L,0,l%L,0,l%L,1,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion((-1)%L,0,l%L,0,l%L,(-1)%L,L,type;μ=μ,Δ=Δ)
    ##
    res -= analytical_four_fermion(0,1,l%L,0,(l+1)%L,0,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion(0,1,l%L,0,(l-1)%L,0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,1,l%L,0,l%L,1,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,1,l%L,0,l%L,(-1)%L,L,type;μ=μ,Δ=Δ)
    ##
    res -= analytical_four_fermion(0,(-1)%L,l%L,0,(l+1)%L,0,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion(0,(-1)%L,l%L,0,(l-1)%L,0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,(-1)%L,l%L,0,l%L,1,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,(-1)%L,l%L,0,l%L,(-1)%L,L,type;μ=μ,Δ=Δ)
    return -res
end

function num_d_pairing_correlation(l::Int64, L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = 0.0+0.0im
    res += numerical_four_fermion(1,0,mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res += numerical_four_fermion(1,0,mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res -= numerical_four_fermion(1,0,mod(l,L),0,mod(l,L),1,L,type;μ=μ,Δ=Δ)
    res -= numerical_four_fermion(1,0,mod(l,L),0,mod(l,L),mod(-1,L),L,type;μ=μ,Δ=Δ)
    ##
    res += numerical_four_fermion(mod(-1,L),0,mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res += numerical_four_fermion(mod(-1,L),0,mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res -= numerical_four_fermion(mod(-1,L),0,mod(l,L),0,mod(l,L),1,L,type;μ=μ,Δ=Δ)
    res -= numerical_four_fermion(mod(-1,L),0,mod(l,L),0,mod(l,L),mod(-1,L),L,type;μ=μ,Δ=Δ)
    ##
    res -= numerical_four_fermion(0,1,mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res -= numerical_four_fermion(0,1,mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res += numerical_four_fermion(0,1,mod(l,L),0,mod(l,L),1,L,type;μ=μ,Δ=Δ)
    res += numerical_four_fermion(0,1,mod(l,L),0,mod(l,L),mod(-1,L),L,type;μ=μ,Δ=Δ)
    ##
    res -= numerical_four_fermion(0,mod(-1,L),mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res -= numerical_four_fermion(0,mod(-1,L),mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res += numerical_four_fermion(0,mod(-1,L),mod(l,L),0,mod(l,L),1,L,type;μ=μ,Δ=Δ)
    res += numerical_four_fermion(0,mod(-1,L),mod(l,L),0,mod(l,L),mod(-1,L),L,type;μ=μ,Δ=Δ)
    return -res
end
mod(-1,L)
L = 17
d_correlation_4_swave = [d_pairing_correlation(l,L,:s) for l in 1:Int((L-1)/2)]
d_correlation_4_dwave = [d_pairing_correlation(l,L,:d) for l in 1:Int((L-1)/2)]
num_d_correlation_4_swave = [num_d_pairing_correlation(l,L,:s) for l in 1:Int((L-1)/2)]
num_d_correlation_4_dwave = [num_d_pairing_correlation(l,L,:d) for l in 1:Int((L-1)/2)]

begin
    plt.figure(figsize=(10,5))
    plt.plot(1:Int((L-1)/2), real(d_correlation_4_swave), label="Analytical d-wave")
    plt.plot(1:Int((L-1)/2), real(d_correlation_4_dwave), label="Analytical d-wave")
    plt.plot(1:Int((L-1)/2), real(num_d_correlation_4_swave), label="Numerical d-wave")
    plt.plot(1:Int((L-1)/2), real(num_d_correlation_4_dwave), label="Numerical d-wave")
    plt.legend()
    plt.show()
end

BCS_G(3, :d;μ=0.5,Δ=5.0)

# Define a square matrix
A = [1.0 0.5; 0.5 1.0]

G = getBCSaij(2, :d;μ=0.5,Δ=5.0)
begin
    plt.imshow(real(G)) 
    plt.show()
end